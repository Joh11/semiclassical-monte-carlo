#=

A script to compute the dynamic structure factor for
square kagome lattice

This is for J1=J2=J3

Each chain is indepentent, and stores the result as a temporary
average. Afterwards, an average is made over all the chains (OK since
they each have the same number of samples). 

=#

using LinearAlgebra
using Base.Threads
using HDF5
using Statistics
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

function collect_samples!(corr, fac; chain=nothing)
    v = randomstate(Ns, L)
    L½ = 1 + L÷2 # like that L=4 => L½= 3

    corr_tmp = zeros(Ns, L, L, Ns, L, L)
    # just a trick to use structurefactor_kpath!
    fac_tmp = zeros(ℂ, 1)
    kpath_tmp = reshape([4π, 0], (2, 1))
    
    println("Starting for chain $chain")
    # thermalization step
    mcstep!(H, v, T, p["thermal"])

    # preallocation for DifferentialEquations.jl
    # (see https://discourse.julialang.org/t/using-differentialequations-jl-for-very-large-system/1030/2 )
    timeseries = []
    ts = []
    ks = []
    
    # sampling step
    for nsample = 1:nsamples_per_chain
        println("Doing sample $nsample / $nsamples_per_chain for chain $chain")
        mcstep!(H, v, T, stride)

        # time evolution
        # we still have to time evolve to improve the sampling
        vs = simulate(H, v, p["dt"], nt, timeseries, ts, ks)
        
        # save measurements of interest
        @views allcorrelations!(vs[:, :, :, 1], vs[:, :, :, 1],
                                Ns, L, corr_tmp)
        
        # now compute the measurements of interest
        corr[nsample] = mean([corr[i, 1, 1, i, L½, L½, 1] for i = 1:Ns])

        # S(4π, 0)
        structurefactor_kpath!(Rs, corr_tmp, kpath_tmp, fac_tmp)
        fac[sample] = abs.(fac_tmp[1])

        # use the last time evolved state
        @views v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples_per_chain (chain $chain)")
        flush(stdout)

        # trigger the garbage collector manually
        if chain == 1 && nsample % 50 == 0
            println("Running GC manually ...")
            GC.gc()
        end
    end
end

function kpath2mat(kpath)
    mat = zeros(2, length(kpath))
    for (n, k) in enumerate(kpath)
        mat[:, n] .= k
    end
    mat
end

# Params
# ======

@assert length(ARGS) == 1
const p = Dict("comment" => "finite size scaling analysis, UUD",
               "J1" => 0,
               "J2" => 1,
               "J3" => 1,
               "L" => parse(Int, ARGS[1]),
               "T" => 0.1,
               "thermal" => 100_000,
               "nchains" => 16, # because 8 cores on my laptop
               "nsamples_per_chain" => 4_000, # so ~ 30k samples
               "stride" => 100,
               # time evolution params
               "dt" => 100,
               "nt" => 2
               )
output = "../data/scaling/scaling_uud_$(ARGS[1]).h5"
const H = loadhamiltonian("../hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const Ns = H.Ns
const L = p["L"]
const nchains = p["nchains"]
const T = p["T"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Precompute the position of each site
const Rs = compute_positions(H, L)

# Running the simulation
# ======================

# storing results per chain should be thread safe

# in this case we don't have much to store, so we store every sample
# to do binning analysis
corr = zeros(nsamples_per_chain, nchains)
fac = zeros(nsamples_per_chain, nchains)

Threads.@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    @views collect_samples!(corr[:, n], fac[:, n]; chain=n)
end

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["corr"] = corr
    f["fac"] = fac
    f["E"] = E
end

println("Finished !")
