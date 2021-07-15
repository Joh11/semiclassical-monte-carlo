#=

TODO this is not working yet

=#

using LinearAlgebra
using Base.Threads
using HDF5
using Statistics
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

function collect_samples!(; chain=nothing)
    v = randomstate(Ns, L)
    L½ = 1 + L÷2 # like that L=4 => L½= 3

    corr = 0
    fac = 0
    
    corr_tmp = zeros(Ns, L, L, Ns, L, L)
    # just a trick to use structurefactor_kpath!
    fac_tmp = zeros(ℂ, 1)
    kpath_tmp = [[4π, 0]]
    
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
        corr += mean([corr_tmp[i, 1, 1, i, L½, L½] for i = 1:Ns])

        # S(4π, 0)
        structurefactor_kpath!(Rs, corr_tmp, kpath_tmp, fac_tmp)
        fac += abs.(fac_tmp[1])

        # use the last time evolved state
        @views v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples_per_chain (chain $chain)")
        flush(stdout)

        # trigger the garbage collector manually
        if chain == 1 && nsample % 1_000 == 0
            println("Running GC manually ...")
            GC.gc()
        end
    end

    corr / nsamples_per_chain, fac / nsamples_per_chain
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

# do a "dumb" error estimation by simply taking the std over the chains
corr = zeros(nchains)
fac = zeros(nchains)

@threads for n in 1:nchains
    corr[n], fac[n] = collect_samples!(; chain=n)
end

Δcorr = std(corr)
Δfac = std(fac)

corr = mean(corr)
fac = mean(fac)

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

    f["Δcorr"] = Δcorr
    f["Δfac"] = Δfac
end

println("Finished !")
