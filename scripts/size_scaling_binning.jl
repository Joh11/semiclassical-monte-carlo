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

function collect_samples!(corr, fac)
    v = randomstate(Ns, L)
    L½ = 1 + L÷2 # like that L=4 => L½= 3

    corr_tmp = zeros(Ns, L, L, Ns, L, L)
    chi_tmp = zeros(ℂ, 1)
    
    # thermalization step
    mcstep!(H, v, T, p["thermal"])

    # preallocation for DifferentialEquations.jl
    # (see https://discourse.julialang.org/t/using-differentialequations-jl-for-very-large-system/1030/2 )
    timeseries = []
    ts = []
    ks = []
    
    # sampling step
    for nsample = 1:nsamples
        println("Doing sample $nsample / $nsamples")
        mcstep!(H, v, T, stride)

        # time evolution
        # we still have to time evolve to improve the sampling
        vs = simulate(H, v, p["dt"], nt, timeseries, ts, ks)
        
        # save measurements of interest
        @views allcorrelations!(vs[:, :, :, 1], vs[:, :, :, 1],
                                Ns, L, corr_tmp)
        
        # now compute the measurements of interest
        corr[nsample] += mean([corr_tmp[i, 1, 1, i, L½, L½] for i = 1:Ns])

        # S(4π, 0)
        fill!(chi_tmp, 0)
        structurefactor_kpath!(Rs, corr_tmp, [[4π, 0]], chi_tmp)
        fac[nsample] = abs(chi_tmp[1])

        # use the last time evolved state
        @views v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples")
        flush(stdout)

        # trigger the garbage collector manually
        # if nsample % 1_000 == 0
        #     println("Running GC manually ...")
        #     GC.gc()
        # end
    end
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
               "nsamples" => 2^15,
               "stride" => 50,
               # time evolution params
               "dt" => 100,
               "nt" => 2
               )
output = "../data/scaling/scaling_uud_$(ARGS[1]).h5"
const H = loadhamiltonian("../hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const Ns = H.Ns
const L = p["L"]
const T = p["T"]
const nsamples = p["nsamples"]
const stride = p["stride"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Precompute the position of each site
const Rs = compute_positions(H, L)

# Running the simulation
# ======================

# use a single chain for simplicity
corr = zeros(nsamples) # the <Si⋅Sj> correlation
fac = zeros(nsamples) # S(4π, 0)

collect_samples!(corr, fac)

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
end

println("Finished !")
