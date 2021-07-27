#=

Similar to size_scaling_binning.jl

=#

using LinearAlgebra
using Base.Threads
using HDF5
using Statistics
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

function collect_samples!(sc!)
    v = randomstate(Ns, L)

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
        @views scalarchirality(vs[:, :, :, 1], sc![:, :, :, nsample])
        
        # use the last time evolved state
        @views v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples")
        flush(stdout)
    end
    nothing
end

# Params
# ======

@assert length(ARGS) == 1
const p = Dict("comment" => "sc, QSL",
               "J1" => 1,
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
output = "../data/scalar_chirality/sc_qsl_$(ARGS[1]).h5"
const H = loadhamiltonian("../hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const Ns = H.Ns
const L = p["L"]
const T = p["T"]
const nsamples = p["nsamples"]
const stride = p["stride"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Running the simulation
# ======================

# use a single chain for simplicity
sc = zeros(2, L, L, nsamples)

collect_samples!(sc)

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["sc"] = sc
end

println("Finished !")
