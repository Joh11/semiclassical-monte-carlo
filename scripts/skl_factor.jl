#=

A script to compute the structure factor (static and dynamical) for
square kagome lattice

Scan over the temperature for J1=J2=J3

Each chain is indepentent, and stores the result as a temporary
average. Afterwards, an average is made over all the chains (OK since
they each have the same number of samples). 

=#

using LinearAlgebra
using Base.Threads
using HDF5
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

# Utility functions
# =================

function logrange(x1, x2, n)
    [10^y for y in range(log10(x1), log10(x2), length=n)]
end

# Params
# ======

const p = Dict("comment" => "Small one to try",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 20,
               "T" => 0.01,
               "thermal" => 100_000,
               "nchains" => 8, # because 72 cores on Piz Daint
               "nsamples_per_chain" => 10000, # 420 => ~ 30K total samples
               "stride" => 100,
               # time evolution stuff
               "dt" => 0.1,
               "nt" => 1)
output = "skl_factor_fast.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const T = p["T"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const dt = p["dt"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Running the simulation
# ======================

total_corr = zeros(ℂ, H.Ns, L, L, H.Ns, L, L)
total_energy = 0

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    
    v = randomstate(H.Ns, L)
    corr = zeros(ℂ, H.Ns, L, L, H.Ns, L, L)
    E = 0

    println("Starting T = $T for chain $n / $nchains")
    # thermalization step
    mcstep!(H, v, T, p["thermal"])
    
    # sampling step
    for nsample = 1:nsamples_per_chain
        println("$nsample / $nsamples_per_chain (chain $n / $nchains)")
        flush(stdout)
        mcstep!(H, v, T, stride)
        
        # save the measurements of interest
        corr += allcorrelations(v)
        E += energy(H, v)
    end
    
    # now that everything is sampled, average
    corr /= nsamples_per_chain
    E /= nsamples_per_chain
    
    # and append to the total results (for all chains)
    global total_corr += corr
    global total_energy += E
end

# average over all chains
total_corr /= nchains
total_energy /= nchains

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["corr"] = total_corr
    f["E"] = total_energy
end

println("Finished !")
