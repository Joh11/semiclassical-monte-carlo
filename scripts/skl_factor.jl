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

const p = Dict("comment" => "First one, hope it's not too big",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 10,
               "Ts" => logrange(5e-3, 0.1, 10),
               "thermal_first" => 100_000,
               "thermal" => 10_000,
               "nchains" => 72, # because 72 cores on Piz Daint
               "nsamples_per_chain" => 420, # ~ 30K total samples
               "stride" => 500,
               # time evolution stuff
               "dt" => 0.1,
               "nt" => 100)
output = "skl_factor.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const Ts = p["Ts"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const dt = p["dt"]
const nt = p["nt"]

const Ns = 6L^2 # number of sites in total

# Running the simulation
# ======================

total_St = zeros(ℂ, L, L, nt, length(Ts))
total_Somega = zeros(ℂ, L, L, nt, length(Ts))
total_energy = zeros(length(Ts))

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    
    v = randomstate(H.Ns, L)
    St = zeros(ℂ, L, L, nt)
    Somega = zeros(ℂ, L, L, nt)
    E = 0

    St_temp = zeros(ℂ, L, L, nt) # tmp var

    for i in eachindex(Ts)
        T = Ts[i]
        println("Starting T = $T for chain $n / $nchains")
        # reset variables
        St .= 0
        Somega .= 0
        E = 0
        # thermalization step
        if i == 1 # hard the first time
            mcstep!(H, v, T, p["thermal_first"])
        else
            mcstep!(H, v, T, p["thermal"])
        end

        # sampling step
        for nsample = 1:nsamples_per_chain
            println("$nsample / $nsamples_per_chain, T = $T (chain $n / $nchains)")
            flush(stdout)
            mcstep!(H, v, T, stride)
            # save the measurements of interest
            vs = simulate(H, v, dt, nt)
            St_temp = structuralfactor(H, vs)
            St += St_temp
            Somega += frequencystructuralfactor(St_temp)            
            E += energy(H, v)
            # use the last time evolved state to improve sampling
            v = vs[:, :, :, end]
        end

        # now that everything is sampled, average
        St /= nsamples_per_chain
        Somega /= nsamples_per_chain
        E /= nsamples_per_chain

        # and append to the total results (for all chains)
        total_St[:, :, :, i] += St
        total_Somega[:, :, :, i] += Somega
        total_energy[i] += E
    end
end

# average over all chains
total_St /= nchains
total_Somega /= nchains
total_energy /= nchains



# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    for n in eachindex(Ts)
        T = Ts[n]
        group = create_group(f, "$n")
        group["T"] = T
        @views group["St"] = total_St[:, :, :, n]
        @views group["Somega"] = total_Somega[:, :, :, n]
        group["E"] = total_energy[n]
    end
end

println("Finished !")
