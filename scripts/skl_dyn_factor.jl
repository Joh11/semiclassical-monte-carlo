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
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

# Params
# ======

const p = Dict("comment" => "Trying with the extended BZ",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 10,
               "T" => 0.01,
               "thermal" => 100_000,
               "nchains" => 8, # because 8 cores on my laptop
               "nsamples_per_chain" => 4_000, # so ~ 30k samples
               "stride" => 100,
               # time evolution params
               "dt" => 1,
               "nt" => 100)
output = "skl_dyn_factor.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const T = p["T"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Running the simulation
# ======================

total_corrs = zeros(H.Ns, L, L, H.Ns, L, L, nt)
total_energy = 0

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    
    v = randomstate(H.Ns, L)
    corrs = zeros(H.Ns, L, L, H.Ns, L, L, nt)
    E = 0

    println("Starting T = $T for chain $n / $nchains")
    # thermalization step
    mcstep!(H, v, T, p["thermal"])
    
    # sampling step
    for nsample = 1:nsamples_per_chain
        time = @elapsed begin
            mcstep!(H, v, T, stride)

            # time evolution
            vs = simulate(H, v, p["dt"], nt)
            
            # save the measurements of interest
            for t = 1:nt
                @views corrs[:, :, :, :, :, :, t] .+= allcorrelations(vs[:, :, :, t])
            end
            E += energy(H, v)

            # use the last time evolved state
            v = vs[:, :, :, end]
        end
        
        println("done $nsample / $nsamples_per_chain (chain $n / $nchains) in $time seconds")
        flush(stdout)
    end
    
    # now that everything is sampled, average
    corrs /= nsamples_per_chain
    E /= nsamples_per_chain
    
    # and append to the total results (for all chains)
    global total_corrs .+= corrs
    global total_energy += E
end

# average over all chains
total_corrs /= nchains
total_energy /= nchains

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["corrs"] = total_corrs
    f["E"] = total_energy
end

println("Finished !")
