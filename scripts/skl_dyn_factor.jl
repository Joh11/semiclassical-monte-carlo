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

function collect_samples!(corr, E; chain=nothing)
    v = randomstate(Ns, L)

    corr_tmp = zeros(Ns, L, L, Ns, L, L, nt)
    
    println("Starting for chain $chain")
    # thermalization step
    mcstep!(H, v, T, p["thermal"])
    
    # sampling step
    for nsample = 1:nsamples_per_chain
        mcstep!(H, v, T, stride)

        # time evolution
        vs = simulate(H, v, p["dt"], nt)

        # save measurements of interest
        for t = 1:nt
            @views allcorrelations!(vs[:, :, :, 1], vs[:, :, :, t],
                                    Ns, L, corr_tmp[:, :, :, :, :, :, t])
        end
        corr .+= corr_tmp
        E += energy(H, v)

        # use the last time evolved state
        v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples_per_chain (chain $n / $nchains)")
        flush(stdout)
    end
end

# Params
# ======

const p = Dict("comment" => "Trying with the extended BZ",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 10,
               "T" => 0.01,
               "thermal" => 1,# 00_000,
               "nchains" => 8, # because 8 cores on my laptop
               "nsamples_per_chain" => 4_000,# 4_000, # so ~ 30k samples
               "stride" => 100,
               # time evolution params
               "dt" => 1,
               "nt" => 100)
output = "skl_dyn_factor.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const Ns = H.Ns
const L = p["L"]
const nchains = p["nchains"]
const T = p["T"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const nt = p["nt"]

const Nsites = H.Ns * L^2 # number of sites in total

# Running the simulation
# ======================

# storing results per chain should be thread safe
corr = zeros(Ns, L, L, Ns, L, L, nt, nchains)
E = zeros(nchains)

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")

end

# normalize for each chain
corr ./= nsamples_per_chain
E ./= nsamples_per_chain

# then compute the means
corr = reshape(mean(corr; dim=8), (Ns, L, L, Ns, L, L, nt))
E = mean(E)

# now compute the structure factor
kxs = 
Sqt = zeros(ℂ, L, L, nt)
for t = 1:nt
    Sqt[:, :, t] = structurefactors(H, corr[:, :, :, :, :, :, t])
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
    f["E"] = E
end

println("Finished !")
