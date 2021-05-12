#=

A script to compute the dimer static correlation for square kagome
lattice

Scan over the temperature for J1=J2=J3

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

# Utility functions
# =================

function logrange(x1, x2, n)
    [10^y for y in range(log10(x1), log10(x2), length=n)]
end

"Inplace version"
function compute_dimers!(v, L, bonds, dimer)
    for i = 1:L
        for j = 1:L
            for (n, bond) in enumerate(bonds)
                dimer[n, i, j] = v[bond.a.index, i, j] ⋅ v[bond.b.index,
                                                           SCMC.wrapindex(i + bond.b.i, L),
                                                           SCMC.wrapindex(j + bond.b.j, L)]
            end
        end
    end
end

"""Compute the dimer operators D_i for each bond of each unit cell. 
Returns a (12, L, L) array of floats. """
function compute_dimers(v, H)
    L = size(v, 2)
    dimer = zeros(12, L, L)
    compute_dimers!(v, L, bonds(H), dimer)
    dimer
end

"In place version"
function compute_dimer2!(dimer2!, dimer, L)
    for i = 1:size(dimer, 1)
        @views @. dimer2![:, i] = dimer[:, 1, 1] * dimer[i]
    end
end

"Takes a (nbonds, L, L) array as input"
function compute_dimer2(dimer, nbonds, L)
    dimer2 = zeros(nbonds, nbonds * L^2)
    compute_dimer2!(dimer2, dimer, L)
    dimer2
end

"Collect `nsamples_per_chain` samples for each temperature"
function collect_samples!(dimer, dimer2, E; chain=nothing)
    println("Starting for chain $chain ...")
    
    v = randomstate(H.Ns, L)
    Di = zeros(nbonds, L, L) # temporary variable for compute_dimers!
    Di2 = zeros(nbonds, nbonds_tot) # temporary variable for compute_dimers2!
    
    for i in eachindex(Ts)
        T = Ts[i]
        println("Starting T = $T for chain $chain")

        # thermalization step
        if i == 1 # hard the first time
            mcstep!(H, v, T, p["thermal_first"])
        else
            mcstep!(H, v, T, p["thermal"])
        end

        # sampling step
        for nsample = 1:nsamples_per_chain
            println("$nsample / $nsamples_per_chain, T = $T (chain $chain)")
            mcstep!(H, v, T, stride)
            # save the measurements of interest
            compute_dimers!(v, L, bs, Di)
            compute_dimer2!(Di2, Di, L)
            @views dimer[:, i] .+= reshape(Di, :)
            @views dimer2[:, :, i] .+= Di2
            E[i] += energy(H, v)
        end
    end
end

# Params
# ======

const p = Dict("comment" => "Only 4x4 to see if I see a dimer order",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 4,
               "Ts" => logrange(5e-3, 0.1, 10),
               "thermal_first" => 100_000,
               "thermal" => 10_000,
               "nchains" => 8, # because I have 8 threads on my laptop
               "nsamples_per_chain" => 5_000,
               "stride" => 100)
output = "skl_dimer_4x4.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])
const bs = bonds(H)

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const Ts = p["Ts"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]

const Ns = H.Ns * L^2 # number of sites in total
const nbonds = length(bs) # number of bonds per unit cell
const nbonds_tot = nbonds * L^2

# Running the simulation
# ======================

# storing results per chain should be thread safe
dimer = zeros(nbonds_tot, length(Ts), nchains) # <D_i> for all i
dimer2 = zeros(nbonds, nbonds_tot, length(Ts), nchains) # <D_i D_j> for i in UC, all j
E = zeros(length(Ts), nchains)

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    # use a function to avoid dealing with global variables
    # should help with performance
    @views collect_samples!(dimer[:, :, n], dimer2[:, :, :, n], E[:, n];
                            chain=n)
end

# normalize for each chain
dimer ./= nsamples_per_chain
dimer2 ./= nsamples_per_chain
E ./= nsamples_per_chain

# then compute the means
dimer = reshape(mean(dimer; dims=3), (nbonds_tot, length(Ts)))
dimer2 = reshape(mean(dimer2; dims=4), (nbonds, nbonds_tot, length(Ts)))
E = reshape(mean(E; dims=2), length(Ts))

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    for i in eachindex(Ts)
        T = Ts[i]
        group = create_group(f, "$i")
        group["T"] = T
        @views group["dimer"] = dimer[:, i]
        @views group["dimer2"] = dimer2[:, :, i]
        group["E"] = E[i]
    end
end

println("Finished !")
