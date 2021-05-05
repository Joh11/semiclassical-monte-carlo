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

"Takes a (12, L, L) array as input"
function compute_dimer2(dimer, L)
    dimer2 = zeros(12, 12L^2)

    for i = 1:12L^2
        @views dimer2[:, i] = dimer[:, 1, 1] .* dimer[i]
    end
    
    dimer2
end

# Params
# ======

const p = Dict("comment" => "A big one to make sure everything has converged",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 10,
               "Ts" => logrange(5e-3, 0.1, 10),
               "thermal_first" => 100_000,
               "thermal" => 10_000,
               "nchains" => 8, # because I have 8 threads on my laptop
               "nsamples_per_chain" => 5000,
               "stride" => 500)
output = "skl_dimer_long2.h5"
const H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])
const bs = bonds(H)

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const Ts = p["Ts"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]

const Ns = 6L^2 # number of sites in total

# Running the simulation
# ======================

total_dimer = zeros(2Ns, length(Ts))
total_dimer2 = zeros(12, 2Ns, length(Ts))
total_energy = zeros(length(Ts))

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    
    v = randomstate(H.Ns, L)
    dimer = zeros(2Ns) # <D_i> for all i
    dimer2 = zeros(12, 2Ns) # <D_i D_j> for i in UC, all j
    E = 0
    Di = zeros(12, L, L) # temporary variable for compute_dimers!
    
    for i in eachindex(Ts)
        T = Ts[i]
        println("Starting T = $T for chain $n / $nchains")
        # reset variables
        dimer .= 0
        dimer2 .= 0
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
            mcstep!(H, v, T, stride)
            # save the measurements of interest
            compute_dimers!(v, L, bs, Di)
            dimer += reshape(Di, :)
            dimer2 += compute_dimer2(Di, L)
            E += energy(H, v)
        end

        # now that everything is sampled, average
        dimer /= nsamples_per_chain
        dimer2 /= nsamples_per_chain
        E /= nsamples_per_chain

        # and append to the total results (for all chains)
        total_dimer[:, i] += dimer
        total_dimer2[:, :, i] += dimer2
        total_energy[i] += E
    end
end

# average over all chains
total_dimer /= nchains
total_dimer2 /= nchains
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
        @views group["dimer"] = total_dimer[:, n]
        @views group["dimer2"] = total_dimer2[:, :, n]
        group["E"] = total_energy[n]
    end
end

println("Finished !")
