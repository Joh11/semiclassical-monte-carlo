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
function compute_dimer2!(dimer2!, dimer, nbonds, L)
    for i = 1:nbonds
        @simd for j = 1:nbonds*L^2
            dimer2![i, j] = dimer[i, 1, 1] * reshape(dimer, :)[j]
        end
    end    
end

"Takes a (nbonds, L, L) array as input"
function compute_dimer2(dimer, nbonds, L)
    dimer2 = zeros(nbonds, nbonds * L^2)
    compute_dimer2!(dimer2, dimer, nbonds, L)
    dimer2
end

"""Takes a dimer array (Nbonds, L, L), for L|2, and returns the order
parameter (a scalar)"""
function skl_order_parameter(dimers)
    L = size(v, 2)
    Nbonds = size(v, 1)
    @assert L % 2 == 0
    @assert Nbonds == 12
    
    # θ for each different bond
    red = 1
    blue = -1
    green = 1
    pink = -1

    # for the two kinds of unit cell (bc symmetry breaking)
    θ1 = [blue, red, 0, 0, red, 0, 0, 0, 0, blue, 0, 0]
    θ2 = [0, 0, green, pink, 0, green, pink, green, pink, 0, green, pink]

    # now compute the order parameter
    order = 0
    for x = 1:L
        for y = 1:L
            if (x + y) % 2 == 0
                @views order += θ1 ⋅ dimers[:, x, y]
            else
                @views order += θ2 ⋅ dimers[:, x, y]
            end
        end
    end

    order
end

"Collect `nsamples_per_chain` samples for each temperature"
function collect_samples!(dimer, dimer2, order_param, E; chain=nothing)
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
            if nsample % 100 == 0
                println(nsample, " / ", nsamples_per_chain,
                        ", T = ", T, " (chain ", chain, ")")
            end
            
            mcstep!(H, v, T, stride)
            # save the measurements of interest
            compute_dimers!(v, L, bs, Di)
            compute_dimer2!(Di2, Di, nbonds, L)
            @views dimer[:, i] .+= reshape(Di, :)
            @views dimer2[:, :, i] .+= Di2
            # take the absolute value, as <O> = 0
            order_param[i] += abs(skl_order_parameter(Di))
            E[i] += energy(H, v)
        end
    end
end

# Params
# ======

const p = Dict("comment" => "Higher T",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 4,
               "Ts" => [0.3, 0.6, 1, 3, 6, 10],
               "thermal_first" => 100_000,
               "thermal" => 10_000,
               "nchains" => 8, # because I have 8 threads on my laptop
               "nsamples_per_chain" => 4_000,
               "stride" => 500)
output = "skl_dimer_4x4_highT.h5"
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
order_param = zeros(length(Ts), nchains) # <|O|>
E = zeros(length(Ts), nchains)

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    # use a function to avoid dealing with global variables
    # should help with performance
    @views collect_samples!(dimer[:, :, n], dimer2[:, :, :, n], order_param[:, n], E[:, n];
                            chain=n)
end

# normalize for each chain
dimer ./= nsamples_per_chain
dimer2 ./= nsamples_per_chain
order_param ./= nsamples_per_chain
E ./= nsamples_per_chain

# then compute the means
dimer = reshape(mean(dimer; dims=3), (nbonds_tot, length(Ts)))
dimer2 = reshape(mean(dimer2; dims=4), (nbonds, nbonds_tot, length(Ts)))
order_param = reshape(mean(order_param; dims=2), length(Ts))
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
