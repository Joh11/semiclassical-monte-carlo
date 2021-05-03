#=

A script to compute the dimer static correlation for square kagome
lattice

Scan over the temperature for J1=J2=J3

Each chain is indepentent, and stores the result as a temporary
average. Afterwards, an average is made over all the chains (OK since
they each have the same number of samples). 

=#

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
function compute_dimers!(v, L, dimer)
    for i = 1:L
        for j = 1:L
            # do this manually for now
            dimer[1, i, j] = v[1, i, j] ⋅ v[6, i, SCMC.wrapindex(j+1, L)]
            dimer[2, i, j] = v[2, i, j] ⋅ v[6, i, SCMC.wrapindex(j+1, L)]
            dimer[3, i, j] = v[1, i, j] ⋅ v[2, i, j]
            dimer[4, i, j] = v[1, i, j] ⋅ v[3, i, j]
            dimer[5, i, j] = v[1, i, j] ⋅ v[4, i, j]
            dimer[6, i, j] = v[2, i, j] ⋅ v[5, i, j]
            dimer[7, i, j] = v[2, i, j] ⋅ v[3, SCMC.wrapindex(i+1, L), j]
            dimer[8, i, j] = v[3, i, j] ⋅ v[4, i, j]
            dimer[9, i, j] = v[5, i, j] ⋅ v[3, SCMC.wrapindex(i+1, L), j]
            dimer[10, i, j] = v[4, i, j] ⋅ v[5, i, j]
            dimer[11, i, j] = v[4, i, j] ⋅ v[6, i, j]
            dimer[12, i, j] = v[5, i, j] ⋅ v[6, i, j]
        end
    end
end

"""Compute the dimer operators D_i for each bond of each unit cell. 
Returns a (12, L, L) array of floats. """
function compute_dimers(v)
    L = size(v)[2]
    dimer = zeros(12, L, L)
    compute_dimers!(v, L, dimer)
    dimer
end

# Params
# ======

const p = Dict("comment" => "first try, mostly to find where is Tc",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 10,
               "Ts" => logrange(5e-3, 0.1, 10),
               "thermal_first" => 10_000,
               "thermal" => 1_000,
               "nchains" => 10,
               "nsamples_per_chain" => 100,
               "stride" => 50)
output = "skl_dimer.h5"
H = loadhamiltonian("hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const L = p["L"]
const nchains = p["nchains"]
const Ts = p["Ts"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]

const Ns = 6 * L^2 # number of sites in total

# Running the simulation
# ======================

total_dimer = zeros(2Ns, length(Ts))
total_dimer2 = zeros(12, 2Ns, length(Ts))

@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    
    v = randomstate(H.Ns, L)
    dimer = zeros(2Ns) # <D_i> for all i
    dimer2 = zeros(12, 2Ns) # <D_i D_j> for i in UC, all j
    Di = zeros(12, L, L) # temporary variable for compute_dimers!
    
    for i in eachindex(Ts)
        T = Ts[i]
        println("Starting T = $T for chain $n / $nchains")
        # reset variables
        dimer .= 0
        dimer2 .= 0
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
            compute_dimers!(v, L, Di)
            dimer += reshape(Di, :)
            dimer2 += reshape(@view Di[:, 1, 1], (12, 1)) .*  reshape(Di, (1, :))
        end

        # now that everything is sampled, average
        dimer /= nsamples_per_chain
        dimer2 /= nsamples_per_chain

        # and append to the total results (for all chains)
        total_dimer[:, i] += dimer
        total_dimer2[:, i] += dimer2
    end
end

# average over all chains
total_dimer /= nchains
total_dimer2 /= nchains

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
        group["dimer"] = @view total_dimer[:, n]
        group["dimer2"] = @view total_dimer2[:, :, n]
    end
end

println("Finished !")
