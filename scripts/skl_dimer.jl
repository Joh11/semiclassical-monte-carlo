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

function logrange(x1, x2, n)
    [10^y for y in range(log10(x1), log10(x2), length=n)]
end

"""Compute the dimer operators D_i for each bond of each unit cell. 

Returns a (12, L, L) array of floats. """
function compute_dimers(v)
    L = size(v)[2]
    dimer = zeros(12, L, L)

    for i = 1:L
        for j = 1:L
            dimer[:, i, j] = 
        end
    end

    dimer
end

# Params
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

# Run the simulation
@threads for n in 1:nchains
    println("Starting for chain $n / $(nchains) ...")
    
    v = randomstate(H.Ns, L)
    dimer = zeros(2Ns) # <D_i> for all i
    dimer2 = zeros(12, 2Ns) # <D_i D_j> for i in UC, all j
    
    for i in eachindex(Ts)
        T = Ts[i]
        # thermalization step
        if i == 1 # hard the first time
            mcstep!(H, v, T, p["thermal_first"])
        else
            mcstep!(H, v, T, p["thermal"])
        end

        # sampling step
        for nsample = 1:nsamples_per_chain
            mcstep!(H, v, T, stride)
            # save the measurements of interest
            Di = compute_dimers(v)
            dimer += Di
            dimer2 += 
        end
    end
    
end

# Save
