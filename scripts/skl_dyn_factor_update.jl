#=
A script to add the whole (extended) BZ structure factor S(q) to a HDF5 dataset
=#

using HDF5
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

# Functions
# =========

function h5overwrite(filename, name, data)
    if haskey(filename, name)
        delete_object(filename, name)
    end
    write(filename, name, data)
end

"Add a whole BZ structure factor for t=0"
function includeBZstructurefactor(f, H)
    # read corr (at zero time) from old (+ all other necessary variables)
    @views corr0 = read(f["corr"])[:, :, :, :, :, :, 1]

    kxs = 4π * Array(-1:0.01:1)
    kys = 4π * Array(-1:0.01:1)
    
    St0 = structurefactors(H, corr0, kxs, kys)

    # save
    h5overwrite(f, "St0", real.(St0))
end

function kpath2mat(kpath)
    mat = zeros(2, length(kpath))
    for (n, k) in enumerate(kpath)
        mat[:, n] .= k
    end
    mat
end



"Update the kpath and kpath specific structure factor (both time and energy resolved)"
function update_kpath_structurefactor(f, H, highsympoints)
    # TODO
    nt, L, nk
    
    kpath = makekpath(highsympoints, nk)
    
    # first the time resolved structure factor
    Sqt = zeros(ℂ, length(kpath), nt)
    Rs = compute_positions(H, L)
    for t = 1:nt
        println("Doing $t / $nt ...")
        @views structurefactor_kpath!(Rs, corr[:, :, :, :, :, :, t], kpath, Sqt[:, t])
    end

    # now compute the frequency resolved structure factor
    println("Now computing (frequency resolved) structure factor ...")
    Sqω = frequencystructuralfactor(Sqt; dim=2)

    println("Saving ...")
    h5overwrite(f, "kpath", kpath2mat(kpath))
    h5overwrite(f, "Sqt", Sqt)
    h5overwrite(f, "Sqω", Sqω)
end

# Script itself
# =============

filename = ARGS[1]
f = h5open(filename, "r+")

J1 = read(attributes(f)["J1"])
J2 = read(attributes(f)["J2"])
J3 = read(attributes(f)["J3"])

const H = loadhamiltonian("hamiltonians/skl.dat", [J1, J2, J3])

# includeBZstructurefactor(f, H)

# highsympoints = [ # UUD
#     [0, 0],
#     [2π, 0],
#     [2π, 2π],
#     [0, 0]
# ]

highsympoints = [ # Neel
                  [0, 0],
                  [2π, 0],
                  [2π, 2π],
                  [0, 0]
                  ]

update_kpath_structurefactor(f, H, highsympoints)

println("Finished !")
