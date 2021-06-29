#=
A script to add the whole (extended) BZ structure factor S(q) to a HDF5 dataset
=#

using HDF5
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

# Functions
# =========

function includeBZstructurefactor(f, H)
    # read corr (at zero time) from old (+ all other necessary variables)
    @views corr0 = read(f["corr"])[:, :, :, :, :, :, 1]

    kxs = 4π * Array(-1:0.01:1)
    kys = 4π * Array(-1:0.01:1)
    
    St0 = structurefactors(H, corr0, kxs, kys)

    # save
    if haskey(f, "St0")
        println("Deleting old content ...")
        delete_object(f, "St0")
    end
    
    write(f, "St0", real.(St0))
end

# Script itself
# =============

filename = ARGS[1]
f = h5open(filename, "r+")

J1 = read(attributes(f)["J1"])
J2 = read(attributes(f)["J2"])
J3 = read(attributes(f)["J3"])

const H = loadhamiltonian("hamiltonians/skl.dat", [J1, J2, J3])

includeBZstructurefactor(f, H)

println("Finished !")
