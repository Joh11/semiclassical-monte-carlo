#=
A script to update only the structure factor from a HDF5 dataset
=#

using HDF5
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

const ℂ = Complex{Float64}

oldfilename = "skl_dyn_factor.h5"
newfilename = "skl_dyn_factor_fixed.h5"

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

const H = loadhamiltonian("../hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])
const nt = p["nt"]

# read from old
f = h5open(oldfilename, "r")
corr = read(f["corr"])
E = read(f["E"])
close(f)

# now compute the structure factor
println("Now computing structure factor ...")
kxs = 8π * (-1:0.05:1)
kys = 8π * (-1:0.05:1)
Sqt = zeros(ℂ, length(kxs), length(kys), nt)
for t = 1:nt
    println("Doing $t / $nt ...")
    @views Sqt[:, :, t] = structurefactors(H, corr[:, :, :, :, :, :, t], kxs, kys)
end

# now compute the frequency resolved structure factor
println("Now computing (frequency resolved) structure factor ...")
Sqω = frequencystructuralfactor(Sqt)

# Saving
# ======

println("Saving to $newfilename ...")
h5open(newfilename, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["corr"] = corr
    f["Sqt"] = Sqt
    f["Sqω"] = Sqω
    f["E"] = E
end

println("Finished !")
