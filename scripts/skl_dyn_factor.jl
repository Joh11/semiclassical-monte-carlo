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

    # preallocation for DifferentialEquations.jl
    # (see https://discourse.julialang.org/t/using-differentialequations-jl-for-very-large-system/1030/2 )
    timeseries = []
    ts = []
    ks = []
    
    # sampling step
    for nsample = 1:nsamples_per_chain
        println("Doing sample $nsample / $nsamples_per_chain for chain $chain")
        mcstep!(H, v, T, stride)

        # time evolution
        vs = simulate(H, v, p["dt"], nt, timeseries, ts, ks)
        
        # save measurements of interest
        for t = 1:nt
            @views allcorrelations!(vs[:, :, :, 1], vs[:, :, :, t],
                                    Ns, L, corr_tmp[:, :, :, :, :, :, t])
        end
        corr .+= corr_tmp
        E += energy(H, v)

        # use the last time evolved state
        @views v .= vs[:, :, :, end]
        
        println("done $nsample / $nsamples_per_chain (chain $chain)")
        flush(stdout)

        # trigger the garbage collector manually
        if chain == 1 && nsample % 50 == 0
            println("Running GC manually ...")
            GC.gc()
        end
    end
end

function kpath2mat(kpath)
    mat = zeros(2, length(kpath))
    for (n, k) in enumerate(kpath)
        mat[:, n] .= k
    end
    mat
end

# Params
# ======

const p = Dict("comment" => "8x8, J1=J2=J3, 4x zoom",
               "J1" => 1,
               "J2" => 1,
               "J3" => 1,
               "L" => 6,
               "T" => 0.1,
               "thermal" => 100_000,
               "nchains" => 8, # because 8 cores on my laptop
               "nsamples_per_chain" => 4_000, # so ~ 30k samples
               "stride" => 100,
               # time evolution params
               "dt" => 2,
               "nt" => 200,
               "nk" => 100
               )
output = "skl_dyn_factor_zoomed_8x8.h5"
const H = loadhamiltonian("../hamiltonians/skl.dat", [p["J1"], p["J2"], p["J3"]])

# variables often used have an alias
const Ns = H.Ns
const L = p["L"]
const nchains = p["nchains"]
const T = p["T"]
const nsamples_per_chain = p["nsamples_per_chain"]
const stride = p["stride"]
const nt = p["nt"]
const nk = p["nk"]

const Nsites = H.Ns * L^2 # number of sites in total

# Running the simulation
# ======================

# storing results per chain should be thread safe
corr = zeros(Ns, L, L, Ns, L, L, nt, nchains)
E = zeros(nchains)

Threads.@threads for n in 1:nchains
    println("Starting for chain $n / $nchains ...")
    @views collect_samples!(corr[:, :, :, :, :, :, :, n], E[n]; chain=n)
end

# normalize for each chain
corr ./= nsamples_per_chain
E ./= nsamples_per_chain

# then compute the means
corr = reshape(mean(corr; dims=8), (Ns, L, L, Ns, L, L, nt))
E = mean(E)

# now compute the structure factor
println("Now computing structure factor ...")
Γ = [0, 0]
X = [π, 0]
twoπ0 = [2π, 0]
M = [π, π]
kpath = makekpath([Γ, X, twoπ0, M, Γ], nk)

Sqt = zeros(ℂ, length(kpath), nt)
const Rs = compute_positions(H, L)
for t = 1:nt
    println("Doing $t / $nt ...")
    @views structurefactor_kpath!(Rs, corr[:, :, :, :, :, :, t], kpath, Sqt[:, t])
end

# now compute the frequency resolved structure factor
println("Now computing (frequency resolved) structure factor ...")
Sqω = frequencystructuralfactor(Sqt; dim=2)

# Saving
# ======

println("Saving to $output ...")
h5open(output, "w") do f
    # save params
    for (name, val) in p
        attributes(f)[name] = val
    end

    f["corr"] = corr
    f["kpath"] = kpath2mat(kpath)
    f["Sqt"] = Sqt
    f["Sqω"] = Sqω
    f["E"] = E
end

println("Finished !")
