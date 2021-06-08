module SemiClassicalMonteCarlo

using LinearAlgebra
using StaticArrays
using Random
using LinearAlgebra
using FFTW
using Statistics
using StaticArrays
using Base.Threads

export loadhamiltonian, randomstate, energy, deltaenergy, magnetization, bonds
export mcstep!, simulate, structuralfactor, frequencystructuralfactor, allcorrelations!, allcorrelations
export structurefactors, compute_dimers!, compute_dimers, compute_dimer2!, compute_dimer2, skl_order_parameter

include("hamiltonian.jl")
include("montecarlo.jl")
include("postprocessing.jl")

end
