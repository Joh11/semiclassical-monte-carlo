module SemiClassicalMonteCarlo

using LinearAlgebra
using StaticArrays
using Random
using LinearAlgebra
using FFTW
using Statistics
using StaticArrays

export loadhamiltonian, randomstate, energy, deltaenergy, magnetization
export mcstep!, simulate, structuralfactor, frequencystructuralfactor

include("hamiltonian.jl")
include("montecarlo.jl")

end
