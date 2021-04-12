module SemiClassicalMonteCarlo

using LinearAlgebra
using StaticArrays
using Random
using LinearAlgebra
using FFTW
using Statistics
using StaticArrays

export loadhamiltonian, energy, deltaenergy, magnetization
export mcstep!, simulate, structuralfactor

include("hamiltonian.jl")
include("montecarlo.jl")

end
