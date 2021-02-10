using Random

include("hamiltonian.jl")

"Generate a random state (normalized) of shape (3, Ns, Nx, Ny)"
function randomstate(Ns, Nx, Ny)
    vec = randn((3, Ns, Nx, Ny))
    mapslices(v, dims=1) do u
        normalize!(u)
    end
    vec
end

"Do a single Monte-Carlo step"
function mcstep(H, v, T)
    
end
