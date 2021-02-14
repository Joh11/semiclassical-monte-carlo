using Random
using LinearAlgebra

push!(LOAD_PATH, pwd())

using HamiltonianMod

"Generate a random state (normalized) of shape (3, Ns, L, L)"
function randomstate(Ns, L)
    vec = randn(3, Ns, L, L)
    mapslices(vec, dims=1) do u
        normalize!(u)
    end
    vec
end

"Returns a random unit vector in 3D"
function randomunitvec()
    u = randn(3)
    normalize(u)
end

"Do several Monte-Carlo steps, modifying the state vector v inplace"
function mcstep!(H, v, T, niter=1; E=nothing)
    L = size(v)[3]
    Ns = H.Ns

    if isnothing(E)
        E = energy(H, v)
    end
    
    for n in 1:niter
        # choose a random spin
        i = rand(1:L)
        j = rand(1:L)
        s = rand(1:L)
        # choose a random orientation
        u = randomunitvec()
        uold = v[:, Ns, L, L]
        
        # compute new energy
        v[:, Ns, L, L] = u
        # TODO make it better
        Enew = energy(H, v)
        ΔE = Enew - E
        # update spin if accepted
        if ΔE < 0 || rand() < exp(-ΔE / T)
            E = Enew # accept !
        else
            # revert to the previous one
            v[:, Ns, L, L] = uold
        end
    end
end
