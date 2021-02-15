using Random
using LinearAlgebra
using Plots

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

"Does several Monte-Carlo steps, modifying the state vector v inplace"
function mcstep!(H, v, T, niter=1; E=nothing)
    L = size(v)[3]
    Ns = H.Ns
    N = Ns * L^2
    niter *= N
        
    if isnothing(E)
        E = energy(H, v)
    end
    
    for n in 1:niter
        # choose a random spin
        i = rand(1:L)
        j = rand(1:L)
        s = rand(1:Ns)
        # choose a random orientation
        u = randomunitvec()
        uold = v[:, Ns, L, L]
        
        # compute new energy
        v[:, s, i, j] = u
        # TODO make it better
        Enew = energy(H, v)
        ΔE = Enew - E
        # update spin if accepted
        if ΔE < 0 || rand() < exp(-ΔE / T)
            E = Enew # accept !
        else
            # revert to the previous one
            v[:, s, i, j] = uold
        end
    end

    E
end

"Computes the magnetization, that is the mean of all spins (thus a 3D vector)"
function magnetization(v)
    sum(v; dims=(2, 3, 4)) / length(v) * 3
end

function main()
    # H = loadhamiltonian("hamiltonians/skl.dat", [1, 0.5, 0.5])
    H = loadhamiltonian("hamiltonians/square.dat", [1, 0])
    L = 5
    T = 1
    N = H.Ns * L^2

    # initial state
    v = randomstate(H.Ns, L)

    nsamples = 100

    E = zeros(nsamples)
    m = zeros(3, nsamples)

    E[1] = energy(H, v)
    m[:, 1] = magnetization(v)

    # thermalization
    E[2] = mcstep!(H, v, T, 1000; E=E[1])
    m[:, 2] = magnetization(v)
    
    for i = 3:nsamples
        E[i] = mcstep!(H, v, T, 100; E=E[i-1])
        m[:, i] = magnetization(v)
        println(i)
    end

    plot(E)
    
    # return
    E, m
end
