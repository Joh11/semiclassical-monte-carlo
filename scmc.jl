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

function correlation(v1, v2)
    sum(v1 .* v2) / length(v1)
end

"Use an 8th order Runge-Kutta scheme to advance dt in time the given state"
function dormandprince(H, v, dt)
    Ns, L = size(v)[2:3]
    
    function f(v)
        ret = zeros(3, Ns, L, L)
        for s in 1:Ns
            for i in 1:L
                for j in 1:L
                    ret[:, s, i, j] = -cross(v[:, s, i, j], localfield(H, v, i, j, s))
                end
            end
        end
        ret
    end

    a21 = 1/5
    a31, a32 = [3/40, 9/40]
    a41, a42, a43 = [44/45, -56/15, 32/9]
    a51, a52, a53, a54 = [19372/6561, -25360/2187, 64448/6561, -212/729]
    a61, a62, a63, a64, a65 = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656]
    a71, a72, a73, a74, a75, a76 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]

    b1, b2, b3, b4, b5, b6, b7 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    
    k1 = f(v)
    k2 = f(v + dt * (a21  * k1))
    k3 = f(v + dt * (a31  * k1 + a32 * k2))
    k4 = f(v + dt * (a41  * k1 + a42 * k2 + a43 * k3))
    k5 = f(v + dt * (a51  * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    k6 = f(v + dt * (a61  * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))
    k7 = f(v + dt * (a71  * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

    v + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7)
end

"""Advances the state v in time using the semiclassical
equations. Returns a (3, Ns, L, L, ndt) vector. """
function simulate(H, v, dt, ndt)
    Ns, L = size(v)[2:3]
    ret = zeros(3, Ns, L, L, ndt)
    ret[:, :, :, :, 1] = v

    for i in 2:ndt
        ret[:, :, :, :, i] = dormandprince(H, ret[:, :, :, :, i - 1], dt)
    end
    
    ret
end

# -----------------------------------------------------------------------------
# Dev stuff
# -----------------------------------------------------------------------------

"""Make sure the energy and magnetization are conserved"""
function checkconservation(H, v)
    dt = 0.1
    ndt = 1000
    vs = simulate(H, v, dt, ndt)

    E = zeros(ndt)
    m = zeros(3, ndt)

    for n in 1:ndt
        E[n] = energy(H, v)
        m[:, n] = magnetization(v)
    end

    display(plot(E))

    E, m
end

"""Use this function to find which stride to put to have a small
enough correlation betwen consecutive samples"""
function findcorrelation()
    H = loadhamiltonian("hamiltonians/square.dat", [1, 0])
    L = 5
    T = 1
    N = H.Ns * L^2
    nsamples = 100

    # initial state
    v = randomstate(H.Ns, L)
    v0 = copy(v)

    # mesurements
    E = zeros(nsamples)
    m = zeros(3, nsamples)
    corr = zeros(nsamples)

    E[1] = energy(H, v)
    m[:, 1] = magnetization(v)
    corr[1] = correlation(v0, v)

    for i = 2:nsamples
        E[i] = mcstep!(H, v, T, 1; E=E[i-1])
        m[:, i] = magnetization(v)
        corr[i] = correlation(v0, v)
        println(i)
    end
    E, m, corr
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

    display(plot(E))
    
    # return
    E, m
end
