using Random
using LinearAlgebra
# TODO find why this throws a warning
# using Plots 
using FFTW
using Printf
using Statistics
using HDF5: h5write, h5read

push!(LOAD_PATH, pwd())

using HamiltonianMod

"Generate a random state (normalized) of shape (3, Ns, L, L)"
function randomstate(Ns, L)
    vec = randn(3, Ns, L, L)
    mapslices(normalize, vec, dims=1)::Array{Float64, 4}
end

"Returns a random unit vector in 3D"
function randomunitvec()
    u = randn(3)
    normalize(u)
end

"Does several Monte-Carlo steps, modifying the state vector v inplace"
function mcstep!(H, v, T, niter=1)
    L = size(v)[3]
    Ns = H.Ns
    N = Ns * L^2
    niter *= N
    
    for n in 1:niter
        # choose a random spin
        i = rand(1:L)
        j = rand(1:L)
        s = rand(1:Ns)
        # choose a random orientation
        u = randomunitvec()
        uold = v[:, s, i, j]
        
        ΔE = deltaenergy(H, v, u, i, j, s)
        # update spin if accepted
        if ΔE < 0 || rand() < exp(-ΔE / T)
            v[:, s, i, j] = u
        else
            # reject (do nothing actually)
        end
    end
end

"Computes the magnetization, that is the mean of all spins (thus a 3D vector)"
function magnetization(v)
    sum(v; dims=(2, 3, 4))
end

function correlation(v1, v2)
    sum(v1 .* v2) / length(v1)
end

"Build the function representing the time derivative of v"
function makef(H)
    function(v)
        Ns = size(v)[2]
        L = size(v)[3]
        
        ret = zeros(3, Ns, L, L)
        for j in 1:L
            for i in 1:L
                for s in 1:Ns
                    ret[:, s, i, j] = - v[:, s, i, j] × localfield(H, v, i, j, s)
                end
            end
        end
        ret
    end
end

"Use an 8th order Runge-Kutta scheme to advance dt in time the given state"
function dormandprince(f, v, dt)
    a21 = 1/5
    a31, a32 = [3/40, 9/40]
    a41, a42, a43 = [44/45, -56/15, 32/9]
    a51, a52, a53, a54 = [19372/6561, -25360/2187, 64448/6561, -212/729]
    a61, a62, a63, a64, a65 = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656]
    a71, a72, a73, a74, a75, a76 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]
    
    k1 = f(v)
    k2 = f(v + dt * (a21  * k1))
    k3 = f(v + dt * (a31  * k1 + a32 * k2))
    k4 = f(v + dt * (a41  * k1 + a42 * k2 + a43 * k3))
    k5 = f(v + dt * (a51  * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    k6 = f(v + dt * (a61  * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))
    k7 = f(v + dt * (a71  * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

    # First solution
    b1, b2, b3, b4, b5, b6, b7 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    ret1 = v + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7)

    # Second solution
    b1, b2, b3, b4, b5, b6, b7 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]
    ret2 = v + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7)

    ret1, ret2
end

"""Advances the state v in time using the semiclassical
    equations. Returns a (3, Ns, L, L, ndt) vector. """
function simulate(H, v, dt, ndt; stride=1)
    Ns, L = size(v)[2:3]
    ret = zeros(3, Ns, L, L, ndt)
    ret[:, :, :, :, 1] = v

    f = makef(H)
    for i in 1:ndt-1
        println("$i / $ndt")
        ret[:, :, :, :, i] = v
        for n = 1:stride
            v, s2 = dormandprince(f, v, dt)
        end
    end
    ret[:, :, :, :, end] = v
    
    ret
end

@doc raw"""Computes the space FT of the given time
evolved state v. It is defined as such: 
``\vec s_{\vec Q}(t) = \sum_{i, j, s}\vec S_{i, j, s}(t) 
e^{-i (\vec R_{ij} + \vec r_s)\cdot \vec Q}``

Practically, takes a (3, Ns, L, L, ndt) array, and returns a (3, L, L,
ndt) array.  
"""
function ftspacespins(H, vs, dt)
    Ns, L = size(vs)[2:3]
    ndt = size(vs)[5]

    kxs = 2π / L * (0:L-1)
    kys = 2π / L * (0:L-1)
    rs = H.rs
    
    sqsublattice = zeros(Complex{Float64}, 3, Ns, L, L, ndt)
    
    for s in 1:Ns
        # kr = [dot([kx, ky], rs[:, s]) for kx in kxs, ky in kys]
        # initialize the phase shift e^(-i k.r)
        phase_shift = zeros(Complex{Float64}, 3, L, L, ndt)
        for i in 1:L
            for j in 1:L
                kx = kxs[i]
                ky = kys[j]
                
                phase_shift[:, i, j, :] .= exp(-1im * [kx, ky] ⋅ rs[:, s])
            end
        end
        
        sqsublattice[:, s, :, :, :] = fft(vs[:, s, :, :, :], [2, 3]) .* phase_shift
    end

    reshape(sum(sqsublattice; dims=[2]), (3, L, L, ndt))
end

@doc raw"""Computes the (dynamical) structural factor of the given time
evolved state v. It is defined as such: 
``S(\vec Q, t) = <\vec s_{-\vec Q}(0) \vec s_{\vec Q}(t)>``, with
``\vec s_{\vec Q}(t) = \sum_{i, j, s}\vec S_{i, j, s}(t) 
e^{-i (\vec R_{ij} + \vec r_s)\cdot \vec Q}``

Practically, takes a (3, Ns, L, L, ndt) array, and returns a (L, L,
ndt) array.
"""
function structuralfactor(H, vs, dt)
    Ns, L = size(vs)[2:3]
    ndt = size(vs)[5]

    sq = ftspacespins(H, vs, dt)

    # now build the structural factor itself
    # s_-Q(0)
    smq0 = conj.(sq[:, :, :, 1])

    reshape(sum(reshape(smq0, (3, L, L, 1)) .* sq; dims=1), (L, L, ndt))
end

@doc raw"""Computes the (dynamical) frequency structural factor
``S(\vec Q, \omega)`` of the given time evolved state v. """
function frequencystructuralfactor(H, vs, dt)
    Sqt = structuralfactor(H, vs, dt)
    Sqω = fft(Sqt, 3)

    Sqω
end

# -----------------------------------------------------------------------------
# Plotting stuff
# -----------------------------------------------------------------------------

function plotfrequencystructuralfactor(Sqω; lognorm=true)
    L, ndt = size(Sqω)[2:3]
    if L % 2 != 0
        throw(DomainError("L should be even"))
    end
    
    # take the mean if necessary
    if ndims(Sqω) == 4
        Sqω = reshape(mean(Sqω; dims=4), (L, L, ndt))
    end
    
    l = Int(L // 2)
    # build kpath
    nkps = 3l+3
    kpath = fill([], 3l+3)
    kpath[1:l+1] = [[1, 1] .* i for i in 1:l+1] # Γ-M
    kpath[l+2:2l+2] = [[l+1, l+1] .+ [0, -1] * i for i in 0:l] # M-X
    kpath[2l+3:3l+3] = [[l+1, 1] .+ [-1, 0] * i for i in 0:l] # X-Γ

    z = zeros(nkps, ndt)
    for nk in 1:nkps
        nx, ny = kpath[nk]
        for nt in 1:ndt
            if lognorm
                z[nk, nt] = log10(1e-8 + abs(Sqω[nx, ny, nt]))
            else
                z[nk, nt] = abs(Sqω[nx, ny, nt])
            end
        end
    end
    # , ["Γ", "M", "X"]
    heatmap(transpose(z); xticks=([1, 7, 13, 19], ["Γ", "M", "X", "Γ"])), z
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

function reproducefig4()
    H = loadhamiltonian("hamiltonians/kagome.dat", [1])
    output = "kagome-fig4-t40.h5"
    
    T = 10
    L = 30
    nsamples = 1000
    dt = 0.1
    nt = 400
    stride = 15
    thermal = 20

    E = 0
    m = zeros(3)
    Sqω = zeros(Complex{Float64}, L, L, nt)
    Sqt = zeros(Complex{Float64}, L, L, nt)
    vs = zeros(3, H.Ns, L, L, nt)
    
    # initial state
    v = randomstate(H.Ns, L)

    # thermalization
    println("Doing thermalization")
    mcstep!(H, v, T, thermal)
    E = energy(H, v)
    m = magnetization(v)

    # time evolution
    println("Doing simulation")
    vs = simulate(H, v, dt, nt)
    println("Computing structural factor")
    Sqω = frequencystructuralfactor(H, vs, dt)
    Sqt = structuralfactor(H, vs, dt)

    function saveh5(output, i, E, m, vs, Sqω, Sqt)
        h5write(output, "$i/E", E)
        h5write(output, "$i/m", m)
        h5write(output, "$i/vs", vs)
        h5write(output, "$i/Sqω", Sqω)
        h5write(output, "$i/Sqt", Sqt)
    end
    
    saveh5(output, 1, E, m, vs, Sqω, Sqt)
    
    for i = 2:nsamples
        println("$i / $nsamples")
        mcstep!(H, v, T, stride)
        E = energy(H, v)
        m = magnetization(v)

        # time evolution
        vs = simulate(H, v, dt, nt)
        Sqω = frequencystructuralfactor(H, vs, dt)
        Sqt = structuralfactor(H, vs, dt)
        saveh5(output, i, E, m, vs, Sqω, Sqt)
    end
end

# reproducefig4()
