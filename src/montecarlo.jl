using DifferentialEquations

const Vec3 = SVector{3, Float64}

"Returns a random unit vector in 3D"
function randomunitvec()
    u = @SVector randn(3)
    normalize(u)
end


"Generate a random state (normalized) of shape (3, Ns, L, L)"
function randomstate(Ns, L)
    vec = zeros(Vec3, Ns, L, L)
    for i in eachindex(vec)
        vec[i] = randomunitvec()
    end
    vec
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
        uold = v[s, i, j]
        
        ΔE = deltaenergy(H, v, u, i, j, s)
        # update spin if accepted
        if ΔE < 0 || rand() < exp(-ΔE / T)
            v[s, i, j] = u
        else
            # reject (do nothing actually)
        end
    end
end

"Computes the magnetization, that is the mean of all spins (thus a 3D vector)"
function magnetization(v)
    sum(v)
end

function correlation(v1, v2)
    sum(v1 .* v2) / length(v1)
end

"Build the function representing the time derivative of v"
function makef(H)
    function f(dv, v, p=nothing, t=nothing)
        Ns = size(v)[1]
        L = size(v)[2]
        for j in 1:L
            for i in 1:L
                for s in 1:Ns
                    @views dv[s, i, j] = localfield(H, v, i, j, s) × v[s, i, j]
                end
            end
        end
    end
end

"""Advances the state v in time using the semiclassical
        equations. Returns a (3, Ns, L, L, ndt) vector. """
function simulate(H, v, tmax)
    # define the problem
    f = makef(H)
    prob = ODEProblem(f, v, (0, tmax))

    # solve it
    # save every second
    sol = solve(prob, saveat=10, progress=true)

    # put it in a format usable by the rest of the code
    nt = length(sol)
    L = size(v)[2]
    ret = zeros(Vec3, H.Ns, L, L, nt)

    for n in eachindex(sol.u)
        @views ret[:, :, :, n] = sol.u[n]
    end
    
    ret
end

@doc raw"""Computes the space FT of the given time
        evolved state v. It is defined as such: 
        ``\vec s_{\vec Q}(t) = \sum_{i, j, s}\vec S_{i, j, s}(t) 
        e^{-i (\vec R_{ij} + \vec r_s)\cdot \vec Q}``

        Practically, takes a (Ns, L, L, ndt) array of Vec3, and
        returns a (3, L, L, ndt) array.  """
function ftspacespins(H, vs)
    Ns, L = size(vs)[1:2]
    ndt = size(vs)[4]

    kxs = 2π / L * (0:L-1)
    kys = 2π / L * (0:L-1)
    rs = H.rs

    # copy everything to a simple array
    vs_temp = zeros(3, Ns, L, L, ndt)
    for t in 1:ndt
        for j in 1:L
            for i in 1:L
                for s in 1:Ns
                    vs_temp[1, s, i, j, t] = vs[s, i, j, t][1]
                    vs_temp[2, s, i, j, t] = vs[s, i, j, t][2]
                    vs_temp[3, s, i, j, t] = vs[s, i, j, t][3]
                end
            end
        end
    end
    
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
        
        sqsublattice[:, s, :, :, :] = fft(vs_temp[:, s, :, :, :], [2, 3]) .* phase_shift
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
    Ns, L = size(vs)[1:2]
    ndt = size(vs)[4]

    sq = ftspacespins(H, vs)

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
