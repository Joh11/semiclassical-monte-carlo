include("scmc.jl")

using Test


@testset "random unit vec" begin
    v = randomunitvec()
    @test size(v) == (3,)
    @test norm(v) ≈ 1
end

@testset "random state" begin
    Ns, L = 10, 100
    v = randomstate(Ns, L)

    @test size(v) == (3, Ns, L, L)
    @test all(mapslices(norm, v; dims=1) .≈ 1)
end

@testset "rec lattice" begin
    lattice = rand(2, 2)
    rec = reciprocallattice(lattice)

    a1, a2 = lattice[:, 1], lattice[:, 2]
    b1, b2 = rec[:, 1], rec[:, 2]

    @test a1⋅b1 ≈ 2π
    @test a1⋅b2 ≈ 0 atol=1e-12
    @test a2⋅b1 ≈ 0 atol=1e-12
    @test a2⋅b2 ≈ 2π
end

@testset "delta energy" begin
    L = 144
    H = loadhamiltonian("hamiltonians/square.dat", [1, 2])
    v = randomstate(H.Ns, L)

    E = energy(H, v)

    # update a single spin
    i, j = rand(1:L), rand(1:L)
    s = rand(1:H.Ns)
    u = randomunitvec()

    ΔE = deltaenergy(H, v, u, i, j, s)

    v[:, s, i, j] = u
    newE = energy(H, v)

    @test (newE - E) ≈ ΔE
end

@testset "two spins system" begin
    H = loadhamiltonian("hamiltonians/two-spins.dat", [1])
    v = randomstate(H.Ns, 1)

    dt = 0.1
    nt = 1000

    # run the simulation
    vs = simulate(H, v, dt, nt)
    
    # make sure the magnetization is conserved
    m0 = v[:, 1, 1, 1] + v[:, 2, 1, 1]
    ms = reshape(mapslices(magnetization, vs; dims=[1, 2, 3, 4]), (3, :))
    @test all(ms .≈ m0)
    
    # make sure the energy is conserved
    E0 = energy(H, v)
    Es = mapslices(vs; dims=[1, 2, 3, 4]) do v
        energy(H, v)
    end
    @test all(Es .≈ E0)

    # make sure it matches the analytical solution
    m0hat = normalize(m0)
    ω = norm(m0) # should be J|M|, but here J = 1
    a = v[:, 1, 1, 1]
    b = m0hat × v[:, 1, 1, 1]
    vs_ana = zeros(3, H.Ns, 1, 1, nt)
    
    for n = 1:nt # construct the analytical solution
        t = (n - 1) * dt
        vs_ana[:, 1, 1, 1, n] = a * cos(ω * t) + b * sin(ω * t) # s1
        vs_ana[:, 2, 1, 1, n] = m0 - vs_ana[:, 1, 1, 1, n] # s2 = M - s1
    end

    @test vs ≈ vs_ana

end
