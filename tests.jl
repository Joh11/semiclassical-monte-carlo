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
    @test isapprox(a1⋅b2, 0; atol=1e-12)
    @test isapprox(a2⋅b1, 0; atol=1e-12)
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
