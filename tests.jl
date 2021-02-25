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
    @test a1⋅b2 ≈ 0
    @test a2⋅b1 ≈ 0
    @test a2⋅b2 ≈ 2π
end
