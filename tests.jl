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
