include("scmc.jl")

using Test
using DifferentialEquations

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

@testset "loading hamiltonians" begin
    function inverse_coupling(s1, coupling)
        (i, j, s2, J) = coupling
        s2, (-i, -j, s1, J)
    end

    "Make sure Jij = Jji"
    function check_symmetryp(H)
        for s1 = 1:H.Ns
            for coupling in H.couplings[s1]
                s2, inv_coupling = inverse_coupling(s1, coupling)
                x = findfirst(x -> x == inv_coupling, H.couplings[s2])
                if isnothing(x)
                    println("no inverse coupling for s1 = $s1, coupling = $coupling")
                    println("excepted: $s2, $inv_coupling")
                    return false
                end
            end
        end
        true
    end

    H = loadhamiltonian("hamiltonians/kagome.dat", [1])
    @test check_symmetryp(H)
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
    s0 = v[:, 1, 1, 1] # initial condition for s_1
    # decompose into parallel and normal to M
    spar = (s0 ⋅ m0hat) * m0hat
    snor = s0 - spar

    a = snor
    b = m0hat × a
    vs_ana = zeros(3, H.Ns, 1, 1, nt)

    for n = 1:nt # construct the analytical solution
        t = (n - 1) * dt
        vs_ana[:, 1, 1, 1, n] = spar +  a * cos(ω * t) + b * sin(ω * t) # s1
        vs_ana[:, 2, 1, 1, n] = m0 - vs_ana[:, 1, 1, 1, n] # s2 = M - s1
    end

    @test vs ≈ vs_ana atol=1e-4
end

@testset "kagome lattice system" begin
    H = loadhamiltonian("hamiltonians/kagome.dat", [1])
    
    T = 0.17
    L = 30
    dt = 0.001
    nt = 100
    thermal = 20

    # thermalization
    v = randomstate(H.Ns, L)
    mcstep!(H, v, T, thermal)

    # simulate
    vs = simulate(H, v, dt, nt)

    # check energy conservation
    E0 = energy(H, v)
    Es = mapslices(vs, dims=[1, 2, 3, 4]) do v
        energy(H, v)
    end
    @test all(Es .≈ E0)
    println("E0 = $E0")
    println("Max of ΔE = $(reduce(max, abs.(Es .- E0)))")

    # check magnetization conservation
    m0 = magnetization(v)
    ms = reshape(mapslices(magnetization, vs; dims=[1, 2, 3, 4]), (3, :))
    @test all(ms .≈ m0)

    # check FT of spins
    sq = ftspacespins(H, vs, dt)
    @test all(sq[:, 1, 1, :] .≈ m0)
end
