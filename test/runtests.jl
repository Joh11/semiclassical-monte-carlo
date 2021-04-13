using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo
using Test
using LinearAlgebra
using Statistics

@testset "SemiClassicalMonteCarlo.jl" begin
    @testset "random unit vec" begin
        v = SCMC.randomunitvec()
        @test size(v) == (3,)
        @test norm(v) ≈ 1
    end

    @testset "random state" begin
        Ns, L = 10, 100
        v = SCMC.randomstate(Ns, L)

        @test size(v) == (Ns, L, L)
        @test all(map(norm, v) .≈ 1)
    end

    
    @testset "rec lattice" begin
        lattice = rand(2, 2)
        rec = SCMC.reciprocallattice(lattice)
        
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
        
        H = loadhamiltonian("../hamiltonians/kagome.dat", [1])
        @test check_symmetryp(H)
    end

    @testset "delta energy" begin
        L = 144
        H = loadhamiltonian("../hamiltonians/square.dat", [1, 2])
        v = SCMC.randomstate(H.Ns, L)

        E = energy(H, v)

        # update a single spin
        i, j = rand(1:L), rand(1:L)
        s = rand(1:H.Ns)
        u = SCMC.randomunitvec()

        ΔE = deltaenergy(H, v, u, i, j, s)

        v[s, i, j] = u
        newE = energy(H, v)

        @test (newE - E) ≈ ΔE
    end

    @testset "two spins system" begin
        H = loadhamiltonian("../hamiltonians/two-spins.dat", [1])
        v = SCMC.randomstate(H.Ns, 1)

        dt = 0.1
        nt = 1000

        # run the simulation
        vs = simulate(H, v, dt, nt)
        
        # make sure the magnetization is conserved
        m0 = v[1, 1, 1] + v[2, 1, 1]
        ms = reshape(mapslices(magnetization, vs; dims=[1, 2, 3]), (3, :))
        @test all(ms .≈ m0)
        
        # make sure the energy is conserved
        E0 = energy(H, v)
        Es = mapslices(vs; dims=[1, 2, 3]) do v
            energy(H, v)
        end
        @test reduce(max, Es .- E0) ≈ 0 atol=1e-5

        # make sure it matches the analytical solution
        m0hat = normalize(m0)
        ω = norm(m0) # should be J|M|, but here J = 1
        s0 = v[1, 1, 1] # initial condition for s_1
        # decompose into parallel and normal to M
        spar = (s0 ⋅ m0hat) * m0hat
        snor = s0 - spar

        a = snor
        b = m0hat × a
        vs_ana = zeros(SCMC.Vec3, H.Ns, 1, 1, nt)

        for n = 1:nt # construct the analytical solution
            t = (n - 1) * dt
            vs_ana[1, 1, 1, n] = spar +  a * cos(ω * t) + b * sin(ω * t) # s1
            vs_ana[2, 1, 1, n] = m0 - vs_ana[1, 1, 1, n] # s2 = M - s1
        end

        @test vs ≈ vs_ana atol=1e-4
    end

    function compute_norm_Sqt_ana(H, L)
        norm_Sqt_ana = zeros(L, L)

        r2 = H.rs[:, 2]
        for nx = 1:L
            for ny = 1:L
                q = 2π * [nx-1, ny-1] / L
                norm_Sqt_ana[nx, ny] = 1/2 * (1 + cos(dot(q, r2)))
            end
        end

        norm_Sqt_ana
    end

    @testset "two spins system - structure factor" begin
        H = loadhamiltonian("../hamiltonians/two-spins.dat", [1])

        # first look at it for the uniform distribution (no MC step)
        nsamples = 800
        L = 2
        dt = 0.05
        nt = 100
        N = L^2

        function runone(n)
            v = SCMC.randomstate(H.Ns, L)
            simulate(H, v, dt, nt)
        end

        samples = map(runone, 1:nsamples)
        
        # 1. check that the decaying term decays
        decay = mapreduce(+, samples) do vs
            v = vs[:, :, :, 1]
            s1, s2 = v[1, 1, 1], v[2, 1, 1]
            M = s1 + s2
            normM = norm(M)
            Mhat = M / normM
            sperp = s1 - dot(s1, Mhat) * Mhat
            
            reshape(map(s -> dot(sperp, s), vs[1, 1, 1, :]), :)
        end
        decay /= nsamples
        # plot(decay)
        # pretty high tolerance because this is random
        @test decay[1] ≈ 1/2 atol=3e-2
        @test decay[end] ≈ 0 atol=4e-2
        
        # 2. check the structure factor
        Sqt = mapreduce(vs -> structuralfactor(H, vs, dt), +, samples)
        Sqt /= nsamples

        # nsamples | max(imag.(Sqt[2, 1, :]))
        # 100      | 0.47
        # 200      | 0.42
        # 400      | 0.25
        # 600      | 0.13
        # 800      | 0.15
        # 1000     | 0.13
        # 2000     | 0.06
        # this is a lame atol but it is slowly vanishing...
        @test reduce(max, imag.(Sqt[2, 1, :]) ./ real.(Sqt[2, 1, :])) ≈ 0 atol=5e-2
        
        # plot(real.(Sqt[2, 1, :]), label="Sqt")
        # plot!(2N * (decay .+ 0.5), label="Sqt from the decaying term")
        @test decay ≈ (real.(Sqt[2, 1, :]) / 2N .- 1/2) atol=0.5

        # now compare with the analytical solution
        norm_Sqt = real.(Sqt[:, :, end] ./ Sqt[:, :, 1])

        norm_Sqt_ana = compute_norm_Sqt_ana(H, L)

        @test norm_Sqt ≈ norm_Sqt_ana atol=5e-2
    end

    @testset "two spins system - structure factor, high L" begin
        H = loadhamiltonian("../hamiltonians/two-spins.dat", [1])
        
        nsamples = 100
        L = 20
        dt = 0.1
        nt = 100
        N = L^2

        function runone(n)
            v = SCMC.randomstate(H.Ns, L)
            simulate(H, v, dt, nt)
        end
        samples = map(runone, 1:nsamples)

        Sqts = map(vs -> structuralfactor(H, vs, dt), samples)

        norm_Sqt_ana = compute_norm_Sqt_ana(H, L)

        ncuterr = map([10, 30, 60, 100]) do ncut
            Sqt = mean(Sqts[1:ncut])
            norm_Sqt = real.(Sqt[:, :, end] ./ Sqt[:, :, 1])
            norm(norm_Sqt - norm_Sqt_ana) / L^2
        end

        @test all(diff(ncuterr) .< 0) # decreasing error
        @test ncuterr[end] < 1e-2
    end

    @testset "kagome lattice system" begin
        H = loadhamiltonian("../hamiltonians/kagome.dat", [1])
        
        T = 0.17
        L = 30
        dt = 0.001
        nt = 100
        thermal = 20

        # thermalization
        v = SCMC.randomstate(H.Ns, L)
        mcstep!(H, v, T, thermal)

        # simulate
        vs = simulate(H, v, dt, nt)

        # check energy conservation
        E0 = energy(H, v)
        Es = mapslices(vs, dims=[1, 2, 3]) do v
            energy(H, v)
        end
        @test all(Es .≈ E0)

        # check magnetization conservation
        m0 = magnetization(v)
        ms = reshape(mapslices(magnetization, vs; dims=[1, 2, 3]), (3, :))
        @test all(ms .≈ m0)

        # check FT of spins
        sq = SCMC.ftspacespins(H, vs)
        @test all(sq[:, 1, 1, :] .≈ m0)
    end
end
