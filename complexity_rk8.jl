include("scmc.jl")

using Plots
using LinearAlgebra

const H = loadhamiltonian("hamiltonians/kagome.dat", [1])

function runone(L)
    dt = 0.1
    nt = 80
    v = randomstate(H.Ns, L)
    @time simulate(H, v, dt, nt)
end

function linear_regression(x, y)
    A = [x ones(length(x))]
    b = y

    inv(transpose(A) * A) * transpose(A) * b
end

function complexity_analysis()
    Ls = [2, 4, 6, 10, 20, 40, 60, 80, 100]
    nsamples = 10
    times = zeros(nsamples, length(Ls))

    # trigger the JIT
    runone(1)

    for iL = 1:length(Ls)
        L = Ls[iL]
        println("Doing L=$L ...")

        for n = 1:nsamples
            println("Sample $n / $nsamples")
            times[n, iL] = @elapsed runone(L)
        end
    end

    avg_time = reshape(mapslices(mean, times; dims=[1]), (:,))

    a, b = linear_regression(log.(Ls), log.(avg_time))
    plot()

    scatter!(Ls, avg_time;
             xaxis=:log,
             yaxis=:log)

    xlabel!("L")
    ylabel!("time [s]")

    # plot the linear regression
    Ls_lin = Array(1:100)
    plot!(Ls_lin, exp.(b) .* Ls_lin.^a;
          xaxis=:log,
          yaxis=:log,
          label="\$\\propto L^{$(round(a, digits=2))}\$")
end

function benchmark()
    # trigger the JIT
    @elapsed runone(1)
    
    L = 40
    nsamples = 10

    time = 0
    for n = 1:nsamples
        println("Doing sample $n / $nsamples ...")
        time += @elapsed runone(L)
    end
    time /= nsamples
    println("L = $L: $time seconds ($nsamples samples) -- ($(round(Int, time / 80 * 1000)) ms per RK8 step)")
end


benchmark()
# complexity_analysis()
