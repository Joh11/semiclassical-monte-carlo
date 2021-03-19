include("scmc.jl")
using Plots
using Base.Threads: @threads, nthreads
using LaTeXStrings
using HDF5: h5write

@show nthreads()

Ls = [1, 20, 50, 100]
nL = length(Ls)
T = 1
dt = 0.1
nt = 1000
thermal = 20
nsamples = 100
stride = 15

H = loadhamiltonian("hamiltonians/kagome.dat", [1])

S00 = zeros(Complex{Float64}, nt, nL)
Spipi = zeros(Complex{Float64}, nt, nL)
S0pi = zeros(Complex{Float64}, nt, nL)
Spi0 = zeros(Complex{Float64}, nt, nL)

# simulation itself
@threads for i = 1:nL
    L = Ls[i]

    # thermalization
    v = randomstate(H.Ns, L)
    mcstep!(H, v, T, thermal)

    println("thermalization for L = $L finished")
    
    for n in 1:nsamples
        mcstep!(H, v, T, stride)
        vs = simulate(H, v, dt, nt)
        Sqt = structuralfactor(H, vs, dt)
        
        S00[:, i] += Sqt[1, 1, :] / Sqt[1, 1, 1]
        Spipi[:, i] += Sqt[end, end, :] / Sqt[end, end, 1]
        Spi0[:, i] += Sqt[end, 1, :] / Sqt[end, 1, 1]
        S0pi[:, i] += Sqt[1, end, :] / Sqt[1, end, 1]
    end
    println("L = $L finished")
end

S00 /= nsamples
Spipi /= nsamples
S0pi /= nsamples
Spi0 /= nsamples

function plotSqt(Sqt, Ls; title="Sqt")
    nL = size(Ls)
    plot()
    for i = 1:nL
        L = Ls[i]
        plot!(real.(Sqt[:, i]), label="L=$L")
    end
end

output = "fig4.h5"
h5write(output, "S00", S00)
h5write(output, "Spipi", Spipi)
h5write(output, "S0pi", S0pi)
h5write(output, "Spi0", Spi0)
# plot(S00, Ls; title=L"S_{\vec 0}(t) / S_{\vec 0}(0)")
