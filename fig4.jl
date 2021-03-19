include("scmc.jl")
using Plots
using Base.Threads: @spawn, threadid, nthreads
using LaTeXStrings
using HDF5: h5write, h5open, attributes

@show nthreads()

# const Ls = [1, 20, 50, 100]
const Ls = [2, 20]
const nL = length(Ls)
const T = 1
const dt = 0.1
const nt = 100
const ts = dt * (0:nt-1)
const thermal = 20
const nsamples = 500
const stride = 15

H = loadhamiltonian("hamiltonians/kagome.dat", [1])

S00 = zeros(Complex{Float64}, nt, nL)
Spipi = zeros(Complex{Float64}, nt, nL)
S0pi = zeros(Complex{Float64}, nt, nL)
Spi0 = zeros(Complex{Float64}, nt, nL)

function runsim(v, i, n)
    npi = 1 + round(Int, Ls[i] / 2)
    vs = simulate(H, v, dt, nt)
    Sqt = structuralfactor(H, vs, dt)
    global S00[:, i] += Sqt[1, 1, :] / Sqt[1, 1, 1]
    global Spipi[:, i] += Sqt[npi, npi, :] / Sqt[npi, npi, 1]
    global Spi0[:, i] += Sqt[npi, 1, :] / Sqt[npi, 1, 1]
    global S0pi[:, i] += Sqt[1, npi, :] / Sqt[1, npi, 1]
    println("Done simulating the $n th iL=$i in thread $(threadid())")
end

# simulation itself
@sync for i = 1:nL
    @spawn begin
        L = Ls[i]

        # thermalization
        v = randomstate(H.Ns, L)
        mcstep!(H, v, T, thermal)
        println("thermalization for L = $L finished")
        
        for n in 1:nsamples
            mcstep!(H, v, T, stride)
            @spawn runsim(v, i, n)
            
            println("L = $L: sample $n done (thread $(threadid()))")
        end
        println("L = $L finished")
    end
end

S00 /= nsamples
Spipi /= nsamples
S0pi /= nsamples
Spi0 /= nsamples

function plotSqt(Sqt; title="Sqt")
    display(plot())
    ylims!(-0.1, 1.1)
    title!(title)
    
    for i = 1:nL
        L = Ls[i]
        display(plot!(ts, real.(Sqt[:, i]), label="L=$L"))
    end
end

output = "fig4.h5"
h5open(output, "w") do f
    attributes(f)["Ls"] = Ls
    attributes(f)["T"] = T
    attributes(f)["dt"] = dt
    attributes(f)["nt"] = nt
    attributes(f)["thermal"] = thermal
    attributes(f)["nsamples"] = nsamples
    attributes(f)["stride"] = stride
    
    write(f, "S00", S00)
    write(f, "Spipi", Spipi)
    write(f, "S0pi", S0pi)
    write(f, "Spi0", Spi0)
end
# plot(S00, Ls; title=L"S_{\vec 0}(t) / S_{\vec 0}(0)")
