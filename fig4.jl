include("scmc.jl")
using Plots
using Base.Threads: @spawn, threadid, nthreads
using LaTeXStrings
using HDF5: h5write, h5open, attributes

@show nthreads()

# const Ls = [144]
const L = 144
const T = 0.17
const dt = 0.1
const nt = 80
const tstride = 100
const ts = dt * (0:nt-1)
const thermal = 20
const nsamples = 1000
const stride = 15
const nk = 7 # take the 7 first kpoints (h, h)

H = loadhamiltonian("hamiltonians/kagome.dat", [1])

Sqt = zeros(Complex{Float64}, L, L, nt)

function runsim(v, n)
    vs = simulate(H, v, dt, nt; stride=tstride)
    S = structuralfactor(H, vs, dt)
    global Sqt += S ./ reshape(S[:, :, 1], (L, L, 1))
    println("Done simulating the $(n)th sample in thread $(threadid())")
end

@sync begin
    # thermalization
    v = randomstate(H.Ns, L)
    mcstep!(H, v, T, thermal)
    println("thermalization finished")
    
    for n in 1:nsamples
        mcstep!(H, v, T, stride)
        @spawn runsim(v, n)
        
        println("sample $n done (thread $(threadid()))")
    end
    println("finished")
end

Sqt /= nsamples

# save data
output = "fig4.h5"
h5open(output, "w") do f
    attributes(f)["L"] = L
    attributes(f)["T"] = T
    attributes(f)["dt"] = dt
    attributes(f)["nt"] = nt
    attributes(f)["tstride"] = tstride
    attributes(f)["thermal"] = thermal
    attributes(f)["nsamples"] = nsamples
    attributes(f)["stride"] = stride
    
    write(f, "Sqt", Sqt)
end

# plot data
plot(xlabel=L"t",
     ylabel=L"S(q, t) / S(q, 0)")
for ik = 1:nk
    q = sqrt(3) / L * (ik - 1)
    plot!(ts, real.(Sqt[ik, ik, :]), label="|q|=$(round(q; digits=3))")
end
ylims!(-0.1, 1.1)
savefig("fig4.png")
