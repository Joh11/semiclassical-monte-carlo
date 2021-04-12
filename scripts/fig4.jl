include("../src/scmc.jl")

using Plots
using Base.Threads: @threads, threadid, nthreads
using LaTeXStrings
using HDF5: h5write, h5open, attributes

@show nthreads()

const L = 144
const T = 0.17
const dt = 0.1
const nt = 80
const tstride = 100
const ts = dt * tstride * (0:nt-1)
const thermal = 20
const nsamples_per_thread = 14
const stride = 15
const nk = 7 # take the 7 first kpoints (h, h)

@assert nthreads() == 72 # so that we get ~ 1,000 samples

const nsamples = nsamples_per_thread * nthreads()

H = loadhamiltonian("hamiltonians/kagome.dat", [1])
Sqt = zeros(Complex{Float64}, L, L, nt)

@threads for n = 1:nthreads()
    v = randomstate(H.Ns, L)
    # thermalization
    mcstep!(H, v, T, thermal)
    println("thermalization $n finished")

    for i = 1:nsamples_per_thread
        time = @elapsed begin
            mcstep!(H, v, T, stride)
            vs = simulate(H, v, dt, nt; stride=tstride)
            S = structuralfactor(H, vs, dt)
            global Sqt += S ./ reshape(S[:, :, 1], (L, L, 1))
        end
        println("Done sample $i / $nsamples_per_thread for thread $n in $time seconds")
    end
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
