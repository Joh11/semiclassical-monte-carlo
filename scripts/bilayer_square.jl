include("../src/scmc.jl")

using Plots
using Base.Threads: @threads, threadid, nthreads
using HDF5: h5write, h5open, attributes

const gs = [0.1, 1, 2, 2.522, 3]

const L = 20
const T = 0.02 # 1/50
const dt = 0.1
const nt = 100
const tstride = 10
const ts = dt * tstride * (0:nt-1)
const thermal = 
const nsamples = 1000
const stride = 15
const nk = 7 # take the 7 first kpoints (h, h)

"Generate samples, compute the dyn. struct. factor, and average it"
function dynfactor(H, name="")
    v = randomstate(H.Ns, L)
    # thermalize
    v = mcstep!(H, v, T)

    for i = 1:nsamples
        time = @elapsed begin
            mcstep!(H, v, T, stride)
            vs = simulate(H, v, dt, nt; stride=tstride)
            S += frequencystructuralfactor(H, vs, dt)
        end
        println("Done sample $i / $nsamples_per_thread for task $name in $(time)s")
    end
    S / nsamples
end

fs = []
for g in gs
    println("Doing g = $g")
    # load the hamiltonian
    # assume J = 1 => J' = g
    H = loadhamiltonian("hamiltonians/bilayer-square.dat", [1, g])

    f = dynfactor(H)
    push!(fs, f)
end

# save everything to a HDF5 file
output = "bilayer-fig1.h5"
h5open(output, "w") do f
    attributes(f)["L"] = L
    attributes(f)["T"] = T
    attributes(f)["dt"] = dt
    attributes(f)["nt"] = nt
    attributes(f)["tstride"] = tstride
    attributes(f)["thermal"] = thermal
    attributes(f)["nsamples"] = nsamples
    attributes(f)["stride"] = stride

    for (g, f) in zip(gs, fs)
        group = create_group(f, "g=$g")
        group["Sqomega"] = f
        attributes(group)["g"] = g
    end
end


