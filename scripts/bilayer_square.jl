using Plots
using Base.Threads
using HDF5
using SemiClassicalMonteCarlo

const gs = [0.1, 1, 2, 2.522, 3]

const comment = "simple bilayer, with only NN, no NNN, for dt=0.6"
const J1 = 1
const J2 = 0
const L = 20
const T = 0.02 # 1/50
const dt = 0.2
const nt = 100
const thermal = 10_000
const nsamples_per_chain = 100
const nchains = 10
const nsamples = nsamples_per_chain * nchains
const stride = 50

"Generate samples, compute the dyn. struct. factor, and average it"
function dynfactor(H, name="")
    Somega = zeros(Complex{Float64}, (L, L, nt))
    St = zeros(Complex{Float64}, (L, L, nt))
    
    @threads for n = 1:nchains
        v = randomstate(H.Ns, L)
        mcstep!(H, v, T, thermal) # thermalize

        for i = 1:nsamples_per_chain
            time = @elapsed begin
                mcstep!(H, v, T, stride)
                vs = simulate(H, v, dt, nt)
                Somega += frequencystructuralfactor(H, vs, dt)
                St += structuralfactor(H, vs)
                # continue the MCMC chain at the end of the time evolution
                v = vs[:, :, :, end]
            end
            println("Done sample $i / $nsamples_per_chain in $(time)s for chain $n / $nchains ($name)")
        end

    end

    Dict("Somega" => Somega / nsamples,
         "St" => St / nsamples)
end

fs = Array{Any}(undef, length(gs))
@threads for n in eachindex(gs)
    g = gs[n]
    println("Doing g = $g")
    # load the hamiltonian
    # abuse of notation: g is actually J3
    H = loadhamiltonian("hamiltonians/bilayer-square.dat", [J1, J2, g])

    fs[n] = dynfactor(H, "g = $g")
end

# save everything to a HDF5 file
output = "bilayer-fig1-02.h5"
println("Saving to $output ...")

h5open(output, "w") do f
    attributes(f)["comment"] = comment
    attributes(f)["J1"] = J1
    attributes(f)["J2"] = J2
    attributes(f)["L"] = L
    attributes(f)["T"] = T
    attributes(f)["dt"] = dt
    attributes(f)["nt"] = nt
    attributes(f)["thermal"] = thermal
    attributes(f)["nsamples_per_chain"] = nsamples_per_chain
    attributes(f)["nchains"] = nchains
    attributes(f)["nsamples"] = nsamples
    attributes(f)["stride"] = stride

    for (g, dict) in zip(gs, fs)
        group = create_group(f, "g=$g")
        for (name, value) in dict
            group[name] = value
        end
        attributes(group)["g"] = g
    end
end

println("Finished !")
