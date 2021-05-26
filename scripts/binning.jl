#=

Binning analysis

=#

using Plots
using LinearAlgebra
using Statistics
using SemiClassicalMonteCarlo

"Bin the 2n samples into n averages"
function bin(Qs)
    N = length(Qs)
    @assert N % 2 == 0
    n = round(Int, N / 2)
    ret = zeros(n)

    for k in eachindex(ret)
        ret[k] = 1/2 * (Qs[2k-1] + Qs[2k])
    end
    
    ret
end

function computeΔQ(Qs)
    N = length(Qs)
    Q0 = mean(Qs)
    
    √(1 / (N * (N-1)) * sum((Qs .- Q0).^2))
end

function computeτQ(ΔQs)
    1/2 * (()^2 - 1)
end

const H = loadhamiltonian("hamiltonians/skl.dat", [1, 1, 1])
const L = 4
const T = 0.01
const nthermal = 100_000
const l = 15 # number of steps
const nsamples = 2^l
const bs = bonds(H)
const nbonds = length(bs)

# thermalize
v = randomstate(H.Ns, L)
mcstep!(H, v, T, nthermal)

Qs = zeros(nsamples) # samples
ΔQs = zeros(l)
Di = zeros(nbonds, L, L) # tmp var for compute_dimers!

naccepted = 0

# generate samples
for n in eachindex(Qs)
    if n % 1024 == 0
        println("Doing $n / $(length(Qs))")
    end
    
    global naccepted += mcstep!(H, v, T)
    # evolve it 10s
    vs = simulate(H, v, 10, 2)
    global v = vs[:, :, :, end]
    
    # compute order parameter
    compute_dimers!(v, L, bs, Di)
    global Qs[n] = abs(skl_order_parameter(Di))
end

println("<Q> = $(mean(Qs))")
println("Acceptance rate r = $(naccepted / (nsamples * H.Ns * L^2))")

histogram(Qs,
          title="T = $T",
          xlabel="Q (abs of order param)")
readline()

# compute each estimate
for n in 1:l
    ΔQs[n] = computeΔQ(Qs)
    global Qs = bin(Qs)
end

# plot to see if there is saturation
display(scatter(ΔQs,
                xlabel="l",
                ylabel="ΔQ^(l)",
                legend=nothing,
                title="T = $T"))
