using LaTeXStrings
using Plots
using LinearAlgebra
using SemiClassicalMonteCarlo
const SCMC = SemiClassicalMonteCarlo

"Plot the given observable to determine the optimal number of thermalization steps"
function thermalization(H, L, T, observable, name="Observable")
    n = 500
    stride = 1
    
    v = randomstate(H.Ns, L)
    Obs = zeros(n)

    println("1 / $n")
    Obs[1] = observable(v)

    for i = 2:n
        println("$i / $n")
        mcstep!(H, v, T, stride)
        Obs[i] = observable(v)
    end

    display(plot(Obs;
                 label = name,
                 xlabel="number of MC steps / $stride"))
end

"Plot the correlation between samples"
function correlation(H, L, T; thermal=0)
    n = 25
    v = randomstate(H.Ns, L)

    # thermalization
    if thermal > 0
        mcstep!(H, v, T, thermal)
    end

    v0 = copy(v)
    corr = zeros(n)
    for i = 1:n
        println("$i / $n")
        mcstep!(H, v, T, 1)
        corr[i] = correlation(v0, v)
    end

    plot(corr;
         xlabel="number of MC steps",
         ylabel="correlation")
end

# Load the hamiltonian
const H = loadhamiltonian("hamiltonians/kagome.dat", [1])
const L = 144
const T = 0.17

observable = function(v)
    # arbitrary points
    i1, j1 = 10, 10, 1
    i2, j2 = 1, 10, 2
    
    vs = reshape(v, (H.Ns, L, L, 1))
    sq = reshape(SCMC.ftspacespins(H, vs), (3, L, L))

    abs(sq[:, i1, j1] ⋅ sq[:, i2, j2])
end

# 1. check the thermalization step
thermalization(H, L, T, observable, "structure factor")
# we can see that 20 steps is enough

# 2. check the decorrelation steps between two samples
# correlation(H, L, T; thermal=20)
# we can see that 15 steps is enough
