using LaTeXStrings
using Plots

include("../src/scmc.jl")

"Plot the energy and magnetization to determine the optimal number of thermalization steps"
function thermalization(H, L, T)
    v = randomstate(H.Ns, L)

    n = 50
    stride = 1
    
    E, m = zeros(n), zeros(3, n)

    println("1 / $n")
    E[1] = energy(H, v)
    m[:, 1] = magnetization(v)
    
    for i = 2:n
        println("$i / $n")
        mcstep!(H, v, T, stride)
        E[i] = energy(H, v)
        m[:, i] = magnetization(v)
    end

    plot(E;
         label = "E",
         xlabel="number of MC steps / $stride")
    plot!(transpose(m), label=[L"m_1" L"m_2" L"m_3"])
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
H = loadhamiltonian("../hamiltonians/bilayer-square.dat", [1])
L = 20
T = 0.02

# 1. check the thermalization step
thermalization(H, L, T)
# we can see that 20 steps is enough

# 2. check the decorrelation steps between two samples
# correlation(H, L, T; thermal=20)
# we can see that 15 steps is enough
