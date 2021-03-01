using LaTeXStrings

include("scmc.jl")

"Plot the energy and magnetization to determine the optimal number of thermalization steps"
function thermalization(H, L, T)
    v = randomstate(H.Ns, L)

    n = 100
    stride = 1
    
    E, m = zeros(n), zeros(3, n)

    E[1] = energy(H, v)
    m[:, 1] = magnetization(v)
    
    for i = 2:n
        println("$i / $n \n")
        mcstep!(H, v, T, stride)
        E[i] = energy(H, v)
        m[:, i] = magnetization(v)
    end

    plot(E;
         label = "E",
         xlabel="number of MC steps / $stride")
    plot!(transpose(m), label=[L"m_1" L"m_2" L"m_3"])
end

H = loadhamiltonian("hamiltonians/kagome.dat", [1])
# H = loadhamiltonian("hamiltonians/square.dat", [1, 0])
thermalization(H, 20, 0.17)

# plot(E, label = "E", xticks = (1:10, 0:10:100), xlabel="number of MC steps")
# plot!(transpose(m), label=[L"m_1" L"m_2" L"m_3"])
