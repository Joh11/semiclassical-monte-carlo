module LoadHamiltonian

export loadparamhamiltonian, loadhamiltonian

function iscomment(line)
    line == "" || line[1] == '#'
end

struct ParamHamiltonian
    # Like Hamiltonian, but the couplings are not replaced yet
    Ns :: Int
    couplings :: Array{Array{Tuple{Int, Int, Int, Int}}}
    rs :: Array{Float64, 2} # shape: (2, Ns)
end

struct Hamiltonian
    Ns :: Int
    couplings :: Array{Array{Tuple{Int, Int, Int, Float64}}}
    rs :: Array{Float64, 2} # shape: (2, Ns)
end

function mkhamiltonian(paramH, couplings)
    Hamiltonian(paramH.Ns, map(paramH.couplings) do (cs)
                map(cs) do (i, j, s, c)
                (i, j, s, couplings[c])
                end
                end, paramH.rs)
end

function loadparamhamiltonian(path)
    open(path, "r") do io
        # compute the number of sites
        line = ""
        while iscomment(line)
            line = readline(io)
        end
        Ns = parse(Int, line)
        
        # parse site coordinates
        rs = zeros(2, Ns)
        
        let n = 1; while n <= Ns
            # skip comments
            line = readline(io)
            if iscomment(line) continue end
            rs[:, n] = [parse(Float64, x) for x in split(line)]
            n += 1
        end end
        
        # now make a NsxNs matrix of arrays
        bonds = fill(Tuple{Int, Int, Int, Int}[], Ns)
        
        for line in eachline(io)
            # skip comments
            if iscomment(line) continue end
            s1, i, j, s2, c = [parse(Int, x) for x in split(line)]
            push!(bonds[s1+1], (i+1, j+1, s2+1, c+1))
        end

        ParamHamiltonian(Ns, bonds, rs)
    end
end

function loadhamiltonian(path, couplings)
    mkhamiltonian(loadparamhamiltonian(path), couplings)
end

end
