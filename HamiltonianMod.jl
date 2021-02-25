module HamiltonianMod

using LinearAlgebra

export loadparamhamiltonian, loadhamiltonian, localfield, energy, reciprocallattice

function wrapindex(i, L)
    1 + mod(i - 1, L)
end

function iscomment(line)
    line == "" || line[1] == '#'
end

function skipcommentreadline(io)
    line = ""
    while iscomment(line)
        line = readline(io)
    end
    line
end

struct ParamHamiltonian
    # Like Hamiltonian, but the couplings are not replaced yet
    Ns :: Int
    couplings :: Array{Array{Tuple{Int, Int, Int, Int}}}
    rs :: Array{Float64, 2} # shape: (2, Ns)
    # (lattice[:, 1] is the first lattice vector)
    lattice :: Array{Float64, 2} # shape: (2, 2)
end

struct Hamiltonian
    Ns :: Int
    couplings :: Array{Array{Tuple{Int, Int, Int, Float64}}}
    rs :: Array{Float64, 2} # shape: (2, Ns)
    # (lattice[:, 1] is the first lattice vector)
    lattice :: Array{Float64, 2} # shape: (2, 2) 
end

"Build an Hamiltonian from a ParamHamiltonian, and the couplings"
function mkhamiltonian(paramH, couplings)
    Hamiltonian(paramH.Ns, map(paramH.couplings) do (cs)
                map(cs) do (i, j, s, c)
                (i, j, s, couplings[c])
                end
                end,
                paramH.rs,
                paramH.lattice)
end

"Load a ParamHamiltonian from the path. "
function loadparamhamiltonian(path)
    open(path, "r") do io
        line = ""

        # parse the lattice vectors
        lattice = zeros(2, 2)
        lattice[:, 1] = [parse(Float64, x) for x in split(skipcommentreadline(io))]
        lattice[:, 2] = [parse(Float64, x) for x in split(skipcommentreadline(io))]
        
        # compute the number of sites
        Ns = parse(Int, skipcommentreadline(io))
        
        # parse site coordinates
        rs = zeros(2, Ns)
        let n = 1; while n <= Ns
            # skip comments
            line = skipcommentreadline(io)
            rs[:, n] = [parse(Float64, x) for x in split(line)]
            n += 1
        end end
        # put them in cartesian coords
        rs = lattice * rs
        
        # now make a NsxNs matrix of arrays
        bonds = fill(Tuple{Int, Int, Int, Int}[], Ns)
        
        for line in eachline(io)
            # skip comments
            if iscomment(line) continue end
            s1, i, j, s2, c = [parse(Int, x) for x in split(line)]
            push!(bonds[s1+1], (i, j, s2+1, c+1))
        end

        ParamHamiltonian(Ns, bonds, rs, lattice)
    end
end

"Load a ParamHamiltonian from the path, and convert it to an Hamiltonian using couplings"
function loadhamiltonian(path, couplings)
    mkhamiltonian(loadparamhamiltonian(path), couplings)
end

function localfield(H, v, i, j, s)
    L = size(v)[3]
    sum(map(H.couplings[s]) do (Δi, Δj, s2, c)
        c * v[:, s2,
              wrapindex(i + Δi, L),
              wrapindex(j + Δj, L)]
        end)
end

"Compute the energy per site of the given state"
function energy(H, v)
    L = size(v)[3]
    # simply the sum of spin . local field
    E = 0
    
    for s in 1:H.Ns
        for j in 1:L
            for i in 1:L
                S = v[:, s, i, j]
                h = localfield(H, v, i, j, s)
                E += dot(S, h)
            end
        end
    end

    E / (H.Ns * L^2)
end

"Takes a 2x2 matrix (lattice[:, 1] is the first vector) and returns
the reciprocal lattice. "
function reciprocallattice(lattice)
    a1, a2 = lattice[:, 1], lattice[:, 2]

    det = a1[1] * a2[2] - a1[2] * a2[1]

    b1 = [a2[2], -a2[1]]
    b2 = [-a1[2], a1[1]]

    ret = zeros(2, 2)
    ret[:, 1] = b1
    ret[:, 2] = b2
    
    (2π / det) * ret
end

end
