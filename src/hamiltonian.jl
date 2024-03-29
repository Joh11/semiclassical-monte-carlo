@inline function wrapindex(i, L)
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
    couplings :: Array{Array{Tuple{Int, Int, Int, Int}, 1}, 1}
    rs :: Array{Float64, 2} # shape: (2, Ns)
    # (lattice[:, 1] is the first lattice vector)
    lattice :: Array{Float64, 2} # shape: (2, 2)
end

struct Hamiltonian
    Ns :: Int
    couplings :: Array{Array{Tuple{Int, Int, Int, Float64}, 1}, 1}
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
        bonds = [Tuple{Int, Int, Int, Int}[] for i = 1:Ns]
        
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

@inline function localfield(H, v, i, j, s)
    L = size(v)[2]
    ret = @MVector zeros(3)

    for (Δi, Δj, s2, c) in H.couplings[s]
        ret .+= c .* v[s2, wrapindex(i + Δi, L), wrapindex(j + Δj, L)]
    end
    
    ret
end

"Compute the energy of the given state"
function energy(H, v)
    L = size(v)[3]
    # simply the sum of spin . local field
    E = 0
    
    for s in 1:H.Ns
        for j in 1:L
            for i in 1:L
                E += v[s, i, j] ⋅ localfield(H, v, i, j, s)
            end
        end
    end
    E
end

"Faster way to compute the variation of energy when updating a single spin"
function deltaenergy(H, v, S, i, j, s)
    ΔS = S - v[s, i, j]
    # WARNING this only holds for symmetric couplings (Jij = Jji)
    2 .* ΔS ⋅ localfield(H, v, i, j, s)
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

struct Site
    index :: Int
    # unit cell indices
    i :: Int
    j :: Int
    # can be outside of [0,1]² (if i, j != 0)
    pos :: Array{Float64, 1}
end

struct Bond
    a :: Site
    b :: Site
    strength :: Float64
end

"Returns a list of the bonds in the hamiltonian"
function bonds(H::Hamiltonian)
    ret = Bond[]

    # only take bonds from a to b > a
    for a in eachindex(H.couplings)
        siteA = Site(a, 0, 0, H.rs[:, a])
        
        for (Δi, Δj, b, strength) in H.couplings[a]
            if b <= a # skip to avoid double counting
                continue
            end
            siteB = Site(b, Δi, Δj, H.rs[:, b] + H.lattice * [Δi, Δj])

            push!(ret, Bond(siteA, siteB, strength))
        end
    end
        
    ret
end
