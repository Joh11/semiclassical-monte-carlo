# Structure factor
# ================

"Computes the structure factor at all the momentum given by the grid kxs, kys"
function structurefactors(H, corr, kxs, kys)
    L = size(corr, 2)
    Ns = H.Ns
    N = Ns * L^2

    # first step: compute the positions
    Rs = zeros(2, N^2)
    n = 1
    for yj = 1:L, xj = 1:L, sj = 1:Ns
        @views rj = H.lattice * [xj, yj] + H.rs[:, sj]
        for yi = 1:L, xi = 1:L, si = 1:Ns
            @views ri = H.lattice * [xi, yi] + H.rs[:, si]
            Rs[:, n] = ri - rj
            n += 1
        end
    end

    # second step: compute the structure factor itself
    Lkx, Lky = length(kxs), length(kys)
    chi = zeros(Complex{Float64}, Lkx, Lky)

    @threads for j = 1:Lky
        println("Doing $j / $Lky ...")
        for i = 1:Lkx
            k = [kxs[i] kxs[j]]
            chi[i, j] = sum(reshape(corr, :) .* reshape(exp.(1im * k * Rs), :))
        end
    end

    chi
end

# Dimer order
# ===========

"Inplace version of compute_dimers"
function compute_dimers!(v, L, bonds, dimer)
    for i = 1:L
        for j = 1:L
            for (n, bond) in enumerate(bonds)
                dimer[n, i, j] = v[bond.a.index, i, j] ⋅ v[bond.b.index,
                                                           wrapindex(i + bond.b.i, L),
                                                           wrapindex(j + bond.b.j, L)]
            end
        end
    end
end

"""Compute the dimer operators D_i for each bond of each unit cell. 
Returns a (12, L, L) array of floats. """
function compute_dimers(v, H)
    L = size(v, 2)
    dimer = zeros(12, L, L)
    compute_dimers!(v, L, bonds(H), dimer)
    dimer
end

"In place version of compute_dimer2"
function compute_dimer2!(dimer2!, dimer, nbonds, L)
    for i = 1:nbonds
        @simd for j = 1:nbonds*L^2
            dimer2![i, j] = dimer[i, 1, 1] * reshape(dimer, :)[j]
        end
    end    
end

"Takes a (nbonds, L, L) array as input"
function compute_dimer2(dimer, nbonds, L)
    dimer2 = zeros(nbonds, nbonds * L^2)
    compute_dimer2!(dimer2, dimer, nbonds, L)
    dimer2
end

"""Takes a dimer array (Nbonds, L, L), for L|2, and returns the order
parameter (a scalar)"""
function skl_order_parameter(dimers)
    L = size(dimers, 2)
    Nbonds = size(dimers, 1)
    @assert L % 2 == 0
    @assert Nbonds == 12
    
    # θ for each different bond
    red = 1
    blue = -1
    green = 1
    pink = -1

    # for the two kinds of unit cell (bc symmetry breaking)
    θ1 = [blue, red, 0, 0, red, 0, 0, green, 0, blue, 0, 0]
    θ2 = [0, 0, green, pink, 0, green, pink, 0, pink, 0, green, pink]

    # now compute the order parameter
    order = 0
    for x = 1:L
        for y = 1:L
            if (x + y) % 2 == 0
                @views order += θ1 ⋅ dimers[:, x, y]
            else
                @views order += θ2 ⋅ dimers[:, x, y]
            end
        end
    end

    order
end
