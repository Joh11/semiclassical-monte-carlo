# just to mess with extended BZ
# to be able to run long computations outside of a notebook

using HDF5
using Plots
using HDF5

f = h5open("skl_factor_fast.h5", "r")
L = read(attributes(f)["L"])
corr = real.(read(f["corr"])) # real. is temporary

H = loadhamiltonian("hamiltonians/skl.dat", [1, 1, 1])

function make_arrays(H, L)
    Ns = H.Ns
    N = Ns * L^2
    Ri = zeros(2, N^2)
    Rj = zeros(2, N^2)

    n = 1
    for yj = 1:L, xj = 1:L, sj = 1:Ns
        @views rj = H.lattice * [xj, yj] + H.rs[:, sj]
        for yi = 1:L, xi = 1:L, si = 1:Ns
            @views ri = H.lattice * [xi, yi] + H.rs[:, si]
            Ri[:, n] = ri
            Rj[:, n] = rj
            n += 1
        end
    end
    
    Ri, Rj
end


Ri, Rj = make_arrays(H, L)

N = H.Ns * L^2
kxs = Array(range(-8π, 8π; length=100))
kys = Array(range(-8π, 8π; length=100))

Lkx, Lky = length(kxs), length(kys)
chi = zeros(Complex{Float64}, Lkx, Lky)

for i = 1:Lkx
    for j = 1:Lky
        k = [kxs[i] kxs[j]]
        println("doing $i $j, $k")
        chi[i, j] = sum(reshape(corr, :) .* reshape(exp.(1im * k * (Ri - Rj)), :))
    end
end

heatmap(real.(chi) / N)
xticks!(1 .+ 99 .* Array(0:1/16:1), ["$(n)π" for n in -8:8])
yticks!(1 .+ 99 .* Array(0:1/16:1), ["$(n)π" for n in -8:8])
title!("\$S(\\vec k, t=0)\$")
