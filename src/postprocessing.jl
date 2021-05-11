
"Computes the structure factor at the momentum k. "
function structurefactor(H, corr, k)
    L = size(corr, 2)
    Ns = H.Ns

    ret = 0
    
    for yj = 1:L, xj = 1:L, sj = 1:Ns
        @views rj = H.lattice * [xj, yj] + H.rs[:, sj]
        for yi = 1:L, xi = 1:L, si = 1:Ns
            @views ri = H.lattice * [xi, yi] + H.rs[:, si]
            ret += corr[si, xi, yi, sj, xj, yj] * exp(1im * k ⋅ (ri - rj))
        end
    end

    ret
end

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

    for j = 1:Lky
        for i = 1:Lkx
            k = [kxs[i] kxs[j]]
            println("doing $i $j, $k")
            chi[i, j] = sum(reshape(corr, :) .* reshape(exp.(1im * k * Rs), :))
        end
    end

    chi
end
