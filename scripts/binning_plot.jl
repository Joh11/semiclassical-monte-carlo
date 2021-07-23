#=

A script to plot both the long range correlation and S(4π, 0) as a
function of the system size

=#

using PyCall
using HDF5
using LinearAlgebra
using Statistics
using LaTeXStrings

@pyimport matplotlib.pyplot as plt

plt.rc("text", usetex=true)
plt.rc("font", family="TeX Gyre Adventor", size=14)
plt.rc("pdf", fonttype=42)

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

function plotΔQ(Qs)
    ΔQs = []
    push!(ΔQs, computeΔQ(Qs))
    
    while length(Qs) > 1
        Qs = bin(Qs)
        push!(ΔQs, computeΔQ(Qs))
    end

    plt.figure()
    plt.scatter(0:length(ΔQs)-1, ΔQs)
    plt.show()
end

function load_corr_and_fac(datapaths)
    corrs = []
    facs = []
    for dp in datapaths
        f = h5open(dp)
        push!(corrs, mean(read(f["corr"])))
        
        L = read(attributes(f)["L"])
        N = 6L^2
        push!(facs, mean(read(f["fac"])) / N^2)
    end
    
    corrs, facs
end

function plot_corr(Ls, corrs, corrs_errors; savefig=nothing, interp=false)
    plt.figure()

    xs = 1 ./ Ls
    plt.errorbar(xs, abs.(corrs), yerr=corrs_errors, linestyle="None", marker=".")
    plt.ylim(bottom=0)

    plt.xlabel(L"1 / L")
    plt.ylabel(L"\left| \langle \mathbf{\hat S}_i \cdot \mathbf{\hat S}_j \rangle \right|")

    if interp # draw the line
        α, β = linfit(xs, corrs)
        println("α=$α, β=$β")
        m = maximum(xs)
        draw_line([0, β], [m, α * m + β])
    end
    
    if !isnothing(savefig)
        plt.savefig(savefig)
    end

    plt.tight_layout()
end

function plot_fac(Ls, fac, facs_errors; savefig=nothing, interp=false)
    plt.figure()

    xs = 1 ./ Ls.^2
    
    plt.errorbar(xs, fac, yerr=facs_errors, linestyle="None", marker=".")
    plt.ylim(bottom=0)

    plt.xlabel(L"1 / L^2")
    plt.ylabel(L"S(4\pi, 0)")
    
    if interp # draw the line
        α, β = linfit(xs, fac)
        println("α=$α, β=$β")
        m = maximum(xs)
        draw_line([0, β], [m, α * m + β])
    end

    if !isnothing(savefig)
        plt.savefig(savefig)
    end

    plt.tight_layout()
end

function draw_line(p1, p2; args...)
    plt.plot([p1[1], p2[1]], [p1[2], p2[2]], "--")
end

function linfit(x, y)
    n = length(x)
    # construct the A matrix
    A = zeros(n, 2)
    A[:, 1] = x
    A[:, 2] .= 1

    # solve Au ≈ y ⇔ AᵀAu = Aᵀy
    u = inv(A' * A) * A' * y
    α, β = u

    α, β
end

Ls = [6, 8, 10, 15, 20]


# for UUD
datapaths = ["../data/scaling/scaling_uud_$(L).h5" for L in Ls]
corrs_errors = [0.0005, 0.0007, 0.00065, 0.00013, 0.0014]
facs_errors = [0.00016, 0.00025, 0.00022, 0.0004, 0.00035]

# load data
corrs, fac = load_corr_and_fac(datapaths)

# plot it
plot_corr(Ls, corrs, corrs_errors;
          interp=true,
          savefig="../plots/binning/corr_uud.pdf")
plot_fac(Ls, fac, facs_errors;
         interp=true,
         savefig="../plots/binning/fac_uud.pdf")




# for QSL
datapaths = ["../data/scaling/scaling_qsl_$(L).h5" for L in Ls]
corrs_errors = [0.0018, 0.0017, 0.0018, 0.0017, 0.0018]
facs_errors = [0.00011, 6.5e-5, 4e-5, 1.8e-5, 1e-5]

# load data
corrs, fac = load_corr_and_fac(datapaths)

# plot it
plot_corr(Ls, corrs, corrs_errors;
          interp=true,
          savefig="../plots/binning/corr_qsl.pdf")
plot_fac(Ls, fac, facs_errors;
         interp=true,
         savefig="../plots/binning/fac_qsl.pdf")

plt.show()

# to estimate the errors by hand
# for (dp, L) in zip(datapaths, Ls)
#     N = 6L^2
#     f = h5open(dp)
#     println("Binning for $dp")
#     # plotΔQ(read(f["corr"]))
#     plotΔQ(read(f["fac"]) / N^2) # normalize
# end
