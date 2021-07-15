"A script to plot various quantities over the system size"

using PyCall
using HDF5
using LinearAlgebra
using Statistics
using LaTeXStrings

@pyimport matplotlib.pyplot as plt

plt.rc("text", usetex=true)
plt.rc("font", family="TeX Gyre Adventor", size=14)
plt.rc("pdf", fonttype=42)

# loop over all files in the path
# check if uud or qsl
# retrieve the S(4π, 0) and correlation
# put it in the correct box
# retrieve also the system sizes
function load_data(dir)
    Ls_uud = []
    corr_uud = []
    fac_uud = []

    Ls_qsl = []
    corr_qsl = []
    fac_qsl = []
    
    for filename in readdir(dir)
        # determine if it's a UUD or QSL
        uudp = contains(filename, "uud")
        
        h5open(joinpath(dir, filename), "r") do f
            L = read(attributes(f)["L"])
            if uudp
                push!(Ls_uud, L)
                push!(corr_uud, compute_correlator(read(f["corr"])))
                push!(fac_uud, compute_fac(read(f["St0"]), L))
            else
                push!(Ls_qsl, L)
                push!(corr_qsl, compute_correlator(read(f["corr"])))
                push!(fac_qsl, compute_fac(read(f["St0"]), L))
            end
        end
    end

    # sort them
    p = sortperm(Ls_uud)
    Ls_uud = Ls_uud[p]
    corr_uud = corr_uud[p]
    fac_uud = fac_uud[p]

    p = sortperm(Ls_qsl)
    Ls_qsl = Ls_qsl[p]
    corr_qsl = corr_qsl[p]
    fac_qsl = fac_qsl[p]
    
    Ls_uud, corr_uud, fac_uud, Ls_qsl, corr_qsl, fac_qsl
end

function compute_correlator(corr)
    Ns = size(corr, 1)
    L = size(corr, 2)
    L½ = 1 + L÷2 # like that L=4 => L½= 3
    
    mean([corr[i, 1, 1, i, L½, L½, 1] for i = 1:Ns])
end

function compute_fac(St0, L)
    number_spins = 6L^2
    abs(St0[1]) / number_spins^2
end

function plot_correlator(Ls, corr; title="", savefig=nothing)
    plt.figure()
    plt.scatter(1 ./ Ls, abs.(corr))

    plt.title(title)
    plt.xlabel(L"1 / L")
    plt.ylabel(L"\left| \langle \mathbf{\hat S}_i \cdot \mathbf{\hat S}_j \rangle \right|")
    plt.tight_layout()
    
    if !isnothing(savefig)
        plt.savefig(savefig)
    end
end

function plot_fac(Ls, fac; title="", savefig=nothing)
    plt.figure()
    plt.scatter(1 ./ Ls .^2, fac)

    plt.title(title)
    plt.xlabel(L"1 / L^2")
    plt.ylabel(L"S(\mathbf{q} = (4\pi, 0))")
    plt.tight_layout()
    
    if !isnothing(savefig)
        plt.savefig(savefig)
    end
end

# main part of the script
dir = "../data/results/scaling/"
Ls_uud, corr_uud, fac_uud, Ls_qsl, corr_qsl, fac_qsl = load_data(dir)

plot_correlator(Ls_uud, corr_uud;
                title="UUD",
                savefig="../plots/scaling/scaling_corr_uud.png")
plot_correlator(Ls_qsl, corr_qsl;
                title="QSL",
                savefig="../plots/scaling/scaling_corr_qsl.png")

plot_fac(Ls_uud, fac_uud;
         title="UUD",
         savefig="../plots/scaling/scaling_fac_uud.png")
plot_fac(Ls_qsl, fac_qsl;
         title="QSL",
         savefig="../plots/scaling/scaling_fac_qsl.png")

plt.show()
