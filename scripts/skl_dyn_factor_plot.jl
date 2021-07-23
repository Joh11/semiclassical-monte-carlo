"A script to plot the dynamical structure factor"

using PyCall
using HDF5
using LinearAlgebra
using LaTeXStrings

@pyimport matplotlib.pyplot as plt

plt.rc("text", usetex=true)
plt.rc("font", family="TeX Gyre Adventor", size=14)
plt.rc("pdf", fonttype=42)

# load hdf5 file
function load_data(datapath)
    f = h5open(datapath, "r")
    # for vars that are used
    L = read(attributes(f)["L"])
    number_spins = 6L^2
    
    Dict(:L => L,
         :nt => read(attributes(f)["nt"]),
         :dt => read(attributes(f)["dt"]),
         :Sω => read(f["Sqω"]) / number_spins^2,
         # little hack to make it optional
         :St0 => haskey(f, "St0") ? read(f["St0"]) / number_spins^2 : nothing,
         :J1 => read(attributes(f)["J1"]))
end

# plot static structure factor
function plot_stat(; St0, title="", savefig=nothing, kwargs...)
    plt.figure()

    plt.imshow(abs.(St0),
               origin=:lower,
               extent=(-4, 4, -4, 4),
               cmap=:magma)
    if title != ""
        plt.title("\$S(\\vec q)\$ " * title)
    end
    
    ticks = Array(-4:2:4)
    plt.xticks(ticks, [string(tick) * L"\pi" for tick in ticks])
    plt.yticks(ticks, [string(tick) * L"\pi" for tick in ticks])

    plt.colorbar()
    plt.tight_layout()

    if !isnothing(savefig)
        plt.savefig(savefig)
    end
end

# plot dynamical structure factor
function plot_dyn(; Sω, dt, kpath, title="", maxlim=nothing, savefig=nothing, kwargs...)
    nt = size(Sω, 2)
    nk = length(kpath)

    plt.figure()
    # go up to ω_½
    # skip the continuous part
    plt.imshow(abs.(@views Sω[:, 2:nt÷2]'),
               origin=:lower,
               extent=(0, 1, 0, 2 / dt),
               aspect=:auto,
               cmap=:magma)
    if !isnothing(maxlim)
        plt.clim(0, maxlim)
    end
    plt.xticks(Array(0:(nk-1)) / (nk - 1), kpath)
    if title != ""
        plt.title("\$S(\\vec q, t)\$ " * title)
    end
    plt.gca().yaxis.set_major_formatter("{x:1.2f}\$\\pi\$")
    plt.ylabel(L"\omega / J")
    plt.colorbar()
    plt.tight_layout()

    if !isnothing(savefig)
        plt.savefig(savefig)
    end
end

uud_kpath = [L"\Gamma", L"(2\pi, 0)", L"(4\pi, 0)", L"(2\pi, 2\pi)", L"\Gamma"]
neel_kpath = [L"\Gamma", L"(2\pi, 0)", L"(2\pi, 2\pi)", L"\Gamma"]

# UUD
# plot_dyn(;kpath=uud_kpath,
#          title="for UUD phase",
#          maxlim=0.005,
#          savefig="../plots/sω_uud.png",
#          load_data("../data/results/uud_20x20.h5")...)
# UUD BZ
# plot_dyn(;kpath=[L"\Gamma", L"X", L"M", L"\Gamma"],
#          title="for UUD no phase shift",
#          maxlim=0.01,
#          savefig="../plots/sω_uud_bz.png",
#          load_data("../data/results/uud_bz.h5")...)

# Neel
# plot_dyn(;kpath=neel_kpath,
#          title="for Néel phase",
#          maxlim=0.5,
#          load_data("../data/results/neel_reverse_T=0.1.h5")...)


# plot UUD -> QSL scan
function plot_uud_qsl_scan()
    path = "../data/results/uud_qsl_scan/"
    for filename in readdir(path)
        data = load_data(path * filename)
        # no title to put them in latex
        # S(q, ω)
        plot_dyn(;kpath=uud_kpath,
                 # title="for \$J_1 = $(data[:J1])\$",
                 maxlim=0.1,
                 savefig="../plots/uud_qsl_scan/sω_$(data[:J1]).png",
                 data...)
        # S(q)
        plot_stat(;
                  # title="for \$J_1 = $(data[:J1])\$",
                  savefig="../plots/uud_qsl_scan/sq_$(data[:J1]).png",
                  data...)
    end
end

plot_uud_qsl_scan()
plt.show()
