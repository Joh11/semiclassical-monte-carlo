#=

A script to plot S(q) over the whole (E)BZ

To get correct xticks, enlarge the window, up until the ±8π ticks
appear, but then enlarge even more until they are on the edges of the
x axis

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

function draw_line(p1, p2; args...)
    plt.plot([p1[1], p2[1]], [p1[2], p2[2]], "--", c="w")
end

function plot_sq(datapath; savefig=nothing, kpath=nothing)
    krange = nothing
    St0 = nothing
    
    h5open(datapath) do f
        # load data
        krange = round(Int, read(f["krange"]) / π) # get the k-range in π units
        L = read(attributes(f)["L"])
        number_spins = 6 * L^2
        St0 = read(f["St0"]) / number_spins^2 # normalized
    end
    
    plt.figure()
    plt.imshow(St0,
               extent=(-krange, krange, -krange, krange),
               cmap=:magma)
    # ticks
    ticks = Array(-krange:2:krange)
    plt.xticks(ticks, ["\$$tick\\pi\$" for tick in ticks])
    plt.yticks(ticks, ["\$$tick\\pi\$" for tick in ticks])

    if !isnothing(kpath) # draw the kpath
        for k in eachindex(kpath)
            # / π since we draw in π units
            draw_line(kpath[k] / π, kpath[k % length(kpath) + 1] / π)
        end
    end
    
    plt.colorbar()
    plt.tight_layout()

    if !isnothing(savefig)
        plt.savefig(savefig)
    end
end

# Néel
datapath = "../data/results/neel2_8pi.h5"
plot_sq(datapath;
        savefig="../plots/sq/sq_neel.pdf",
        kpath=[[0, 0], [2π, 0], [2π, 2π]])

# UUD
datapath = "../data/results/uud2_8pi.h5"
plot_sq(datapath;
        savefig="../plots/sq/sq_uud.pdf",
        kpath=[[0, 0], [4π, 0], [2π, 2π]])

# QSL
datapath = "../data/results/qsl_8pi.h5"
plot_sq(datapath;
        savefig="../plots/sq/sq_qsl.pdf")

plt.show()
