#=

A script to plot the scalar chirality as a
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

function load_sc(datapaths)
    sc = []
    for dp in datapaths
        h5open(dp) do f
            push!(sc, read(f["sc"]))
        end
    end
    sc
end

Ls = [6]

# for QSL
datapaths = ["../data/scalar_chirality/sc$(L).h5" for L in Ls]

# load data
sc = load_sc(datapaths)

# plot it
plt.figure()

plt.scatter(Ls, sc)
    
plt.show()
