### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b3aad6f6-a9f1-11eb-307e-d1f21cd78c7b
begin
	using PlutoUI
	using Plots
	using HDF5
end

# ╔═╡ d0ead2c6-25fc-4274-bffc-894378e05c2b
md"# Replicating the paper \"Dynamical structure factors and excitation modes of the bilayer Heisenberg model\""

# ╔═╡ ae4116a7-5594-4300-be63-69c967d6247c
md"## Replicating figure 1, top row"

# ╔═╡ 9922da8b-6883-4e48-add2-7da92a9ab1e5
f = h5open("../data/bilayer-fig1.h5", "r") # load data file

# ╔═╡ 0103eb95-d5cf-4988-a8d0-55154dd1d5aa
# load all the variables
begin
	const L = read(attributes(f)["L"])
	const Ns = 2 * L^2
	const nt = read(attributes(f)["nt"])
	const dt = read(attributes(f)["dt"])
	const T = nt * dt
	Dict("L" => L,
	"Ns" => Ns,
	"nt" => nt,
	"dt" => dt,
	"T" => T)
end

# ╔═╡ 39fde87d-17f6-4757-8a9f-2c9b2a9f5d48
"Select only the values of S along a standard Γ-M-X-Γ path"
function plot_diagram(S, L, nt)
    # build kpath
    @assert L % 2 == 0
    Lhalf = floor(Int, L / 2)
    nkps = 3 * Lhalf + 1
    kpath = zeros(Int, 2, nkps)

    for i = 0:Lhalf
        # Gamma - (pi, 0)
        kpath[:, i + 1] = 1 .+ [i, 0]
        # (pi, 0) - (pi, pi)
        kpath[:, Lhalf + 1 + i] = 1 .+ [Lhalf, i]
        # (pi, 0) - Gamma
        kpath[:, 2*Lhalf + 1 + i] = 1 .+ (Lhalf - i) * [1, 1]
    end

    ret = zeros(Complex{Float64}, nt, nkps)
    for nk in 1:nkps
        kx, ky = kpath[:, nk]
        ret[:, nk] = S[kx, ky, :]
    end

    ret
end

# ╔═╡ 7b46eff8-fe7f-46da-afc3-d160b887c308
md"Use the following to cycle between the possible interlayer coupling strengths: "

# ╔═╡ 95134fbc-dd03-460a-8cca-0f03a8dba318
@bind g Radio(keys(f))

# ╔═╡ 160d9e55-bed0-440b-9b57-5af5aad16416
# structure factor over the full BZ
Sqomega_bz = read(f["$g/Somega"]) / Ns; nothing

# ╔═╡ 8ae451f1-4d21-4cc6-9093-e05c563bdc2f
St = read(f["$g/St"]) / Ns; nothing

# ╔═╡ 23d58720-25b3-46e4-92d5-7a8c5bc03348
Sqomega = plot_diagram(Sqomega_bz, L, nt); nothing

# ╔═╡ fbebc16e-b58e-4d55-b28c-175591ce6987
heatmap(real.(Sqomega[1:30, :]),
	xticks=([1, 11, 21, 31], ["Γ", "(π, 0)", "(π, π)", "Γ"]),
	yticks=(T / 2π * Array(0:9), 0:9),
	ylabel="ω",
	aspect_ratio=31/30,
	clims=(0, 1), # can be uncommented to fix the colorbar scale for all plots
	title="symmetrical dyn. structure factor, $g"
)

# ╔═╡ e9d4e665-f8a0-4c7c-9554-47614f90a4e4
md"## How did I fix the ugly double bands I was obtaining

Contrary to what I was thinking, the frequencies do not take values in \$[0, 2\pi[\$, but it \$[0, 2\pi / dt[\$ (of course, dimensions...). This means that with a smaller timestep, I can reach greater frequencies (sounds obvious now...). 

Probably due to the symmetry of the correlation, the negative and positive frequencies have the same dynamical structure factor values. And we know that the excitations we should get have bounded energy. So if the frequency range is big enough, the positive and negative frequency bands should be totally disconnected. 

I first tried with \$dt = 0.1\$, which worked, but the frequency range was too large compared to the effective spectrum, so I doubled the timestep to 0.2. "

# ╔═╡ d08064e9-3c98-4ec4-9d57-fa94cb1aff71
md"## Note about the ±ω convention

In my code, since it is done by FFTW, both the spatial and temporal FTs are done with a minus sign in the exponent, that is \$e^{-i\omega t}\$ or \$e^{-i\vec q \cdot \vec r}\$. However, using [this paper](https://arxiv.org/pdf/1508.07816.pdf)'s definition of the dynamical structure factor, we need to reverse it for the frequencies to \$e^{i\omega t}\$.

I could have fixed this discrepancy by reordering the components of my array correctly, however this is error prone. If we look at what we are FTing, the correlation between two spins, since it is symmetrical, the structure factor \$S(\vec q, \omega)\$ computed with my -ω convention is equivalent to the structure factor \$S(-\vec q, \omega)\$ in the normal convention. Since we have inversion symmetry it should be equivalent. "

# ╔═╡ 2e287fd1-48df-4664-acb3-7e96c9b87e39
md"## Note about scale

Even if the heatmaps seem qualitatively correct, we do not have the same magnitude for the dynamical structure factor. Couple things could cause this:
- classical / quantum correspondance: this should add a \$3/4\$ factor to \$\vec S^2\$ right ?
- could it be influenced by our choice of time steps \$dt\$, \$nt\$ ?
- global normalization constant in the definition: I think I did this correctly, dividing by the number of sites \$Ns\$
- I'm taking the real value of each component. Perhaps the modulus could do better"

# ╔═╡ Cell order:
# ╟─d0ead2c6-25fc-4274-bffc-894378e05c2b
# ╟─ae4116a7-5594-4300-be63-69c967d6247c
# ╠═b3aad6f6-a9f1-11eb-307e-d1f21cd78c7b
# ╠═9922da8b-6883-4e48-add2-7da92a9ab1e5
# ╠═0103eb95-d5cf-4988-a8d0-55154dd1d5aa
# ╟─39fde87d-17f6-4757-8a9f-2c9b2a9f5d48
# ╟─7b46eff8-fe7f-46da-afc3-d160b887c308
# ╠═95134fbc-dd03-460a-8cca-0f03a8dba318
# ╠═160d9e55-bed0-440b-9b57-5af5aad16416
# ╠═8ae451f1-4d21-4cc6-9093-e05c563bdc2f
# ╠═23d58720-25b3-46e4-92d5-7a8c5bc03348
# ╠═fbebc16e-b58e-4d55-b28c-175591ce6987
# ╟─e9d4e665-f8a0-4c7c-9554-47614f90a4e4
# ╟─d08064e9-3c98-4ec4-9d57-fa94cb1aff71
# ╟─2e287fd1-48df-4664-acb3-7e96c9b87e39
