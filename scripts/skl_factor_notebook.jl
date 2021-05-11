### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 76541f94-ae5e-11eb-0a69-19c54f8171fc
begin
	using Plots
	using PlutoUI
	using HDF5
	using LaTeXStrings
	using FFTW
	
	push!(LOAD_PATH, "..")
	using SemiClassicalMonteCarlo
	
	const ℂ = Complex{Float64}
end

# ╔═╡ 701dd7b9-423d-49ab-846a-53fef921cc48
f = h5open("../skl_factor_fast.h5", "r")

# ╔═╡ f4a8548c-1aea-4a15-be67-96354f927aaf
begin
	L = read(attributes(f)["L"])
	Ns = 6L^2
	nsamples_per_chain = read(attributes(f)["nsamples_per_chain"])
	nchains = read(attributes(f)["nchains"])
	nsamples = nsamples_per_chain * nchains
end

# ╔═╡ 6222058c-a9ad-476c-b59f-114a82dd5b33
# static structure factor
corr = real.(read(f["corr"])) # real. is temporary

# ╔═╡ 155c9cea-516e-4f91-8c70-7cf53d3dbab1
H = loadhamiltonian("../hamiltonians/skl.dat", [1, 1, 1])

# ╔═╡ e5d2707a-7d74-4b98-9863-4eec4bdab86f
# construct the structure factor
begin
	Sq = zeros(ℂ, 8L, 8L) # goes from -8π to 8π
	# fill it
	for i = 1:8L, j = 1:8L
		k = 2π * [i-1, j-1] .- 8π
		Sq[i, j] = structurefactor(H, corr, k)
		println("doing $i $j")
	end
end

# ╔═╡ 2cb689eb-4d63-4b82-84ec-f405bf5b715d
heatmap(static,
	xlabel=L"k_x",
	xticks=nothing,
	ylabel=L"k_y",
	yticks=nothing)

# ╔═╡ 6c5c58ca-e1da-4930-b5f3-fe713c78b649
static[1, 1]

# ╔═╡ Cell order:
# ╠═76541f94-ae5e-11eb-0a69-19c54f8171fc
# ╠═701dd7b9-423d-49ab-846a-53fef921cc48
# ╠═f4a8548c-1aea-4a15-be67-96354f927aaf
# ╠═6222058c-a9ad-476c-b59f-114a82dd5b33
# ╠═155c9cea-516e-4f91-8c70-7cf53d3dbab1
# ╠═e5d2707a-7d74-4b98-9863-4eec4bdab86f
# ╠═2cb689eb-4d63-4b82-84ec-f405bf5b715d
# ╠═6c5c58ca-e1da-4930-b5f3-fe713c78b649
