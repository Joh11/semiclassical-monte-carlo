### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ ce493ae8-c850-11eb-0ff8-153f2e04f562
begin
	using HDF5
	using Plots
	using Statistics
	
	plotly()
end

# ╔═╡ 2aabfed3-cc32-4928-a5af-a0e13e0ab24c
f = h5open("../data/skl_dyn_factor_fixed.h5", "r")

# ╔═╡ e33e3de8-7c80-4b1a-823d-36a970aa6c6d
begin
	Sqω = read(f["Sqω"])
	Sqt = read(f["Sqt"])
	T = read(attributes(f)["T"])
end

# ╔═╡ 317462f9-68df-451f-9209-9dfe4bb94a4a
heatmap(abs.(Sqt[:, :, 1]))

# ╔═╡ 0ea8db17-7a34-4d5b-bc7a-b16eec157499
begin
	kxs = 8π * Array(-1:0.01:1)
	kys = 8π * Array(-1:0.01:1)
end

# ╔═╡ b30037f1-fe93-41d3-84b2-4157137b5893
function band_diagram_standard(Sqω)
	N = size(Sqω, 1)
	n = round(Int, (N-1)/2)
	nt = size(Sqω, 3)

	ret = zeros(3n, nt)
	ret[1:n, :] = abs.(Sqω[n+1:end-1, n+1, :])
	ret[n+1:2n, :] = abs.(Sqω[end, n+1:end-1, :])
	for t = 1:nt
		ret[2n+1:end, t] = abs.([Sqω[k, k, t] for k = N:-1:n+2])
	end
	
	ret
end

# ╔═╡ c84c7955-ccd0-4840-8fab-0a06cd8ddcbe
function band_diagram_large(Sqω)
	N = size(Sqω, 1)
	n = round(Int, (N-1)/2)
	nhalf = round(Int, (N-1)/4)
	nquart = round(Int, (N-1)/8)
	nt = size(Sqω, 3)

	ret = zeros(4nquart, nt)
	ret[1:nhalf, :] = abs.(Sqω[n+1:n+nhalf, n+1, :])
	for t = 1:nt
		ret[nhalf+1:nhalf+nquart, t] = abs.([Sqω[n+nhalf+2-k, n+k, t] for k = 1:nquart])
	end
	for t = 1:nt
		ret[nhalf+nquart+1:end, t] = abs.([Sqω[n+nquart+2-k, n+nquart+2-k, t] for k = 1:nquart])
	end
	
	ret
end

# ╔═╡ 76efa531-8055-4253-a0ba-a367821cb579
40//2

# ╔═╡ 705d5b29-6ae7-44c9-9000-874c672aba45
band_diagram_large(Sqω)

# ╔═╡ f44d135c-d0c9-4008-8fe5-c9a79d0486ad
#=
heatmap((band_diagram_standard(Sqω)'[2:50, :]), 
	xticks=([1, 101, 201, 300], ["Γ", "X", "M", "Γ"]),
	yticks=(0:10:100, ["2π×$(n/100)" for n = 0:10:100]),
	ylabel="ω",
	title="dynamical structure factor"
)
=#

# ╔═╡ 01908ab5-69a5-45d3-bddf-f2c20adf402f
heatmap((band_diagram_large(Sqω)'[2:50, :]), 
	xticks=([1, 6, 11, 16, 20], ["Γ", "X", "(2π, 2π)", "M", "Γ"]),
	yticks=(0:10:100, ["2π×$(n/100)" for n = 0:10:100]),
	ylabel="ω",
	title="dynamical structure factor"
)

# ╔═╡ 621a2a26-8c4b-4625-a4df-04330af035e6
Sqω

# ╔═╡ 41304844-7ed6-4687-81b4-ece512c29d36
heatmap(reshape(real.(mean(Sqω; dims=3)), 41, 41),
	xticks=([1, 21, 41], ["-8π", "0", "8π"]),
	yticks=([1, 21, 41], ["-8π", "0", "8π"]))

# ╔═╡ Cell order:
# ╠═ce493ae8-c850-11eb-0ff8-153f2e04f562
# ╠═2aabfed3-cc32-4928-a5af-a0e13e0ab24c
# ╠═e33e3de8-7c80-4b1a-823d-36a970aa6c6d
# ╠═317462f9-68df-451f-9209-9dfe4bb94a4a
# ╠═0ea8db17-7a34-4d5b-bc7a-b16eec157499
# ╠═b30037f1-fe93-41d3-84b2-4157137b5893
# ╠═c84c7955-ccd0-4840-8fab-0a06cd8ddcbe
# ╠═76efa531-8055-4253-a0ba-a367821cb579
# ╠═705d5b29-6ae7-44c9-9000-874c672aba45
# ╠═f44d135c-d0c9-4008-8fe5-c9a79d0486ad
# ╠═01908ab5-69a5-45d3-bddf-f2c20adf402f
# ╠═621a2a26-8c4b-4625-a4df-04330af035e6
# ╠═41304844-7ed6-4687-81b4-ece512c29d36
