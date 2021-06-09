### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ ce493ae8-c850-11eb-0ff8-153f2e04f562
begin
	using HDF5
	using Plots
	
	plotly()
end

# ╔═╡ 2aabfed3-cc32-4928-a5af-a0e13e0ab24c
f = h5open("../data/skl_dyn_factor.h5", "r")

# ╔═╡ e33e3de8-7c80-4b1a-823d-36a970aa6c6d
begin
	Sqω = read(f["Sqω"])
end

# ╔═╡ 317462f9-68df-451f-9209-9dfe4bb94a4a
heatmap(abs.(Sqω[:, :, 10]))

# ╔═╡ 0ea8db17-7a34-4d5b-bc7a-b16eec157499
begin
	kxs = 8π * -1:0.05:1
	kys = 8π * -1:0.05:1
end

# ╔═╡ b30037f1-fe93-41d3-84b2-4157137b5893
function band_diagram(Sqω)
	size(Sqω, 1)
end

# ╔═╡ 705d5b29-6ae7-44c9-9000-874c672aba45
band_diagram(Sqω)

# ╔═╡ f44d135c-d0c9-4008-8fe5-c9a79d0486ad
Array(kxs)

# ╔═╡ Cell order:
# ╠═ce493ae8-c850-11eb-0ff8-153f2e04f562
# ╠═2aabfed3-cc32-4928-a5af-a0e13e0ab24c
# ╠═e33e3de8-7c80-4b1a-823d-36a970aa6c6d
# ╠═317462f9-68df-451f-9209-9dfe4bb94a4a
# ╠═0ea8db17-7a34-4d5b-bc7a-b16eec157499
# ╠═b30037f1-fe93-41d3-84b2-4157137b5893
# ╠═705d5b29-6ae7-44c9-9000-874c672aba45
# ╠═f44d135c-d0c9-4008-8fe5-c9a79d0486ad
