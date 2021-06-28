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

# ╔═╡ ce493ae8-c850-11eb-0ff8-153f2e04f562
begin
	using PlutoUI
	using HDF5
	using Plots
	using Statistics
	
	plotly()
end

# ╔═╡ d9f5224d-913a-42b9-ae3d-204ccad8223a
pathnames = ["../data/results/skl_dyn_factor_high_res.h5" => "default one",
		     "../data/results/skl_dyn_factor_zoomed.h5" => "zoomed 4x",
			 "../data/results/skl_dyn_factor_unzoomed.h5" => "unzoomed 4x",
			 "../data/results/skl_dyn_factor_uud.h5" => "UUD",
             "../data/results/skl_dyn_factor_neel.h5" => "Néel"]

# ╔═╡ eeed684e-e53f-4b93-ba45-567d372c9e39
@bind path Radio(pathnames, default=first(pathnames[1]))

# ╔═╡ c30ed7e4-6726-4de3-abed-a9eecbdb6067
path

# ╔═╡ 2aabfed3-cc32-4928-a5af-a0e13e0ab24c
f = h5open(path, "r")

# ╔═╡ e33e3de8-7c80-4b1a-823d-36a970aa6c6d
begin
	Sqω = read(f["Sqω"])
	Sqt = read(f["Sqt"])
	T = read(attributes(f)["T"])
	kpath = read(f["kpath"])
	Dict("dt" => read(attributes(f)["dt"]),
		 "nt" => read(attributes(f)["nt"]))
end

# ╔═╡ 317462f9-68df-451f-9209-9dfe4bb94a4a
Sqω

# ╔═╡ 85e0cc07-58a6-43b7-a394-e574055da765
kpath

# ╔═╡ a9fc0eb5-1cca-4706-b176-b67392b0d91d
heatmap(abs.(Sqω)'[1:100, 1:end],
	xlabel="kpath",
	ylabel="ω",
	xticks=([1, 101, 201, 301, 400], ["Γ", "X", "(2π, 0)", "M", "Γ"]))

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

# ╔═╡ 705d5b29-6ae7-44c9-9000-874c672aba45
band_diagram_large(Sqω)

# ╔═╡ f44d135c-d0c9-4008-8fe5-c9a79d0486ad
scatter(kpath[1, 1:250], kpath[2, 1:250])

# ╔═╡ Cell order:
# ╠═ce493ae8-c850-11eb-0ff8-153f2e04f562
# ╠═d9f5224d-913a-42b9-ae3d-204ccad8223a
# ╠═eeed684e-e53f-4b93-ba45-567d372c9e39
# ╠═c30ed7e4-6726-4de3-abed-a9eecbdb6067
# ╠═2aabfed3-cc32-4928-a5af-a0e13e0ab24c
# ╠═e33e3de8-7c80-4b1a-823d-36a970aa6c6d
# ╠═317462f9-68df-451f-9209-9dfe4bb94a4a
# ╠═85e0cc07-58a6-43b7-a394-e574055da765
# ╠═a9fc0eb5-1cca-4706-b176-b67392b0d91d
# ╠═b30037f1-fe93-41d3-84b2-4157137b5893
# ╠═c84c7955-ccd0-4840-8fab-0a06cd8ddcbe
# ╠═705d5b29-6ae7-44c9-9000-874c672aba45
# ╠═f44d135c-d0c9-4008-8fe5-c9a79d0486ad
