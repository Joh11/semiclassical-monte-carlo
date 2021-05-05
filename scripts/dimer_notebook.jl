### A Pluto.jl notebook ###
# v0.14.4

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

# ╔═╡ f217d8a6-ac0f-11eb-3614-e1b3fd065e68
begin
	using HDF5
	using PlutoUI
	using Plots
	using LaTeXStrings
	using Statistics
end

# ╔═╡ 79fb766c-3537-4c4b-aacc-355ab293d14c
f = h5open("../skl_dimer_long.h5", "r")

# ╔═╡ 552ccc21-bf70-40b6-be18-c219b01d48c3
# load all the variables
begin
	const L = read(attributes(f)["L"])
	const Ns = 12L^2
	const Ts = read(attributes(f)["Ts"])
end

# ╔═╡ b1283663-56c9-4c12-a6d0-821ff7f50cc0
md"## Correlation sign plot"

# ╔═╡ e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
@bind Tindex Radio(["$n" => "T = $(round(Ts[n]; digits=3))" for n in eachindex(Ts)])

# ╔═╡ 03dce486-d258-4eb8-bdb1-b3361fd7007c
const T = Ts[parse(Int, Tindex)]

# ╔═╡ 592219b3-2a6a-4863-ad72-97f6ab8bfc0c
const dimer = read(f["$Tindex/dimer"])

# ╔═╡ d1322808-afce-4611-b223-d5d964243960
const dimer2 = read(f["$Tindex/dimer2"])

# ╔═╡ 23367950-201e-468c-a6ab-e9402a487c10
@bind ref_bond_str Radio(["$n" for n in 1:12])

# ╔═╡ 02f9308b-a02b-48ca-88fd-092ccc750c0b
const ref_bond = parse(Int, ref_bond_str)

# ╔═╡ 314899e2-75a2-42a8-864c-76d5a3f561c9
# correlation wrt the ref bond
const corr = dimer2[ref_bond, :] - dimer[ref_bond] * dimer

# ╔═╡ 3e6bb03a-723a-457c-8f78-ef5923dc9a5e
const E = read(f["$Tindex/E"])

# ╔═╡ bc3e3d14-040e-4b77-b376-0c84afe378d5
"Takes a (12, L, L) array of bond strengths"
function plot_diagram(bonds)
	L = size(bonds, 2)
	plot()
	
	function draw(x1, y1, x2, y2, val; color=:auto)
		plot!([x1, x2], [y1, y2],
			legend=nothing,
			color=color != :auto ? color : (val > 0 ? :red : :blue),
			linewidth=abs(val))
	end
	
	for x0 = 0:L-1
		for y0 = 0:L-1
			draw(x0+.25, y0+.75, x0+.5, y0+1, bonds[1, x0+1, y0+1])
			draw(x0+.5, y0+1, x0+.75, y0+.75, bonds[2, x0+1, y0+1])
			draw(x0+.25, y0+.75, x0+.75, y0+.75, bonds[3, x0+1, y0+1])
			draw(x0, y0+.5, x0+.25, y0+.75, bonds[4, x0+1, y0+1])
			draw(x0+.25, y0+.75, x0+.25, y0+.25, bonds[5, x0+1, y0+1])
			draw(x0+.75, y0+.75, x0+.75, y0+.25, bonds[6, x0+1, y0+1])
			draw(x0+.75, y0+.75, x0+1, y0+.5, bonds[7, x0+1, y0+1])
			draw(x0, y0+.5, x0+.25, y0+.25, bonds[8, x0+1, y0+1])
			draw(x0+.75, y0+.25, x0+1, y0+.5, bonds[9, x0+1, y0+1])
			draw(x0+.25, y0+.25, x0+.75, y0+.25, bonds[10, x0+1, y0+1])
			draw(x0+.25, y0+.25, x0+.5, y0, bonds[11, x0+1, y0+1])
			draw(x0+.5, y0, x0+.75, y0+.25, bonds[12, x0+1, y0+1])
		end
	end
	
	# reference bond
	# TODO draw it automatically
	draw(.25, .75, .5, 1, bonds[ref_bond, 1, 1]; color=:black)
	
	plot!(aspect_ratio=:equal,
			title="correlation for T=$(round(T; digits=3))") # to make sure it renders
end

# ╔═╡ 59b79a21-a592-46b5-9719-9a8a39f4a76e


# ╔═╡ 162831b0-73d0-46c6-95fc-2317fd961cf2
plot_diagram(reshape(corr, (12, L, L)) / corr[ref_bond])

# ╔═╡ a31d8d78-baa8-469d-9d5e-7c9401ca320f
const self_corr = [dimer2[i, i] for i=1:12] - dimer[1:12] .^2

# ╔═╡ 73a122b0-0246-4a89-9e54-b01d34b962eb
scatter(self_corr)

# ╔═╡ 2fbabd80-2613-4d57-a4b7-a2606da28d1c
[dimer2[i, i] for i=1:12], dimer[1:12]

# ╔═╡ fc4c1e7f-8e5d-444a-8dfc-b634ebd9681f
mean(corr), std(corr), findmin(corr), findmax(corr)

# ╔═╡ a70fe460-be56-451d-bb61-5b008f73d6af
# try to see if the correlations decay
heatmap(reshape(corr, (12, L, L))[1, :, :],
	)

# ╔═╡ Cell order:
# ╠═f217d8a6-ac0f-11eb-3614-e1b3fd065e68
# ╠═79fb766c-3537-4c4b-aacc-355ab293d14c
# ╠═552ccc21-bf70-40b6-be18-c219b01d48c3
# ╟─b1283663-56c9-4c12-a6d0-821ff7f50cc0
# ╠═e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
# ╠═03dce486-d258-4eb8-bdb1-b3361fd7007c
# ╠═592219b3-2a6a-4863-ad72-97f6ab8bfc0c
# ╠═d1322808-afce-4611-b223-d5d964243960
# ╠═23367950-201e-468c-a6ab-e9402a487c10
# ╠═02f9308b-a02b-48ca-88fd-092ccc750c0b
# ╠═314899e2-75a2-42a8-864c-76d5a3f561c9
# ╠═3e6bb03a-723a-457c-8f78-ef5923dc9a5e
# ╠═bc3e3d14-040e-4b77-b376-0c84afe378d5
# ╠═59b79a21-a592-46b5-9719-9a8a39f4a76e
# ╠═162831b0-73d0-46c6-95fc-2317fd961cf2
# ╠═a31d8d78-baa8-469d-9d5e-7c9401ca320f
# ╠═73a122b0-0246-4a89-9e54-b01d34b962eb
# ╠═2fbabd80-2613-4d57-a4b7-a2606da28d1c
# ╠═fc4c1e7f-8e5d-444a-8dfc-b634ebd9681f
# ╠═a70fe460-be56-451d-bb61-5b008f73d6af
