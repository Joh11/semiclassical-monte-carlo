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
	using Revise
	
	push!(LOAD_PATH, "..")
	using SemiClassicalMonteCarlo
	const SCMC = SemiClassicalMonteCarlo
end

# ╔═╡ 79fb766c-3537-4c4b-aacc-355ab293d14c
f = h5open("../skl_dimer_4x4_highT2.h5", "r")

# ╔═╡ 552ccc21-bf70-40b6-be18-c219b01d48c3
# load all the variables
begin
	const L = read(attributes(f)["L"])
	const Ns = 6L^2
	const Ts = read(attributes(f)["Ts"])
	const J1 = read(attributes(f)["J1"])
	const J2 = read(attributes(f)["J2"])
	const J3 = read(attributes(f)["J3"])
	
	const H = loadhamiltonian("../hamiltonians/skl.dat", [J1, J2, J3])
	
	const nsamples_per_chain = read(attributes(f)["nsamples_per_chain"])
	const nchains = read(attributes(f)["nchains"])
	const nsamples = nsamples_per_chain * nchains
	
	const stride = read(attributes(f)["stride"])
	const comment = read(attributes(f)["comment"])
	
	Dict("L" => L, "Ns" => Ns, "Ts" => Ts, "nsamples_per_chain" => nsamples_per_chain,
		"nchains" => nchains, "nsamples" => nsamples, "stride" => stride, 
		"comment" => comment)
end

# ╔═╡ b1283663-56c9-4c12-a6d0-821ff7f50cc0
md"## Correlation sign plot

Choose the temperature:"

# ╔═╡ e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
@bind Tindex Radio(["$n" => "T = $(round(Ts[n]; digits=3))" for n in eachindex(Ts)])

# ╔═╡ 03dce486-d258-4eb8-bdb1-b3361fd7007c
const T = Ts[parse(Int, Tindex)]

# ╔═╡ 592219b3-2a6a-4863-ad72-97f6ab8bfc0c
const dimer = read(f["$Tindex/dimer"]);

# ╔═╡ d1322808-afce-4611-b223-d5d964243960
const dimer2 = read(f["$Tindex/dimer2"]);

# ╔═╡ ab926e03-f10b-474c-b267-c54ca9f3cd0c
md"Choose the reference bond (shown in black on the next plot): "

# ╔═╡ 23367950-201e-468c-a6ab-e9402a487c10
@bind ref_bond_str Radio(["$n" for n in 1:12])

# ╔═╡ 02f9308b-a02b-48ca-88fd-092ccc750c0b
const ref_bond = parse(Int, ref_bond_str)

# ╔═╡ 314899e2-75a2-42a8-864c-76d5a3f561c9
# correlation wrt the ref bond
const corr = dimer2[ref_bond, :] - dimer[ref_bond] * dimer;

# ╔═╡ 3e6bb03a-723a-457c-8f78-ef5923dc9a5e
const E = read(f["$Tindex/E"]);

# ╔═╡ bc3e3d14-040e-4b77-b376-0c84afe378d5
"Takes an Hamiltonian, and a (12, L, L) array of bond strengths"
function plot_diagram(H::SCMC.Hamiltonian, values)
	bs = bonds(H)
	L = size(values, 2)
	plot()
	
	function draw(r1, r2, val; color=:auto)
		x1, y1 = r1
		x2, y2 = r2
		
		plot!([x1, x2], [y1, y2],
			legend=nothing,
			color=color != :auto ? color : (val > 0 ? :red : :blue),
			linewidth=100abs(val))
	end
	
	for x0 = 0:L-1
		for y0 = 0:L-1
			for (n, bond) in enumerate(bs)
				draw(bond.a.pos + [x0, y0], 
					bond.b.pos + [x0, y0],
					values[n, x0+1, y0+1])
			end
		end
	end
	
	# reference bond
	draw(bs[ref_bond].a.pos, bs[ref_bond].b.pos, values[ref_bond, 1, 1]; color=:black)
	
	plot!(aspect_ratio=:equal,
		xticks=nothing,
		yticks=nothing) # to make sure it renders
end

# ╔═╡ 162831b0-73d0-46c6-95fc-2317fd961cf2
begin
	# plot_diagram(H, reshape(corr, (12, L, L)) / corr[ref_bond])
	plot_diagram(H, reshape(corr, (12, L, L)))
	plot!(title="correlation for T=$(round(T; digits=3))")
end

# ╔═╡ fc665d1f-638c-4055-8440-6ea8e2a8553d
findmax(corr)

# ╔═╡ 5d7e180c-1fa6-4670-a49d-ca0c8a327186
histogram(corr,
	title="Distribution of the correlation",
	legend=nothing,
	xlabel=L"<D_iD_j> - <D_i><D_j>")

# ╔═╡ 0186d3a7-ba9d-4ed1-9853-fdf8d05b3804
md"## Order parameter"

# ╔═╡ 3f1cfc3c-aa15-47c3-b436-042e435a91e0
const order_param = [read(f["$n/order_param"]) for n in eachindex(Ts)]

# ╔═╡ 878e6432-55af-46f6-9b9c-6acd8321ad35
scatter(Ts, order_param / 2Ns,
	title="order parameter for 4x4",
	xlabel="T",
	ylabel=L"<|O|> / N_{bonds}",
	xaxis=:log,
	legend=nothing)

# ╔═╡ 8f06d92c-2229-4578-8cff-c774d1d0b4a6
function sign_structure(L)
    Nbonds = 12
    @assert L % 2 == 0
    
    # θ for each different bond
    red = :red
    blue = :blue
    green = :green
    pink = :pink

    # for the two kinds of unit cell (bc symmetry breaking)
    θ1 = [blue, red, 0, 0, red, 0, 0, green, 0, blue, 0, 0]
    θ2 = [0, 0, green, pink, 0, green, pink, 0, pink, 0, green, pink]

    # now compute the order parameter
    ret = Array{Any}(undef, (Nbonds, L, L))
    for x = 1:L
        for y = 1:L
            if (x + y) % 2 == 0
                ret[:, x, y] = θ1
            else
				ret[:, x, y] = θ2
            end
        end
    end

    ret
end

# ╔═╡ bd3dbe56-edfa-44e4-a17d-47b6e1cf020a
sign_structure(4)

# ╔═╡ 58674407-d1ec-4ac6-8df0-9948e3889b20
"Takes an Hamiltonian, and a (12, L, L) array of colored bonds"
function plot_bonds(H::SCMC.Hamiltonian, colors)
	bs = bonds(H)
	L = size(colors, 2)
	plot()
	
	function draw(r1, r2, color)
		x1, y1 = r1
		x2, y2 = r2
		
		if color == 0
			return nothing
		end
		
		plot!([x1, x2], [y1, y2],
			legend=nothing,
			color=color,
			linewidth=2)
	end
	
	for x0 = 0:L-1
		for y0 = 0:L-1
			for (n, bond) in enumerate(bs)
				draw(bond.a.pos + [x0, y0], 
					bond.b.pos + [x0, y0],
					colors[n, x0+1, y0+1])
			end
		end
	end
	
	plot!(aspect_ratio=:equal,
		xticks=nothing,
		yticks=nothing) # to make sure it renders
end

# ╔═╡ 193d3495-c394-4f9e-aac6-86c8fdd7edec
plot_bonds(H, sign_structure(4))

# ╔═╡ Cell order:
# ╠═f217d8a6-ac0f-11eb-3614-e1b3fd065e68
# ╠═79fb766c-3537-4c4b-aacc-355ab293d14c
# ╠═552ccc21-bf70-40b6-be18-c219b01d48c3
# ╟─b1283663-56c9-4c12-a6d0-821ff7f50cc0
# ╠═e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
# ╟─03dce486-d258-4eb8-bdb1-b3361fd7007c
# ╠═592219b3-2a6a-4863-ad72-97f6ab8bfc0c
# ╠═d1322808-afce-4611-b223-d5d964243960
# ╟─ab926e03-f10b-474c-b267-c54ca9f3cd0c
# ╠═23367950-201e-468c-a6ab-e9402a487c10
# ╟─02f9308b-a02b-48ca-88fd-092ccc750c0b
# ╠═314899e2-75a2-42a8-864c-76d5a3f561c9
# ╠═3e6bb03a-723a-457c-8f78-ef5923dc9a5e
# ╠═bc3e3d14-040e-4b77-b376-0c84afe378d5
# ╠═162831b0-73d0-46c6-95fc-2317fd961cf2
# ╠═fc665d1f-638c-4055-8440-6ea8e2a8553d
# ╠═5d7e180c-1fa6-4670-a49d-ca0c8a327186
# ╟─0186d3a7-ba9d-4ed1-9853-fdf8d05b3804
# ╠═3f1cfc3c-aa15-47c3-b436-042e435a91e0
# ╠═878e6432-55af-46f6-9b9c-6acd8321ad35
# ╠═8f06d92c-2229-4578-8cff-c774d1d0b4a6
# ╠═bd3dbe56-edfa-44e4-a17d-47b6e1cf020a
# ╠═58674407-d1ec-4ac6-8df0-9948e3889b20
# ╠═193d3495-c394-4f9e-aac6-86c8fdd7edec
