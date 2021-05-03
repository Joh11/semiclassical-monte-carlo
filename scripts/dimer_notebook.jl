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
end

# ╔═╡ 05606f63-e841-42a6-94b6-ad5bb68d1b47
f = h5open("../skl_dimer.h5", "r")

# ╔═╡ 552ccc21-bf70-40b6-be18-c219b01d48c3
# load all the variables
begin
	const L = read(attributes(f)["L"])
	const Ns = 12L^2
	const Ts = read(attributes(f)["Ts"])
end

# ╔═╡ b1283663-56c9-4c12-a6d0-821ff7f50cc0
md"## Checking that the dimer values make sense"

# ╔═╡ e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
@bind Tindex Radio(["$n" => "T = $(round(Ts[n]; digits=3))" for n in eachindex(Ts)])

# ╔═╡ d1322808-afce-4611-b223-d5d964243960
const dimer2 = read(f["$Tindex/dimer2"])

# ╔═╡ 314899e2-75a2-42a8-864c-76d5a3f561c9
# show only <D_i D_i> for i in UC
const DiDi = [dimer2[i, i] for i = 1:12]

# ╔═╡ 4e4f2a77-127b-488f-b760-34886e2c3946
scatter(DiDi,
	xlabel="i",
	ylabel=L"<D_i D_i>",
	xticks=1:12,
	legend=nothing)

# ╔═╡ Cell order:
# ╠═f217d8a6-ac0f-11eb-3614-e1b3fd065e68
# ╠═05606f63-e841-42a6-94b6-ad5bb68d1b47
# ╠═552ccc21-bf70-40b6-be18-c219b01d48c3
# ╠═b1283663-56c9-4c12-a6d0-821ff7f50cc0
# ╠═e7022954-f8ad-4f52-bcd3-bcde7d4cfb1b
# ╠═d1322808-afce-4611-b223-d5d964243960
# ╠═314899e2-75a2-42a8-864c-76d5a3f561c9
# ╠═4e4f2a77-127b-488f-b760-34886e2c3946
