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

# ╔═╡ 2a04c69e-ad7e-11eb-238a-0b946693580e
begin
	using Plots
	using PlutoUI
	using HDF5
	using LaTeXStrings
end

# ╔═╡ 790a8e3f-1b56-47fc-8c9d-478fb4c82ec4
fs_paths = ["../data/skl-dimer-convergence/skl_dimer_1k.h5", "../data/skl-dimer-convergence/skl_dimer_10k.h5","../data/skl-dimer-convergence/skl_dimer_20k.h5", "../data/skl-dimer-convergence/skl_dimer_30k.h5", "../data/skl-dimer-convergence/skl_dimer_40k.h5", "../data/skl-dimer-convergence/skl_dimer_50k.h5", "../data/skl-dimer-convergence/skl_dimer_100k.h5"];

# ╔═╡ 9204ea5e-ef8e-4e9c-8074-062e4a0df2ac
md"Choose at which dataset to start:"

# ╔═╡ b716e78a-874b-4faa-8b63-eace88d16e3d
@bind first_fs Slider(1:length(fs_paths); show_value=true)

# ╔═╡ 5989056c-186a-4fbb-841a-e7b38d25a4fe
md"Starting at $(fs_paths[first_fs]). "

# ╔═╡ 64281da3-615f-4a7c-af46-2d13b03aebd3
# datasets
fs = map(fs_paths[first_fs:end]) do path
	h5open(path, "r")
end; fs[1]

# ╔═╡ d557ee10-d1a4-4d0e-aa2d-838e71dd507a
# retrieve number of sample
nsamples = [read(attributes(f)["nsamples_per_chain"]) * read(attributes(f)["nchains"]) for f in fs]

# ╔═╡ a12ca816-cd5a-4bc8-85e0-a0c65a26504e
# retrieve the temperatures (assume all the same)
Ts = read(attributes(fs[1])["Ts"]);

# ╔═╡ b5ab96a6-2a83-44ca-a96f-eedc48f3b92a
md"Selected temperatures to show:"

# ╔═╡ 7d25dd3f-03ef-4856-85a4-1c94bb55e6e1
@bind selected_Ts_str MultiCheckBox(["$n" => "$(round(T;digits=3))" for (n, T) in enumerate(Ts)]; select_all=true)

# ╔═╡ aaa1adb1-13ec-4390-adf3-0c12282e181f
selected_Ts_indices = map(x -> parse(Int, x), selected_Ts_str);

# ╔═╡ 1ba984c2-6779-4040-b0be-578a548f30e6
selected_Ts = [Ts[n] for n in selected_Ts_indices];

# ╔═╡ c941e2b0-4acd-43d5-90f4-fc2c479204c2
@bind obs_name Radio(["Energy", "Dimer", "Dimer²"])

# ╔═╡ c401fa51-9154-4bc5-9d73-c898fa56d389
# compute our observable
begin
	# obs = [read(f["$n/E"]) for f in fs, n in selected_Ts_indices]
	obs = obs_name == "Dimer" ? [read(f["$(n)/dimer"])[1] for f in fs, n in selected_Ts_indices] : 
	obs_name == "Energy" ? [read(f["$n/E"]) for f in fs, n in selected_Ts_indices] :
	obs_name == "Dimer²" ? [read(f["$n/dimer2"])[1, 100] for f in fs, n in selected_Ts_indices] : missing
end

# ╔═╡ 3a371291-24a9-461d-8328-2990240fe479
plot(nsamples, obs, 
	label=reshape(["T = $(round(T; digits=3))" for T in selected_Ts], (1, :)),
	xlabel="Number of samples",
	xaxis=:log,
	legend=:right,
	ylabel=obs_name, m=:diamond)

# ╔═╡ Cell order:
# ╠═2a04c69e-ad7e-11eb-238a-0b946693580e
# ╠═790a8e3f-1b56-47fc-8c9d-478fb4c82ec4
# ╟─9204ea5e-ef8e-4e9c-8074-062e4a0df2ac
# ╟─b716e78a-874b-4faa-8b63-eace88d16e3d
# ╟─5989056c-186a-4fbb-841a-e7b38d25a4fe
# ╠═64281da3-615f-4a7c-af46-2d13b03aebd3
# ╠═d557ee10-d1a4-4d0e-aa2d-838e71dd507a
# ╠═a12ca816-cd5a-4bc8-85e0-a0c65a26504e
# ╟─b5ab96a6-2a83-44ca-a96f-eedc48f3b92a
# ╠═7d25dd3f-03ef-4856-85a4-1c94bb55e6e1
# ╠═aaa1adb1-13ec-4390-adf3-0c12282e181f
# ╠═1ba984c2-6779-4040-b0be-578a548f30e6
# ╠═c941e2b0-4acd-43d5-90f4-fc2c479204c2
# ╠═c401fa51-9154-4bc5-9d73-c898fa56d389
# ╠═3a371291-24a9-461d-8328-2990240fe479
