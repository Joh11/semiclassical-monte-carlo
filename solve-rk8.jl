include("scmc.jl")

using Plots

H = loadhamiltonian("hamiltonians/kagome.dat", [1])

function compute_diff(H, dt, t)
    T = 0.17
    L = 20
    nt = round(Int, t / dt)
    thermal = 20

    # thermalization
    v = randomstate(H.Ns, L)
    mcstep!(H, v, T, thermal)

    # simulate
    vs = simulate(H, v, dt, nt)

    # check energy conservation
    E0 = energy(H, v)
    Es = mapslices(vs, dims=[1, 2, 3, 4]) do v
        energy(H, v)
    end
    # ΔE = reduce(max, abs.(Es .- E0))
    ΔE = reshape(Es .- E0, nt) / E0 # relative error as a function of time
    
    # check magnetization conservation
    m0 = magnetization(v)
    ms = reshape(mapslices(magnetization, vs; dims=[1, 2, 3, 4]), (3, :))
    # Δm = reduce(max, mapslices(norm, ms .- m0, dims=[1]))
    Δm = mapslices(norm, ms .- m0, dims=[1]) ./ m0

    ΔE, Δm
end

# plot the differences

# dts = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
dts = [0.1, 0.05, 0.01]
diffs = [compute_diff(H, dt, 10) for dt in dts]
ΔEs = [ΔE for (ΔE, Δm) in diffs]
Δms = [Δm for (ΔE, Δm) in diffs]

display(plot())
for n in 1:length(dts)
    ΔE = ΔEs[n]
    dt = dts[n]
    
    display(plot!(range(0, 10, step=dt)[1:end-1], ΔE, label="dt=$dt"))
end

# scatter(ΔEs, label="ΔE",
#         xticks=(1:length(dts), dts),
#         xlabel="time step")
        
# scatter!(Δms, label="Δm")
