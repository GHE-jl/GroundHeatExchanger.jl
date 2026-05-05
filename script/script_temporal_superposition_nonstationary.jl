"""
Test for the non-stationary convolution from temporal_superposition.jl.
This script does not correspond to real application of the non-stationary convolution, as it is more
for operating conditions. Here, the test allows comparison with the Matlab implementation of the
non-stationary convolution.
"""

using BenchmarkTools
using CairoMakie

includet("../src/temporal_superposition.jl")

# Define paremeters
includet("../src/Utils.jl")
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = range(60.0, 3600.0*24*6, step=60) # Time (linear)

# Define operating conditions
Q = 10000.0 * ones(60*24*6)     # Constant heating power signal for every minutes in 6 days
q = Q / H

ks = [3.0, 2.0, 3.0, 1.0]       # Varying thermal conductivity for the 4 states
H = [125.0, 150.0, 125.0]       # Varying borehole depth for the 3 states
V1 = repeat(ks, inner = Integer(60 * 24 * 1.5))
V2 = repeat(H, inner = Integer(60 * 24 * 2))

# Convert operating conditions in indices, states and state vector
ind, state, ind_unique = state_indices(V1, V2)
state_vec = state_vector(ind, state, length(t))

# Create transfer functions
g_fls = Matrix{Float64}(undef, length(t), length(unique(state)))
for (i, j) in enumerate(ind_unique)
    g_fls[:, i] = fls(t, V1[j], Cs, rb, V2[j], D)
end

# Non-stationary convolution
dT1 = similar(t)
f = impulse_func_ns(q, state_vec)
@time convolution_ns!(dT1, f, g_fls)
@time dT2 = convolution_ns(f, g_fls)
@time dT3 = convolution_ns(q, g_fls, state_vec)
@time dT4 = convolution_ns(q, g_fls, ind, state)

T_verif = dT1[[1000, 2800, 3000, 6000, 8000, 8600]]

fig = Figure()

ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = "g-function (-)", xscale = log10)
for (i, gi) in enumerate(eachcol(g_fls))
    lines!(ax, t, gi, linewidth = 3, label="State $i")
end
axislegend(ax, position = :lt)

ax = Axis(fig[2, 1], xlabel = L"$t$ (s)", ylabel = "dT (°C)")
lines!(ax, t, dT1, color=:black, linestyle=:solid, linewidth=3, label="dT1")
lines!(ax, t, dT2, color=:grey, linestyle=:dash, linewidth=3, label="dT2")
lines!(ax, t, dT3, color=:blue, linestyle=:dashdot, linewidth=3, label="dT3")
lines!(ax, t, dT4, color=:green, linestyle=:dot, linewidth=3, label="dT4")
axislegend(ax, position=:lt, nbanks=2)
display(fig);