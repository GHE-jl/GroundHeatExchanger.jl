"""
Test for the non-stationary convolution from Convolutions.jl. This script does not
correspond to real application of the non-stationary convolution, as it is more for operating
conditions. Here, the test allows comparison with the Matlab implementation of the non-stationary
convolution.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GHEModels.jl")
using .GHEModels

# Define paremeters
t = range(60.0, 3600.0*24*6, step=60) # Time (linear)
H = 150.0                       # Borehole depth
D = 2.0                         # Borehole buried depth
s = 0.05                        # Shank spacing (s/2 is the half-shank spacing)
r = (b = 0.08,                  # Borehole radius
    o = 0.022,                  # Pipe outlet radius
    i = 0.017)                  # Pipe inlet radius
k = (s = 3.0,                   # Ground thermal conductivity
    g = 1.6,                    # Grout thermal conductivity
    p = 0.4,                    # Pipe thermal conductivity
    f = 0.6)                    # Fluid thermal conductivity
C = (s = 2.11e6,                # Ground volumetric specific heat
    g = 2.25e6,                 # Grout volumetric specific heat
    p = 1.9e6,                  # Pipe volumetric specific heat
    f = 4.2e6)                  # Fluid volumetric specific heat
ρ = (s = 1000.0,                # Groud density
    g = 1000.0,                 # Grout density
    p = 1000.0,                 # Pipe density
    f = 1000.0)                 # Fluid density

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
    g_fls[:, i] = fls(t, V1[j], C.s, r.b, V2[j], D)
end

# Non-stationary convolution
f = impulse_func_ns(q, state_vec)
@time dT1 = convolution_ns(f, g_fls, state_vec)
@time dT2 = convolution_ns(f, g_fls, ind, state)
@time dT3 = convolution_ns(q, g_fls, state_vec)
@time dT4 = convolution_ns(q, g_fls, ind, state)

T_verif = dT1[[1000, 2800, 3000, 6000, 8000, 8600]]

fig = Figure()

ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = "g-function (-)", xscale = log10)
for i in eachcol(g_fls)
    lines!(ax, t, i, linewidth = 1.5)
end

ax = Axis(fig[2, 1], xlabel = L"$t$ (s)", ylabel = "dT (°C)")
lines!(ax, t, dT1, color = :black, linestyle = :solid, linewidth = 1.5)
lines!(ax, t, dT2, color = :grey, linestyle = :dash, linewidth = 1.5)
lines!(ax, t, dT3, color = :blue, linestyle = :dashdot, linewidth = 1.5)
lines!(ax, t, dT4, color = :green, linestyle = :dot, linewidth = 1.5)
display(fig);