"""
Script showcasing the stationary convolution from Convolutions.jl.
"""

using Random
using BenchmarkTools
using CairoMakie

includet("../src/Convolutions.jl")
includet("../src/GroundModels.jl")

includet("../src/FigOptions.jl")
update_fig_theme()
col = fig_color()

# Define paremeters
t = range(60.0, 3600.0*24*6, step=60) # Time (lin)
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

# Define the transfer function (here, a g-function for impulse in W/m from a finite line source)
g = fls(t, k.s, C.s, r.b, H, D)

# Define the heat load profile with added noise
q = [12000.0 * ones(60*24*2); 8000.0 * ones(60*24*2); 10000.0 * ones(60*24*2)] / H
rnd = rand(length(t))
q_rnd = q .+ rnd

# Convolution to have a temperature vector
ΔT = convolution(diff([0; q_rnd]), g)

# Plot the results
f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = L"$ΔT$ (°C)")
lines!(ax, t, ΔT)
display(f)