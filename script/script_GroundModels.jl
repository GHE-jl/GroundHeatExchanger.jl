"""
Script showcasing the analytical models from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GHEModels.jl")
using .GHEModels

# Define paremeters
#t = range(3600.0, 3600.0*24*365*100, step=3600)                                     # Time (lin)
#s = set_nodes(length(t), 150)                                                       # Nodes
t = exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500))          # Time (log)
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
V = 30.0 / 60000                # Circulating flow rate
vD = 1e-7                       # Groundwater flow
#Q = 10000.0 * ones(length(t)) + (300 * randn(length(t))) # Constant heating power signal

# Run models
@time g_ils = ils(t, k.s, C.s, r.b)
@time g_ics = ics(t, k.s, C.s, r.b)
@time g_fls = fls(t, k.s, C.s, r.b, H, D)
@time g_mfls = mfls(t, k.s, C.s, r.b, H, D, vD)
@time g_ilsβ = βils(t, k.s, C.s, r.b, H, V, 0.01)

t̃ = t / (3600 * 24 * 365)
f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = "Time (yr)", ylabel = "g-function (-)", xscale = log10)
lines!(ax, t̃, g_ils, linewidth = 4, label = "ILS")
lines!(ax, t̃, g_ics, linewidth = 2, label = "ICS")
lines!(ax, t̃, g_fls, linewidth = 2, label = "FLS")
lines!(ax, t̃, g_mfls, linewidth = 2, label = "MFLS")
lines!(ax, t̃, g_ilsβ, linewidth = 2, label = "βILS")
axislegend(ax, position = :lt)
display(f);