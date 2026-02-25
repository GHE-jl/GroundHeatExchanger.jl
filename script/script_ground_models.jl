"""
Script showcasing the analytical models from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
#t = range(3600.0, 3600.0*24*365*100, step=3600)                                     # Time (lin)
#s = set_nodes(length(t), 150)                                                       # Nodes
t = exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500))          # Time (log)
H = 150.0                       # Borehole depth
D = 2.0                         # Borehole buried depth
r = 0.08                        # Radius at which to compute the models
rb = 0.08                       # Borehole radius
ks = 3.0                        # Ground thermal conductivity
Cs = 2.11e6                     # Ground volumetric specific heat
Cf = 4.2e6                      # Fluid volumetric specific heat
vD = 5e-7                       # Groundwater flow

# Run models
@time g_ils = ils(t, ks, Cs, r)
@time g_ics = ics(t, ks, Cs, r)
@time g_fls = fls(t, ks, Cs, r, H, D)
@time g_mils = mils(t, ks, Cs, Cf, r, vD)
@time g_mfls = mfls(t, ks, Cs, Cf, [0, 0], rb, H, D, vD)

t̃ = t / (3600 * 24 * 365)
f = Figure()
ax = Axis(f[1, 1], xlabel = "Time (yr)", ylabel = "g-function (°Cm/W)", xscale = log10)
lines!(ax, t̃, g_ils, linewidth = 5, linestyle = :solid, label = "ILS")
lines!(ax, t̃, g_ics, linewidth = 4, linestyle = :dash, label = "ICS")
lines!(ax, t̃, g_fls, linewidth = 3, linestyle = :dot, label = "FLS")
lines!(ax, t̃, g_mils, linewidth = 2, linestyle = :dashdot, label = "MILS")
lines!(ax, t̃, g_mfls, linewidth = 1, linestyle = :dashdotdot, label = "MFLS")
axislegend(ax, position = :lt)
display(f);