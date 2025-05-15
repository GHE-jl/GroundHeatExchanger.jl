"""
Script showing the analytical models from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using Plots

includet("../src/GHEModels.jl")
using .GHEModels

# Paremeters
t2 = range(3600.0, 3600.0*24*365*100, step=3600)                       # Time (lin)
#s = set_nodes(length(t), 150)                                                  # Nodes
t = collect(exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500))) # Time (log)
H = 150.0                       # Borehole depth
D = 2.0                         # Borehole buried depth
s = 0.05                        # Shank spacing (s/2 is the half-shank spacing)
r = (b = 0.08,                  # Borehole radius
    i = 0.017,                  # Pipe inlet radius
    o = 0.022)                  # Pipe outlet radius
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
vD = 1e-9                       # Groundwater flow
#Q = 10000.0 * ones(length(t)) + (300 * randn(length(t))) # Constant heating power signal

# Run models
@time g_ils = ils(t, k.s, C.s, r.b)
@time g_fls = fls(t, k.s, C.s, r.b, H, D)
@time g_mfls = mfls_single_borehole(t, k.s, C.s, r.b, H, D, vD)
@time g_scwm = scwm(t, k.s, C.s, r.b, H, V, 0.01)

ghe_model(t, k.s, C.s, r.b, H, D, [0.0 0.0])

col = palette(:tab10)
plot(t, g_ils, lw = 4, lc = col[1], label = "IlS")
plot!(t, g_fls, lw = 3, lc = col[2], label = "FlS")
plot!(t, g_mfls, lw = 2, lc = col[3], label = "MFlS")
plot!(t, g_scwm, lw = 1, lc = col[4], label = "SCWM")
plot!(xaxis = "Time (s)",
    yaxis = "g (-)",
    xscale = :log10,
    framestyle = :box,
    grid = false,
    xlabelfontsize = 8,
    ylabelfontsize = 8,
    xtickfontsize = 8,
    ytickfontsize = 8,
    legendfontsize = 8,
    legend = :topleft
)