"""
Script testing the models from GHEModels.jl.
"""

using BenchmarkTools
using Plots

includet("../src/GHEModels.jl")
using .GHEModels

# Paremeters
#t = collect(range(3600.0, 3600.0*24*365*100, step=3600))                       # Time (lin)
#ᵢ = set_nodes(length(t), 150)                                                   # Nodes
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

nx, ny, B = 1, 1, 5.
xy = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_mfls_1 = mfls_borefield_I(t, k.s, C.s, r.b, H, D, vD, xy)

nx, ny, B = 2, 2, 5.
xy = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_mfls_2 = mfls_borefield_I(t, k.s, C.s, r.b, H, D, vD, xy)

nx, ny, B = 10, 10, 5.
xy = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_mfls_3 = mfls_borefield_I(t, k.s, C.s, r.b, H, D, vD, xy)

#@time g_fls = fls(t, ks, Cs, rb, H, D)

# Figure
col = fig_color()

plot(t, g_mfls_1, lw=3, lc=col[1], label="MFLS - 1x1")
plot!(t, g_mfls_2, lw=2, lc=col[2], label="MFLS - 2x2")
plot!(t, g_mfls_3, lw=1, lc=col[3], label="MFLS - 10x10")
plot!(xaxis="Time (s)",
    yaxis="g (-)",
    xscale=:log10,
    framestyle=:box,
    grid=false,
    xlabelfontsize=8,
    ylabelfontsize=8,
    xtickfontsize=8,
    ytickfontsize=8,
    legendfontsize=8,
    legend=:topleft)