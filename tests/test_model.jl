"""
Script testing the models from GHEModels.jl.
"""

using BenchmarkTools
using Plots

includet("../src/GHEModels.jl")
using .GHEModels

# Paremeters
#t = collect(range(3600., 3600*24*365*100, step=3600))                          # Time range (lin)
#s = set_nodes(length(t), 150)                                                  # Nodes
t = collect(exp10.(range(log10(60.), log10(3600*24*365*100), length=500)))      # Time range (log)
H = 150.                        # Borehole depth
D = 2.                          # Borehole buried depth
s = 0.05                        # Shank spacing (s/2 is the half-shank spacing)
rb = 0.08                       # Borehole radius
ri = 0.017                      # Pipe inlet radius
ro = 0.022                      # Pipe outlet radius
ks = 3.                         # Ground thermal conductivity
kg = 1.6                        # Grout thermal conductivity
kp = 0.4                        # Pipe thermal conductivity
kf = 0.6                        # Fluid thermal conductivity
Cs = 2.11e6                     # Ground volumetric specific heat
Cg = 2.25e6                     # Grout volumetric specific heat
Cp = 1.9e6                      # Pipe volumetric specific heat
Cf = 4.2e6                      # Fluid volumetric specific heat
ρs = 1000.                      # Groud density
ρg = 1000.                      # Grout density
ρp = 1000.                      # Pipe density
ρf = 1000.                      # Fluid density
V = 30. / 60000                 # Circulating flow rate
vD = 1e-12                      # Groundwater flow
Q = 10000*ones(length(t))       # Constant heating power signal

xy = [0., 0.]

# Run models
@time g_fls = fls(t, ks, Cs, rb, H, D)
@time g_mfls = mfls(t, ks, Cs, rb, H, D, vD, xy)

col = fig_color()
plot(t, g_fls, lw=3, lc=col[1], ls=:dash, label="FlS")
plot!(t, g_mfls, lw=3, lc=col[2], label="MFlS")
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