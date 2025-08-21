"""
Script testing the models from GHEModels.jl.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GHEModels.jl")
using .GHEModels
includet("../src/FigOptions.jl")
update_fig_theme()
col = fig_color()

# Define paremeters
#t = range(3600., 3600*24*365*100, step=3600)                          # Time range (lin)
#s = set_nodes(length(t), 150)                                                  # Nodes
t = exp10.(range(log10(1.), log10(3600*24*365*100), length=500))      # Time range (log)
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
vD = 1e-9                       # Groundwater flow

#Q = 10000*ones(length(t))       # Constant heating power signal
xy = [0., 0.]

# Run models
#@time g_fls = fls(t, k.s, C.s, r.b, H, D)
#@time g_ics = ics(t, k.s, C.s, r.b)
#@time g_mfls = mfls(t, k.s, C.s, r.b, H, D, vD, xy)

#=f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = "g (-)", xscale = log10)
lines!(ax, t, g_ics, color = col[1], label = "ICS")
lines!(ax, t, g_fls, color = col[2], label = "FLS")
display(f);=#

# Validate the Figure 1 from Nguyen et Pasquier (2021) - Successive flux
t2 = range(3600, 3600*24*365*40, length = 100)
r2 = range(0.08, 100, length = 100)
g_fls_mat = Matrix{Float64}(undef, 100, 100)
for (j1, j2) in enumerate(r2)
    g_fls_mat[:, j1] = fls(t2, 2.5, 2e6, j2, 150., 2.)
end

f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = "r (m)", zlabel = "g (-)")
surface!(ax, t2, r2, g_fls_mat)
display(f);