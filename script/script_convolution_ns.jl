"""
Script showcasing the non-stationary convolution from Convolutions.jl.
"""

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

# Define operating conditions
q = 10000.0 * ones(60*24*6) / H # Constant heating power signal for every minutes in 6 days
V = repeat([15.0, 30.0, 45.0] / 60000, inner = 60 * 24 * 2)                # Circulating flow rate
n_state = 3
ind, s = state_transitions(V)

# Create transfer functions
ks = [1., 2., 3.]
g_fls = hcat([fls(t, k, C.s, r.b, H, D) for k in ks]...)

# Non-stationary convolution
dT = convolution_ns(Q, g_fls, ind, s)

f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = "g-function (-)", xscale = log10)
for i in 1:3
    lines!(ax, t, g_fls[:, i], color = col[i], linewidth = 1.5)
end
ax = Axis(f[2, 1], xlabel = L"$t$ (s)", ylabel = "dT (°C)", xscale = log10)
lines!(ax, t, dT, color = :black, linewidth = 1.5)
display(f)