"""
Script showcasing the non-stationary convolution from Convolutions.jl.
"""

using BenchmarkTools
using CairoMakie

includet("../src/GHEModels.jl")
using .GHEModels

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
Q = 10000.0 * ones(60*24*6)
q = 10000.0 * ones(60*24*6) / H # Constant heating power signal for every minutes in 6 days

#ks = [1., 2., 3., 2., 3., 1.]
#V1 = repeat(ks, inner = Integer(60 * 24 * 1))
#@time ind, s, ind_unique = state_transitions(V1)

ks = [3., 2., 3., 1]
H = [125.0, 150.0, 125.0]
V1 = repeat(ks, inner = Integer(60 * 24 * 1.5))
V2 = repeat(H, inner = Integer(60 * 24 * 2))
@time ind, s, ind_unique = state_transitions(V1, V2)

# Create transfer functions
#g_fls = hcat([fls(t, k, C.s, r.b, H, D) for k in ks[1:3]]...)
g_fls = Matrix{Float64}(undef, length(t), length(unique(s)))
for (i, j) in enumerate(ind_unique)
    g_fls[:, i] = fls(t, V1[j], C.s, r.b, V2[j], D)
end

# Non-stationary convolution
@time dT3 = convolution_ns(q, g_fls, ind, s)

f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = "g-function (-)", xscale = log10)
for i in 1:length(unique(s))
    lines!(ax, t, g_fls[:, i], color = col[i], linewidth = 1.5)
end
ax = Axis(f[2, 1], xlabel = L"$t$ (s)", ylabel = "dT (°C)")
lines!(ax, t, dT3, color = :black, linewidth = 1.5)
display(f)

#x1 = [8,8,8,8,8,4,4,4,4,4]
#x2 = [3,3,6,6,7,7,9,9,1,1]