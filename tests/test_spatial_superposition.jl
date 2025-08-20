"""
Script testing the finite line source (FLS) model from GHEModels.jl
"""

using BenchmarkTools
using Plots

includet("../src/GHEModels.jl")
using .GHEModels

# Define paremeters
#t = collect(range(3600.0, 3600.0*24*365*100, step=3600))                       # Time (lin)
#ᵢ = set_nodes(length(t), 150)                                                   # Nodes
t = collect(exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500))) # Time (log)
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
#Q = 10000.0 * ones(length(t)) + (300 * randn(length(t))) # Constant heating power signal

# Spatial superposition with bloc matrix
nx, ny, B = 1, 1, 5.
xy1 = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_fls_1 = bloc_matrix(t, k.s, C.s, r.b, H, D, xy1)

nx, ny, B = 2, 2, 5.
xy2 = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_fls_2 = bloc_matrix(t, k.s, C.s, r.b, H, D, xy2)
#@time g_fls_2 = successive_flux()

nx, ny, B = 5, 5, 5.
xy3 = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'.-B
@time g_fls_3 = bloc_matrix(t, k.s, C.s, r.b, H, D, xy3)

# Figure