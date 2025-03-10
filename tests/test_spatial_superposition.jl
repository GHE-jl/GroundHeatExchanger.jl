"""
Script testing the finite line source (FLS) model from GHEModels.jl.
"""

using BenchmarkTools
using Plots

includet("../src/GroundTSimulations.jl")
using .GroundTSimulations

Revise.revise()

# Paremeters
params = GHE_param(
    collect(exp10.(range(log10(60.), log10(3600*24*365*100), length=500))), # Time range (log)
    150.,                       # Borehole depth
    2.,                         # Borehole buried depth
    0.05,                       # Shank spacing (s/2 is the half-shank spacing)
    0.08,                       # Borehole radius
    0.017,                      # Pipe inlet radius
    0.022,                      # Pipe outlet radius
    3.,                         # Ground thermal conductivity
    1.6,                        # Grout thermal conductivity
    0.4,                        # Pipe thermal conductivity
    0.6,                        # Fluid thermal conductivity
    2.11e6,                     # Ground volumetric specific heat
    2.25e6,                     # Grout volumetric specific heat
    1.9e6,                      # Pipe volumetric specific heat
    4.2e6,                      # Fluid volumetric specific heat
    1000.,                      # Groud density
    1000.,                      # Grout density
    1000.,                      # Pipe density
    1000.,                      # Fluid density
    30. / 60000,                # Circulating flow rate
    1e-12,                      # Groundwater flow
)

# Positions of borehole(s) in a borefield
# xy = [0. 0.]

nx, ny, B = 2, 2, 8.
xy = B * hcat([[i, j] for i in 1:nx for j in 1:ny]...)'
# xy = [0. 0.; xy]

# Create matrix of g-functions
@time gg = gfunc_matrix(params, xy, "fls")

# Spatial superposition with bloc matrix
@time g_bloc = bloc_matrix(params, gg, xy)

plot(params.t, g_bloc, lw=3, lc="black", label="Bloc matrix")
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