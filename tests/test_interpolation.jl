"""
Script testing the nodes generation to apply to analytical models.
"""

using BenchmarkTools
using Plots

includet("../src/GHEModels.jl")
using .GHEModels

# Paremeters
#t = collect(range(3600.0, 3600.0*24*365*100, step=3600))                       # Time (lin)
t = collect(exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500))) # Time (log)
ᵢ = set_nodes(length(t), 150)                                                   # Nodes
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

function ils_interp(t, tᵢ, kg, Cg, rb)
    gᵢ = ils(tᵢ, kg, Cg, rb)
    g = similar(t)
    @time g = pchip_interpolation(tᵢ, gᵢ, t)
    return g
end

# Run ILS model
#@time g = ils(t, ks, Cs, rb)
gᵢ = ils_interp(t, t[ᵢ], k.s, C.s, r.b)

# Plot
#plot(t, g; label="full", lw=1.5, lc=:black)
#plot!(t[ᵢ], gᵢ; label="nodes", lw=1.5, mc=:blue, seriestype=:scatter)

#=plot!(;
    xaxis="Time (s)",
    yaxis="g (-)",
    xscale=:log10,
    framestyle=:box,
    grid=false,
    xlabelfontsize=8,
    ylabelfontsize=8,
    xtickfontsize=8,
    ytickfontsize=8,
    legendfontsize=8
)=#
