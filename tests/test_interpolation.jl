"""
Script testing the nodes generation to apply to analytical models.
"""

includet("../src/GHEModels.jl")

using Plots
using BenchmarkTools
using .GHEModels

function ils_interp(t, tᵢ, kg, Cg, rb)
    gᵢ = ils(tᵢ, kg, Cg, rb)
    g = similar(t)
    @time g = pchip_interpolation(tᵢ, gᵢ, t)
    return g
end

# Set parameters
t = collect(range(60., 3600*24*365*100, step=60))
ᵢ = set_nodes(length(t), 100)
ks = 2.5
Cs = 2.2e6
H = 250.0
rb = 0.08

# Run ILS model
#@time g = ils(t, ks, Cs, rb)
gᵢ = ils_interp(t, t[ᵢ], ks, Cs, rb)

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
