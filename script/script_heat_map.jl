"""
Script to generate heat maps using the moving finite line source model.
"""

using CairoMakie
includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()

t = [3600.0 * 24, 3600.0 * 24 * 30, 3600.0 * 24 * 365, 
         3600.0 * 24 * 365 * 5, 3600.0 * 24 * 365 * 10, 3600.0 * 24 * 365 * 25]
labels = ["1 Day", "1 Month", "1 Year", "5 Years", "10 Years", "25 Years"]

# Grid setup
n = 200
x = range(-5.0, 15.0, length=n)
y = range(-10.0, 10.0, length=n)

# Compute grid data and store results
results = [_mfls(t_, ks, Cs, Cf, xi, yi, rb, H, D, vD) for t_ in t, xi in x, yi in y]

# Global min/max for consistent colorbar
g_max = maximum(maximum.(results))

# Plot
fig = Figure()
axes = [Axis(fig[i, j], title = labels[(i-1)*3 + j], xlabel = L"$x$ (m)", ylabel = L"$y$ (m)") 
        for i in 1:2, j in 1:3]

for i in 1:2, j in 1:3
    local g = results[(i-1)*3 + j, :, :]
    # Heatmap, Contours and Borehole
    heatmap!(axes[i, j], x, y, g, colormap = :dense, colorrange = (0, g_max))
    contour!(axes[i, j], x, y, g, color = :white, linewidth = 0.5)
    scatter!(axes[i, j], [0.0], [0.0], color = :grey, markersize = 10, label = "Borehole")
end

# Add a shared colorbar
Colorbar(fig[:, 4], limits = (0, g_max), colormap = :dense, label = "g-function (°Cm/W)")

display(fig)