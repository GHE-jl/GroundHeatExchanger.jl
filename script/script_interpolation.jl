# Validation of logarithmic node generation and PCHIP interpolation on the ILS model.
# Demonstrates the compression used by ground_response.

using CairoMakie
using GroundHeatExchanger

# Parameters
ti, H, D, _, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = collect(3600.0:3600:ti[end])   # hourly steps up to 100 years [s]

# Compute ILS at reduced logarithmic node set (using grout parameters at borehole wall)
id = set_nodes(length(t), 100)
@time gᵢ = ils(t[id], rb, kg, Cg)

# Compute ILS at full time vector
@time g_full = ils(t, rb, kg, Cg)

# PCHIP interpolation from nodes to full vector
@time g̃ = pchip_interpolation(t[id], gᵢ, t)

# Error metric
rmse = sqrt(sum((g_full .- g̃).^2) / length(t))
println("RMSE = $(round(rmse, sigdigits=4)) °C·m/W  ($(length(id)) nodes / $(length(t)) steps)")

# Plot
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1],
    xlabel="Time (s)", ylabel="ILS g-function (°C·m/W)",
    title="PCHIP compression: $(length(id)) nodes → $(length(t)) steps",
    xscale=log10)
lines!(ax, t, g_full, label="Full ILS ($(length(t)) pts)")
lines!(ax, t, g̃, label="PCHIP ($(length(id)) nodes)", linestyle=:dash)
scatter!(ax, t[id], gᵢ, markersize=4, label="Nodes")
text!(ax, 0.05, 0.5; text="RMSE = $(round(rmse, sigdigits=4))", space=:relative, align=(:left, :center))
axislegend(ax, position=:lt)
display(fig)
