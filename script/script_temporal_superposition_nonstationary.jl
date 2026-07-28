# Non-stationary temporal superposition — a minimal, side-by-side companion to
# `script_temporal_superposition_stationary.jl` (same 6-day, 1-minute time base).
#
# What makes this "non-stationary"?
#   The flow rate V changes during operation. The effective borehole resistance Rbₑ depends on V,
#   so the system's transfer function is no longer a single fixed curve — it switches between
#   states. We capture this with one combined transfer function per flow state, for an impulse of
#   1 W/m:
#       h_s(t) = Rbₑ(V_s) + g(t)
#   and the mean fluid temperature follows from a single non-stationary convolution:
#       Tf(t) = T0 + (f * h)(t)
#   where f segregates the incremental heat load by the state active at each time step.
#
# To keep the effect unambiguous, the HEAT LOAD IS HELD CONSTANT. The only thing that varies is
# the flow rate, so every difference between the non-stationary and stationary curves below is due
# to non-stationarity alone.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per project, to install dependencies

using CairoMakie
using GroundHeatExchanger

# Fixed GHE parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, vD, V_nom = GHE()
t  = 60.0:60:3600*24*6   # 6 days, 1-minute steps [s] (same as the stationary script)
th = t ./ 3600            # hours (for plotting)

# Constant heat load — isolates the flow-rate (non-stationary) effect
Q = fill(12e3, length(t))   # 12 kW, constant [W]
q = Q ./ H                  # heat load per unit length [W/m]

# Ground response: one g-function for all states (soil and geometry are fixed)
g = ground_response(collect(t), rb, [0.0 0.0], FLSModel(H, D, ks, Cs))   # FLS g-function [°C·m/W]

# Three flow-rate states, one per 2-day block (mirrors the stationary script's 3 blocks).
# Block 1 is the nominal flow, so the non-stationary and stationary curves start identical and
# only split once the flow rate changes.
V_states = [V_nom, V_nom / 2, V_nom * 3/2]   # nominal → low → high [m³/s]
n_states = length(V_states)

nb = length(t) ÷ 3
state_vec = vcat(fill(1, nb), fill(2, nb), fill(3, length(t) - 2nb))
V_t = V_states[state_vec]   # flow rate at each time step [m³/s]

println("Flow-rate states [L/min]: ", round.(V_states .* 6e4, digits=1))

# Effective borehole resistance for each flow state (BoreholeResistance.jl).
# Lower flow → higher Rbₑ (more thermal short-circuiting), and vice versa.
Rbₑ = [resistance_ULoop_effective(Vi, H, s, rb, ro, ri, ks, kg, kp, kf, Cf/ρf, ρf, μf)
       for Vi in V_states]
println("Rbₑ states [m·K/W]:       ", round.(Rbₑ, digits=4))

# Combined transfer function per state: h_s(t) = Rbₑ(V_s) + g(t)
h = hcat([Rbₑ[i] .+ g for i in 1:n_states]...)

# Non-stationary mean fluid temperature: Tf(t) = T0 + (f * h)(t)
Tf = T0 .+ convolution_ns(q, h, state_vec)

# Stationary reference: same load, but the flow rate stays at nominal for the whole simulation
h_nom  = Rbₑ[1] .+ g
Tf_nom = T0 .+ convolution(q, h_nom)

# Sanity check: the h-formulation must reproduce fluid_temperature() in the stationary case
@assert Tf_nom ≈ fluid_temperature(t, q, g, T0, Rbₑ[1])  "h-formulation must match fluid_temperature"

# Figure
fig = Figure(size = (1000, 700))

ax1 = Axis(fig[1, 1],
    xlabel = "Time (h)", ylabel = "Flow rate (L/min)",
    title  = "Driver of non-stationarity — flow rate V(t)")
lines!(ax1, th, V_t .* 6e4, linewidth = 2)

ax2 = Axis(fig[1, 2],
    xlabel = "Time (h)", ylabel = "h-function (°C·m/W)",
    title  = "Transfer function per flow state")
for i in 1:n_states
    lines!(ax2, th, h[:, i], linewidth = 2,
           label = "V = $(round(V_states[i] * 6e4, digits=1)) L/min")
end
axislegend(ax2, position = :rb)

ax3 = Axis(fig[2, 1],
    xlabel = "Time (h)", ylabel = "Temperature (°C)",
    title  = "Mean fluid temperature")
lines!(ax3, th, Tf,     color = :black,        linewidth = 2, label = "non-stationary")
lines!(ax3, th, Tf_nom, color = (:green, 0.7), linewidth = 2, label = "stationary (V nominal)")
axislegend(ax3, position = :rb)

ax4 = Axis(fig[2, 2],
    xlabel = "Time (h)", ylabel = "ΔTf (°C)",
    title  = "Non-stationarity effect: Tf − Tf(stationary)")
lines!(ax4, th, Tf .- Tf_nom, color = :purple, linewidth = 2)

display(fig)
