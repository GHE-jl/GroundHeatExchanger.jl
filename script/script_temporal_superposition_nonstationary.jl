# Non-stationary GHE simulation with time-varying flow rate.
#
# Non-stationarity enters through Rbₑ(V(t)).  Instead of separating the borehole and ground
# contributions, we define a per-state combined transfer function for an impulse of 1 W/m:
#
#   h(t) = Rbₑ(V) + g(t) / (2π·ks)
#
# so that the fluid temperature reduces to a single NS convolution:
#
#   Tf(t) = T0 + (f * h)(t)
#
# where f = impulse_func_ns(q, state_vec) is the incremental heat load segregated by state.
# This formulation is equivalent to fluid_temperature() for the stationary case (verified below).

using CairoMakie
using GroundHeatExchanger

# Fixed GHE parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V_nom = GHE()
t  = collect(1.0:3600:3600*24*365)   # 1 year, hourly [s]
th = t ./ 3600                        # hours (for plotting)

# Heat load
Q = heat_load_profile(th)   # [W]
q = Q ./ H                  # [W/m]

# Ground response (same g for all flow states: soil and geometry are fixed)
model = FLSModel(H, D, ks, Cs)
g = ground_response(t, rb, [0.0 0.0], model)   # FLS g-function, PCHIP-compressed [°C·m/W]

# Three flow rate states (low / nominal / high)
V_states = [V_nom / 2, V_nom, V_nom * 3/2]   # [m³/s]
n_states = length(V_states)

n1 = length(t) ÷ 3
n2 = 2 * n1
state_vec = vcat(fill(1, n1), fill(2, n2 - n1), fill(3, length(t) - n2))

println("Flow rate states [L/min]: ", round.(V_states .* 6e4, digits=1))

# Borehole resistance per state (BoreholeResistance.jl)
Rbₑ = [resistance_borehole_effective(Vi, H, s, rb, ro, ri, ks, kg, kp, kf, Cf/ρf, ρf, μf)
       for Vi in V_states]
println("Rbₑ states [m·K/W]:       ", round.(Rbₑ, digits=4))

# Combined transfer function h_s(t) = Rbₑ_s + g(t)/(2π·ks)
h = Matrix{Float64}(undef, length(t), n_states)
for i in 1:n_states
    h[:, i] = Rbₑ[i] .+ g ./ (2π * ks)
end

# Non-stationary convolution: Tf(t) = T0 + (f * h)(t)
f_ns = impulse_func_ns(q, state_vec)
Tf = T0 .+ convolution_ns(f_ns, h)

# Outlet and inlet temperatures (state-dependent V)
V_t  = [V_states[state_vec[i]] for i in eachindex(t)]
Tout = Tf .- Q ./ (2 .* V_t .* Cf)
Tin  = Tf .+ Q ./ (2 .* V_t .* Cf)

# Stationary reference at nominal flow rate
h_nom  = Rbₑ[2] .+ g ./ (2π * ks)
Tf_nom = T0 .+ convolution(q, h_nom)

# Verify: h-formulation must match fluid_temperature for the stationary case
Tf_check = fluid_temperature(t, q, g, T0, ks, Rbₑ[2])
@assert Tf_nom ≈ Tf_check  "h-formulation must match fluid_temperature"

Tout_nom = Tf_nom .- Q ./ (2 * V_nom * Cf)
Tin_nom  = Tf_nom .+ Q ./ (2 * V_nom * Cf)

# Figure
fig = Figure(size=(1100, 750))

ax1 = Axis(fig[1, 1],
    xlabel = "Time (h)", ylabel = "Flow rate (L/min)",
    title  = "Operating states — time-varying flow rate")
lines!(ax1, th, V_t .* 6e4)

ax2 = Axis(fig[1, 2],
    xlabel = "Time (h)", ylabel = "h-function (°C·m/W)",
    title  = "Ground heat exchanger response — FLS h-function")
for i in 1:n_states
    lines!(ax2, th, h[:, i], linewidth = 2,
           label = "V = $(round(V_states[i] * 6e4, digits=1)) L/min")
end
axislegend(ax2, position = :cb, nbanks = 2)

ax3 = Axis(fig[2, 1],
    xlabel = "Time (h)", ylabel = "Temperature (°C)",
    title  = "Non-stationary temperatures")
lines!(ax3, th, Tin, color = :red, linewidth = 2, label = "Tin non-stationary")
lines!(ax3, th, Tout, color = :blue, linewidth = 2, label = "Tout non-stationary")
axislegend(ax3, position = :cb, nbanks = 2)

ax4 = Axis(fig[2, 2],
    xlabel = "Time (h)", ylabel = "Temperature difference (°C)",
    title  = "Non-stationary vs. stationary fluid temperatures")
lines!(ax4, th, Tf, color = :black, linewidth = 2, label = "Tf non-stationary")
lines!(ax4, th, Tf_nom, color = (:green, 0.5), linewidth = 2, label = "Tf stationary (Vnom)")

axislegend(ax4, position = :cb, nbanks = 2)

display(fig)