# End-to-end GHE temperature simulation: g-function → Rb → Tf → Tin/Tout.
# Demonstrates the full three-package workflow for a single U-tube borehole over 1 year.

using CairoMakie
using GroundHeatExchanger

# Parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t  = collect(3600.0:3600:3600*24*365)   # 1 year, hourly [s]
th = t ./ 3600                        # hours (for plotting)

# Step 1: Ground response (GroundResponse.jl)
model = FLSModel(H, D, ks, Cs)
g = ground_response(t, rb, [0.0 0.0], model)   # FLS g-function, PCHIP-compressed to 150 nodes
g_fls = fls(t, rb, H, D, ks, Cs)                        # FLS g-function, full time vector

# Step 2: Borehole thermal resistance (BoreholeResistance.jl)
Rbₑ = resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, Cf/ρf, ρf, μf)

# Step 3: Heat load profile
Q = heat_load_profile(th)   # [W]
q = Q ./ H                  # heat load per unit length [W/m]

# Step 4: Mean fluid temperature (GroundHeatExchanger.jl)
Tf = fluid_temperature(t, q, g, T0, ks, Rbₑ)

# Step 5: Outlet and inlet temperatures
Tout = outlet_temperature(Tf, Q, V, Cf)
Tin  = inlet_temperature(Tf, Q, V, Cf)

# Workflow figure
fig = Figure(size=(1100, 700))

ax_g = Axis(fig[1, 1],
    xlabel = "Time (h)", ylabel = "g-function (°C·m/W)",
    title  = "Step 1 — FLS g-function")
lines!(ax_g, th, g, linewidth = 2, label = "PCHIP-compressed")
lines!(ax_g, th, g_fls, linewidth = 2, linestyle = :dash, color = :black, label = "Full FLS")
axislegend(ax_g, position = :rb)

ax_q = Axis(fig[1, 2],
    xlabel = "Time (h)", ylabel = "Heat load (W)",
    title  = "Step 2 — Heat load profile")
lines!(ax_q, th, Q)

ax_T = Axis(fig[2, 1:2],
    xlabel = "Time (h)", ylabel = "Temperature (°C)",
    title  = "Steps 3–5 — Mean, inlet and outlet fluid temperatures")
    lines!(ax_T, th, Tin,  label = "Tin (inlet)",   linewidth = 1.5, color = :red)
lines!(ax_T, th, Tf,   label = "Tf  (mean)",    linewidth = 2,   color = :black)
lines!(ax_T, th, Tout, label = "Tout (outlet)", linewidth = 1.5, color = :blue)
axislegend(ax_T, position = :lb)

display(fig)
