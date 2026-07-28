# End-to-end GHE temperature simulation: g-function → Rb → Tf → Tin/Tout.
# Demonstrates the full three-package workflow for a single U-tube borehole over 1 year.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per project, to install dependencies

using CairoMakie
using GroundHeatExchanger

# Parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, vD, V = GHE()
t  = collect(3600.0:3600:3600*24*365)   # 1 year, hourly [s]
th = t ./ 3600                        # hours (for plotting)

# Step 1: Ground response (GroundResponse.jl)
g = ground_response(t, rb, [0.0 0.0], FLSModel(H, D, ks, Cs))   # FLS g-function, [°C·m/W]
g_fls = fls(t, rb, H, D, ks, Cs)                        # FLS g-function, full time vector

# Step 2: Borehole thermal resistance (BoreholeResistance.jl)
Rbₑ = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, Cf/ρf, ρf, μf)

# Step 3: Heat load profile
Q = ground_load_profile(th)   # [W]
q = Q ./ H                  # heat load per unit length [W/m]

# Step 4: Mean fluid temperature (GroundHeatExchanger.jl)
Tf = fluid_temperature(t, q, g, T0, Rbₑ)

# Step 5: Outlet and inlet temperatures
Tout = outlet_temperature(Tf, Q, V, Cf)
Tin  = inlet_temperature(Tf, Q, V, Cf)

# Workflow figure
fig = Figure()

ax_g = Axis(fig[1, 1],
    xlabel = "Time (h)", ylabel = "g-function (°C·m/W)",
    title  = "Step 1 — FLS g-function")
lines!(ax_g, th, g, linewidth = 2, label = "PCHIP-compressed")
lines!(ax_g, th, g_fls, linewidth = 2, linestyle = :dash, color = :black, label = "Full FLS")
axislegend(ax_g, position = :rb)

ax_q = Axis(fig[1, 2],
    xlabel = "Time (h)", ylabel = "Ground heat load (W)",
    title  = "Step 2 — Ground heat load profile")
lines!(ax_q, th, Q)

ax_T = Axis(fig[2, 1:2],
    xlabel = "Time (h)", ylabel = "Temperature (°C)",
    title  = "Steps 3–5 — Mean, inlet and outlet fluid temperatures")
    lines!(ax_T, th, Tin,  label = "Tin (inlet)",   linewidth = 1.5, color = :red)
lines!(ax_T, th, Tf,   label = "Tf  (mean)",    linewidth = 2,   color = :black)
lines!(ax_T, th, Tout, label = "Tout (outlet)", linewidth = 1.5, color = :blue)
axislegend(ax_T, position = :lb)

display(fig)
