"""
Script showcasing the workflow for the whole GHE simulation
"""

using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = 1.0:3600:3600*24*365  # Time vector for 1 year with 1 hour time steps [s]

# Compute g-function for FLS
g_fls = fls(t, H, rb, D, ks, Cs)

# Compute the thermal resistance for a single U-tube
Rb = resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, Cf, ρf, μf)

# Create a heat load profile
Q = heat_load_profile(t / 3600)  # Heat load profile with daily variation [W]
dT = Q / (V * Cf)  # Temperature difference corresponding to the heat load [K]

# Compute the average fluid temperature at the borehole outlet
Tf = fluid_temperature(t, Q / H, g_fls, T0, ks, Rb)
Tf2 = fluid_temperature(t, Q, g_fls, H, T0, ks, Rb)
@assert all(Tf .≈ Tf2) "The two methods for computing fluid temperature should give the same result"

# Compute the outlet temperature
Tout = outlet_temperature(Tf, Q, V, Cf)

# Compute the inlet temperature
Tin = inlet_temperature(Tf, Q, V, Cf)

# Plot results
fig = Figure()
ax1 = Axis(fig[1, 1:2], xlabel="Time (s)", ylabel="Temperature (°C)", title="Fluid Temperatures")
lines!(ax1, t, Tin, label="Inlet Temperature")
lines!(ax1, t, Tout, label="Outlet Temperature")
axislegend(ax1, position=:cb)
ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="g-function (°Cm/W)", title="g-function")
lines!(ax2, t, g_fls, label="FLS g-function")
axislegend(ax2, position=:rb)
ax3 = Axis(fig[2, 2], xlabel="Time (s)", ylabel="Heat Load (W)", title="Heat Load Profile")
lines!(ax3, t, Q, label="Heat Load")
axislegend(ax3, position=:cb)
display(fig)