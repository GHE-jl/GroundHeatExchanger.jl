"""
Script showing the thermal resistance models from ThermalResistances.jl using either coaxial, 
single or double U-loop ground heat exchanger.
"""

using BenchmarkTools

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
V̇ = V / (π * ri^2)
cf = Cf / ρf

# Fluid thermal resistance
# Calculate Reynolds and Prandtl numbers
Re = Reynold(V̇, ri, ρf, μf)
Pr = Prandtl(kf, cf, μf)

# Calculate Nusselt number
Nu1 = Nusselt(Re, Pr, ri, ϵ)
Nu2 = Nusselt(V̇, ri, kf, cf, ρf, μf, ϵ)

# Calculate fluid resistance
Rf_1 = resistance_fluid(Nu1, ri, kf)
Rf_2 = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ)

# Pipe thermal resistance (single U-loop)
Rp = resistance_pipe(ro, ri, kp)

# Borehole thermal resistance
Rb_0_1 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_1, 0)
Rb_0_2 = resistance_borehole_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ, 0)
Rb_1_1 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_1, 1)
Rb_1_2 = resistance_borehole_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ, 1)

# Total internal thermal resistance
Ra_0_1 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_1, 0)
Ra_0_2 = resistance_total_internal_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ, 0)
Ra_1_1 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_1, 1)
Ra_1_2 = resistance_total_internal_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ, 1)

Rbe_1 = resistance_borehole_effective(V, H, cf, ρf, Rb_1_1, Ra_1_1)
Rbe_2 = resistance_borehole_effective(V, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_1)
Rbe_3 = resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)