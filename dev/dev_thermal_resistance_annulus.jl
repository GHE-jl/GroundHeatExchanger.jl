"""
Script that computes the thermal resistance ground heat exchanger in an annulus region.
"""

includet("../src/borehole_thermal_resistance/resistance_fluid.jl")
includet("../src/borehole_thermal_resistance/resistance_pipe.jl")
includet("../src/borehole_thermal_resistance/resistance_borehole.jl")

# Define the parameters of the borehole
includet("../src/utils.jl")
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
V̇ = V / (π * ri^2)  # Volumetric flow rate of the fluid

# Compute Nusselt numbers depending on the flow regime
Re = Reynold(V̇, rb - ro, ρf, μf)
Pr = Prandtl(kf, Cf/ρf, μf)

Nu1 = Nusselt(Re, Pr, rb - ro, ϵ)
Nu2 = Nusselt_annulus(Re, Pr, rb, ro, ϵ)

println("Nusselt number for single pipe and equivalent diameter: ", Nu1)
println("Nusselt number for annulus: ", Nu2)