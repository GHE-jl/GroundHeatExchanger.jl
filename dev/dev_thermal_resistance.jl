"""
Script to compute the thermal resistance of a double loop ground heat exchanger.
"""

includet("../src/borehole_thermal_resistance/resistance_fluid.jl")
includet("../src/borehole_thermal_resistance/resistance_pipe.jl")
includet("../src/borehole_thermal_resistance/resistance_borehole.jl")

# Define the parameters of the borehole
includet("../src/utils.jl")
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
V̇ = V / (π * ri^2)  # Volumetric flow rate of the fluid

# Compute thermal resistance of the fluid and pipe
Rf = resistance_fluid(V̇, ri, kf, Cf/ρf, ρf, μf)
Rp = resistance_pipe(ri, ro, kp)

# Compute the thermal resistance of the borehole for a double loop configuration
Rb0 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 0)
Rb1 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 1)

println("`Rb` (order 0): ", Rb0, " mK/W")
println("`Rb` (order 1): ", Rb1, " mK/W")

Ra0_1 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 0, "diagonal")
Ra0_2 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 0, "adjacent")
Ra1_1 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 1, "diagonal")
Ra1_2 = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 2, 1, "adjacent")

println("`Ra` (order 0, diagonal): ", Ra0_1, " mK/W")
println("`Ra` (order 0, adjacent): ", Ra0_2, " mK/W")
println("`Ra` (order 1, diagonal): ", Ra1_1, " mK/W")
println("`Ra` (order 1, adjacent): ", Ra1_2, " mK/W")