"""
Script showing the thermal resistance models from ThermalResistances.jl using either coaxial, 
single or double U-loop ground heat exchanger.
"""

using BenchmarkTools

includet("../src/ThermalResistances.jl")

# Define parameters
r = (b = 0.08,                  # Borehole radius
    o = 0.022,                  # Pipe outlet radius
    i = 0.017)                  # Pipe inlet radius
k = (s = 3.0,                   # Ground thermal conductivity
    g = 1.6,                    # Grout thermal conductivity
    p = 0.4,                    # Pipe thermal conductivity
    f = 0.6)                    # Fluid thermal conductivity
c = (s = 2.11e3,                # Ground specific heat
    g = 2.25e3,                 # Grout specific heat
    p = 1.9e3,                  # Pipe specific heat
    f = 4.2e3)                  # Fluid specific heat
ρ = (s = 1000.0,                # Groud density
    g = 1000.0,                 # Grout density
    p = 1000.0,                 # Pipe density
    f = 1000.0)                 # Fluid density
s = 0.05                        # Shank spacing (s/2 is the half-shank spacing)
V = 30.0 / 60000                # Circulating flow rate
V̇ = V / (π * r.i^2)             # Fluid speed in pipe
μf = 1.3e-3                     # Fluid viscosity
H = 150.0                       # Borehole length
n = 2                           # Number of pipes in the borehole

# Fluid thermal resistance
@show Rf = R_f(V̇, k.f, r.i, c.f, ρ.f, μf)
V̇ᵥ = range(5, 60, 1000) / (60000 * π * r.i^2)
RF = zeros(1000)
for (i, j) in enumerate(V̇ᵥ)
    RF[i] = R_f(j, k.f, r.i, c.f, ρ.f, μf)
end

# Pipe thermal resistance
@show Rp = R_p(k.p, r.o, r.i)

# Grout thermal resistance
@show Rg_0 = R_b_zeroth_order_multipole(k.s, k.g, r.b, r.o, s, 0.0, 0.0)
@show Rg_1 = R_b_first_order_multipole(k.s, k.g, r.b, r.o, s, 0.0, 0.0)

# Single U-loop
@show Rb_0 = R_b_zeroth_order_multipole(k.s, k.g, r.b, r.o, s, Rp, Rf)
@show Rb_1 = R_b_first_order_multipole(k.s, k.g, r.b, r.o, s, Rp, Rf)

@show Rb = R_b(V̇, k.s, k.g, k.p, k.f, r.b, r.o, r.i, s, n, c.f, ρ.f, μf)

@show Ra = R_a_first_order_multipole(k.s, k.g, r.b, r.o, s, Rp, Rf)

@show Rbₑ1 = R_bₑ(V, c.f, ρ.f, 150.0, Rb, Ra)
@show Rbₑ2 = R_bₑ(V, k.s, k.g, r.b, r.o, s, c.f, ρ.f, 150.0, Rp, Rf)
@show Rbₑ3 = R_bₑ(V, k.s, k.g, k.p, k.f, r.b, r.o, r.i, s, c.f, ρ.f, μf, H)

# Validation with data from Javed and Spitler 2017
#Rg = R_b_first_order_multipole(3.0, 0.6, 0.192, 0.032, 2*0.032, 0.00, 0.0) # see Figure 6. (left)