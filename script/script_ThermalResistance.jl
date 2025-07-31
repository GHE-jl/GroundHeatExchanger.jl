"""
Script showing the thermal resistance models from ThermalResistance.jl using either single or double
or coaxial U-loop.
"""

includet("../src/ThermalResistance.jl")

# Parameters
r = (b = 0.08,                  # Borehole radius
    i = 0.017,                  # Pipe inlet radius
    o = 0.022)                  # Pipe outlet radius
k = (s = 3.0,                   # Ground thermal conductivity
    g = 1.6,                    # Grout thermal conductivity
    p = 0.4,                    # Pipe thermal conductivity
    f = 0.6)                    # Fluid thermal conductivity
C = (s = 2.11e6,                # Ground volumetric specific heat
    g = 2.25e6,                 # Grout volumetric specific heat
    p = 1.9e6,                  # Pipe volumetric specific heat
    f = 4.2e6)                  # Fluid volumetric specific heat
ρ = (s = 1000.0,                # Groud density
    g = 1000.0,                 # Grout density
    p = 1000.0,                 # Pipe density
    f = 1000.0)                 # Fluid density
V = 30.0 / 60000                # Circulating flow rate
V̇ = V / (π * r.i^2)

# Single U-loop
@time Rb.Hellstrom = Rb_Hellstrom(V̇, k.s, k.g, k.p, k.f, r.b, r.i, r.o, s, 
    n=2, cf = 4200.0, ρf = 1000., μf = 1.3e-3)

@time Rb.FOM = Rb_first_order_multipole(V̇, k.s, k.g, k.p, k.f, r.b, r.i, r.o, s, 
    n=2, cf = 4200.0, ρf = 1000., μf = 1.3e-3)

@time Ra.FOM = Ra_first_order_multipole(V̇, k.s, k.g, k.p, k.f, r.b, r.i, r.o, s, 
    n=2, cf = 4200.0, ρf = 1000., μf = 1.3e-3)