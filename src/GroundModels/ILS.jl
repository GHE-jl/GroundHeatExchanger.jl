using SpecialFunctions: expinti

"""
    ils(t, ks, Cs, r)

Compute the infinite line source (ILS) model based on Ingersol (1954). If calculated at the borehole
radius, the output is a g-function that requires a heat load per unit of borehole length [W/m] to 
provide the borehole wall temperature.
# Arguments
    - t: Time vector [s]
    - ks: Ground thermal conductivity [W/mK]
    - Cs: Ground volumetric specific heat [J/m³K]
    - r: Radius at which to compute the ILS (typically the borehole radius) [m]
# Output
    - g: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. 
    ASHVE Journal Section, Heating, Piping and Air Conditioning.
"""
function ils(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, r::Real)
    
    g = -expinti.(-r^2 ./ (4 * (ks / Cs) * t)) ./ (4 * π * ks)
    return g
end