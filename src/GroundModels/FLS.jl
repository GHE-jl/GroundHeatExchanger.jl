using SpecialFunctions: erf
using QuadGK: quadgk

"""
    fls(t, ks, Cs, rb, H, D)

Computes the finite line source (FLS) model based on Claesson and Javed (2011). The output is a 
g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
wall temperature.
# Arguments
    - t: Time vector [s]
    - ks: Ground thermal conductivity [W/mK]
    - Cs: Ground volumetric specific heat [J/m³K]
    - rb: Borehole radius [m]
    - H: Borehole depth [m]
    - D: Buried depth [m]
# Output
    - g: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid 
    temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 
    279–288.
"""
function fls(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real, D::Real)
    # Set initial parameters
    const_π = 1 / sqrt(π)
    nt = length(t)              # Number of element in the time vector
    g = zeros(nt)               # Preallocation of the borehole wall temperature
    α = ks / Cs
    lim_int = 1 ./ sqrt.(4 * α * t)

    """
        integrand_fls(s, r, H, D)

    Integrand of the FLS model. Assumes constant heat flux boundary condition.
    """
    function integrand_fls(s::Real, r::Real, H::Real, D::Real)
        """
            ierf(x)

        Inverse "erf" function used in the FLS model.
        """
        function ierf(x::Real)
            return x * erf(x) - const_π * (1 - exp(-x^2))
        end

        return exp(-r^2 * s^2) * (2 * ierf(H * s) + 2 * ierf(H * s + 2 * D * s) -
                ierf(2 * H * s + 2 * D * s) - ierf(2 * D * s)) / (H * s^2)
    end
    # Compute, in a loop, each value of the fls
    for i in 1:nt
        integral, _ = quadgk(s -> integrand_fls(s, rb, H, D), lim_int[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (4 * π * ks)
end