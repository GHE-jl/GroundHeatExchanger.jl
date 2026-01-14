using SpecialFunctions: besselj0, besselj1, bessely0, bessely1
using QuadGK: quadgk

"""
    ics(t, ks, Cs, rb, rm=rb)

Computes the infinite cylindre source (ICS) model based on Carlsaw and Jaeger (1959). The output
is a g-function that requires a heat load per unit of borehole length [W/m] to provide the
borehole wall temperature.
# Arguments
    - t: Time vector [s]
    - ks: Ground thermal conductivity [W/mK]
    - Cs: Ground volumetric specific heat [J/m³K]
    - rb: Borehole radius (or any radius at which to compute the ICS) [m]
    - rm: (optional) Radius where the model is evaluated (usually equal to rb) [m]
# Output
    - g: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference:
    Ingersoll, L. R., Zabel, O. J., Ingersoll, A. C., & others. (1954). Heat conduction with 
    engineering, geological, and other applications. University of Wisconsin Press.
"""
function ics(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, rm::Real = rb)
    # Set initial parameters
    nt = length(t)                      # Number of time step
    r̃ = float(rm / rb)                  # Ratio of location to the temperature and cylinder radius
    t̃ = t .* ks ./ (Cs * rb^2)          # Fourier number
    g = Vector{Float64}(undef, nt)      # Preallocation

    function integrand_ics(s::Real, r̃::Real, tᵢ::Real)
        """
            integrand_ics(s, r̃, tᵢ)

        Integrand of the ICS model for scalar call.
        """
        if s < 1e-12
            return 0.0
        end
        return (exp(-s^2 * tᵢ) - 1) * (((besselj0(r̃ * s) * bessely1(s)) -
            (bessely0(r̃ * s) * besselj1(s))) / (s^2 * (besselj1(s)^2 + bessely1(s)^2)))
    end

    for (i, tᵢ) in enumerate(t̃)
        integral, _ = quadgk(s -> integrand_ics(s, r̃, tᵢ), 1e-6, Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (π^2 * ks)
end