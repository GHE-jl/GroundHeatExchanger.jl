using SpecialFunctions: besselj0, besselj1, bessely0, bessely1
using QuadGK: quadgk

"""
    _ics_integrand(s, r̃, t̃)

Computes the integrand of the infinite cylindrical source model.
"""
function _ics_integrand(s::T, r̃::T, t̃::T) where {T<:AbstractFloat}
    if s < 1e-12
        return zero(T)
    end
    # Pre-calculate Bessel terms to keep the expression readable
    j0_rs = besselj0(r̃ * s)
    y1_s  = bessely1(s)
    y0_rs = bessely0(r̃ * s)
    j1_s  = besselj1(s)
    
    return ((exp(-s^2 * t̃) - 1) * (j0_rs * y1_s - y0_rs * j1_s)) / (s^2 * (j1_s^2 + y1_s^2))
end

"""
    _ics(t, ks, Cs, r, rc)

Kernel function for the infinite cylindrical source model based on Carlsaw and Jaeger (1959). The 
response function is based on an impulse of 1 W/m.
"""
function _ics(t::T, ks::T, Cs::T, r::T, rc::T) where {T<:AbstractFloat}
    # Define parameters
    r̃ = r / rc
    t̃ = (t * ks) / (Cs * rc^2)
    
    # Perform numerical integration. T(1e-8) ensures the lower bound matches input precision.
    integral, _ = quadgk(s -> _ics_integrand(s, r̃, t̃), T(1e-8), T(Inf), rtol = T(1e-6))
    
    return integral / (T(π)^2 * ks)
end

"""
    ics(t, ks, Cs, r, rc)

Computes the infinite cylindre source (ICS) model based on Carlsaw and Jaeger (1959). The output is
a g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole 
wall temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
    - `rc`: Radius of the cylinder [m]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference:
    - Carslaw, H. S., & Jaeger, J. C. (1959). Conduction of heat in solids. Oxford: Clarendon Press,
        1959, 2nd Ed.
# Example
    ics(60:60:3600, 3.0, 2e6, 0.076, 0.076)
"""
function ics(t::Real, ks::Real, Cs::Real, r::Real, rm::Real = r)
    T = float(promote_type(typeof(t), typeof(ks), typeof(Cs), typeof(r), typeof(rm)))
    return _ics(T(t), T(ks), T(Cs), T(r), T(rm))
end

function ics(t::AbstractVector{<:Real}, ks::Real, Cs::Real, r::Real, rm::Real = r)
    # Check type
    T = float(promote_type(eltype(t), typeof(ks), typeof(Cs), typeof(r), typeof(rm)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and ICS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _ics(t_T[i], T(ks), T(Cs), T(r), T(rm))
    end
    return g
end

function ics!(g::AbstractVector{T}, t::AbstractVector, ks::Real, Cs::Real, r::Real, rm::Real = r
    ) where {T<:AbstractFloat}
    # Check for same vector length
    @assert length(g) == length(t)

    # Convert parameters to T once
    ks_T, Cs_T, rb_T, rm_T = T(ks), T(Cs), T(r), T(rm)
    
    # ICS
    @inbounds @simd for i in eachindex(g, t)
        g[i] = _ics(T(t[i]), ks_T, Cs_T, rb_T, rm_T)
    end
    return g
end

function ics_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, r::Real, rm::Real = r)
    # Set initial parameters
    nt = length(t)                      # Number of time step
    r̃ = float(rm / r)                  # Ratio of location to the temperature and cylinder radius
    t̃ = t .* ks ./ (Cs * r^2)          # Fourier number
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