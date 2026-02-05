using SpecialFunctions: erf
using QuadGK: quadgk

"""
    _ierf(x)

Inverse "erf" function used in the finite line source model.
"""
function _ierf(x::T) where {T<:AbstractFloat}
    return x * erf(x) - inv(sqrt(T(π))) * (one(T) - exp(-x^2))
end

"""
    _fls_integrand(s, r, H, D)

Computes the integrand of the finite line source model.
"""
function _fls_integrand(s::T, r::T, H::T, D::T) where {T<:AbstractFloat}
    # Calculate terms
    term1 = 2 * _ierf(H * s)
    term2 = 2 * _ierf(H * s + 2 * D * s)
    term3 = _ierf(2 * H * s + 2 * D * s)
    term4 = _ierf(2 * D * s)
    
    return (exp(-r^2 * s^2) * (term1 + term2 - term3 - term4)) / (H * s^2)
end

"""
    _fls(t, ks, Cs, r, H, D)

Kernel function for the finite line source model based on Claesson and Javed (2011). The response
function is based on an impulse of 1 W/m.
"""
function _fls(t::T, ks::T, Cs::T, r::T, H::T, D::T) where {T<:AbstractFloat}
    # The lower limit of the integral depends on time
    α = ks / Cs
    lower_lim = inv(sqrt(4 * α * t))
    
    # Perform numerical integration
    integral, _ = quadgk(s -> _fls_integrand(s, r, H, D), lower_lim, T(Inf), rtol = T(1e-6))
    
    return integral / (4 * T(π) * ks)
end

"""
    fls(t, ks, Cs, r, H, D)

Computes the finite line source (FLS) model based on Claesson and Javed (2011). The output is a 
g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
wall temperature.
# Arguments
    - `t`: Time vector [s]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
    - `H`: Borehole depth [m]
    - `D`: Buried depth [m]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid 
        temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 
        279–288.
# Example
    fls(60:60:3600, 3.0, 2e6, 0.076, 150, 4)
"""
function fls(t::Real, ks::Real, Cs::Real, r::Real, H::Real, D::Real)
    T = float(promote_type(typeof(t), typeof(ks), typeof(Cs), typeof(r), typeof(H), typeof(D)))
    return _fls(T(t), T(ks), T(Cs), T(r), T(H), T(D))
end

function fls(t::AbstractVector{<:Real}, ks::Real, Cs::Real, r::Real, H::Real, D::Real)
    # Check type
    T = float(promote_type(eltype(t), typeof(ks), typeof(Cs), typeof(r), typeof(H), typeof(D)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and FLS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _fls(t_T[i], T(ks), T(Cs), T(r), T(H), T(D))
    end
    return g
end

function fls!(g::AbstractVector{T}, t::AbstractVector, ks::Real, Cs::Real, r::Real, H::Real, 
    D::Real) where {T<:AbstractFloat}
    # Check for same vector length
    @assert length(g) == length(t)

    # Convert parameters to T once
    ks_T, Cs_T, rb_T, H_T, D_T = T(ks), T(Cs), T(r), T(H), T(D)
    
    # FLS
    @inbounds @simd for i in eachindex(g, t)
        g[i] = _fls(T(t[i]), ks_T, Cs_T, rb_T, H_T, D_T)
    end
    return g
end

function fls_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, r::Real, H::Real,
    D::Real)
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
        integral, _ = quadgk(s -> integrand_fls(s, r, H, D), lim_int[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (4 * π * ks)
end