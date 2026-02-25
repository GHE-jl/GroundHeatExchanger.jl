using SpecialFunctions: expinti

"""
    _ils(t, ks, Cs, r)

Kernel function for the infinite line source model based on Ingersol (1954). The response function 
is based on an impulse of 1 W/m.
"""
@inline function _ils(t::T, ks::T, Cs::T, r::T) where {T<:AbstractFloat}
    α = ks / Cs
    x = -r^2 / (4 * α * t)
    return -expinti(x) / (4 * T(π) * ks)
end

"""
    ils(t, ks, Cs, r)

Compute the infinite line source (ILS) model based on Ingersol (1954). The output is a g-function
that requires a heat load per unit of borehole length [W/m] to provide the borehole wall
temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. 
        ASHVE Journal Section, Heating, Piping and Air Conditioning.
# Example
    ils(60:60:3600, 3.0, 2e6, 0.076)
"""
function ils(t::Real, ks::Real, Cs::Real, r::Real)
    T = float(promote_type(typeof(t), typeof(ks), typeof(Cs), typeof(r)))
    return _ils(T(t), T(ks), T(Cs), T(r))
end

function ils(t::AbstractVector{<:Real}, ks::Real, Cs::Real, r::Real)
    # Check type
    T = float(promote_type(eltype(t), typeof(ks), typeof(Cs), typeof(r)))
    t_T  = convert(Vector{T}, t)

    # Preallocate and ILS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _ils(t_T[i], T(ks), T(Cs), T(r))
    end
    return g
end

function ils!(g::AbstractVector{T}, t::AbstractVector, ks::Real, Cs::Real, r::Real
    ) where {T<:AbstractFloat}
    # Check for same vector length
    @assert length(g) == length(t)

    # Convert parameters once to the target type
    ks_T, Cs_T, r_T = T(ks), T(Cs), T(r)

    # ILS
    @inbounds @simd for i in eachindex(t)
        g[i] = _ils(T(t[i]), ks_T, Cs_T, r_T)
    end
    return g
end