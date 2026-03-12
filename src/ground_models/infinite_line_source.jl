using SpecialFunctions: expinti

"""
    _ils(t, r, ks, Cs)

Kernel function for the infinite line source model based on Ingersol (1954). The response function 
is based on an impulse of 1 W/m.
"""
@inline function _ils(t::T, r::T, ks::T, Cs::T) where {T<:AbstractFloat}
    α = ks / Cs
    x = -r^2 / (4 * α * t)
    return -expinti(x) / (4 * T(π) * ks)
end

"""
    ils(t, r, ks, Cs)

Compute the infinite line source (ILS) model based on Ingersol (1954). The output is a g-function
that requires a heat load per unit of borehole length [W/m] to provide the borehole wall
temperature.
# Arguments
    - `t`: Time value or vector [s]
    - `r`: Radius at which to computed (typically the borehole radius) [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. 
        ASHVE Journal Section, Heating, Piping and Air Conditioning.
# Example
    ils(60:60:3600, 3.0, 2e6, 0.076)
"""
function ils(t::Real, r::Real, ks::Real, Cs::Real)
    T = float(promote_type(typeof(t), typeof(r), typeof(ks), typeof(Cs)))
    return _ils(T(t), T(r), T(ks), T(Cs))
end
function ils(t::AbstractVector{<:Real}, r::Real, ks::Real, Cs::Real)
    # Check type
    T = float(promote_type(eltype(t), typeof(r), typeof(ks), typeof(Cs)))
    t_T  = convert(Vector{T}, t)

    # Preallocate and ILS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _ils(t_T[i], T(r), T(ks), T(Cs))
    end
    return g
end