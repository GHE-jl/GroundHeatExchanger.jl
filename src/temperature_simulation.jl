include("temporal_superposition.jl")

#TODO: Validate that this file is complete. Check to put the function g_model here?

"""
    fluid_temperature(t, q, g, T0, ks, Rb)
    fluid_temperature(t, Q, g, H, T0, ks, Rb)

Function that computes the average temperature of a ground heat exchanger based on the response
function and the heat load profile. The general function to solve is
`T_out(t) = T0 + q(t)⋅Rb + (f*g)(t)` (see Eq. 8 of Pasquier et al. 2018)
# Arguments
    - `t`: Time vector [s]
    - `q`: Heat load vector [W/m]
    - `Q`: Total heat load [W]
    - `g`: Response function for an impulse of 1 W/m (output of ground_models)[°Cm/W]
    - `H`: Borehole depth [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Rb`: Borehole thermal resistance (either single, double or coaxial U-loop) [°Cm/W]
# Output
    - `T_out`: Outlet temperature vector [°C]
# Reference
    - Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks to
        near-instant construction of short-term g-functions. Applied Thermal Engineering, 143,
        910–921. https://doi.org/10.1016/j.applthermaleng.2018.07.137
"""
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    g::AbstractVector{<:Real}, T0::Real, ks::Real, Rb::Real)  
    # Ensure that q and g are vectors of the same length as t
    if length(q) != length(t) || length(g) != length(t)
        throw(ArgumentError("Length of q and g must match length of t"))
    end    
    return T0 .+ q .* Rb .+ convolution(q, g) / (2 * π * ks)
end
function fluid_temperature(t::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    g::AbstractVector{<:Real}, H::Real, T0::Real, ks::Real, Rb::Real)
    # Ensure that Q and g are vectors of the same length as t
    if length(Q) != length(t) || length(g) != length(t)
        throw(ArgumentError("Length of Q and g must match length of t"))
    end    
    q = Q / H
    return T0 .+ q .* Rb .+ convolution(q, g) / (2 * π * ks)
end

"""
    outlet_temperature(Tf, Q, V, C)
    outlet_temperature(Tf, q, H, V, C)

Function that computes the outlet temperature of a ground heat exchanger based on the average fluid
temperature at the borehole outlet (Tin+Tout)/2. The outlet is solved with `Tout = 2*Tf - Tin`.
# Arguments
    - `Tf`: Average fluid temperature at the borehole outlet (Tin+Tout)/2 [°C]
    - `Q`: Total heat load [W]
    - `q`: Heat load vector [W/m]
    - `H`: Borehole depth [m]
    - `V`: Fluid flow rate [m³/s]
    - `C`: Fluid volumetric specific heat [J/m³K]
# Output
    - `Tout`: Outlet temperature vector [°C]
"""
function outlet_temperature(Tf::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    V::Real, C::Real)
    # Ensure that Q and Tf are vectors of the same length
    if length(Q) != length(Tf)
        throw(ArgumentError("Length of Q and Tf must match"))
    end
    return Tf .- Q ./ (2 * V * C)
end
function outlet_temperature(Tf::AbstractVector{<:Real}, q::AbstractVector{<:Real}, H::Real,
    V::Real, C::Real)
    # Ensure that q and Tf are vectors of the same length
    if length(q) != length(Tf)
        throw(ArgumentError("Length of q and Tf must match"))
    end
    return Tf .- (q .* H) ./ (2 * V * C)
end

"""
    inlet_temperature(Tf, Q, V, C)
    inlet_temperature(Tf, q, H, V, C)

Function that computes the inlet temperature of a ground heat exchanger based on the average fluid
temperature at the borehole outlet (Tin+Tout)/2. The inlet is solved with `Tin = 2*Tf - Tout`.
# Arguments
    - `Tf`: Average fluid temperature at the borehole outlet (Tin+Tout)/2 [°C]
    - `Q`: Total heat load [W]
    - `q`: Heat load vector [W/m]
    - `H`: Borehole depth [m]
    - `V`: Fluid flow rate [m³/s]
    - `C`: Fluid volumetric specific heat [J/m³K]
"""
function inlet_temperature(Tf::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    V::Real, C::Real)
    # Ensure that Q and Tf are vectors of the same length
    if length(Q) != length(Tf)
        throw(ArgumentError("Length of Q and Tf must match"))
    end
    return Tf .+ Q ./ (2 * V * C)
end
function inlet_temperature(Tf::AbstractVector{<:Real}, q::AbstractVector{<:Real}, H::Real,
    V::Real, C::Real)
    # Ensure that q and Tf are vectors of the same length
    if length(q) != length(Tf)
        throw(ArgumentError("Length of q and Tf must match"))
    end
    return Tf .+ (q .* H) ./ (2 * V * C)
end