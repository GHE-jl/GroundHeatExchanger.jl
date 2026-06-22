"""
    ground_response(t, rb, xy, m; n_nodes=150)

GHE-level g-function interface. Wraps `GroundResponse.ground_response` with automatic PCHIP
interpolation: when `length(t) > n_nodes`, evaluates on a logarithmically-spaced subset and
reconstructs the full time vector via PCHIP interpolation. Pass `n_nodes=0` to skip interpolation.
# Arguments
    - `t`: Time vector [s]
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — `[0.0 0.0]` for a single borehole
    - `m`: Ground model (`FLSModel`, `ILSModel`, `ICSModel`, `MILSModel`, `MFLSModel`)
    - `n_nodes`: Max evaluation nodes before PCHIP interpolation; `0` disables (default 150)
# Output
    - `g`: g-function [°C·m/W], same length as `t`
"""
function ground_response(t::AbstractVector{<:Real}, rb::Real,
    xy::AbstractMatrix{<:Real}, m::AbstractGroundModel; n_nodes::Int=150)
    tv = collect(Float64, t)
    if n_nodes > 0 && length(tv) > n_nodes
        id = set_nodes(length(tv), n_nodes)
        gᵢ = GroundResponse.ground_response(tv[id], rb, xy, m)
        return pchip_interpolation(tv[id], gᵢ, tv)
    end
    return GroundResponse.ground_response(tv, rb, xy, m)
end

"""
    fluid_temperature(t, q, g, T0, ks, Rbₑ)
    fluid_temperature(t, Q, g, H, T0, ks, Rbₑ)
    fluid_temperature(t, q, m, rb, T0, ks, Rbₑ; n_nodes=150)
    fluid_temperature(t, q, m, rb, xy, T0, ks, Rbₑ; n_nodes=150)
    fluid_temperature(t, Q, H, m, rb, T0, ks, Rbₑ; n_nodes=150)

Compute the mean fluid temperature at the borehole mid-depth:
`Tf(t) = T0 + q(t)·Rbₑ + (q*g)(t) / (2π·ks)`
where q*g is the temporal superposition (FFT convolution) of the incremental heat load with the
g-function. When called with a ground model, `ground_response` is invoked internally with
`n_nodes` — pass `n_nodes=0` to skip PCHIP compression.
# Arguments
    - `t`: Time vector [s]
    - `q`: Heat load per unit length [W/m]
    - `Q`: Total heat load [W]
    - `g`: Pre-computed g-function [°C·m/W] — from `ground_response`
    - `H`: Borehole depth [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `ks`: Ground thermal conductivity [W/m·K]
    - `Rbₑ`: Effective borehole thermal resistance [m·K/W]
    - `m`: Ground model (`FLSModel`, `ILSModel`, `MFLSModel`, etc.)
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — omit for a single borehole
    - `n_nodes`: Nodes for PCHIP compression; `0` disables (default 150)
# Output
    - `Tf`: Mean fluid temperature vector [°C]
# Reference
    - Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks to
        near-instant construction of short-term g-functions. Applied Thermal Engineering, 143,
        910–921.
"""
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    g::AbstractVector{<:Real}, T0::Real, ks::Real, Rbₑ::Real)
    length(q) == length(t) && length(g) == length(t) ||
        throw(ArgumentError("Length of q and g must match length of t"))
    return T0 .+ q .* Rbₑ .+ convolution(q, g) / (2 * π * ks)
end
function fluid_temperature(t::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    g::AbstractVector{<:Real}, H::Real, T0::Real, ks::Real, Rbₑ::Real)
    length(Q) == length(t) && length(g) == length(t) ||
        throw(ArgumentError("Length of Q and g must match length of t"))
    q = Q ./ H
    return T0 .+ q .* Rbₑ .+ convolution(q, g) / (2 * π * ks)
end
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    m::AbstractGroundModel, rb::Real, T0::Real, ks::Real, Rbₑ::Real; n_nodes::Int=150)
    xy = reshape([0.0, 0.0], 1, 2)
    g = ground_response(t, rb, xy, m; n_nodes=n_nodes)
    return fluid_temperature(t, q, g, T0, ks, Rbₑ)
end
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    m::AbstractGroundModel, rb::Real, xy::AbstractMatrix{<:Real},
    T0::Real, ks::Real, Rbₑ::Real; n_nodes::Int=150)
    g = ground_response(t, rb, xy, m; n_nodes=n_nodes)
    return fluid_temperature(t, q, g, T0, ks, Rbₑ)
end
function fluid_temperature(t::AbstractVector{<:Real}, Q::AbstractVector{<:Real}, H::Real,
    m::AbstractGroundModel, rb::Real, T0::Real, ks::Real, Rbₑ::Real; n_nodes::Int=150)
    xy = reshape([0.0, 0.0], 1, 2)
    g = ground_response(t, rb, xy, m; n_nodes=n_nodes)
    return fluid_temperature(t, Q, g, H, T0, ks, Rbₑ)
end

"""
    outlet_temperature(Tf, Q, V, C)
    outlet_temperature(Tf, q, H, V, C)

Borehole outlet temperature from the mean fluid temperature: `Tout = Tf − Q/(2·V·C)`.
# Arguments
    - `Tf`: Mean fluid temperature [°C]
    - `Q`: Total heat load [W]
    - `q`: Heat load per unit length [W/m]
    - `H`: Borehole depth [m]
    - `V`: Volumetric flow rate [m³/s]
    - `C`: Fluid volumetric specific heat [J/m³·K] (`water_cp(T) * water_ρ(T)`)
# Output
    - `Tout`: Outlet fluid temperature [°C]
"""
function outlet_temperature(Tf::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    V::Real, C::Real)
    length(Q) == length(Tf) || throw(ArgumentError("Length of Q and Tf must match"))
    return Tf .- Q ./ (2 * V * C)
end
function outlet_temperature(Tf::AbstractVector{<:Real}, q::AbstractVector{<:Real}, H::Real,
    V::Real, C::Real)
    length(q) == length(Tf) || throw(ArgumentError("Length of q and Tf must match"))
    return Tf .- (q .* H) ./ (2 * V * C)
end

"""
    inlet_temperature(Tf, Q, V, C)
    inlet_temperature(Tf, q, H, V, C)

Borehole inlet temperature from the mean fluid temperature: `Tin = Tf + Q/(2·V·C)`.
# Arguments
    - `Tf`: Mean fluid temperature [°C]
    - `Q`: Total heat load [W]
    - `q`: Heat load per unit length [W/m]
    - `H`: Borehole depth [m]
    - `V`: Volumetric flow rate [m³/s]
    - `C`: Fluid volumetric specific heat [J/m³·K] (`water_cp(T) * water_ρ(T)`)
# Output
    - `Tin`: Inlet fluid temperature [°C]
"""
function inlet_temperature(Tf::AbstractVector{<:Real}, Q::AbstractVector{<:Real},
    V::Real, C::Real)
    length(Q) == length(Tf) || throw(ArgumentError("Length of Q and Tf must match"))
    return Tf .+ Q ./ (2 * V * C)
end
function inlet_temperature(Tf::AbstractVector{<:Real}, q::AbstractVector{<:Real}, H::Real,
    V::Real, C::Real)
    length(q) == length(Tf) || throw(ArgumentError("Length of q and Tf must match"))
    return Tf .+ (q .* H) ./ (2 * V * C)
end
