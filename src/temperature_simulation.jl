"""
    fluid_temperature(t, q, g, T0, ks, Rbₑ)
    fluid_temperature(t, q, m, rb, T0, ks, Rbₑ; interp=true)
    fluid_temperature(t, q, m, rb, xy, T0, ks, Rbₑ; interp=true)

Compute the mean fluid temperature at the borehole mid-depth:
`Tf(t) = T0 + q(t)·Rbₑ + (q*g)(t) / (2π·ks)`
where q*g is the temporal superposition (FFT convolution) of the incremental heat load with the
g-function. When called with a ground model, `GroundResponse.ground_response` is invoked internally
with `interp`, pass `interp=false` to evaluate the g-function exactly at every `t` (valid only for
a uniformly spaced `t` with the temporal solvers).
The load is always supplied as a heat load per unit length `q` [W/m]. If you have a total heat
load `Q` [W], divide it by the **entire** borehole length (`L = H · nb. boreholes`) beforehand:
`q = Q / L`. This is intentionally left to the caller to avoid the common confusion of dividing
`Q` by a single borehole depth `H` rather than by the full field length.
# Arguments
    - `t`: Time vector [s]
    - `q`: Heat load per unit length q = Q/L [W/m] (`L` is the total borehole length)
    - `g`: Pre-computed g-function [°C·m/W] — from `ground_response`
    - `T0`: Undisturbed ground temperature [°C]
    - `ks`: Ground thermal conductivity [W/m·K]
    - `Rbₑ`: Effective borehole thermal resistance [m·K/W]
    - `m`: Ground model (`FLSModel`, `ILSModel`, `MFLSModel`, etc.)
    - `rb`: Borehole radius [m]
    - `xy`: Borehole coordinates (nb × 2) [m] — omit for a single borehole
    - `interp`: enable GroundResponse's constant-step sub-sampling + PCHIP interpolation (default
      `true`); `false` evaluates exactly at every `t` (valid only for a uniformly spaced `t`)
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
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    m::AbstractGroundModel, rb::Real, T0::Real, ks::Real, Rbₑ::Real; interp::Bool=true)
    xy = reshape([0.0, 0.0], 1, 2)
    g = ground_response(collect(Float64, t), rb, xy, m; interp=interp)
    return fluid_temperature(t, q, g, T0, ks, Rbₑ)
end
function fluid_temperature(t::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    m::AbstractGroundModel, rb::Real, xy::AbstractMatrix{<:Real},
    T0::Real, ks::Real, Rbₑ::Real; interp::Bool=true)
    g = ground_response(collect(Float64, t), rb, xy, m; interp=interp)
    return fluid_temperature(t, q, g, T0, ks, Rbₑ)
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

"""
    outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s, Rbₑ, m;
                             xy=[0.0 0.0], clamp=true, interp=false)

Complete dimensionless transfer function of the borehole-**outlet** fluid temperature (`T_out`,
i.e. the entering water temperature on the heat-pump source side), spanning minutes to decades.
It joins the short-term ANN transfer function of Pasquier et al. (2018) to a long-term ground
model `m`:
  1. the short-term is `short_term_response` on its 7-day validity horizon;
  2. the long-term is the borehole-wall temperature response of `m` (already in [°C·m/W] per
     1 W/m — the `1/2πkₛ` conductivity scaling is baked into `ground_response`), converted to a
     dimensionless outlet transfer function through the effective borehole resistance `Rbₑ`
     (Eq. 9 of the paper), `ḡ = (nb·g_wall + Rbₑ) · (V̇·Cf / H)`;
  3. the long-term curve is shifted vertically to meet the short-term curve at the 7-day
     contact time, then the two halves are spliced.
`m` may be any `AbstractGroundModel`; spatial superposition (`successive_flux`) is applied
automatically when `xy` describes more than one borehole (in which case the result is the
field-average outlet transfer function under the equal-inlet-temperature assumption).
For a borefield, `ground_response` returns the field-average wall response normalised to 1 W/m of
*total* field heat rate (flux constraint `Σqᵢ = 1`), i.e. `1/nb` of the per-borehole response,
whereas the short-term ANN transfer function is per single borehole. The `nb` factor above rescales
the ground response back to a per-borehole basis so both halves join with a matching slope; without
it the long-term slope is `1/nb` too shallow and the junction kinks (`nb = 1` makes this a no-op).
The impulse of the outlet transfer function is of 1 °C (as opposed to 1 W/m for the models in
GroundResponse.jl).
# Arguments
    - `t`: Time vector [s] (minutes to decades)
    - `ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s`: physical inputs of the short-term ANN
        (see `short_term_response`; `s` is the shank spacing)
    - `Rbₑ`: Effective borehole thermal resistance `Rb*` [m·K/W]. The paper's *equivalent*
        borehole resistance, e.g. `resistance_ULoop_effective(...)`.
    - `m`: Long-term ground model (typically `FLSModel(H, D, ks, Cs)`; note its buried
        depth `D` is a distinct quantity from the ANN shank spacing `s`)
    - `xy`: Borehole coordinates (nb × 2) [m] — `[0.0 0.0]` for a single borehole
    - `clamp`: out-of-range handling for the ANN inputs (see `short_term_response`)
    - `interp`: enable GroundResponse's sub-sampling + PCHIP interpolation for the long-term
      evaluation; `false` evaluates exactly at every `t` (default `false`)
# Output
    - `g`: Dimensionless outlet transfer function at `t` [-] (add `1` for the inlet transfer
        function)
# Reference
    - Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks
        to near-instant construction of short-term g-functions. Applied Thermal Engineering,
        143, 910–921.
"""
function outlet_transfer_function(t::AbstractVector{<:Real}, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb,
    H, V̇, s, Rbₑ, m::AbstractGroundModel;
    xy::AbstractMatrix{<:Real} = reshape([0.0, 0.0], 1, 2), clamp::Bool = true,
    interp::Bool = false)

    tv = collect(Float64, t)

    # Short-term: ANN transfer function on its 85 native nodes; the last node is the
    # short-/long-term contact time (7 days).
    ts_nodes, g_nodes = short_term_nodes(ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro,
        rb, H, V̇, s; clamp = clamp)
    ts_nodes = Float64.(ts_nodes)
    tc = ts_nodes[end]                            # short-/long-term contact time (7 days)
    short = tv .<= tc
    g_short = any(short) ? pchip_interpolation(ts_nodes, g_nodes, tv[short]) : Float64[]
    g_st_contact = g_nodes[end]                   # short-term value at tc (== last node)

    # Long-term: field-average wall response → dimensionless outlet g (Eq. 9). `ground_response`
    # already carries the 1/2πkₛ scaling, but for a borefield it is normalised to 1 W/m of *total*
    # field heat rate (1/nb of the per-borehole response); the short-term ANN half is per single
    # borehole, so we rescale by nb to join them on the same per-borehole basis (nb = 1 is a no-op).
    nb = size(xy, 1)
    to_outlet(gw) = (nb .* gw .+ Rbₑ) .* (V̇ * Cf / H)
    # Prepend tc so the first evaluated value provides the continuity shift.
    tq = vcat(tc, tv[.!short])
    g_long_full = to_outlet(ground_response(tq, rb, xy, m; interp = interp))
    shift = g_long_full[1] - g_st_contact
    g_long = g_long_full[2:end] .- shift

    return vcat(g_short, g_long)
end
