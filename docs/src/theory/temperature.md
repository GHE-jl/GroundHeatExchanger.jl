# Fluid temperature

This page documents the three temperature outputs of the simulation: the **mean** fluid temperature
``T_f`` and the **outlet** and **inlet** temperatures that bracket it. Together they are what a
sizing or control calculation actually consumes.

## Mean fluid temperature

[`fluid_temperature`](@ref) evaluates the master equation (see [Simulation pipeline](@ref)):

```math
T_f(t) = T_0 + q(t)\,R_b^* + (q \star g)(t).
```

The first term is the undisturbed ground temperature, the second the instantaneous resistance drop
across the borehole, and the third the cumulative ground temperature rise from the
[Temporal superposition](@ref) of the load with the g-function.

The function comes in several overloads of increasing convenience:

```julia
# With a pre-computed g-function vector
fluid_temperature(t, q, g, T0, Rb)          # q in W/m

# With a ground model — g computed and PCHIP-compressed automatically
fluid_temperature(t, q, m, rb, T0, Rb;        interp=true)   # single borehole
fluid_temperature(t, q, m, rb, xy, T0, Rb;    interp=true)   # borefield (xy: nb×2)
```

The load is always supplied as the heat load **per unit length** `q` [W/m]. To use a total heat
load `Q` [W], divide it yourself by the **entire** borehole length `L = H · (number of boreholes)`
before calling: `q = Q / L`. Overloads taking `Q` and `H` directly were removed because dividing
`Q` by a single borehole depth `H` instead of the full field length is an easy mistake to make.

The model overloads call `ground_response` (re-exported from `GroundResponse.jl`) for you, so the
g-function is evaluated on a logarithmic node subset and reconstructed by PCHIP interpolation — pass
`interp=false` to disable that and evaluate exactly. See [g-function compression](@ref).

## Outlet and inlet temperatures

The mean temperature ``T_f`` is the depth-average of the fluid. The fluid actually *enters* the
borehole warmer (in injection) and *leaves* it cooler, split symmetrically about the mean by the
sensible heat it exchanges:

```math
T_{out} = T_f - \frac{Q}{2\,V\,C_f}, \qquad
T_{in}  = T_f + \frac{Q}{2\,V\,C_f},
```

where ``Q`` [W] is the total load, ``V`` [m³/s] the volumetric flow rate and ``C_f`` [J/m³·K] the
**volumetric** fluid heat capacity. [`outlet_temperature`](@ref) and [`inlet_temperature`](@ref)
each accept either the total load `Q` or the per-length load `q` together with the depth `H`:

```julia
Tout = outlet_temperature(Tf, Q, V, Cf)        # total load Q [W]
Tout = outlet_temperature(Tf, q, H, V, Cf)     # per-length load q [W/m]
Tin  = inlet_temperature(Tf, Q, V, Cf)
```

!!! warning "Use the volumetric heat capacity here"
    ``C_f`` in these formulas is the **volumetric** specific heat ``C_f = c_f\,\rho_f`` [J/m³·K],
    *not* the mass-specific ``c_f`` [J/kg·K] used by `resistance_ULoop_effective`. Compute it as
    `Cf = water_cp(T) * water_ρ(T)`.

## Short-term outlet transfer function

The models above describe the long-term (Eskilson-type) response. For the first hours to days the
borehole geometry, thermal capacity and vertical fluid advection dominate, and a dedicated
short-term transfer function is required. [`outlet_transfer_function`](@ref) builds the complete
dimensionless transfer function of the borehole **outlet** temperature (``T_{out}``, i.e. the
source-side entering water temperature) from minutes to decades by joining a trained short-term
ANN — [`short_term_response`](@ref) / [`short_term_nodes`](@ref), either the wider-range
[`DeepANN`](@ref) (Pasquier & Marcotte, 2020; the **default**, 21-day horizon) or the original
[`PublishedANN`](@ref) (Pasquier, Zarrella & Labib, 2018; 7-day horizon) — to any long-term
`AbstractGroundModel`. The long-term borehole-wall response ``\tilde g`` [°C·m/W] is converted to an
outlet transfer function through the effective borehole resistance ``R_b^*`` and shifted to meet
the short-term curve at the model's contact time (7 or 21 days):

```math
\bar g = \left(N_b\,\tilde g + R_b^*\right)\frac{\dot V\,C_f}{H},
```

where ``N_b`` rescales a borefield's field-average response (normalised by `GroundResponse.jl` to
1 W/m of *total* field load) back to the per-borehole basis the ANN uses (``N_b = 1`` for a single
borehole).

```julia
Rbₑ   = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
m     = FLSModel(H, D_buried, ks, Cs)
g_out = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s, Rbₑ, m)
```

## Functions on this page

```@docs
fluid_temperature
outlet_temperature
inlet_temperature
outlet_transfer_function
AbstractANN
DeepANN
PublishedANN
short_term_response
short_term_nodes
```
