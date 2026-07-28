# Tutorial

This tutorial runs a complete one-year, hourly simulation of a single borehole, from geometry to
inlet/outlet temperatures, then shows the borefield and non-stationary variants. Every function used
here is documented in the [API reference](@ref).

## 1. Geometry, ground and a time vector

```julia
using GroundHeatExchanger

H, D, rb, s, ro, ri = 150.0, 2.0, 0.08, 0.05, 0.022, 0.017   # geometry [m]
ks, Cs, kg, kp      = 3.0, 2.11e6, 1.6, 0.4                   # ground + materials
T0, V               = 10.0, 30/6e4                            # 10 °C, 30 L/min

t = collect(3600.0:3600:3600*24*365)    # 1 year, hourly [s]
```

## 2. Fluid properties and borehole resistance

The water properties and the resistance functions are re-exported from `BoreholeResistance.jl`.
Note the two heat-capacity conventions: ``c_f`` (mass-specific) for the resistance, ``C_f =
c_f\,\rho_f`` (volumetric) for the inlet/outlet split.

```julia
kf = water_k(T0);  cf = water_cp(T0);  ρf = water_ρ(T0);  μf = water_μ(T0)
Cf = cf * ρf                              # volumetric specific heat [J/m³·K]

Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
```

## 3. A ground model and a load

Pick a ground response model (here the finite line source) and generate an annual load profile with
[`ground_load_profile`](@ref) (which takes time in **hours**):

```julia
model = FLSModel(H, D, ks, Cs)            # re-exported from GroundResponse
Q     = ground_load_profile(t ./ 3600)      # [W]
```

## 4. Mean fluid temperature

[`fluid_temperature`](@ref) ties it together. Pass the per-length load `Q ./ H`, the model, the
geometry and the resistance — the g-function is computed and PCHIP-compressed automatically:

```julia
Tf = fluid_temperature(t, Q ./ H, model, rb, T0, Rb)
```

Under the hood this evaluates ``T_f = T_0 + q\,R_b^* + (q \star g)`` — see
[Simulation pipeline](@ref).

## 5. Outlet and inlet temperatures

The two ends split symmetrically about the mean. Use the **volumetric** heat capacity `Cf` here:

```julia
Tout = outlet_temperature(Tf, Q, V, Cf)
Tin  = inlet_temperature(Tf, Q, V, Cf)
```

`Tin` and `Tout` straddle `Tf` by ``\pm Q/(2 V C_f)``. See [Fluid temperature](@ref).

## 6. A borehole field

To simulate a field instead of a single borehole, build a coordinate matrix with `borefield`
(re-exported from `GroundResponse.jl`) and pass it to the matching `fluid_temperature` overload —
spatial superposition is applied inside `ground_response`:

```julia
xy = borefield(:rectangle, 3, 4, 6.0)     # 3×4 grid, 6 m spacing
Tf = fluid_temperature(t, Q ./ (H * size(xy, 1)), model, rb, xy, T0, Rb)
```

## 7. Controlling g-function compression

The `interp` keyword controls PCHIP compression. The default (`true`) is accurate for most cases;
set `interp=false` to evaluate the g-function exactly at every step (slower, for reference):

```julia
Tf_exact = fluid_temperature(t, Q ./ H, model, rb, T0, Rb; interp=false)
```

See [g-function compression](@ref).

## 8. Time-varying operating conditions

If the flow rate or fluid properties change during the year, the response is no longer described by
a single g-function. Build one g-function column per operating state and use the non-stationary
convolution ([`convolution_ns`](@ref)):

```julia
ks_states = [3.0, 2.0, 1.5]
g_matrix  = hcat([fls(t, rb, H, D, ks_i, Cs) for ks_i in ks_states]...)

ind, state, _ = state_indices(operating_condition_vec)
state_vec     = state_vector(ind, state, length(t))
ΔT            = convolution_ns(Q ./ H, g_matrix, state_vec)
```

See [Temporal superposition](@ref).

## Validation scripts

The `script/` directory contains runnable scripts. Run them from the package root (the temperature
script's filename contains a space — quote it):

```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=script/ "script/script_temperature simulation.jl"
```

| Script | What it shows |
|---|---|
| `script_temperature simulation.jl` | Full pipeline: FLS g → ``R_b^*`` → ``T_f`` / ``T_{in}`` / ``T_{out}``, 1-year hourly. |
| `script_temporal_superposition_stationary.jl` | FFT convolution with an FLS g-function over a 6-day load. |
| `script_temporal_superposition_nonstationary.jl` | Non-stationary convolution across three operating states. |
| `script_outlet_transfer_function.jl` | Short- + long-term outlet transfer function: Pasquier et al. (2018) Figs. 2a/3, plus a `PublishedANN` vs `DeepANN` borefield comparison. |
| `script_ann_validation.jl` | `PublishedANN` and `DeepANN` plotted against their MATLAB references at native nodes, with the interpolated `short_term_response` overlaid. |
| `script_validation_Lamarche2023.jl` | `@test`-based validation of g-functions, resistances and fluid properties against Lamarche's textbook examples. |
