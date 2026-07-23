# API reference

This page indexes every symbol **defined in** `GroundHeatExchanger.jl`. Symbols re-exported from
`BoreholeResistance.jl` and `GroundResponse.jl` are documented in their own sites (see
[Re-exported APIs](@ref)).

```@index
Modules = [GroundHeatExchanger]
```

## By topic

### Temporal superposition functions

- [`convolution`](@ref), [`convolutionf`](@ref) — stationary FFT convolution
- [`impulse_func`](@ref) — incremental load function
- [`convolution_ns`](@ref), [`convolution_ns!`](@ref) — non-stationary convolution
- [`impulse_func_ns`](@ref) — per-state incremental load matrix
- [`state_indices`](@ref), [`state_vector`](@ref) — operating-state bookkeeping

### Fluid temperatures

- [`fluid_temperature`](@ref) — mean fluid temperature
- [`outlet_temperature`](@ref), [`inlet_temperature`](@ref) — ends of the loop
- [`outlet_transfer_function`](@ref) — complete short-/long-term outlet (EWT) transfer function

### Short-term ANN model

- [`short_term_response`](@ref) — ANN-based short-term outlet transfer function
- [`short_term_nodes`](@ref) — raw ANN transfer function on its native time nodes (the building
  block behind `short_term_response` and `outlet_transfer_function`)
- [`AbstractANN`](@ref) — model tag selecting which trained network to evaluate
  - [`DeepANN`](@ref) — Pasquier & Marcotte (2020); wider validity ranges, 95 nodes, 21-day
    horizon; **the default** for `outlet_transfer_function`
  - [`PublishedANN`](@ref) — Pasquier, Zarrella & Labib (2018); 85 nodes, 7-day horizon

### Ground response and compression

- `ground_response` — g-function interface, re-exported from
  [GroundResponse.jl](https://GHE-jl.github.io/GroundResponse.jl); compression via its `interp` keyword
- [`pchip_interpolation`](@ref) — shape-preserving reconstruction

### Utility functions

- [`ground_load_profile`](@ref) — synthetic annual load
- [`GHE`](@ref) — example parameter set

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Temporal superposition](@ref), [Fluid temperature](@ref), [g-function compression](@ref)
    and [Utilities](@ref).

!!! note "Re-exported symbols are not indexed here"
    The index above lists only symbols defined in this package. The water-property, resistance,
    ground-model and borefield functions reached through `using GroundHeatExchanger` are documented
    in the [BoreholeResistance.jl](https://GHE-jl.github.io/BoreholeResistance.jl) and
    [GroundResponse.jl](https://GHE-jl.github.io/GroundResponse.jl) sites.
