# Simulation pipeline

This page describes how the three packages combine into a temperature simulation and states the
master equation that the rest of the theory pages expand on.

## The three layers

A ground heat exchanger simulation needs three physically distinct pieces, supplied by three
packages:

| Layer | Question it answers | Package |
|---|---|---|
| **Fluid → borehole wall** | How large is the temperature drop across the grout and pipes? | `BoreholeResistance.jl` → ``R_b^*`` |
| **Borehole wall → ground** | How does the ground temperature evolve under a unit load? | `GroundResponse.jl` → ``g(t)`` |
| **Load history → temperature** | How do these combine over a time-varying load? | `GroundHeatExchanger.jl` |

`GroundHeatExchanger.jl` owns the third layer: it convolves the load history with the g-function and
adds the resistance drop to recover the fluid temperature.

## The master equation

The mean fluid temperature at the borehole mid-depth is

```math
T_f(t) = T_0 + q(t)\,R_b^* + \frac{(q \star g)(t)}{2\pi k_s},
```

where

- ``T_0`` is the undisturbed ground temperature,
- ``q(t)`` is the heat load per unit borehole length [W/m],
- ``R_b^*`` is the effective borehole resistance from `BoreholeResistance.jl`,
- ``g`` is the ground thermal response from `GroundResponse.jl`, and
- ``q \star g`` is the **temporal superposition** (discrete convolution) of the load with the
  response, evaluated by [`fluid_temperature`](@ref).

The middle term is the *instantaneous* drop across the borehole; the last term is the *cumulative*
ground temperature rise built up by the entire load history. Their split is what makes the model
both fast and accurate: the resistance is a steady-state constant, while all the transient memory
lives in the convolution.

## Outlet and inlet temperatures

The mean fluid temperature ``T_f`` is the average of the down-going and up-coming legs. With a load
``Q`` [W] carried by a flow ``V`` [m³/s] of fluid with volumetric heat capacity ``C_f``, the two
ends split symmetrically about the mean:

```math
T_{out} = T_f - \frac{Q}{2\,V\,C_f}, \qquad
T_{in}  = T_f + \frac{Q}{2\,V\,C_f}.
```

These are [`outlet_temperature`](@ref) and [`inlet_temperature`](@ref); see
[Fluid temperature](@ref).

## Where each step is documented

| Step | Function(s) | Page |
|---|---|---|
| Borehole resistance ``R_b^*`` | `resistance_ULoop_effective` (re-exported) | [BoreholeResistance docs](https://GHE-jl.github.io/BoreholeResistance.jl) |
| Ground response ``g(t)`` | `ground_response`, `FLSModel`, … | [g-function compression](@ref) + [GroundResponse docs](https://GHE-jl.github.io/GroundResponse.jl) |
| Convolution ``q \star g`` | [`convolution`](@ref), [`convolution_ns`](@ref) | [Temporal superposition](@ref) |
| Fluid temperatures | [`fluid_temperature`](@ref), [`outlet_temperature`](@ref), [`inlet_temperature`](@ref) | [Fluid temperature](@ref) |

## A note on g-function compression

A one-year hourly simulation has 8760 time steps; a multi-decade one has many times more. Evaluating
the analytical g-function at *every* step is the bottleneck. `ground_response` (re-exported from
`GroundResponse.jl`) evaluates the g-function only at a logarithmically spaced subset of time nodes
and reconstructs the full vector by PCHIP interpolation — transparently, inside the
`fluid_temperature` overloads that take a ground model. See [g-function compression](@ref).
