# GroundHeatExchanger.jl

*Full thermal simulation of ground heat exchangers, in pure Julia.*

`GroundHeatExchanger.jl` is the **integration layer** of a three-package geothermal stack. It takes
the borehole thermal resistance from
[`BoreholeResistance.jl`](https://github.com/GHE-jl/BoreholeResistance.jl) and the ground
*g*-function from
[`GroundResponse.jl`](https://github.com/GHE-jl/GroundResponse.jl), and assembles them into a
time-domain simulation of the **mean fluid, inlet and outlet temperatures** of a ground heat
exchanger under a time-varying load.

```
 BoreholeResistance.jl                 GroundResponse.jl
   water_k / cp / ρ / μ                  ILSModel, FLSModel, …
   resistance_ULoop_effective        ground_response, successive_flux
              ↘                                  ↙
                   GroundHeatExchanger.jl
   temporal superposition (FFT) · fluid_temperature · outlet/inlet_temperature
   pchip_interpolation · ground_load_profile
```

A single `using GroundHeatExchanger` re-exports the whole ecosystem, so the entire pipeline — water
properties, resistances, ground models, superposition and temperatures — is available from one
import.

## What this package adds

On top of the two upstream packages, `GroundHeatExchanger.jl` provides:

- **Temporal superposition** — [`convolution`](@ref) and friends, an FFT-based ``O(n\log n)``
  convolution of the load history with the g-function, including a **non-stationary** variant
  ([`convolution_ns`](@ref)) for time-varying operating conditions.
- **Fluid temperatures** — [`fluid_temperature`](@ref), [`outlet_temperature`](@ref) and
  [`inlet_temperature`](@ref) from the load, resistance and g-function.
- **Short-term outlet ANN** — [`short_term_response`](@ref) / [`short_term_nodes`](@ref), the
  artificial neural network of Pasquier, Zarrella & Labib (2018) that emulates a combined
  borehole + ground short-term response, joined to any long-term ground model by
  [`outlet_transfer_function`](@ref).
- **g-function compression helper** — [`pchip_interpolation`](@ref) for reconstructing an expensive
  signal from a node subset with shape-preserving PCHIP interpolation (the same kind of compression
  `ground_response` applies internally via `interp`).
- **Utilities** — a synthetic annual [`ground_load_profile`](@ref) and a [`GHE`](@ref) parameter set
  for the examples.

## Installation

The package and its two siblings are not yet registered. The most reliable setup is to clone all
three side by side and develop them locally:

```julia
using Pkg
Pkg.develop(path = "../BoreholeResistance.jl")
Pkg.develop(path = "../GroundResponse.jl")
Pkg.develop(path = ".")           # GroundHeatExchanger.jl
Pkg.instantiate()
```

(`GroundHeatExchanger.jl` declares the two siblings as path `[sources]`, so developing it pulls in
the matching local versions.)

## Quick start

```julia
using GroundHeatExchanger

# Geometry and ground properties
H, D, rb, s, ro, ri = 150.0, 2.0, 0.08, 0.05, 0.022, 0.017
ks, Cs, kg, kp      = 3.0, 2.11e6, 1.6, 0.4
T0, V               = 10.0, 30/6e4        # 10 °C undisturbed, 30 L/min flow

# Time vector: 1 year, hourly
t = collect(3600.0:3600:3600*24*365)

# Fluid properties (re-exported from BoreholeResistance)
kf = water_k(T0);  cf = water_cp(T0);  ρf = water_ρ(T0);  μf = water_μ(T0)
Cf = cf * ρf                              # volumetric specific heat [J/m³·K] for Tin/Tout

# Effective borehole resistance (multipole method)
Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

# Ground model and heat load
model = FLSModel(H, D, ks, Cs)
Q     = ground_load_profile(t ./ 3600)      # [W], sinusoidal annual profile

# Full simulation — g-function computed and PCHIP-compressed automatically
Tf   = fluid_temperature(t, Q ./ H, model, rb, T0, ks, Rb)
Tout = outlet_temperature(Tf, Q, V, Cf)
Tin  = inlet_temperature(Tf, Q, V, Cf)
```

## Manual outline

- **[Tutorial](@ref)** — a complete one-year hourly simulation, built up step by step.
- **Modeling theory** — the equations behind the pipeline:
  - [Simulation pipeline](@ref) — how the three packages fit together and the master equation.
  - [Temporal superposition](@ref) — the FFT convolution, stationary and non-stationary.
  - [Fluid temperature](@ref) — mean, outlet and inlet temperatures.
  - [g-function compression](@ref) — the PCHIP node-subset reconstruction.
- **[Utilities](@ref)** — the load profile, head-loss and parameter helpers.
- **[API reference](@ref)** — the complete docstring reference for every symbol defined here.
- **[References](@ref)** — the bibliography underpinning the methods.

## Conventions used throughout

| Symbol | Meaning | Unit |
|---|---|---|
| ``T_f`` | Mean fluid temperature | °C |
| ``T_{in}, T_{out}`` | Inlet / outlet fluid temperature | °C |
| ``T_0`` | Undisturbed ground temperature | °C |
| ``q`` | Heat load per unit borehole length | W/m |
| ``Q`` | Total heat load | W |
| ``g`` | Ground thermal response | °C·m/W |
| ``R_b^*`` | Effective borehole resistance | m·K/W |
| ``k_s`` | Ground thermal conductivity | W/m·K |
| ``V`` | Volumetric flow rate per borehole | m³/s |
| ``C_f = c_f\,\rho_f`` | Fluid **volumetric** specific heat | J/m³·K |

!!! note "Two heat-capacity conventions"
    `resistance_ULoop_effective` takes the **mass-specific** heat ``c_f`` [J/kg·K] and
    the density ``\rho_f`` separately, while [`outlet_temperature`](@ref) and
    [`inlet_temperature`](@ref) take the **volumetric** specific heat
    ``C_f = c_f\,\rho_f`` [J/m³·K]. Compute `Cf = water_cp(T) * water_ρ(T)` explicitly when crossing
    that boundary.

## Re-exported APIs

`GroundHeatExchanger.jl` re-exports every public symbol of its two upstream packages. Their full
documentation lives in their own sites:

| From | Re-exported symbols | Documentation |
|---|---|---|
| **BoreholeResistance.jl** | `water_k`/`cp`/`ρ`/`μ`, `Reynolds`, `Prandtl`, `Nusselt`, `resistance_*` | [docs](https://GHE-jl.github.io/BoreholeResistance.jl) |
| **GroundResponse.jl** | `ILSModel` … `MFLSModel`, `ils` … `mfls`, `ground_response`, `successive_flux`, `bloc_matrix`, `borefield*` | [docs](https://GHE-jl.github.io/GroundResponse.jl) |

The pages here document only the symbols **defined in** `GroundHeatExchanger.jl`. The g-function
interface `ground_response` is re-exported from `GroundResponse.jl` (its PCHIP compression is
toggled by the `interp` keyword) and documented in that package's site.
