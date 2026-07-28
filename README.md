# GroundHeatExchanger.jl

Julia package for the complete thermal simulation of ground heat exchangers (GHEs). Acts as
the integration layer for a three-package ecosystem:

```
BoreholeResistance.jl          GroundResponse.jl
  water_k/cp/ρ/μ                ILSModel, FLSModel, ...
  resistance_ULoop_effective  ground_response, successive_flux
           ↘                            ↙
           GroundHeatExchanger.jl
           convolution (FFT), fluid_temperature, outlet_temperature, inlet_temperature
           ground_response (with PCHIP), ground_load_profile
```

A single `using GroundHeatExchanger` exposes the full pipeline.

## Quick start

```julia
using GroundHeatExchanger

# Geometry and ground properties
H, D, rb, s, ro, ri = 150.0, 2.0, 0.08, 0.05, 0.022, 0.017
ks, Cs, kg, kp      = 3.0, 2.11e6, 1.6, 0.4
T0, V               = 10.0, 30/6e4    # 10°C undisturbed, 30 L/min flow

# Time vector: 1 year, hourly
t = collect(3600.0:3600:3600*24*365)

# Fluid properties (from BoreholeResistance)
kf = water_k(T0);  cf = water_cp(T0);  ρf = water_ρ(T0);  μf = water_μ(T0)
Cf = cf * ρf       # volumetric specific heat [J/m³·K] for Tin/Tout

# Borehole thermal resistance (Rb* from multipole method)
Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

# Ground model and heat load
model = FLSModel(H, D, ks, Cs)
Q = ground_load_profile(t ./ 3600)      # [W], sinusoidal annual profile

# Full simulation — g-function computed and PCHIP-compressed automatically
Tf   = fluid_temperature(t, Q ./ H, model, rb, T0, Rb)
Tout = outlet_temperature(Tf, Q, V, Cf)
Tin  = inlet_temperature(Tf, Q, V, Cf)
```

## Temperature simulation

The mean fluid temperature follows:
```
Tf(t) = T0 + q(t)·Rb + (f★g)(t)
```
where `f★g` is the discrete convolution (temporal superposition) solved in O(n log n) via FFT.

Outlet and inlet split symmetrically around the mean:
```
Tout = Tf − Q / (2·V·Cf)
Tin  = Tf + Q / (2·V·Cf)
```

### Units

| Quantity | Symbol | Unit |
|---|---|---|
| Borehole resistance | `Rb` | m·K/W |
| g-function | `g` | °C·m/W |
| Fluid specific heat | `cf` | J/kg·K — for `resistance_ULoop_effective` |
| Fluid volumetric heat | `Cf = cf·ρf` | J/m³·K — for `outlet_temperature`, `inlet_temperature` |

### `fluid_temperature` overloads

```julia
# With pre-computed g-function vector
fluid_temperature(t, q, g, T0, Rb)              # q in W/m

# With ground model (g computed + PCHIP-compressed automatically)
fluid_temperature(t, q, m, rb, T0, Rb; interp=true)       # single borehole
fluid_temperature(t, q, m, rb, xy, T0, Rb; interp=true)   # borefield (xy: nb×2 matrix)
```

The load is always the heat load **per unit length** `q` [W/m]. Convert a total load `Q` [W]
yourself by dividing by the **entire** borehole length `L = H · (number of boreholes)`:
`q = Q / L`. This is left to the caller on purpose — dividing `Q` by a single borehole depth `H`
instead of the full field length is a common mistake, so the `Q`/`H` overloads were removed.

## PCHIP g-function compression

For long hourly simulations (8760 steps/year × many years), evaluating g-functions at every
time step is expensive. `ground_response` (re-exported from GroundResponse.jl) evaluates g on a
sub-sampled set of logarithmically-spaced time nodes and uses PCHIP (Piecewise Cubic Hermite
Interpolating Polynomial) to reconstruct the full vector. This is toggled by the `interp` keyword:

```julia
g = ground_response(t, rb, xy, model)               # interp=true (default): PCHIP compression
g = ground_response(t, rb, xy, model; interp=false) # exact, no compression
```

The `fluid_temperature` overloads apply this automatically via the same `interp` kwarg. For
reconstructing your own signals from a node subset, GHE also exposes `pchip_interpolation(tᵢ, gᵢ, t)`.

## Temporal superposition

```julia
convolution(q, g)                        # stationary (single state), FFT O(n log n)
convolutionf(f, g)                       # stationary: pre-computed impulse function
convolution_ns(q, g_matrix, state_vec)   # non-stationary: varying flow / properties
```

The non-stationary path handles time-varying operating conditions (e.g., variable flow rate
or temperature-dependent fluid properties). Each state has its own column in the g-function
matrix `g_matrix` (nt × n_states), and `state_vec` maps each time step to a column:

```julia
# Build per-state g-function matrix
ks_states = [3.0, 2.0, 1.5]
g_matrix = hcat([fls(t, rb, H, D, ks_i, Cs) for ks_i in ks_states]...)

# State assignment and non-stationary temperature response
ind, state, _ = state_indices(operating_condition_vec)
state_vec = state_vector(ind, state, length(t))
ΔT = convolution_ns(q, g_matrix, state_vec)
```

## Ground models (from GroundResponse.jl)

| Struct | Parameters | Reference |
|---|---|---|
| `ILSModel(ks, Cs)` | Ground conductivity, heat capacity | Ingersol (1948) |
| `ICSModel(rc, ks, Cs)` | Cylinder radius, ground properties | Carslaw & Jaeger (1959) |
| `FLSModel(H, D, ks, Cs)` | Borehole depth and burial depth | Claesson & Javed (2011) |
| `MILSModel(rb, ks, Cs, Cf, vD)` | Borehole radius, ground + fluid properties, Darcy velocity | Pasquier & Lamarche (2022) |
| `MFLSModel(H, rb, D, ks, Cs, Cf, vD)` | Full FLS + groundwater flow | Guo et al. (2020) |

All models are dispatched through `ground_response(t, rb, xy, m::AbstractGroundModel)`:
- Single borehole (`size(xy,1) == 1`): direct g-function evaluation at the borehole wall.
- Borefield: `successive_flux` spatial superposition is applied automatically.

```julia
# Borefield spatial superposition
xy = borefield(:rectangle, 3, 4, 6.0)    # 3×4 rectangular grid, 6 m spacing
g  = ground_response(t, rb, xy, model)   # successive_flux applied internally
```

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`.
First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```
Note: the temperature simulation script filename contains a space — quote it on the command line.

| Script | Description |
|---|---|
| `script_temperature simulation.jl` | Full pipeline: FLS g → Rb → Tf/Tin/Tout, 1-year hourly |
| `script_temporal_superposition_stationary.jl` | FFT convolution, FLS g-function, 6-day load |
| `script_temporal_superposition_nonstationary.jl` | Non-stationary convolution, 3 operating states |
| `script_outlet_transfer_function.jl` | Short- + long-term outlet (EWT) transfer function: Pasquier et al. (2018) Figs. 2a/3 (`PublishedANN`), plus a `PublishedANN` vs `DeepANN` borefield comparison |
| `script_ann_validation.jl` | Short-term `PublishedANN` and `DeepANN` (5 cases each) plotted against their MATLAB references (`ANN_gfunction.m` / `SULoop_TRCMz_ANN.m`) at native nodes, with the interpolated `short_term_response` overlaid |
| `script_validation_Lamarche2023.jl` | `@test`-based validation of Reynolds number, fluid properties, ILS/ICS/FLS g-functions and borehole/internal resistance against worked examples from Lamarche's textbook |

## Installation

These packages are not yet registered. Install as local dev dependencies from the parent folder:

```julia
using Pkg
Pkg.develop(path="../BoreholeResistance.jl")
Pkg.develop(path="../GroundResponse.jl")
Pkg.instantiate()
```

## Dependencies

### Library

| Package | Used for |
|---------|----------|
| [BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl) | Borehole and fluid thermal resistances, water properties |
| [GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) | Ground thermal response models and borefield layouts |
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl) | `conv` / `conv!` for FFT-based temporal superposition |
| [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) | PCHIP g-function compression |

### Scripts only

| Package | Used in |
|---------|---------|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |

## References

- Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for GHE systems.
  Geothermics, 37(6), 651–665.
- Beaudry, G., Pasquier, P., & Nguyen, A. (2024). Non-stationary convolutions for time-variant
  flowrates in GHEs. Science and Technology for the Built Environment, 30(3), 208–219.
- Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function
  construction. Renewable Energy, 165, 359–368.
- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for g-function
  construction. Renewable Energy, 121, 249–260.
- Claesson, J., & Javed, S. (2011). An analytical method for borehole fluid temperatures.
  ASHRAE Transactions, 117(PART 2), 279–288.
- Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance methods. Applied Energy, 187.
- Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of ANNs to short-term g-functions
  (`PublishedANN`). Applied Thermal Engineering, 143, 910–921.
- Pasquier, P., & Marcotte, D. (2020). Robust identification of volumetric heat capacity and
  analysis of thermal response tests by Bayesian inference with correlated residuals (`DeepANN`,
  the default short-term model). Applied Energy, 261, 114394.
