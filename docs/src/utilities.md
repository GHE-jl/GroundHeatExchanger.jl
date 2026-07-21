# Utilities

Beyond the simulation core, `GroundHeatExchanger.jl` ships a few helpers used throughout the
examples.

## Synthetic load profile

[`ground_load_profile`](@ref) generates a realistic **annual heat-load signal** for testing the
convolution — a smooth seasonal sine modulated by daily harmonics and a building on/off pattern,
following Eqs. 7–8 of Bernier et al. (2004). It takes time **in hours** and returns the load in
watts:

```julia
t = collect(3600.0:3600:3600*24*365)   # 1 year, hourly [s]
Q = ground_load_profile(t ./ 3600)        # [W]
```

All seven shape parameters (amplitude, phase, harmonics, …) have sensible defaults; override them to
change the amplitude or the heating/cooling balance.

## Example parameter set

[`GHE`](@ref) returns a complete, named tuple of typical geometry, ground, grout, pipe and fluid
properties used by the package's validation scripts. It is a convenience for reproducing the
examples without re-typing two dozen parameters:

```julia
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, vD, V = GHE()
```

## Functions on this page

```@docs
ground_load_profile
GHE
```
