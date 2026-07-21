# g-function compression

Evaluating the analytical g-function (an integral or series, possibly with spatial superposition
over a field) at *every* time step is the dominant cost of a long hourly simulation. Because the
g-function is **smooth and slowly varying on a logarithmic time axis**, it can be evaluated at a
small subset of nodes and reconstructed almost exactly by interpolation. This page documents where
that compression happens and the interpolation helpers `GroundHeatExchanger.jl` exposes for reuse.

## Where compression happens

Compression lives in `GroundResponse.ground_response` (re-exported by `GroundHeatExchanger.jl`,
so `using GroundHeatExchanger` brings it into scope). It is toggled by the `interp` keyword: with
`interp=true` the g-function is evaluated on an internal constant-step sub-sampling grid and
PCHIP-interpolated to `t`; with `interp=false` it is evaluated exactly at every `t`.

```julia
g = ground_response(t, rb, xy, model)                # interp=true (default): compressed
g = ground_response(t, rb, xy, model; interp=false)  # exact — no compression
```

For the temporal solvers on a non-uniform `t`, `interp=true` is a *correctness* requirement (their
convolution assumes a constant step) as well as a performance win; `interp=false` is valid only for
a uniformly spaced `t`. See the
[GroundResponse.jl docs](https://GHE-jl.github.io/GroundResponse.jl) for the full details. The same
`interp` keyword flows through the [`fluid_temperature`](@ref) overloads that take a ground model,
so compression is applied transparently in a full simulation.

## Why PCHIP

The reconstruction uses a **Piecewise Cubic Hermite Interpolating Polynomial**
([`PCHIPInterpolation.jl`](https://github.com/gerlero/PCHIPInterpolation.jl)). PCHIP is
**shape-preserving**: it does not overshoot or introduce spurious oscillations between nodes, which
matters because a g-function is monotone and an overshoot would inject unphysical wiggles into the
temperature response. On a logarithmic time grid a few hundred nodes reproduce the full g-function to
well within other modelling errors, at a fraction of the evaluation cost.

## Custom interpolation

[`pchip_interpolation`](@ref) is exposed directly as a reusable `GroundHeatExchanger.jl` helper —
for example, to reconstruct any other slowly-varying signal from a node subset:

```julia
id = round.(Int, exp10.(range(0, log10(length(t)); length=150)))  # log-spaced node indices
gᵢ = expensive_function(t[id])         # evaluate only at the nodes
g  = pchip_interpolation(t[id], gᵢ, t) # reconstruct the full vector
```

## Functions on this page

```@docs
pchip_interpolation
```
