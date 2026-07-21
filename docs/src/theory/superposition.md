# Temporal superposition

A real ground heat exchanger carries a load that changes every hour. **Temporal superposition** is
the principle that the temperature response to a varying load equals the superposition of the
responses to a sequence of step changes — each step convolved with the unit g-function. This page
covers the FFT-based convolution that evaluates it, in both its stationary and non-stationary forms.

## From load steps to a convolution

Decompose the load history into increments. The **incremental load function**

```math
f_1 = q_1, \qquad f_k = q_k - q_{k-1} \quad (k \ge 2)
```

is the per-step change in load ([`impulse_func`](@ref)). The temperature response is then the
discrete convolution of these increments with the g-function:

```math
(q \star g)_n = \sum_{k=1}^{n} f_k \, g_{\,n-k+1}.
```

This is a *Duhamel* superposition: each load step starts a fresh g-function response, and the
responses add. The package evaluates it in the spectral domain via the fast Fourier transform
([`DSP.jl`](https://github.com/JuliaDSP/DSP.jl)), giving ``O(n\log n)`` cost instead of the
``O(n^2)`` of the direct sum. Zero-padding is used so the FFT computes a *linear*, not circular,
convolution. Time steps are assumed equally spaced.

## Stationary convolution

When the operating conditions (flow rate, fluid properties) are constant, a single g-function
applies for the whole simulation. [`convolution`](@ref) takes the raw load vector `q`, forms the
increments internally, and returns the convolved response:

```julia
ΔT = convolution(q, g)          # q in W/m (or W, °C); g the unit response
```

If the same increments are reused against several g-functions, [`convolutionf`](@ref) skips the
incremental step and takes a pre-computed `f` from [`impulse_func`](@ref) directly.

This is exactly the ``q \star g`` term in the master equation (see [Simulation pipeline](@ref)); it
is what [`fluid_temperature`](@ref) calls internally.

## Non-stationary convolution

When the operating conditions **change during the simulation** — a variable-speed circulation pump,
or temperature-dependent fluid properties — a single g-function no longer applies: each operating
*state* has its own response. [`convolution_ns`](@ref) (Beaudry et al., 2024) handles this by
carrying **one g-function column per state**.

The bookkeeping is done by three helpers:

- [`state_indices`](@ref) — given one or more parameter vectors (flow, temperature, …), find where
  the operating state changes and label each segment;
- [`state_vector`](@ref) — expand those change indices into a per-time-step state index;
- [`impulse_func_ns`](@ref) — build the per-state incremental load matrix by masking the load with
  each state in turn.

[`convolution_ns!`](@ref) is the in-place core: it convolves each state's incremental column with the
matching g-function column and accumulates the result. A typical call:

```julia
# One g-function column per operating state
ks_states = [3.0, 2.0, 1.5]
g_matrix  = hcat([fls(t, rb, H, D, ks_i, Cs) for ks_i in ks_states]...)

# Map each time step to a state, then superpose
ind, state, _ = state_indices(operating_condition_vec)
state_vec     = state_vector(ind, state, length(t))
ΔT            = convolution_ns(q, g_matrix, state_vec)
```

The result is the same temperature response as the stationary `convolution`, but with the correct
g-function applied over each interval of constant operating conditions.

## Functions on this page

```@docs
convolution
convolutionf
impulse_func
convolution_ns!
convolution_ns
impulse_func_ns
state_vector
state_indices
```
