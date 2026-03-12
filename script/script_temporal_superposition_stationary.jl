"""
Script showcasing the stationary convolution from Convolutions.jl.
"""

using Random
using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = 60.0:60:3600*24*6

# Define the transfer function (here, a g-function for impulse in W/m from a finite line source)
g = fls(t, ks, Cs, rb, H, D)

# Define the heat load profile with added noise
q = [12000.0 * ones(60*24*2); 8000.0 * ones(60*24*2); 10000.0 * ones(60*24*2)] / H
rnd = rand(length(t))
q_rnd = q .+ rnd

# Convolution to have a temperature vector
ΔT = convolution(q_rnd, g)

# Plot the results
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = L"$q$ (W/m)")
lines!(ax, t, q_rnd)
ax = Axis(fig[1, 2], xlabel = L"$t$ (s)", ylabel = L"$ΔT$ (°C)")
lines!(ax, t, ΔT)
display(fig)