"""
Script showcasing the workflow for the whole GHE simulation
"""

using BenchmarkTools
using CairoMakie

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define paremeters
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()

#TODO provide a realistic workflow to simulate the whole GHE system, including the computation of
# the g-function, spatial superposition, the convolution and the temperature response.