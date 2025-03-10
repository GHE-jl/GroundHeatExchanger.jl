"""
Script testing the nodes generation to apply to analytical models.
"""

using Plots
using BenchmarkTools
includet("../src/GHEModels.jl")
using .GHEModels

@btime id = set_nodes(1000, 20)