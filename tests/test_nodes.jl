"""
Script testing the nodes generation to apply to analytical models.
"""

includet("../src/GHEModels.jl")

using Plots
using BenchmarkTools
using .GHEModels

@btime id = set_nodes(1000, 20)