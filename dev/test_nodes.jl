"""
Script testing the nodes generation to apply to analytical models.
"""

includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

@btime id = set_nodes(1000, 20)