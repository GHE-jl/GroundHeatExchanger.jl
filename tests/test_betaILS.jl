"""
Script testing the enhanced \beta-ILS model from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using CairoMakie
includet("../src/GHEModels.jl")
using .GHEModels

