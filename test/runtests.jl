"""
    runtests.jl

Main test entry point for the GroundHeatExchanger package.
This file runs all test modules for the package.
"""

using Test
using GroundHeatExchanger

# Include test modules
include("test_finite_line_source.jl")
include("test_infinite_line_source.jl")
include("test_infinite_cylindrical_source.jl")
include("test_moving_finite_line_source.jl")
include("test_moving_infinite_line_source.jl")
