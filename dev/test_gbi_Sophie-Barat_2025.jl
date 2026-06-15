"""
Test for generation of transfer functions for Sophie-Barat 2025.
"""

using CairoMakie
using Revise
includet("../src/GroundHeatExchanger.jl")
using .GroundHeatExchanger

# Define the parameters for the test
t = (1:3600:3600*24*365*10) # Time vector (1 hour steps for 10 years)
H = 537-5.29 # Borehole depth (m)
rb = 0.0895 # Borehole radius (m)
D = 5.29 # Borehole buried depth (m)
ks = 3.9 # Soil thermal conductivity (W/mK)
Cs = 2.7149e6 # Soil volumetric specific heat (J/m³K)
XY = [0 0;
5.6 -12.5;
5.6 -22.5;
5.0 -35.5;
-24.1 -12.2;
-34.1 -12.2;
-20.4 -26.2;
-30.5 -26.1;
-29.0 -36.0;
-57.3 -15.5;
-67.3 -15.5;
-77.3 -15.5]

# Finite line source
g = successive_flux(t, H, rb, D, ks, Cs, XY)

plot(t, g, label="g-function")