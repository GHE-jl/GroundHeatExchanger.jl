"""
Collection of thermal transfer models to simulate a ground heat exchanger (GHE). The models
available are:
    - Infinite line source of Ingersol (1948) (ils)
    - Infinite cylindrical source of Ingersol (1959) (ics)
    - Finite line source of Claesson and Javed (2011) (fls)
    - The standing column well analytical model of Nguyen et al. (2025) (ilsβ)
"""

# Include all ground model files
include("GroundModels/ILS.jl")
include("GroundModels/ICS.jl")
include("GroundModels/FLS.jl")
include("GroundModels/MFLS.jl")
include("GroundModels/betaILS.jl")