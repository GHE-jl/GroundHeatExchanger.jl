# Validation script containig the validations for the functions of GroundHeatExchanger.jl

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per project, to install dependencies

using Revise

# TODO : add multiple dispatch to the accuracy function (ex: array)
function accuracy(val_ref::Real, val_GHE::Real)
    """
    accuracy(val_ref, val_GHE)
Compute the ratio between the two given values, to estimate if the value computed with the
function is accurate.
# Arguments
    - `val_ref`: value of reference obtained from Lamarche's Book
    - `val_GHE`: value computed using the GroundHeatExchanger package
# Output
    - `accuracy` : percent of accuracy of the value computed using the GHE's function
    """
    return min(val_ref, val_GHE)/max(val_ref, val_GHE) *100   
end

########## Reynolds Number - Chap2: p.20, ex.2.1 ##########
# Finding the Reynold Number of water at 300k

includet("../src/borehole_thermal_resistance/resistance_fluid.jl")

# References value of the Reynold Number
Re_ref = 12428

# Reynold Number calculated with the function
ρf = 997
μf = 8.54e-4
r = 0.01
Ac = π/4*(0.04^2 - 0.02^2)
V̇ = 0.5/(ρf*Ac)

Re = Reynold(V̇, r, ρf, μf)

# accuracy verification
accuracy(Re_ref, Re) # 99.97%


########## CoolProp.jl - Chap3: p.45, ex.3.1 ##########
# Finding the properties of the R134 fluid of the four different points of the ideal cycle

# References values of the temperature and the entropy
T1_ref = 260.44 # [°K]
T2_ref = 319.35 # [°K]
T3_ref = 312.54 # [°K]
T4_ref = 260.44 # [°K]
S1_ref = 1735 # [J/kg-K]
S2_ref = 1735 # [J/kg-K]

# Temperature and entropy calculated with the CoolProp package
using CoolProp
T1 = PropsSI("T", "P", 1.8e5, "H", 391020, "R134a")
T2 = PropsSI("T", "P", 1e6, "H", 426760, "R134a")
T3 = PropsSI("T", "P", 1e6, "H", 255500, "R134a")
T4 = PropsSI("T", "P", 1.8e5, "H", 255500, "R134a")
S1 = PropsSI("S", "P", 1.8e5, "H", 391020, "R134a")
S2 = PropsSI("S", "P", 1e6, "H", 426760, "R134a")

# accuracy verification
accuracy(T1_ref, T1) # 99.99 %
accuracy(T2_ref, T2) # 99.99 %
accuracy(T3_ref, T3) # 99.99 %
accuracy(T4_ref, T4) # 99.99 %
accuracy(S1_ref, S1) # 99.98 %
accuracy(S2_ref, S2) # 99.98 %

########## ILS and ICS - Chap4: p ex.4.1 ##########


