# Validation script containig the validations for the functions of GroundHeatExchanger.jl

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

###################################################################################################
# Reynolds Number - Chap2: p.20, ex. 2.1 
# Finding the ReynoldsNumber of water at 300k using the Reynoldsfunction
# Accuracy : 99.97%
###################################################################################################

#includet("../src/borehole_thermal_resistance/resistance_fluid.jl")

using BoreholeResistance
# References value of the ReynoldsNumber
Re_ref = 12428

# Reynolds Number calculated with the function
ρf = 997
μf = 8.54e-4
r = 0.01
Ac = π/4*(0.04^2 - 0.02^2)
V̇ = 0.5/(ρf*Ac)

Re = Reynolds(V̇, r, ρf, μf)

# accuracy verification
accuracy(Re_ref, Re) # 99.97%

###################################################################################################
# CoolProp.jl - Chap3: p.45, ex. 3.1 
# Finding the properties of the R134 fluid of the four different points of the ideal cycle
# Accuracy : 99.99%
###################################################################################################

using CoolProp

# References values of the temperature and the entropy
T1_ref = 260.44 # [°K]
T2_ref = 319.35 # [°K]
T3_ref = 312.54 # [°K]
T4_ref = 260.44 # [°K]
S1_ref = 1735 # [J/kg-K]
S2_ref = 1735 # [J/kg-K]

# Temperature and entropy calculated with the CoolProp package
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

###################################################################################################
# ILS and ICS - Chap4: p.60, ex. 4.1 
# Compute the g function
# Accuracy : 40 % (with ils function)
#            99.92% (with ils_test_lamarche)
###################################################################################################

using GroundResponse
using SpecialFunctions # to remove when GroundResponse/src/infinite_line_source will be modified

# Reference value of the g function
g_ils_ref = 0.2132

# g function calculated with the ils function from GroundResponse

# variables
α = 0.1/(24*3600)
t = 250*24*3600
r = 2
ks = 2.5
Cs = ks/α

# g function calculated with the ils_test_lamarche function
g_ils = ils(t, r, ks, Cs)
# accuracy verification
accuracy(g_ils_ref, g_ils) # 40%

# test of a different ILS function computed with the Lamarche method
function ils_test_lamarche(t, r, ks, Cs) 
    α = ks / Cs
    x = -r^2 / (4 * α * t)
    return -expinti(x) / (4*π)
end

g_ils_test = ils_test_lamarche(t, r, ks, Cs)

# accuracy verification
accuracy(g_ils_ref, g_ils_test ) # 99.92%


###################################################################################################
# ILS, ICS and FLS - Chap4: p.66, ex. 4.2 
# Valider les résultats pour la ICS, ILS et FLS. L'objectif serait de reproduire la figure à la fin de l'exercice. Pour ça, interpoler 10 points par heures et utiliser la fonction de "convolution".
###################################################################################################


###################################################################################################
# Resistance Borehole - Chap5: p.75, ex. 5.1 
# Démontrer les réponses pour les sous questions (c) et (d). La fonction à utiliser est resistance_borehole_multipole avec nLoop = 1 et order = 0 et 1.
###################################################################################################



###################################################################################################
# Line-Source Method - Chap5: p.76, ex. 5.2
# Démontrer les résultats avec la "line-source method" (mêmes paramètres qu'en 5.1, mais en gardant order = 0). pour les cas A, B, et C.
###################################################################################################


###################################################################################################
# Multipole - Chap5: p.77, ex. 5.3
# Obtenir les résultats pour la line source et la first-order multipole methods. Cette fois, la méthode est dans la fonction resistance_total_internal_multipole avec nLoop = 1 et order = 0 et 1.
###################################################################################################



###################################################################################################
# Double U-loop chap5: p.80, ex. 5.4
# Obtenir les résultats demandés pour un double U-loop (c'est le paramètre network dans la méthode resistance_total_internal_multipole qui permet de jouer entre les cas 13 et 24).
###################################################################################################



###################################################################################################
# Linear approximation chap5: p.84, ex. 5.5
# Je l'ai regardé rapidement et on devrait pouvoir avoir les résultats pour la linear approximation.
###################################################################################################
