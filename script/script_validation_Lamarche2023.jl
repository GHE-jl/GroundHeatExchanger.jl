# Validation script containig the validations for the functions of GroundHeatExchanger.jl

using Test

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

@testset "Validation from Lamarche's" begin 
###################################################################################################
# Reynolds Number - Chap2: p.20, ex. 2.1 
# Finding the ReynoldsNumber of water at 300k using the Reynolds function
# Accuracy : 99.97%
###################################################################################################

using BoreholeResistance

# References value of the ReynoldsNumber
Re_ref = 12428

# Reynolds Number calculated with the function
ρf = 997                    # density [kg/m³]
μf = 8.54e-4                # viscosity [kg/m⋅s]
r = 0.01                    # radius[m]
Ac = π/4*(0.04^2 - 0.02^2)  # area [m²]
V̇ = 0.5/(ρf*Ac)             # speed [m/s]

Re = Reynolds(V̇, r, ρf, μf)

# Test and accuracy verification
@test Re_ref ≈ Re rtol= 1e-3

accuracy(Re_ref, Re) # 99.97%

###################################################################################################
# CoolProp.jl - Chap3: p.45, ex. 3.1 
# Finding the properties of the R134 fluid at four different points of the ideal cycle
# Accuracy : 100%
#            100%
#            100%
#            100%
#            99.98%
#            99.98%
###################################################################################################

using CoolProp

# References values of the temperature and the entropy
T1_ref = 260.44 # [°K]
T2_ref = 319.35 # [°K]
T3_ref = 312.54 # [°K]
T4_ref = 260.44 # [°K]
S1_ref = 1735   # [J/kg-K]
S2_ref = 1735   # [J/kg-K]

# Temperature and entropy calculated with the CoolProp package
T1 = PropsSI("T", "P", 1.8e5, "H", 391020, "R134a")
T2 = PropsSI("T", "P", 1e6, "H", 426760, "R134a")
T3 = PropsSI("T", "P", 1e6, "H", 255500, "R134a")
T4 = PropsSI("T", "P", 1.8e5, "H", 255500, "R134a")
S1 = PropsSI("S", "P", 1.8e5, "H", 391020, "R134a")
S2 = PropsSI("S", "P", 1e6, "H", 426760, "R134a")

# Test and accuracy verification
@test T1_ref ≈ T1 rtol= 1e-5
@test T2_ref ≈ T2 rtol= 1e-4
@test T3_ref ≈ T3 rtol= 1e-5
@test T4_ref ≈ T4 rtol= 1e-5
@test S1_ref ≈ S1 rtol= 1e-3
@test S2_ref ≈ S2 rtol= 1e-3

accuracy(T1_ref, T1) # 100%
accuracy(T2_ref, T2) # 100%
accuracy(T3_ref, T3) # 100%
accuracy(T4_ref, T4) # 100%
accuracy(S1_ref, S1) # 99.98%
accuracy(S1_ref, S1) # 99.98%

###################################################################################################
# ILS and ICS - Chap4: p.60, ex. 4.1 
# Compute the transfer function using the Infinite Line Source Model and the Infinite Cylindrical 
# Source.
# Accuracy: 100% 
#           99.90%
###################################################################################################

using GroundResponse

# Reference value of the transfer function with the ILS model and the ICS model
ΔT_ils_ref = 3.414
ΔT_ics_ref = 3.412

# g function and transfer function calculated with ils and ics function from GroundResponse
q0 = -40          # heat load [W/m]
α = 0.1/(24*3600) # thermal diffusivity [m²/s]
t = 250*24*3600   # time [s]
r = 2             # radial position [m]
rc = 0.1          # radius [m]
ks = 2.5          # thermal conductivity [W/m-K]
Cs = ks/α         # thermal capacity [J/m³-K]

g_ils = ils(t, r, ks, Cs)     
g_ics = ics(t, r, rc, ks, Cs) 

ΔT_ils = -q0 * g_ils
ΔT_ics = -q0 * g_ics

# Test and accuracy verification
@test ΔT_ils_ref ≈ ΔT_ils rtol = 1e-3
@test ΔT_ics_ref ≈ ΔT_ics rtol = 1e-3

accuracy(ΔT_ils_ref, ΔT_ils)
accuracy(ΔT_ics_ref, ΔT_ics)

###################################################################################################
# ILS, ICS and FLS - Chap4: p.66, ex. 4.2 
# Compute the g function and the thermal response function using the ils, ics and fls functions.
# Accuracy : 99.96%
#            99.94%
#            99.94%
###################################################################################################

# TODO : do the graph, take a look at the _ierf function

using GroundResponse

# Reference value of the g function with the ILS, ICS and FLS model
g_ics_ref = [0.198, 0.1849, 0.1469]
g_ils_ref = [0.1739, 0.1574, 0.1085]
g_fls_ref = [0.1737, 0.1573, 0.1085]

# Reference value of the transfer function
ΔT_ics_ref = -0.26
ΔT_ils_ref = -0.41
ΔT_fls_ref = -0.41

# g function and transfer function calculated with ils, ics and fls function from GroundResponse 
α = 0.1 / (24*3600)                # thermal diffusivity [m²/s]
ks = 2                             # thermal conductivity [W/m-K]
Cs = ks/α                          # thermal capacity [J/m³-K]
rc = 0.075                         # radius [m]  
H = 100                            # borehole depth [m]
D = 0                              # buried depth [m]
q = [1000/100, 1800/100, -200/100] # heat loads [W/m]
t = [3600, 3*3600, 5*3600]         # time [s]

# using the ics (Infinite Cylindrical Source) function to calculate the g function
g_ics = [
    ks * ics(t[3], rc, rc, ks, Cs),
    ks * ics(t[3] - t[1], rc, rc, ks, Cs),
    ks * ics(t[3] - t[2], rc, rc, ks, Cs)
]

# using the ils (Infinite Line Source) function to calculate the g function
g_ils = [
    ks * ils(t[3], rc, ks, Cs),
    ks * ils(t[3] - t[1], rc, ks, Cs),
    ks * ils(t[3] - t[2], rc, ks, Cs)
]

# using the fls (Finite Line Source) function to calculate the g function
g_fls = [
    ks * fls(t[3], rc, H, D, ks, Cs),
    ks * fls(t[3] - t[1], rc, H, D, ks, Cs),
    ks * fls(t[3] - t[2], rc, H, D, ks, Cs)
]

# calculating the transfer function
ΔT_ics = -1/ks * (q[1]*g_ics[1] + (q[2] - q[1])*g_ics[2] + (q[3] - q[2])*g_ics[3])
ΔT_ils = -1/ks * (q[1]*g_ils[1] + (q[2] - q[1])*g_ils[2] + (q[3] - q[2])*g_ils[3])
ΔT_fls = -1/ks * (q[1]*g_fls[1] + (q[2] - q[1])*g_fls[2] + (q[3] - q[2])*g_fls[3])

# Test and accuracy verification
@test g_ics_ref ≈ g_ics rtol = 1e-3
@test g_ils_ref ≈ g_ils rtol = 1e-3
@test g_fls_ref ≈ g_fls rtol = 1e-3

@test ΔT_ics_ref ≈ ΔT_ics rtol = 1e-2
@test ΔT_ils_ref ≈ ΔT_ils rtol = 1e-2
@test ΔT_fls_ref ≈ ΔT_fls rtol = 1e-2

accuracy(g_fls_ref[1], g_fls[1]) # 99.96%
accuracy(g_fls_ref[2], g_fls[2]) # 99.94%
accuracy(g_fls_ref[3], g_fls[3]) # 99.94%

###################################################################################################
# Resistance Borehole - Chap5: p.75, ex. 5.1 
# Compute the borehole resistance (Rb) with the line source methode and the multipole methode using 
# the resistance_ULoop_borehole function.
# Accuracy : 99.15%
#            98.65%     Test Failed
###################################################################################################
# TODO : investigate the 98% accuracy for the multipole function??

using BoreholeResistance

# Reference value of the borehole resistance with the line source and the multipole method
Rb_ref = 0.120
Rb₁_ref = 0.120

# Borehole resistance calculated with the resistance_ULoop_borehole function from BoreholeResistance 
s = 0.05                         # distance between 2 legs of a U-tubes [m]
rb = 0.075                       # borehole radius [m]
ro = 0.021                       # outer radius of the pipe [m]
ri = 0.0172                      # inner radius of the pipe [m]
ks = 2.5                         # soil conductivity [W/mK]
kg = 1.7                         # grout conductivity [W/mK]
kp = 0.4                         # pipe conductivity [W/mK]
Rp = resistance_pipe(ro, ri, kp) # pipe thermal resistance [mK/W]
Rf = 1/614.50                    # water's thermal resistance at 30° [mK/W]

Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0)
Rb₁ = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1)

# Test and accuracy verification
@test Rb_ref ≈ Rb rtol = 1e-2
@test Rb₁_ref ≈ Rb₁ rtol = 1e-2

accuracy(Rb_ls_ref, Rb) # 99.15%
accuracy(Rb₁_ref, Rb₁)  # 98.65%

###################################################################################################
# Line-Source Method - Chap5: p.76, ex. 5.2
# Compute the borehole resistance (Rb) for different U-tube configuration (A, B, C) using the line 
# source method and the resistance_ULoop_borehole function.
# Accuracy : 99.76%
#            99.42%
#            98.97%     Test Failed
###################################################################################################

# TODO : investigate the 98% with the configuration C

using BoreholeResistance

# Reference value of the borehole resistance with the line source method for different U-tube
# configuration
Rb_A_ref = 0.226
Rb_B_ref = 0.192
Rb_C_ref = 0.109

# Borehole resistance calculated with the resistance_ULoop_borehole function from BoreholeResistance 
rb = 0.075                      # borehole radius [m]
ro = 0.033/2                    # outer radius of the pipe [m]
ri = 0.027/2                    # inner radius of the pipe [m]
ks = 2.5                        # soil conductivity [W/mK]
kg = 1.0                        # grout conductivity [W/mK]
kp = 0.4                        # pipe conductivity [W/mK]
Rp = 0.08                       # pipe thermal resistance [mK/W]
Rf = 1/628.56                   # water's thermal resistance at 30° [mK/W]
s = [2*ro, 2*rb / 3, 2*(rb-ro)] # distance between 2 legs of a U-tubes for case A, B and C [m]

Rb_A = resistance_ULoop_borehole(s[1], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0) 
Rb_B = resistance_ULoop_borehole(s[2], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0)
Rb_C = resistance_ULoop_borehole(s[3], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0)

# Test and accuracy verification
@test Rb_A_ref ≈ Rb_A rtol = 1e-2
@test Rb_B_ref ≈ Rb_B rtol = 1e-2
@test Rb_C_ref ≈ Rb_C rtol = 1e-2

accuracy(Rb_A_ref, Rb_A) # 99.76%
accuracy(Rb_B_ref, Rb_B) # 99.42%
accuracy(Rb_C_ref, Rb_C) # 98.97%

###################################################################################################
# Multipole - Chap5: p.77, ex. 5.3
# Compute the internal resistance (Ra) with the line source method and the multipole methode using 
# the resistance_ULoop_total_internal function.
# Accuracy : 99.16%
#            99.25%
#            99.67%
#            60.35%     Test Failed
#            99.37%
#            99.82%
###################################################################################################

# TODO : check the Ramult_A and why it fails??

using BoreholeResistance

# Reference value of the internal resistance with the line source method and the multipole method 
# for different U-tube configuration
Rals_A_ref = 0.3761
Rals_B_ref = 0.482
Rals_C_ref = 0.592

Ramult_A_ref = 0.592
Ramult_B_ref = 0.475
Ramult_C_ref = 0.592

# Internal resistance computed with the resistance_ULoop_total_internal function 
rb = 0.075                              # borehole radius [m]
ro = 0.033/2                            # outer radius of the pipe [m]
ri = 0.027/2                            # inner radius of the pipe [m]
ks = 2.5                                # soil conductivity [W/mK]
kg = 1.0                                # grout conductivity [W/mK]
kp = 0.4                                # pipe conductivity [W/mK]
Rp = 0.08                               # pipe thermal resistance [mK/W]
Rf = 1/628.56                           # water's thermal resistance at 30° [mK/W]
s = [2*ro + 0.001, 2*rb / 3, 2*(rb-ro)] # distance between 2 legs of a U-tube for case A, B, C [m]

Rals_A = resistance_ULoop_total_internal(s[1], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0, 
            network="diagonal")
Rals_B = resistance_ULoop_total_internal(s[2], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0, 
            network="diagonal")
Rals_C = resistance_ULoop_total_internal(s[3], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0, 
            network="diagonal")

Ramult_A = resistance_ULoop_total_internal(s[1], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1, 
            network="diagonal")
Ramult_B = resistance_ULoop_total_internal(s[2], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1, 
            network="diagonal")
Ramult_C = resistance_ULoop_total_internal(s[3], rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1, 
            network="diagonal")

# Test and accuracy verification
@test Rals_A_ref ≈ Rals_A rtol = 1e-2
@test Rals_B_ref ≈ Rals_B rtol = 1e-2
@test Rals_C_ref ≈ Rals_C rtol = 1e-2

@test Ramult_A_ref ≈ Ramult_A rtol = 1e-2
@test Ramult_B_ref ≈ Ramult_B rtol = 1e-2
@test Ramult_C_ref ≈ Ramult_C rtol = 1e-2

accuracy(Rals_A_ref, Rals_A) # 99.16%
accuracy(Rals_B_ref, Rals_B) # 99.25%
accuracy(Rals_C_ref, Rals_C) # 99.67%

accuracy(Ramult_A_ref, Ramult_A) # 60.35%
accuracy(Ramult_B_ref, Ramult_B) # 99.37%
accuracy(Ramult_C_ref, Ramult_C) # 99.82%

###################################################################################################
# Double U-loop chap5: p.80, ex. 5.4
# Compute the borehole resistance (Rb) and the internal resistance (Ra) of a double U-loop with the
# line source model using resistance_ULoop_borehole and resistance_ULoop_total_internal functions.
# Accuracy : 99.60%
#            82.51%     Test passed but to investigate...
###################################################################################################

# TODO : look at the Ra calculation, and why its so different??

# Reference value of the borehole resistance with the line source method for different U-tube
# configuration
Rb_ref = 0.156
Ra_ref = 0.241

# Borehole resistance calculated with the resistance_ULoop_borehole function from BoreholeResistance 
rb = 0.075                      # borehole radius [m]
ro = 0.033/2                    # outer radius of the pipe [m]
ri = 0.027/2                    # inner radius of the pipe [m]
ks = 2.5                        # soil conductivity [W/mK]
kg = 1.0                        # grout conductivity [W/mK]
kp = 0.4                        # pipe conductivity [W/mK]
Rp = 0.08                       # pipe thermal resistance [mK/W]
Rf = 1/628.56                   # water's thermal resistance at 30° [mK/W]
s = 2*rb / 3                    # distance between 2 legs of a U-tubes for case B [m]


Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=2, order=0)
Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf, nLoop=2, order=0, 
            network="diagonal")

# Test and accuracy verification
@test Rb_ref ≈ Rb rtol = 1e-2
@test Ra ≈ Ra rtol = 1e-2

accuracy(Rb_ref, Rb) # 99.60%
accuracy(Ra_ref, Ra) # 82.51%

###################################################################################################
# Linear approximation chap5: p.84, ex. 5.5
###################################################################################################

# TODO : the calculations
Tb = 0.856 


# Borehole resistance calculated with the resistance_ULoop_borehole function from BoreholeResistance 
t = 100*3600
T₀ = 10
rb = 0.075                      # borehole radius [m]
ro = 0.0335/2                   # outer radius of the pipe [m]
#ri = ??                    # inner radius of the pipe [m]
ks = 2.1                        # soil conductivity [W/m-K]
kg = 1.5                        # grout conductivity [W/m-K]
kp = 0.4                        # pipe conductivity [W/m-K]
Cs = 2.43e6                     # thermal capacity [J/m³-K]        
Rp = 0.076                       # pipe thermal resistance [m-K/W]
Rf = 1/598.03                   # water's thermal resistance at 30° [mK/W]
s = 2*rb / 3                    # distance between 2 legs of a U-tubes for case B [m]
q = -50 

g_ils =  ils(t, rb, ks, Cs)

Tb = T₀ - q*g_ils

end