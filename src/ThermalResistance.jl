"""
Collection of functions that computes the borehole thermal resistance of a ground heat exchanger.
As mentionnend in Javed and Spitler (2017), the first-order multipole method is quite the best and
simplest approach.
"""

function Rb_first_order_multipole(V::T, ks::T, kg::T, kp::T, kf::T, cf::T, ρf::T, μf::T,
    rb::T, ri::T, ro::T, s::T) where T<:Real
    """
    Computes the first-order multipole method for the borehole thermal resistance
    (Eq. 13 of Javed and Spitler 2017) for a single U-tube ground heat exchanger.
    
    Inputs:
        - ̇V: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Ground, grout, pipe, fluid thermal condcutvitiy [W/mK]
        - cf: Fluid specific heat [J/kgK]
        - ρf: Fluid density [kg/m^3]
        - μf: Fluid viscosity [kg/sm]
        - rb, ri, ro: Borehole, inlet pipe and outlet pipe radius [m]
        - s: Shank spacing (distance between the two U-tubes) [m]
    
    Output:
        - Rb: Borehole thermal resistance [mK/W]
    
    Reference:
    Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods for
    grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
    https://doi.org/10.1016/j.apenergy.2016.11.079
    
    Author: Gabriel Dion
    Date: 2025-02
    """

    # Compute β
    Pr = cf * μf / kf                                               # Prandtl number
    Re = 2 * ρf * V / μf                                            # Reynold number
    Nu = 0.023 * Re^0.8 * Pr^0.4                                    # Nusselt number
    h = Nu * kf / (2 * ri)                                          # Convection coef.
    Rₚ = (log(ro / ri) / (2 * π * kp)) + (1 / (2 * π * ri * h))     # Pipe/fluid resistance
    β = 2 * π * kg * Rₚ                                             # β coef.

    # Compute σ
    σ = (kg - ks) / (kg + ks)                                       # σ coef.

    # Compute θ₁ to θ₃
    θ₁ = s / (2 * rb)
    θ₂ = rb / ro
    θ₃ = ro / s

    # Compute Rb with Eq. 13 from Javed and Spitler 2017
    Rb = (1 / (4 * π * kg)) * ((β + log(θ₂ / (2 * θ₁ * (1 - θ₁^4)^σ))) -
    (θ₃^2 * (1 - (4 * σ * θ₁^4) / (1 - θ₁^4))^2)/
    (((1 + β) / (1 - β)) + θ₃^2 * (1 + ((16 * σ * θ₁^4) / (1 - θ₁^4)^2))))
    return Rb
end

function Ra_first_order_multipole(V::T, ks::T, kg::T, kp::T, kf::T, cf::T, ρf::T, μf::T,
    rb::T, ri::T, ro::T, s::T) where T<:Real
    """
    Computes the first-order multipole method for the total internal resistance
    (Eq. 26 of Javed and Spitler 2017) for a single U-tube ground heat exchanger.
    
    Inputs:
        - ̇V: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Soil, grout, pipe, fluid thermal condcutvitiy [W/mK]
        - cf: Fluid specific heat [J/kgK]
        - ρf: Fluid density [kg/m^3]
        - μf: Fluid viscosity [kg/sm]
        - rb, ri, ro: Borehole, inlet pipe and outlet pipe radius [m]
        - s: Shank spacing (distance between the two U-tubes) [m]
    
    Output:
        - Rₐ: Borehole thermal resistance [mK/W]
    
    Reference:
    Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods for
    grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
    https://doi.org/10.1016/j.apenergy.2016.11.079
    
    Author: Gabriel Dion
    Date: 2025-02
    """

    # Compute β
    Pr = cf * μf / kf                                               # Prandtl number
    Re = 2 * ρf * V / μf                                            # Reynold number
    Nu = 0.023 * Re^0.8 * Pr^0.4                                    # Nusselt number
    h = Nu * kf / (2 * ri)                                          # Convection coef.
    Rₚ = (log(ro / ri) / (2 * π * kp)) + (1 / (2 * π * ri * h))     # Pipe/fluid resistance
    β = 2 * π * kg * Rₚ                                             # β coef.

    # Compute σ
    σ = (kg - ks) / (kg + ks)                                       # σ coef.

    # Compute θ₁ and θ₃
    θ₁ = s / (2 * rb)
    θ₃ = ro / s

    # Compute Rₐ with Eq. 26 from Javed and Spitler 2017
    Rₐ = (1 / (π * kg)) * ((β + log((1+θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ))) -
    (θ₃^2 * (1 - θ₁^4 + 4 * σ * θ₁^2)^2) /
    (((1 + β) / (1 - β)) * (1 - θ₁^4)^2 - θ₃^2 * (1 - θ₁^4)^2 + 8 * σ * θ₁^2 * θ₃^2 * (1 + θ₁^4)))
    return Rₐ
end