"""
Module that allows to compute the borehole thermal resistance of a ground heat exchanger. The basic
equation is:
T_f = T_b - q_b * R_b
where T_f is the average fluid temperature, T_b is the average borehole wall temperature, q_b is the
heat flux at the borehole wall and R_b is the borehole thermal resistance.
R_b is defined in 3 components:
R_b = R_g + (R_p + R_f)/n_p
where R_g is the grout thermal thermal resistance, R_p is the pipe thermal resistance, R_f is the
fluid convection thermal resistance and n_p is the number of pipes in the borehole (2, 4, 6).

Author: Gabriel Dion
Date: 2025-07
"""

module ThermalResistance

export R_p,
       R_f,
       Rb_first_order_multipole

function R_f(V̇::T, kf::T, ri::T, cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3, 
    mode::String = "Heating") where {T <: Real}
    """
        R_f()

    Function that computes the convective thermal resistance of a fluid in contact with a surface.
    The fluid is assumed to be flowing in a cylinder pipe.
    Inputs:
        - V̇: Fluid *speed* in pipe [m/s]
        - kf: Fluid thermal conductivity [W/mK]
        - ri: Pipe inside radius [m]
        - cf (default 4200.0): Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf (default 1000.0): Fluid density [kg/m³]. Default for water at 10 °C.
        - μf (default 1.3e-3): Fluid viscosity [kg/m⋅s]. Default for water at 10 °C.
        - mode (default "Heating"): Specify if the heat load is in Heating or Cooling in Nusslet
            number evaluation.
    Output:
        - Rf: Fluid convective thermal resistance [mK/W]
    Reference:
    Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
    Springer Nature Switzerland.
    """
    # Reynold number
    Re = 2 * ri * ρf * V̇ / μf                   # Eq. 2.33 of Lamarche 2023

    # Prandtl number
    Pr = cf * μf / kf                           # Eq. 2.34 of Lamarche 2023

    # Nusselt number
    if Re < 2300                                # Laminar phase
        Nu = 4.364                              # Constant heat flux (Eq. 2.42 of Lamarche 2023)
    elseif Re >= 2300 && Re < 4000              # Transition phase
        γ = (Re - 2300) / (4000 - 2300)         # Eq. 2.49 of Lamarche 2023
        if mode == "Heating" || "heating"
            Nu_4k = 0.023 * 4000^0.8 * Pr^0.4
        elseif mode == "Cooling" || "cooling"
            Nu_4k = 0.023 * 4000^0.8 * Pr^0.3
        end
        Nu = (1 - γ) * 4.364 + γ * Nu_4k        # Eq. 2.48 of Lamarche 2023
    elseif Re >= 4000                           # Turbulent flow in a pipe
         if mode == "Heating" || "heating"
            Nu = 0.023 * Re^0.8 * Pr^0.4        # Dittus-Boelter heating (Eq. 2.43a Lamarche 2023)
         elseif mode == "Cooling" || "cooling"
            Nu = 0.023 * Re^0.8 * Pr^0.3        # Dittus-Boelter cooling (Eq. 2.43a Lamarche 2023)
         end
    end
    
    # Convection coefficient
    h = Nu * kf / (2 * ri)                      # Eq. 2.32 of Lamarche 2023

    # Fluid convective resistance
    return 1 / (2 * pi * ri * h * n)            # Eq. 5.6 of Lamarche 2023
end

function R_p(kp::T, ri::T, ro::T) where {T <: Real}
    """
        R_p(kp, ri, ro)

    Function that computes the pipe thermal resistance (radial conduction resistance around a
    cylinder pipe).
    Inputs:
        - kp: Pipe thermal conductivity [W/mK]
        - ri, ro: Pipe inlet and pipe outlet radius [m]
        - np: Number of pipe in a U-loop configuration (2, 4, 6) (-)
    Output:
        - Rp: Pipe conductive thermal resistance [mK/W]
    Reference:
    Bergman, T.L., Incropera, F.P.: Fundamentals of Heat and Mass Transfer, 7th edn. Wiley, New York
    (2011)
    Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
    Springer Nature Switzerland.
    """
    return log(ro / ri) / (2 * π * kp)
end

function Rb_Hellstrom(
        V̇::T, ks::T, kg::T, kp::T, kf::T, rb::T, ri::T, ro::T, s::T,
        cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3, n::Int = 2) where {T <: Real}
    """
        Rb_Hellstrom(V̇, ks, kg, kp, kf, rb, ri, ro, s, cf = 4200.0, ρf = 1000.0, μf = 1.3e-3, n=2)

    Function that computes the thermal borehole resistance of a ground heat exchanger based on the
    method of Hellström (1991). This method is based on the combination of fluid, pipe and grout 
    thermal resistance.
    Inputs:
        - V̇: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Ground, grout, pipe and fluid thermal conductivity [W/mK]
        - rb, ri, ro: Borehole, pipe inlet and pipe outlet radius [m]
        - s: Shank spacing (distance between the two U-tubes) [m]
        - n (default 2): Number of pipes (=2 for single U-tube and =4 for double U-tube) [-]
        - cf (default 4200.0): Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf (default 1000.0): Fluid density [kg/m³]. Default for water at 10 °C.
        - μf (default 1.3e-3): Fluid viscosity [kg/m⋅s]. Default for water at 10 °C.
    Output:
        - Rb: Thermal borehole resistance [mK/W]
    Reference:
        - Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    """

    # Rf - Fluid resistance
    Pr = cf * μf / kf                                               # Prandtl number
    Re = 2 * ri * ρf * V̇ / μf                                       # Reynold number
    Nu = Nusselt_number(Re, Pr)                                     # Nusselt number
    h = Nu * kf / (2 * ri)                                          # Convective heat transfer coef.
    Rf = 1 / (2 * pi * ri * h * n)

    # Rp - Pipe resistance
    Rp = R_p(kp, ri, ro) / n

    # Rg - Grout resistance
    σ = (kg - ks) / (kg + ks)
    Rg = (1 / (4 * pi * kg)) *
         (log(rb / ro) + log(rb / s) + σ * log(((rb / (s / 2))^4) / ((rb / (s / 2))^4 - 1)))

    # Rb - Borehole thermal resistance
    return Rf + Rp + Rg
end

function Rb_first_order_multipole(
        V̇::T, ks::T, kg::T, kp::T, kf::T, rb::T, ri::T, ro::T, s::T,
        cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3) where {T <: Real}
    """
        Rb_first_order_multipole(V̇, ks, kg, kp, kf, rb, ri, ro, s,
            cf = 4200.0, ρf = 1000.0, μf = 1.3e-3)

    Computes the first-order multipole method for the borehole thermal resistance (Eq. 13 of Javed
    and Spitler 2017) for a single U-tube ground heat exchanger.
    Inputs:
        - ̇V: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Ground, grout, pipe, fluid thermal condcutvitiy [W/mK]
        - cf: Fluid specific heat [J/kgK]
        - ρf: Fluid density [kg/m^3]
        - μf: Fluid viscosity [kg/sm]
        - rb, ri, ro: Borehole, pipe inlet and pipe outlet radius [m]
        - s: Shank spacing (distance between the two U-tubes) [m]
    Output:
        - Rb: Borehole thermal resistance [mK/W]
    Reference:
        Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    """

    # Compute β
    Pr = cf * μf / kf                                               # Prandtl number
    Re = 2 * ri * ρf * V̇ / μf                                       # Reynold number
    Nu = Nusselt_number(Re, Pr)                                     # Nusselt number
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
          (θ₃^2 * (1 - (4 * σ * θ₁^4) / (1 - θ₁^4))^2) /
          (((1 + β) / (1 - β)) + θ₃^2 * (1 + ((16 * σ * θ₁^4) / (1 - θ₁^4)^2))))
    return Rb
end

function Ra_first_order_multipole(
        V̇::T, ks::T, kg::T, kp::T, kf::T, rb::T, ri::T, ro::T, s::T,
        cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3) where {T <: Real}
    """
        Ra_first_order_multipole(V̇, ks, kg, kp, kf, rb, ri, ro, s,
            cf = 4200.0, ρf = 1000.0, μf = 1.3e-3)

    Computes the first-order multipole method for the total internal resistance (Eq. 26 of Javed 
    and Spitler 2017) for a single U-tube ground heat exchanger.
    Inputs:
        - V̇: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Soil, grout, pipe, fluid thermal condcutvitiy [W/mK]
        - cf: Fluid specific heat [J/kgK]
        - ρf: Fluid density [kg/m^3]
        - μf: Fluid viscosity [kg/sm]
        - rb, ri, ro: Borehole, inlet pipe and outlet pipe radius [m]
        - s: Shank spacing (distance between the two U-tubes) [m]
    Output:
        - Rₐ: Borehole thermal resistance [mK/W]
    Reference:
        Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    """

    # Compute β
    Pr = cf * μf / kf                                               # Prandtl number
    Re = 2 * ri * ρf * V̇ / μf                                       # Reynold number
    Nu = Nusselt_number(Re, Pr)                                     # Nusselt number
    h = Nu * kf / (2 * ri)                                          # Convection coef.
    Rₚ = (log(ro / ri) / (2 * π * kp)) + (1 / (2 * π * ri * h))     # Pipe/fluid resistance
    β = 2 * π * kg * Rₚ                                             # β coef.

    # Compute σ
    σ = (kg - ks) / (kg + ks)                                       # σ coef.

    # Compute θ₁ and θ₃
    θ₁ = s / (2 * rb)
    θ₃ = ro / s

    # Compute Rₐ with Eq. 26 from Javed and Spitler 2017
    Rₐ = (1 / (π * kg)) * ((β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ))) -
          (θ₃^2 * (1 - θ₁^4 + 4 * σ * θ₁^2)^2) /
          (((1 + β) / (1 - β)) * (1 - θ₁^4)^2 - θ₃^2 * (1 - θ₁^4)^2 +
           8 * σ * θ₁^2 * θ₃^2 * (1 + θ₁^4)))
    return Rₐ
end

end