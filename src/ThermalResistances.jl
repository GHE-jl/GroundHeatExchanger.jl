"""
Fuctions that allows to compute the borehole thermal resistance Rb of a ground heat exchanger.
The basic local equation is:
    Rb = (Tf - Tb) / qb
where Tf is the local average fluid temperature, Tb is the borehole wall temperature, q_b is the
heat flux at the borehole wall.
Rb is defined in 3 components:
    Rb = Rg + (Rp + Rf)/n
where Rg is the grout thermal thermal resistance, Rp is the pipe thermal resistance, Rf is the
fluid convection thermal resistance and n is the number of pipes in the borehole.
Furthermore, to have whole borehole thermal resistance, the local one (Rb) have to be transfered in
effective borehole thermal resistance Rbₑ (or Rb*).

Author: Gabriel Dion
Date: 2025-08
"""

function R_f(V̇::T, kf::T, ri::T, cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3) where {T <: Real}
    """
        R_f(V̇, kf, ri, cf, ρf, μf)

    Function that computes the convective thermal resistance of a fluid in contact with a surface.
    The fluid is assumed to be flowing in a cylinder pipe.
    Inputs:
        - V̇: Fluid *speed* in pipe [m/s]
        - kf: Fluid thermal conductivity [W/mK]
        - ri: Pipe inside radius [m]
        - cf (default 4200.0): Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf (default 1000.0): Fluid density [kg/m³]. Default for water at 10 °C.
        - μf (default 1.3e-3): Fluid viscosity [kg/m⋅s]. Default for water at 10 °C.
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
        # Average between 4.364 for UHF and 3.657 for UBW boundary conditions
        Nu = 4.0                                # Eq. 2.42 of Lamarche 2023
    elseif Re >= 2300 && Re < 4000              # Transition phase
        γ = (Re - 2300) / (4000 - 2300)         # Eq. 2.49 of Lamarche 2023
        # Dittus-Boelter (Eq. 2.43a of Lamarche 2023)
        #=mode = "Heating"
        if mode == "Heating" || "heating"
            Nu_4k = 0.023 * 4000^0.8 * Pr^0.4
        elseif mode == "Cooling" || "cooling"
            Nu_4k = 0.023 * 4000^0.8 * Pr^0.3
        end=#
        # Gnielinski (Eq. 2.43b of Lamarche 2023)
        f = (0.79 * log(4000) - 1.64)^(-2)
        Nu_4k = (f / 8) * (4000 - 1000) * Pr / (1 + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
        Nu = (1 - γ) * 4.0 + γ * Nu_4k        # Eq. 2.48 of Lamarche 2023
    elseif Re >= 4000                           # Turbulent flow in a pipe
        # Dittus-Boelter (Eq. 2.43a Lamarche 2023)
        #=mode = "Heating"
         if mode == "Heating" || "heating"
            Nu = 0.023 * Re^0.8 * Pr^0.4        
         elseif mode == "Cooling" || "cooling"
            Nu = 0.023 * Re^0.8 * Pr^0.3
         end=#
         # Gnielinski (Eq. 2.43b of Lamarche 2023)
        f = (0.79 * log(Re) - 1.64)^(-2)
        Nu = (f / 8) * (Re - 1000) * Pr / (1 + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
    end
    
    # Convection coefficient
    h = Nu * kf / (2 * ri)                      # Eq. 2.32 of Lamarche 2023

    # Fluid convective resistance
    return 1 / (2 * pi * ri * h)            # Eq. 5.6 of Lamarche 2023
end

function R_p(kp::T, ro::T, ri::T) where {T <: Real}
    """
        R_p(kp, ro, ri)

    Function that computes the pipe thermal resistance (radial conduction resistance around a
    cylinder pipe).
    Inputs:
        - kp: Pipe thermal conductivity [W/mK]
        - ri, ro: Pipe inlet and pipe outlet radius [m]
        - np: Number of pipe in a U-loop configuration (2, 4, 6) (-)
    Output:
        - Rp: Pipe conductive thermal resistance [mK/W]
    Reference:
        Bergman, T.L., Incropera, F.P.: Fundamentals of Heat and Mass Transfer, 7th edn. Wiley, New 
        York (2011)
        Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
        Springer Nature Switzerland.
    """
    return log(ro / ri) / (2 * π * kp)
end

function R_b_zeroth_order_multipole(ks::T, kg::T, rb::T, ro::T, s::T, Rp::T, Rf::T
    ) where {T <: Real}
    """
        R_b_zeroth_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    Function that computes the thermal borehole resistance of a ground heat exchanger based on the
    zeroth-order multipole formula for a single U-loop of Hellström (1991) for grout thermal 
    resistance. It is the equivalent of the line source theory.
    See Eq. 8.36 of the reference.
    Note: To obtain only the grout thermal resistance, write Rp and Rf as 0.0.
    Inputs:
        - ks, kg: Ground and grout thermal conductivity [W/mK]
        - rb, ro: Borehole and pipe outlet radius [m]
        - s: Shank spacing (distance between 2 legs of a U-tubes) [m]
        - Rp: Pipe thermal resistance [mK/W]
        - Rf: Fluid thermal resistance [mK/W]
    Output:
        - Rb: Borehole thermal resistance [mK/W]
    Reference:
        Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
        Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    """
    # Compute β
    Rₚ = Rp + Rf            # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ

    # Compute σ
    σ = (kg - ks) / (kg + ks)

    # Compute θ₁ to θ₃
    θ₁ = s / (2 * rb)   # equal to D/rb
    θ₂ = rb / ro
    #@show [θ₁, θ₂]
    
    # Compute Rb with Eq. 12 from Javed and Spitler 2017
    Rb = (1 / (4 * π * kg)) * (β + log(θ₂ / (2 * θ₁ * (1 - θ₁^4)^σ)))

    #Compute Rb with Eq. 8.36 from Hellström (1991) with division by 4, since there is 2 pipe
    #Rb = (1 / (4 * pi * kg)) * (log(rb^2 / (ro * s)) + σ * log(rb^4 / (rb^4  - (s / 2)^4))) + Rp
    return Rb
end

function R_b_first_order_multipole(ks::T, kg::T, rb::T, ro::T, s::T, Rp::T, Rf::T) where {T <: Real}
    """
        R_b_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    Computes the first-order multipole method for the borehole thermal resistance (Eq. 13 of Javed
    and Spitler 2017) by Hellström 1991 and valid for a single U-tube ground heat exchanger.
    See Eq. 8.69 of the reference.
    Note: To obtain only the grout thermal resistance, write Rp and Rf as 0.0.
    Inputs:
        - ks, kg: Ground and grout thermal condcutvitiy [W/mK]
        - rb, ro: Borehole and pipe outlet radius [m]
        - s: Shank spacing (distance between 2 legs of a U-tubes) [m]
        - Rp: Pipe thermal resistance [mK/W]
        - Rf: Fluid thermal resistance [mK/W]
    Output:
        - Rb: Borehole thermal resistance [mK/W]
    Reference:
        Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
        Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    """
    # Compute β
    Rₚ = Rp + Rf            # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ

    # Compute σ
    σ = (kg - ks) / (kg + ks)

    # Compute θ₁ to θ₃
    θ₁ = s / (2 * rb)   # equal to D/rb
    θ₂ = rb / ro
    θ₃ = ro / s
    #@show [θ₁, θ₂, θ₃]

    # Compute Rb with Eq. 13 from Javed and Spitler 2017
    tmp = (1 - θ₁^4)
    a = log(θ₂ / (2 * θ₁ * tmp^σ))
    b = θ₃^2 * (1 - (4 * σ * θ₁^4 / tmp))^2
    c = ((1 + β) / (1 - β)) + (θ₃^2 * (1 + (16 * σ * θ₁^4 / (tmp^2))))
    Rb = (1 / (4 * π * kg)) * (β + a - (b / c))

    # Compute Rb with Eq. 8.69 from Hellström (1991) with division by 4, since there is 2 pipe
    #=D = s / 2
    tmp = rb^4 - D^4
    tmp2 = ro^2 / s^2
    Rb = (1 / (4 * π * kg)) * ((β + log(rb / ro) + log(rb / s) + (σ*log(rb^4 / tmp))) - 
        ((tmp2 * (1 - σ * (4 * D^4 / tmp))^2) / (((1 + β) /
        (1 - β)) + tmp2 * (1 + (σ * 16 * D^4 * rb^4 / tmp^2)))))=#
    return Rb
end

function R_b(V̇::T, ks::T, kg::T, kp::T, kf::T, rb::T, ro::T, ri::T, s::T, n::Integer = 2,
    cf::T = 4200.0, ρf::T = 1000.0, μf::T = 1.3e-3) where {T <:Real}
    """
        R_b(V̇, ks, kg, kp, kf, rb, ro, ri, s, n, cf, ρf, μf)
    
    Function that sums the thermal resistance components of a borehole as a function of its
    construction parameter and operating conditions.
    Inputs:
        - V̇: Fluid *speed* in pipe [m/s]
        - ks, kg, kp, kf: Ground, grout, pipe and fluid thermal conductivity [W/mK]
        - rb, ro, ri: Borehole, pipe outside and pipe inside radius [m]
        - s: Shank spacing (distance between 2 legs of a U-loop) [m]
        - n (default 2): Number of pipes: Pairs (2, 4 or 6) for U-loops, and 1 for coaxial [-]
        - cf (default 4200.0): Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf (default 1000.0): Fluid density [kg/m³]. Default for water at 10 °C.
        - μf (default 1.3e-3): Fluid viscosity [kg/m⋅s]. Default for water at 10 °C.
    Output:
        - Rb: Borehole thermal resistance [mK/W]
    """
    # Fluid thermal resistance
    Rf = R_f(V̇, kf, ri, cf, ρf, μf)

    # Pipe thermal resistance
    Rp = R_p(kp, ro, ri)

    # Borehole thermal resistance
    if n == 1
        # Coaxial configuration
        # TODO: Implement a method
    elseif n == 2
        # Single U-loop configuration
        Rb = R_b_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)
    elseif n == 4
        # Double U-loop configuration
        # TODO: Implement a valid method
    elseif n % 2 != 0
        throw(ArgumentError("'n' must be 1 (coaxial), 2 (single U-loop) or 4 (double U-loop)"))
    end
    return Rb
end

function R_a_first_order_multipole(ks::T, kg::T, rb::T, ro::T, s::T, Rp::T, Rf::T) where {T <: Real}
    """
        R_a_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    Computes the first-order multipole method for the total internal resistance (Eq. 26 of Javed 
    and Spitler 2017) by Hellström 1991 for a single U-tube ground heat exchanger. This is used to 
    convert the borehole thermal resistance (Rb) in effective borehole thermal resistance 
    (Rb*, or Rbₑ).
    Inputs:
        - ks, kg: Ground and grout thermal condcutvitiy [W/mK]
        - rb, ro: Borehole and pipe outlet radius [m]
        - s: Shank spacing (distance between 2 legs of a U-tubes) [m]
        - Rp: Pipe thermal resistance [mK/W]
        - Rf: Fluid thermal resistance [mK/W]
    Output:
        - Ra: Total internal thermal resistance [mK/W]
    Reference:
        Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
        Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    """
    # Compute β
    Rₚ = Rf + Rp            # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ

    # Compute σ
    σ = (kg - ks) / (kg + ks)

    # Compute θ₁ and θ₃
    θ₁ = s / (2 * rb)
    θ₃ = ro / s

    # Compute Rₐ with Eq. 26 from Javed and Spitler 2017
    Ra = (1 / (π * kg)) * (β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ)) -
          ((θ₃^2 * (1 - (θ₁^4) + (4 * σ * θ₁^2))^2) /
          (((1 + β) / (1 - β)) * (1 - θ₁^4)^2 - (θ₃^2 * (1 - θ₁^4)^2) +
           (8 * σ * θ₁^2 * θ₃^2 * (1 + θ₁^4)))))
    return Ra
end

function R_bₑ(V::T, cf::T, ρf::T, H::T, Rb::T, Ra::T) where {T <:Real}
    """
        R_bₑ(V, cf, ρf, H, Rb, Ra)
    
    Function that computes the effective thermal borehole resistance (also named Rb*). Effective Rb
    allows considering the thermal short-circuiting along the borehole. Two types of boundary 
    conditions are commonly used: (1) uniform borehole wall temperature (UBW) or (2) uniform heat 
    flux (UHF). The most practical approach is to use an average of both approach.
    Note: This application is valid for single U-loop (n = 2).
    Inputs:
        - V: Fluid flow rate [m³/s]
        - cf: Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf: Fluid density [kg/m³]. Default for water at 10 °C.
        - H: Borehole length [m]
        - Rb: Borehole thermal resistance [mK/W]
        - Ra: Total internal thermal resistance [mK/W]
    Outputs:
        - Rbₑ: Effective borehole thermal resistance [mK/W]
    Reference:
        Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal 
        resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
        Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. 
        Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. 
        https://doi.org/10.1016/B978-0-08-100311-4.00003-0
    """
    # Compute Rbₑ
    return _Rbₑ(V, cf, ρf, H, Rb, Ra)
end

function R_bₑ(V::T, ks::T, kg::T, rb::T, ro::T, s::T, cf::T, ρf::T, H::T, Rp::T, Rf::T
    ) where {T <:Real}
    """
        R_bₑ(V, ks, kg, rb, ro, s, cf, ρf, H, Rp, Rf)
    
    Function that computes the effective thermal borehole resistance (also named Rb*). Effective Rb
    allows considering the thermal short-circuiting along the borehole. Two types of boundary 
    conditions are commonly used: (1) uniform borehole wall temperature (UBW) or (2) uniform heat 
    flux (UHF). The most practical approach is to use an average of both approach.
    Note: This application is valid for single U-loop (n = 2).
    Inputs:
        - V: Fluid flow rate [m³/s]
        - ks, kg: Ground and grout thermal conductivity [W/mK]
        - rb, ro: Borehole and pipe ouside radius [m]
        - s: Shank spacing (distance between 2 legs of a U-loop) [m]
        - cf: Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf: Fluid density [kg/m³]. Default for water at 10 °C.
        - H: Borehole length [m]
        - Rp: Pipe thermal resistance [mK/W]
        - Rf: Fluid thermal resistance [mK/W]
    Outputs:
        - Rbₑ: Effective borehole thermal resistance [mK/W]
    Reference:
        Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal 
        resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
        Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. 
        Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. 
        https://doi.org/10.1016/B978-0-08-100311-4.00003-0
    """

    # Compute Rb and Ra with the first-order multipole method
    Rb = R_b_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    Ra = R_a_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    # Compute Rbₑ
    return _Rbₑ(V, cf, ρf, H, Rb, Ra)
end

function R_bₑ(V::T, ks::T, kg::T, kp::T, kf::T, rb::T, ro::T, ri::T, s::T, cf::T, ρf::T, μf::T, 
    H::T) where {T <:Real}
    """
        R_bₑ(V, ks, kg, kp, kf, rb, ro, ri, s, cf, ρf, μf, H)
    
    Function that computes the effective thermal borehole resistance (also named Rb*). Effective Rb
    allows considering the thermal short-circuiting along the borehole. Two types of boundary 
    conditions are commonly used: (1) uniform borehole wall temperature (UBW) or (2) uniform heat 
    flux (UHF). The most practical approach is to use an average of both approach.
    Note: This application is valid for single U-loop (n = 2).
    Inputs:
        - V: Fluid flow rate [m³/s]
        - ks, kg, kp, kf: Ground, grout, pipe and fluid thermal conductivity [W/mK]
        - rb, ro, ri: Borehole, pipe outside and pipe inside radius [m]
        - s: Shank spacing (distance between 2 legs of a U-loop) [m]
        - cf (default 4200.0, opt.): Fluid specific heat [J/kgK]. Default for water at 10 °C.
        - ρf (default 1000.0, opt.): Fluid density [kg/m³]. Default for water at 10 °C.
        - μf (default 1.3e-3, opt.): Fluid viscosity [kg/m⋅s]. Default for water at 10 °C.
        - H: Borehole length [m]
    Outputs:
        - Rbₑ: Effective borehole thermal resistance [mK/W]
    Reference:
        Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal 
        resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
        Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. 
        Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. 
        https://doi.org/10.1016/B978-0-08-100311-4.00003-0
    """
    V̇ = V / (π * ri^2)
    
    # Compute Rf
    Rf = R_f(V̇, kf, ri, cf, ρf, μf)

    # Compute Rp
    Rp = R_p(kp, ro, ri)

    # Compute Rb and Ra with the first-order multipole method
    Rb = R_b_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    Ra = R_a_first_order_multipole(ks, kg, rb, ro, s, Rp, Rf)

    # Compute Rbₑ
    return _Rbₑ(V, cf, ρf, H, Rb, Ra)

end

function _Rbₑ(V::T, cf::T, ρf::T, H::T, Rb::T, Ra::T) where {T <: Real}
    """
        _Rbₑ(V, cf, ρf, H, Rb, Ra)
    """
    # UBW - See Eq. 3.68-3.70 of Javec et Spitler (2016)
    R1b = 2 * Rb                                    # Eq. 3.12
    R12 = (2 * Ra * R1b) / (2 * R1b - Ra)           # Eq. 3.14
    tmp = H / (V * cf * ρf)
    η = tmp / (2 * Rb) * sqrt(1 + (4 * Rb / R12))   # Eq. 3.69
    if η <= 1
        Rbₑ1 = Rb + (1 / (3 * R12)) * tmp^2 + (1 / (12 * Rb)) * tmp^2
    else
        Rbₑ1 = Rb * η * coth(η)
    end

    # UHF - See Eq. 3.67 of Javec et Spitler (2016)
    Rbₑ2 = Rb + (1 / (3 * Ra)) * tmp^2

    # Final calculation of the  effective borehole thermal resistance
    return (Rbₑ1 + Rbₑ2) / 2
end