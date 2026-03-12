includet("resistance_fluid.jl")
includet("resistance_pipe.jl")

"""
    resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, order=1)
    resistance_borehole_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0, order=1)

Function that computes the thermal borehole resistance of a ground heat exchanger based on the
zeroth- or first-order multipole formula for a single U-loop of Hellström (1991) for grout thermal 
resistance. It is the equivalent of the line source theory. See Eq. 8.36 of the Hellström (1991).
Note: To obtain only the grout thermal resistance, set `Rp` and `Rf` as 0.0.
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cf(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_rho(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_mu(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `Rp`: Pipe thermal resistance [mK/W]
    - `Rf`: Fluid thermal resistance [mK/W]
    - `order`: Order of the multipole method (default 1)
# Output
    - `Rb`: Borehole thermal resistance [mK/W]
# Reference
    - Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    - Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
"""
function resistance_borehole_multipole(s::Real, rb::Real, ro::Real, ks::Real, kg::Real, 
    Rp::Real, Rf::Real, order::Int=1)
    # Compute β
    Rₚ = Rp + Rf            # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ

    # Compute σ
    σ = (kg - ks) / (kg + ks)

    # Compute θ₁ to θ₃
    θ₁ = s / (2 * rb)   # equal to D/rb
    θ₂ = rb / ro
    if order == 1
        θ₃ = ro / s
    end

    if order == 0
        # Compute Rb with Eq. 12 from Javed and Spitler 2017
        Rb = (1 / (4 * π * kg)) * (β + log(θ₂ / (2 * θ₁ * (1 - θ₁^4)^σ)))
    elseif order == 1
        # Compute Rb with Eq. 13 from Javed and Spitler 2017
        tmp = (1 - θ₁^4)
        a = log(θ₂ / (2 * θ₁ * tmp^σ))
        b = θ₃^2 * (1 - (4 * σ * θ₁^4 / tmp))^2
        c = ((1 + β) / (1 - β)) + (θ₃^2 * (1 + (16 * σ * θ₁^4 / (tmp^2))))
        Rb = (1 / (4 * π * kg)) * (β + a - (b / c))
    else
        error("Only order 0 and 1 are implemented for the multipole method.")
    end
    return Rb
end
function resistance_borehole_multipole(V::Real, s::Real, rb::Real, ro::Real, ri::Real, ks::Real,
    kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0, order::Int=1)
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp, 2)

    # Compute Rb with the first-order multipole method
    return resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, order)
end

"""
    resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, order=1)
    resistance_total_internal_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0,
        order=1)

Computes the zeroth- or first-order multipole method for the total internal resistance (Eq. 25 or 26
of Javed and Spitler 2017) by Hellström 1991 for a single U-tube ground heat exchanger. This is used
to convert the borehole thermal resistance (Rb) in effective borehole thermal resistance (named Rb*
or Rbₑ).
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cf(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_rho(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_mu(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `Rp`: Pipe thermal resistance [mK/W]
    - `Rf`: Fluid thermal resistance [mK/W]
    - `order`: Order of the multipole method (default 1)
# Output
    - `Ra`: Total internal thermal resistance [mK/W]
# Reference
    - Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    - Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
"""
function resistance_total_internal_multipole(s::Real, rb::Real, ro::Real, ks::Real, kg::Real,
    Rp::Real, Rf::Real, order::Int=1)
    # Compute β
    Rₚ = Rf + Rp            # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ

    # Compute σ
    σ = (kg - ks) / (kg + ks)

    # Compute θ₁ and θ₃
    θ₁ = s / (2 * rb)
    θ₃ = ro / s

    # Compute Rₐ with Eq. 26 from Javed and Spitler 2017
    if order == 0
        Ra = (1 / (π * kg)) * (β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ)))
    elseif order == 1
        Ra = (1 / (π * kg)) * (β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ)) -
            ((θ₃^2 * (1 - (θ₁^4) + (4 * σ * θ₁^2))^2) /
            (((1 + β) / (1 - β)) * (1 - θ₁^4)^2 - (θ₃^2 * (1 - θ₁^4)^2) +
            (8 * σ * θ₁^2 * θ₃^2 * (1 + θ₁^4)))))
    else
        error("Only order 0 and 1 are implemented for the multipole method.")
    end
    return Ra
end
function resistance_total_internal_multipole(V::Real, s::Real, rb::Real, ro::Real, ri::Real, 
    ks::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0, order::Int=1)
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp, 2)

    # Compute Ra with the first-order multipole method
    return resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, order)
end

"""
    resistance_borehole_effective(V, H, cf, ρf, Rb, Ra)
    resistance_borehole_effective(V, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf)
    resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

Function that computes the effective thermal borehole resistance (also named Rb*). Effective Rb
allows considering the thermal short-circuiting along the borehole. Two types of boundary 
conditions are commonly used: (1) uniform borehole wall temperature (UBW) or (2) uniform heat 
flux (UHF). The most practical approach is to use an average of both approach.
Note: This application is valid for single U-loop (n = 2).
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `H`: Borehole length [m]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cf(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_rho(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_mu(T)`)
    - `Rb`: Borehole thermal resistance [mK/W]
    - `Ra`: Total internal thermal resistance [mK/W]
    - `ϵ`: Pipe roughness [m] (default 0.0)
# Outputs
    - `Rbₑ`: Effective borehole thermal resistance [mK/W]
# Reference
    - Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal 
        resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
    - Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. 
        Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. 
        https://doi.org/10.1016/B978-0-08-100311-4.00003-0
"""
function resistance_borehole_effective(V::Real, H::Real, cf::Real, ρf::Real, Rb::Real, Ra::Real)
    # UBW - See Eq. 3.68-3.70 of Javed et Spitler (2016)
    R1b = 2 * Rb                                    # Eq. 3.12
    R12 = (2 * Ra * R1b) / (2 * R1b - Ra)           # Eq. 3.14
    tmp = H / (V * cf * ρf)
    η = tmp / (2 * Rb) * sqrt(1 + 4 * Rb / R12)     # Eq. 3.69
    if η <= 1
        Rbₑ1 = Rb + (1 / (3 * R12)) * tmp^2 + (1 / (12 * Rb)) * tmp^2
    else
        Rbₑ1 = Rb * η * coth(η)
    end

    # UHF - See Eq. 3.67 of Javed et Spitler (2016)
    Rbₑ2 = Rb + (1 / (3 * Ra)) * tmp^2

    # Final calculation of the effective borehole thermal resistance
    return 0.5 * (Rbₑ1 + Rbₑ2)
end
function resistance_borehole_effective(V::Real, H::Real, s::Real, rb::Real, ro::Real, ks::Real,
    kg::Real, cf::Real, ρf::Real, Rp::Real, Rf::Real)
    # Compute Rb and Ra with the first-order multipole methods
    Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, 1)
    Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 1)

    # Compute Rbₑ
    return resistance_borehole_effective(V, H, cf, ρf, Rb, Ra)
end
function resistance_borehole_effective(V::Real, H::Real, s::Real, rb::Real, ro::Real, ri::Real,
    ks::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0)    
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp, 2)

    # Compute Rb and Ra with the first-order multipole method
    Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf, 1)
    Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf, 1)

    # Compute Rbₑ
    return resistance_borehole_effective(V, H, cf, ρf, Rb, Ra)
end