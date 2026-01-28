
using SpecialFunctions: erf
using Roots: find_zero
includet("ILS.jl")
includet("../Convolutions.jl")

"""
    βils_outlet(t, k, kp, rb, ro, ri, H, Hp, V, β, T0; K = K)

Compute the borehole outlet transfer function from the β-ILS model and the thermal resistance
of the SCW. The function first computes the borehole wall (g-function) using the β-ILS model, then
adjusts it to obtain the outlet temperature transfer function using the development from Jacques et 
al. (2025).
# Arguments
    - t: Time vector (nₜ x 1) [s]
    - k: A (ix3) matrix of "i" geological layers and their thermal properties:
        - Column 1: Thickness [m]
        - Column 2: Thermal conductivity [W/mK]
        - Column 3: Volumetric specific heat [J/m³K]
        - Example: [20 k₁ C₁; 80 K₂ C₂] for 2 layers of 20m and 80m thicknesses
        - Note: For a single layer ground, the first value is equal to `H`.
    - kp: Pipe thermal conductivity [W/mK]
    - rb: Borehole radius [m]
    - ro: Outer pipe radius [m]
    - ri: Inner pipe radius [m]
    - H: Borehole depth [m]
    - Hp: Pipe length [m]
    - V: Circulating fluid flow rate ([m³/s])
    - β: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - Tf: Initial ground temperature [°C] (used for water properties)
    - K: (Optional) A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
        - Note: If only one hydrogeological layer is present, this argument can be omitted.
# Output
    - g_outlet: A g-function corresponding to the borehole outlet temperature of the SCW [°Cm/W]
# Reference
    - Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
        standing column wells operating with bleed. Applied Thermal Engineering.
        https://doi.org/10.1016/j.applthermaleng.2024.124543
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function βils_outlet(t, k::AbstractMatrix{T}, kp::T, rb::T,
        ro::T, ri::T, H::T, Hp::T, V::T, β::T, Tf::T; K::Union{Nothing, AbstractMatrix{T}}=nothing
        ) where {T<:AbstractFloat}
    # 1. Compute the thermal resistance of the SCW
    Rfb, Rv, f = Rb_SCW(V, kp, rb, ro, ri, H, Hp, Tf)

    # 2. Compute the borehole wall g-function using the β-ILS model
    gb = βils(t, k, rb, H, V, β, f; K = K)

    # 3. Adjust the borehole wall g-function to obtain the outlet g-function for 1 W/m
    go = gb .- (1.0 - β) * (Rv / (2.0 - β) - Rfb)
    return go
end

"""
    βils(t, k, rb, H, V, β, f = 0.05; K = K)

Compute the analytical model developped by Nguyen et al. (2025) for standing column
wells (SCW). The output is a g-function that requires a heat load per unit of borehole 
length [W/m] to provide the borehole wall temperature.
# Arguments
    - t: Time vector (nₜ x 1) [s]
    - k: A (ix3) matrix of "i" geological layers and their thermal properties:
        - Column 1: Thickness [m]
        - Column 2: Thermal conductivity [W/mK]
        - Column 3: Volumetric specific heat [J/m³K]
        - Example: [20 k₁ C₁; 80 K₂ C₂] for 2 layers of 20m and 80m thicknesses
        - Note: For a single layer ground, the first value is equal to `H`.
    - V: Circulating fluid flow rate ([m³/s])
    - β: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - H: Borehole depth [m]
    - rb: Borehole radius [m]
    - K: (Optional) A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
        - Note: If only one hydrogeological layer is present, this argument can be omitted.
    - f: (Optional) Friction factor when computing the thermal resistance of the SCW (Rb_SCW()) [-]
# Output
    - g: A g-function corresponding to the borehole wall temperature of the SCW [°Cm/W]
# Reference
    - Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
        standing column wells operating with bleed. Applied Thermal Engineering.
        https://doi.org/10.1016/j.applthermaleng.2024.124543
"""
function βils(t, k::AbstractMatrix{T}, rb::T, H::T, V::T, β::T, f::T = T(0.05);
    K::Union{Nothing, AbstractMatrix{T}} = nothing) where {T<:AbstractFloat}

    ts = t isa AbstractVector{T} ? t : Ref(t)
    out = similar(ts, length(ts))

    @inbounds for i in eachindex(ts)
        out[i] = _βils(ts[i], k, rb, H, V, β, f; K = K)
    end

    return t isa AbstractVector ? out : out[1]
end

function _βils(t, k::AbstractMatrix{T}, rb::T, H::T, V::T, β::T, f::T = T(0.05);
    K::Union{Nothing, AbstractMatrix{T}}=nothing) where {T<:AbstractFloat}
    # 0. Check inputs
    if size(k, 2) != 3
        error("k must be a (ix3) matrix of geological layers and their thermal properties.")
    end
    if !isnothing(K) && size(K, 2) != 2
        error("K must be a (jx2) matrix of hydrogeological layers and their hydraulic properties.")
    end

    # 1. Effective thermal properties (arithmetic weighted average)
    ks = sum(k[:, 2] .* k[:, 1]) / sum(k[:, 1])
    Cs = sum(k[:, 3] .* k[:, 1]) / sum(k[:, 1])

    # 2. Induced convergent flow consideration
    Vb = V * β                              # Bleed flow rate [m³/s]
    H₁, K̃₁, K̃₂, Kₕ₁, Kₕ₂, Kᵥ₁, Kᵥ₂ = effective_K(K, H) # Compute the effective hydraulic properties
    H₂ = H - H₁
    println("$H₁, $H₂, $K̃₁, $K̃₂, $Kₕ₁, $Kₕ₂, $Kᵥ₁, $Kᵥ₂")

    !isnothing(K) ? H̃e = H₁ : H̃e = H/2
    ΔHa = 0.22                              # Assumption: drawdown caused by recirculation
    ΔH̃ = ΔHa - (f * H * V^2 / (4 * π^2 * 9.81 * (rb - ro) * (rb^2 - ro^2)^2)) # Compute ΔH̃
    Ṽ = 2 * π * K̃₁* ΔH̃ * H̃e / (log(100 / rb)) # Assumption: 100 of drawdown horizontal radius
    println(Vb)
    println(ΔH̃)
    println(Ṽ)

    # 2. Layered heterogeneity - Effective length (if present)
    if !isnothing(K)
        Kₕ = sum(K[:, 2] .* K[:, 1]) / H    # Arithmetic weighted average
        Kᵥ = H / sum(K[:, 1] ./ K[:, 2])    # Harmonic weighted average
        K̃ = sqrt(Kₕ * Kᵥ)                   # Effective hydraulic conductivity
        He = H * K̃ / Kₕ                     # Effective length only considering heterogeneity
    else                                    # If no layered heterogeneity is present
        He = H
    end

    # 3. Compute the Peclet number and dimensionless time
    ta = @. (Ṽ / H̃e + Vb / He) * t / (π * rb^2) # Dimensionless time
    Pe = (Ṽ / H̃e + Vb / He) / (2 * π * (ks / Cs)) # Effective Peclet number
    println(Pe)
    
    # 4. Range of validity of the inputs based on scenarios in Nguyen et al. (2025)
    values_validity(t, ks, Cs, rb, H, V, β, Vb, Pe)

    # 5. Compute the scaling function h
    coefs₁ = [-0.4478, 0.1288, 0.6458, 0.1483]
    coefs₂ = [0.7110, -0.7595, 0.7838, 0.5114]
    coefs₃ = [1.0011, 0.7339, 0.2383, -0.0100]
    ϵ = @. coefs₁ * Pe^coefs₂ + coefs₃
    h = @. ϵ[1] * (1 - ϵ[2] * β) *
           (1 + erf((sqrt(Pe) / ϵ[3]) * ((1 - ta^ϵ[4]) / (ta^ϵ[4]))))
    
    # 6. Initial transfer function g_0 based on the SLI
    g₀ = ils(t, ks, Cs, rb)

    # 7. Combine initial and scaling function
    gb = convolution(diff([0; g₀]), h)
    return gb
end

"""
    Rb_SCW(V, kp, rb, ro, ri, H, Hp, Tf, dT)

Computation of theoretical thermal resistance in a SCW based on Jacques et al. 2025. This 
formulation assumes uniform flux along the borehole wall.
# Arguments
    - V: Circulating fluid flow rate [m³/s]
    - kp: Pipe thermal conductivity [W/mK]
    - rb: Borehole radius [m]
    - ro: Outer pipe radius [m]
    - ri: Inner pipe radius [m]
    - H: Borehole depth [m]
    - Hp: Pipe length [m]
    - Tf: Fluid temperature [°C]
    - dT: (Optional) Temperature difference at the heat exchanger 
        - Positive value: Ground heating -> the fluid is being cooled (Tb<Tf)
        - Negative value: Ground cooling -> the fluid is being heating (Tb>Tf)
# Output
    - Rfb: Subsurface effective borehole thermal resistance in a SCW [K·m/W]
    - Rv: Hydraulic resistance due to the fluid flow [K·m/W]
# Reference
    - [1] Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
    - [2] Seol, H., Jeong, S., Cho, C., & You, K. (2008). Shear load transfer for rock-
        socketed drilled shafts based on borehole roughness and geological strength index (GSI). 
        International Journal of Rock Mechanics and Mining Sciences, 45(6), 848–861. 
        https://doi.org/10.1016/j.ijrmms.2007.09.008
    - [3] Todorov, O., Alanne, K., Virtanen, M., & Kosonen, R. (2021). Different Approaches for 
        Evaluation and Modeling of the Effective Thermal Resistance of Groundwater-Filled Boreholes.
"""
function Rb_SCW(V::T, kp::T, rb::T, ro::T, ri::T, H::T, Hp::T, Tf::T, dT::T = 2.0, 
    ) where {T<:AbstractFloat}
    # Initialize parameters
    ρf = water_ρ(Tf)        # Water density at working fluid temperature
    cpf = water_cp(Tf)      # Water specific heat at working fluid temperature
    kf = water_k(Tf)        # Water thermal conductivity at working fluid temperature
    μf = water_μ(Tf)        # Water viscosity at working fluid temperature
    μf_dT = water_μ(Tf + dT) # Assumption: 2°C temperature difference between inlet and outlet
    np = 1                  # number of pipes (1 in a SCW)
    L = [Hp, H]             # Lengths for [1] pipe and [2] borehole for [m]
    ϵb = 5e-4               # Assumption: Borehole roughness [m]
    ϵp = 0.0                # Assumption: Pipe roughness [m]

    # Preallocation
    Pr = zeros(T, 2);
    rₑ = zeros(T, 2)
    Re = zeros(T, 2)
    Nu = zeros(T, 2)
    h = zeros(T, 2)
    Rf = zeros(T, 2)
    f = zero(T)
    
    # Prandtl numbers
    Pr[1] = μf * cpf / kf
    Pr[2] = μf * cpf / kf

    # Reynold numbers in pipe (i == 1) or borehole (i == 2)
    rₑ[1] = ri              # Equivalent radius in the pipe
    rₑ[2] = (rb - ro)       # Equivalent radius in the borehole
    Re[1] = ρf * (V / (π * ri^2)) * 2.0 * rₑ[1] / μf
    Re[2] = ρf * (V / ((π * rb^2) - (π * ro^2))) * 2.0 * rₑ[2] / μf

    # Nusselt number, loop for regions i = 1 (pipe), i = 2 (borehole)
    for i in 1:2
        if Re[i] < 2300                 # Laminar flow in pipe or borehole
            println("   -> Nusselt calculation for laminar flow")
            term1 = (Re[i] * Pr[i] * 2 * rₑ[i] / L[i])^0.333 * (μf / μf_dT)^0.14
            if i == 1 && term1 >= 2
                Nu[i] = 4.36
            elseif i == 1
                Nu[i] = 1.86 * term1
            elseif i == 2
                f = 64 / Re[i]
                Nu[i] = 3.66
            else
                error("Unexpected case in Nusselt number calculation. Laminar flow.")
            end

        elseif Re[i] > 10000 || i == 2  # Turbulent flow (Re>2300 in borehole, Ref[3])
            if i == 1                   # Turbulent flow in pipe
                println("   -> Nusselt calculation for turbulent flow in pipe")
                # Dittus-Boelter equation to solve for Nu
                dT > 0 ? n = 0.3 : n = 0.4 # Define fluid heating or cooling
                Nu[i] = 0.023 * Re[i]^0.8 * Pr[i]^n            
            elseif i == 2               # Turbulent flow in borehole
                println("   -> Nusselt calculation for turbulent flow in borehole")
                # Colebrook-White equation to solve for f in borehole
                ϵ = (ϵb * rb + ϵp * ro) / (rb + ro) # Equivalent roughness [m]
                obj_fun(F) = 1 / sqrt(F) + 2 * log10(ϵ / (3.7 * (2 * (rb - ro))) + 2.51 / 
                    (Re[i] * sqrt(F)))
                f = find_zero(obj_fun, 1e-2)
                # Nusselt number in borehole (see Ref[3])
                Nu[i] = ((f / 8) * (Re[i] - 1000) * Pr[i]) / (1 + 12.7 * (f / 8)^0.5 * 
                    (Pr[i]^(2 / 3) - 1))
            else
                error("Unexpected case in Nusselt number calculation. Turbulent flow.")
            end

        else    # Transitional flow in pipe only (i = 1 and 2300 < Re < 10000)
            println("   -> Nusselt calculation for transition flow in pipe.")
            f = (1.58 * log(Re[i]) - 3.28)^(-2) # Jacques et al. 2025
            # f = (1.82 * log(Re[i]) - 3.28)^(-2) # Petukhov 1970
            # f = (0.79 * log(Re[i]) - 1.64)^(-2) # Lamarche 2023
            Nu[i] = ((f / 2) * (Re[i] - 1000) * Pr[i]) / (1 + 12.7 * (f / 2)^0.5 * 
                (Pr[i]^(2 / 3) - 1))
        end

        # Film coefficient (h) and fluid resistance (Rf)
        h[i] = Nu[i] * kf / (2 * rₑ[i])
        Rf[i] = 1 / (2 * np * π * rₑ[i] * h[i]) # Fluid thermal resistance
    end

    # Convective thermal resistances for the fluid in pipe and borehole
    RCi = 1 / (Nu[1] * kf * π)          # In pipe at ri
    RCo = 1 / (Nu[2] * kf * π) * rₑ[2] / ro # In borehole at ro
    RCb = 1 / (Nu[2] * kf * π) * rₑ[2] / rb # In borehole at rb (Nu equivalent in ro and rb)

    # Pipe thermal resistance
    Rp = log(ro / ri) / (2 * π * np * kp * Hp)

    # Total borehole thermal resistance in SCW
    Rv = H / (cpf * ρf * V)             # Hydraulic resistance due to the fluid flow
    R1 = RCb
    R12 = RCi + RCo + Rp
    Rsb = R1 * (1 + (Rv^2 / (3 * R1 * R12))) # Subsurface borehole resistance in a SCW
    Rfb = Rsb - RCb                     # Fluid thermal resistance of the borehole in a SCW
    # Rfb = R12 - RCb
    return Rfb, Rv, f
end

"""
    effective_K(K, H)

Evaluate the balance point "m" and the hydraulic conductivities over and under it in a SCW caused by
recirculation flow. This is used to compute the effective length in the β-ILS model that considers
induced convergent flow in a SCW, as defined by Jacques et al. (2025).
# Arguments
    - K: A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
    - H: Total borehole depth [m]
# Output
    - K̃₁: Hydraulic conductivity above the balance point [m/s]
    - K̃₂: Hydraulic conductivity above the balance point [m/s]
    - m: Depth positon of the balance point [m]
# Reference
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function effective_K(K::AbstractMatrix{T}, H::T) where {T<:AbstractFloat}
    # 0. Constants and initialization
    nlayer = size(K, 1)                 # Number of hydrogeologic layers
    ΔH = sum(@view K[:, 1])             # Total zone thickness where the property are calculated
    Hᵢ = H - ΔH                         # Initial depth where the property are calculated
    mᵢ = ΔH / 2                         # Initial guess of balance point m
    i = 1                               # Iterator
    iₘ = 1000                           # Max number of iteration to convergence
    Δm = 1                              # Middle point convergence
    tol = 1e-6                          # Tolerance criteria before convergence in meters
    
    # Preallocation
    m = zero(T)
    Kh1 = zero(T)
    Kv1 = zero(T)
    Kh2 = zero(T)
    Kv2 = zero(T)
    K̃₁ = zero(T)
    K̃₂ = zero(T)
    zbot = zero(T)

    # 1. Iteration loop
    while i ≤ iₘ && Δm ≥ tol
        # Depth of the middle point
        zm = ΔH / mᵢ + Hᵢ # Weighted average depth of middle point

        # Initialize variables for accumulation
        H_up  = zero(T)
        H_low = zero(T)
        ΣKh1 = zero(T)
        ΣKv1 = zero(T)
        ΣKh2 = zero(T)
        ΣKv2 = zero(T)
        ztop = zero(T)

        # Loop for each layer
        @inbounds for j in 1:nlayer
            thickness = K[j, 1]
            Kval      = K[j, 2]
            zbot = ztop + thickness

            # Upper zone
            if ztop < zm
                dz = min(zbot, zm) - ztop
                if dz > 0
                    H_up  += dz
                    ΣKh1 += dz * Kval
                    ΣKv1 += dz / Kval
                end
            end

            # Lower zone
            if zbot > zm
                dz = zbot - max(ztop, zm)
                if dz > 0
                    H_low  += dz
                    ΣKh2 += dz * Kval
                    ΣKv2 += dz / Kval
                end
            end

            ztop = zbot
        end

        # Effective conductivities
        Kh1 = ΣKh1 / H_up
        Kv1 = H_up / ΣKv1
        K̃₁ = sqrt(Kh1 * Kv1)

        Kh2 = ΣKh2 / H_low
        Kv2 = H_low / ΣKv2
        K̃₂ = sqrt(Kh2 * Kv2)

        # Update balance point
        zmp = Kh2 * ΔH / (Kh1 + Kh2)
        m   = ΔH / (zmp + Hᵢ)

        Δm = abs(m - mᵢ)
        mᵢ = m
        i += 1
    end
    return m, K̃₁, K̃₂, Kh1, Kh2, Kv1, Kv2
end

"""
Validate the values of the parameters used to generate the transfer function. See 
Table 1 in Nguyen et al. (2025).
"""
function values_validity(t, ks, Cs, rb, H, V, β, B, Pe)
    if any(t -> t < 0 || t > 50 * 365 * 24 * 3600, t)
        @warn "t is out of bounds"
    end
    if ks < 1 || ks > 5
        @warn "ks is out of bounds"
    end
    if Cs < 2.1e6 || Cs > 2.6e6
        @warn "Cs is out of bounds"
    end
    if rb < 0.076 || rb > 0.102
        @warn "rb is out of bounds"
    end
    if H < 100 || H > 500
        @warn "L is out of bounds"
    end
    if any(V -> V < 54 / 60000 || V > 600 / 60000, V)
        @warn "V is out of bounds"
    end
    if any(β -> β < 0.0096 || β > 1, β)
        @warn "β is out of bounds"
    end
    if any(B -> B < 0.72 / 60000 || B > 449 / 60000, B)
        @warn "B is out of bounds"
    end
    if Pe < 0.05 || Pe > 1
        @warn "Pe is out of bounds"
    end
end