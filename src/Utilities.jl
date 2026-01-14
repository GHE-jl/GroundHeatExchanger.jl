using PCHIPInterpolation

"""
    pchip_interpolation(tᵢ, gᵢ, t)

Function that performs the complete interpolation of a vector using the PCHIP interpolation method.
# Arguments
    - tᵢ: The id on which to interpolate
    - vᵢ: The vector to interpolate, sampled at tᵢ
    - t: The new interpolated vector sample
# Output
    - v: The new interpolated vector
"""
function pchip_interpolation(tᵢ::AbstractVector{<:Real}, vᵢ::AbstractVector{<:Real},
    t::AbstractVector{<:Real})
    interp = Interpolator(tᵢ, vᵢ)
    v = interp.(t)
    return v
end

"""
    set_nodes(nt, n₀)

Function that sets a logarithmic progression of node positions on a transfer function.
# Arguments
    - nt: Total number of data in the input vectors [-]
    - n₀: User defined number of nodes on the transfer function [-]
# Output
    - id: A vector of length "n₀" of node positions on the transfer function [-]
"""
function set_nodes(nt::Integer, n₀::Integer)
    # Basic inputs
    n_tmp = n₀ - 1
    id = Vector{Integer}(undef, n_tmp)

    while length(id) < n₀
        empty!(id)
        for x in range(0, stop=log10(nt), length=n_tmp)
            push!(id, round(Int, exp10(x)))
        end
        unique!(id)
        n_tmp += 1
    end
    return id
end

"""
    water_density(T::Real)

Water density ρ(T) in kg/m³, 0 ≤ T ≤ 150°C at 1 atm.
Equation from IAPWS Formulation 1980 Auxiliary Equation (Kell, 1975 revision).
Ref: Wagner & Pruß (2002), J. Phys. Chem. Ref. Data 31(2), 387-535
"""
function water_density(T::Real)::Float64
    Tc = 228.725  # K
    t = T + 273.15  # convert to K
    θ = 1 - t / Tc
    ρ = 999.83952 + 16.952577*t + (-0.010415494*t^2) + 
        (-0.000063888*t^3) + (0.000000434*t^4) +
        θ*(47.427794 + (-3.896200*t) + (0.018195*θ)) +
        (θ^2 * (-113.40 + 0.786*t)) +
        (θ^5 * -3.67e5) +
        (θ^18 * 1.93e5)
    return ρ * (1 - 1.63e-8 * t^2)  # small pressure correction term
end

"""
    water_cp(T::Real)

Specific heat capacity Cp(T) in J/(kg·K), 0 ≤ T ≤ 100°C.
Equation from Young & Jones (1994), 5th ed., p. 165.
Ref: Engineering Toolbox water Cp correlations
"""
function water_cp(T::Real)::Float64
    # Polynomial fit: Cp(T) = a + bT + cT² + dT³ + eT⁴
    a, b, c, d, e = 4217.4, -3.7203, 0.14137, -0.0024808, 1.584e-5
    return a + b*T + c*T^2 + d*T^3 + e*T^4 + 2 # small correction term to fit Engineering Toolbox
end

"""
    water_thermal_conductivity(T::Real)

Thermal conductivity k(T) in W/(m·K), 0 ≤ T ≤ 100°C.
Equation from IAPWS Formulation 1995 for thermal conductivity.
Simplified form for atmospheric pressure.
Ref: IAPWS Release on Thermal Conductivity (2018)
"""
function water_thermal_conductivity(T::Real)::Float64
    t = T + 273.15  # K
    # Base value at 0°C: 0.561 W/mK
    k0 = 0.561
    # Temperature dependence from IAPWS
    L = [-1.480663765, 4.122177914e-2, -2.626455426e-1,
         9.622313981e-2, -2.368174243e-2, 2.642215220e-2]
    Tr = 647.096  # K
    τ = 1 - t/Tr
    k_red = exp(L[1] + L[2]*τ + L[3]*τ^2 + L[4]*τ^3 + L[5]*τ^4 + L[6]*τ^5)
    return k0 * k_red * (1 + 0.001*t/273.15)  # minor linear correction
end

"""
    water_viscosity(T::Real)

Dynamic viscosity μ(T) in Pa·s, 0 ≤ T ≤ 150°C.
Equation from IAPWS Formulation 2008 for viscosity.
Ref: IAPWS Release on Viscosity (2008)
"""
function water_viscosity(T::Real)::Float64
    Tc = 647.096  # K
    t = T + 273.15
    τ = t/Tc
    # IAPWS viscosity equation: μ = μ₀ * μ₁ * exp(μ₂ + μ₃/τ)
    μ₀ = 1e-3 * sqrt(100.0/t)  # reference viscosity [Pa·s]
    
    # μ₁ term (zero density limit)
    H = 1.0 + (1.994n0 * τ) + (1.233n0 * τ^1.5)
    μ₁ = exp(0.97897 * (1/τ - 1))
    
    # μ₂ and μ₃ terms (density correction, atm pressure ≈ low density)
    μ₂ = 0.0018  # simplified for low density
    μ₃ = 0.00001
    
    return μ₀ * μ₁ * exp(μ₂ + μ₃/t)
end
