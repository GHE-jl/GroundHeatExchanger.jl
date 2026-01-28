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
    water_ρ(T::Real)

Water density ρ(T), 0 ≤ T ≤ 100°C at 1 atm. The polynomial equation is a fit to the data from the 
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - ρ: Density [kg/m³]
# Reference
    - The Engineering ToolBox (2003). Water Density, Specific Weight and Thermal Expansion 
        Coefficients - Temperature and Pressure Dependence. [online] Available at: 
        https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html 
        [Accessed 2026-01-14].
"""
function water_ρ(T::Real)::Float64
    return 999.8475436930158 + 0.06180756931966996*T - 0.008309049138917115*T^2 + 
        6.35713412865478e-5*T^3 - 3.8404497053894326e-7*T^4 + 1.0249871031879443e-9*T^5
end

"""
    water_cp(T::Real)

Water specific heat capacity cp(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the isobaric 
specific heat capacity data from the Engineering Toolbox. The values can be validates with a 
temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - cp: Specific heat capacity [J/(kg·K)]
# Reference
    - The Engineering ToolBox (2004). Specific Heat Capacity of Water: Temperature-Dependent Data 
        and Calculator. [online] Available at: 
        https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html 
        [Accessed 2026-01-14].
"""
function water_cp(T::Real)
    return 4219.849078078278 - 3.266686616602623*T + 0.09969277880041719*T^2 - 
        0.0014860911377001344*T^3 + 1.161963811563561e-5*T^4 - 3.5034316470844105e-8*T^5
end

"""
    water_k(T::Real)

Water thermal conductivity k(T), 0 ≤ T ≤ 99.6°C at 1 bar. The polynomial equation is a fit to the 
data from the Engineering Toolbox. The values can be validates with a temperature vector 
`T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - k: Thermal conductivity [W/mK]
# Reference
    - The Engineering ToolBox (2018). Thermal Conductivity of Water: Temperature and Pressure Data. 
        [online] Available at: https://www.engineeringtoolbox.com/water-liquid-gas-thermal-
        conductivity-temperature-pressure-d_2012.html [Accessed 2026-01-14].
"""
function water_k(T::Real)
    return 0.5557250521318174 + 0.002490814640452007*T - 2.1170044416971473e-5*T^2 + 
        1.285515973680875e-7*T^3 - 4.546428806458628e-10*T^4 + 9.750314739837196e-14*T^5
end

"""
    water_μ(T::Real)

Dynamic viscosity μ(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the data from the 
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - μ: Dynamic viscosity [Pa·s or kg/(m·s)]
# Reference
    - The Engineering ToolBox (2004). Water - Dynamic and Kinematic Viscosity at Various 
        Temperatures and Pressures. [online] Available at: https://www.engineeringtoolbox.com/water-
        dynamic-kinematic-viscosity-d_596.html [Accessed 2026-01-14].
"""
function water_μ(x::Real)
    return 0.001790966556989398 - 5.965082369793418e-5*x + 1.3185191782991122e-6*x^2 - 
        1.8236868027209892e-8*x^3 + 1.3644271817518522e-10*x^4 - 4.137645533574321e-13*x^5
end