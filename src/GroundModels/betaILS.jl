using SpecialFunctions: erf
includet("ILS.jl")
includet("../Convolutions.jl")

"""
    βils(t, ks, Cs, rb, H, V, β)

Compute the analytical model developped by Nguyen et al. (2025) for standing column
wells (SCW). The output is a g-function that requires a heat load per unit of borehole 
length [W/m] to provide the borehole wall temperature.
# Arguments
    - t: Time vector [s]
    - V: Circulating fluid flow rate ([m³/s])
    - β: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - ks: Ground thermal conductivity [W/mK]
    - Cs: Ground volumetric specific heat [J/m³K]
    - H: Borehole depth [m]
    - rb: Borehole radius [m]
# Output
    - g: A g-function corresponding to the borehole wall temperature of the SCW [°Cm/W]
# Reference
    Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
    standing column wells operating with bleed. Applied Thermal Engineering.
    https://doi.org/10.1016/j.applthermaleng.2024.124543
"""
function βils(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real,
    V::Real, β::Real)
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
    #1. Initial transfer function g_0 based on the SLI
    g₀ = ils(t, ks, Cs, rb)

    # 2. Compute the scaling function h
    B = V * β
    ta = @. B * t / (π * H * rb^2)
    Pe = B / (2 * π * H) / (ks / Cs)

    # 3. Range of validity of the inputs based on scenarios in Nguyen et al. (2025)
    values_validity(t, ks, Cs, rb, H, V, β, B, Pe)

    # 4. Compute the correction function
    coefs₁ = [-0.4478, 0.1288, 0.6458, 0.1483]
    coefs₂ = [0.7110, -0.7595, 0.7838, 0.5114]
    coefs₃ = [1.0011, 0.7339, 0.2383, -0.0100]
    ϵ = @. coefs₁ * Pe^coefs₂ + coefs₃
    h = @. ϵ[1] * (1 - ϵ[2] * β) * (1 + erf((sqrt(Pe) / ϵ[3]) * ((1 - ta^ϵ[4]) / (ta^ϵ[4]))))

    # 5. Combine initial and scaling function
    g = convolution(diff([0; g₀]), h)
    return g
end

"""
    effective_Peclet(k, K, L, b)

Effective Peclet number calculation based on thw roek of Ref[1] and Ref[2].
# Arguments
    - k: A (jx3) matrix of "j" geological layers and their thermal properties:
        - Column 1: Thickness [m]
        - Column 2: Thermal conductivity [W/mK]
        - Column 3: Volumetric specific heat [J/m³K]
        - Example: [20 k₁ C₁; 80 K₂ C₂] for 2 layers of 20m and 80m thicknesses
    - K: A (ix2) matrix of "i" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
    - L: SCW active length [m]
    - b: Radial convergent flow towards the well (bleed flow rate) [m³/s]
# Output
    - Hydrogeological anisotropiy factor (unique k and K combination, 1) [-]
# References
    - [1] Todd, D. K., & Mays, L. W. (2005). Groundwater Hydrology (Third). John Wiley & Sons.
    - [2] Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
            Thermal Engineering, 279, 127620.
"""
function effective_Peclet(k::Matrix{Real}, K::Matrix{Real}, L::Real, b::Vector{Real})
# Verify that the thicknesses match
if sum(k[:,1]) != L || sum(K[:,1]) != L
    error("The sum of the geological layers thicknesses must be equal to L.")
end

# Computation of average units and layers properties
kₐ = k[:,1]' * k[:,2] / L
Cₐ = k[:,1]' * k[:,3] / L

# Compute horizontal and vertical hydraulic conductivities
Kₕ = K[:, 2]' * K[:, 1] / L
Kᵥ = L / sum(K[:, 1] ./ K[:, 2])

# Computation of effective properties
Kₑ = sqrt(Kₕ * Kᵥ)
Lₑ = Kₑ / Kₕ * L

# Computation of Péclet numbers
Pe = B ./ (2 * pi * L) ./ (kₐ / Cₐ)       
Peₑ = B ./ (2 * pi * Lₑ) ./ (kₐ / Cₐ)

return Peₑ ./ Pe
end

"""
    Rb_SCW(r, L_mat, Vdot, k_p, T, dT, Options, Para, keps)

Computation of theoretical thermal resistance in a SCW based on Ref[1].
# Arguments
    - r: Borehole radius [m]
    - L_mat: Borehole length [m]
    - Vdot: Circulating fluid flow rate [m³/s]
    - k_p: Ground thermal conductivity [W/mK]
    - T: Time vector [s]
    - dT: Time step [s]
    - Options: A dictionary containing model options
    - Para: A dictionary containing model parameters
    - keps: Hydrogeological anisotropy factor [-]
"""
function Rb_SCW(r::Real, L_mat::Real, Vdot::Real, k_p::Real, T::Vector{Real}, dT::Real,
    Options::Dict, Para::Dict, keps::Real)

    return Rb_SCW
end