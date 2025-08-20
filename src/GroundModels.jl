"""
Collection of thermal transfer models to simulate a ground heat exchanger (GHE). The models
available are:
    - Infinite line source of Ingersol (1948) (ils)
    - Infinite cylindrical source of Ingersol (1959) (ics)
    - Finite line source of Claesson and Javed (2011) (fls)
    - The standing column well analytical model of Nguyen et al. (2025) (ilsβ)
"""

using SpecialFunctions
using QuadGK

function ils(t::Union{T, AbstractVector{T}}, ks::T, Cs::T, rb::T) where {T <: AbstractFloat}
    """
        ils(t, ks, Cs, rb)

    Compute the infinite line source (ILS) model based on Ingersol (1954). The output is a 
    g-function that requires a heat load per unit of borehole length [W/m] to provide the 
    borehole wall temperature.
    Inputs:
        - t: Time vector [s]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - rb: Borehole radius [m]
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the borehole [-]
    Reference:
        Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. 
        ASHVE Journal Section, Heating, Piping and Air Conditioning.
    """

    g = -expinti.(-rb^2 ./ (4 * (ks / Cs) * t)) ./ (4 * π * ks)
    return g
end

function ics(t::Union{T, AbstractVector{T}}, ks::T, Cs::T, rb::T,
     rm::T = rb) where {T <: AbstractFloat}
    """
        ics(t, ks, Cs, rb, rm=rb)

    Computes the infinite cylindre source (ICS) model based on Carlsaw and Jaeger (1959). The output
    is a g-function that requires a heat load per unit of borehole length [W/m] to provide the
    borehole wall temperature.
    Inputs:
        - t: Time vector [s]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - rb: Borehole radius [m]
        - rm: (optional) Radius where the model is evaluated (usually equal to rb) [m]
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the borehole [-]
    Reference:
        Ingersoll, L. R., Zabel, O. J., Ingersoll, A. C., & others. (1954). Heat conduction with 
        engineering, geological, and other applications. University of Wisconsin Press.
    """

    # Set initial parameters
    nt = length(t)                      # Number of time step
    r̃ = float(rm / rb)                  # Ratio of location to the temperature and cylinder radius
    t̃ = t .* ks ./ (Cs * rb^2)          # Fourier number
    g = Vector{Float64}(undef, nt)      # Preallocation

    function integrand_ics(s::T, r̃::T, tᵢ::T)
        """
            integrand_ics(s, r̃, tᵢ)

        Integrand of the ICS model for scalar call.
        """
        if s < 1e-12
            return 0.0
        end
        return (exp(-s^2 * tᵢ) - 1) * (((besselj0(r̃ * s) * bessely1(s)) -
            (bessely0(r̃ * s) * besselj1(s))) / (s^2 * (besselj1(s)^2 + bessely1(s)^2)))
    end

    for (i, tᵢ) in enumerate(t̃)
        integral, _ = quadgk(s -> integrand_ics(s, r̃, tᵢ), 1e-6, Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (π^2 * ks)
end

function fls(t::Union{T, AbstractVector{T}}, ks::T, Cs::T, rb::T, H::T, D::T
    ) where {T <: AbstractFloat}
    """
        fls(t, ks, Cs, rb. H, D)

    Computes the finite line source (FLS) model based on Claesson and Javed (2011). The output is a 
    g-function that requires a heat load per unit of borehole length [W/m] to provide the borehole
    wall temperature.
    Inputs:
        - t: Time vector [s]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - rb: Borehole radius [m]
        - H: Borehole depth [m]
        - D: Buried depth [m]
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the borehole [-]
    Reference:
        Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid 
        temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 
        279–288.
    """

    # Set initial parameters
    const_π = 1 / sqrt(π)
    nt = length(t)              # Number of element in the time vector
    g = zeros(nt)               # Preallocation of the borehole wall temperature
    α = ks / Cs
    lim_int = 1 ./ sqrt.(4 * α * t)

    function integrand_fls(s::T, r::T, H::T, D::T) where {T <: AbstractFloat}
        """
            integrand_fls(s, r, H, D)

        Integrand of the FLS model. Assumes constant heat flux boundary condition.
        """
        function ierf(x::T) where {T <: AbstractFloat}
            """
                ierf(x)

            Inverse "erf" function used in the FLS model.
            """
            return x * erf(x) - const_π * (1 - exp(-x^2))
        end

        return exp(-r^2 * s^2) * (2 * ierf(H * s) + 2 * ierf(H * s + 2 * D * s) -
                ierf(2 * H * s + 2 * D * s) - ierf(2 * D * s)) / (H * s^2)
    end

    # Compute, in a loop, each value of the fls
    for i in 1:nt
        integral, _ = quadgk(s -> integrand_fls(s, rb, H, D), lim_int[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g / (4 * π * ks)
end

function ilsβ(t::Union{T, AbstractVector{T}}, ks::T, Cs::T, rb::T, H::T, V::T, β::T
    ) where {T <: AbstractFloat}
    """
        ilsβ(t, ks, Cs, rb, H, V, β)

    Compute the analytical model developped by Nguyen et al. (2025) for standing column
    wells (SCW). The output is a g-function that requires a heat load per unit of borehole 
    length [W/m] to provide the borehole wall temperature.
    Inputs:
        - t: Time vector [s]
        - V: Circulating fluid flow rate ([m³/s])
        - β: Fluid bleed rate ratio (between 0 and 1) [-]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - H: Borehole depth [m]
        - rb: Borehole radius [m]
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the SCW [-]
    Reference:
        Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
        standing column wells operating with bleed. Applied Thermal Engineering.
        https://doi.org/10.1016/j.applthermaleng.2024.124543
    """

    function corr_coef(Pe, dim)
        """
        Function to compute ϵ 1 to 4 proposed by Nguyen et al. (2025). The input "dim"
        indicates if the coefficients are taken from the 1D or 2D validation.
        """
        if dim == 1
            coefs₁ = [-0.2262, 0.4955, 0.4950, 0.1284]
            coefs₂ = [0.6390, -0.0001, 0.7417, 0.5631]
            coefs₃ = [0.9920, 0.5046, 0.3947, -0.0009]
        elseif dim == 2
            coefs₁ = [-0.4478, 0.1288, 0.6458, 0.1483]
            coefs₂ = [0.7110, -0.7595, 0.7838, 0.5114]
            coefs₃ = [1.0011, 0.7339, 0.2383, -0.0100]
        else
            error("Input 'dim' must be either 1 or 2")
        end
        ϵ = @. coefs₁ * Pe^coefs₂ + coefs₃
        return ϵ
    end

    function values_validity(t, ks, Cs, rb, H, V, β, B, Pe)
        """
        Validate the values of the parameters used to generate the transfer function. See 
        Table 1 in Nguyen et al. (2025).
        """
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
    ϵ = corr_coef(Pe, 2)
    h = @. ϵ[1] * (1 - ϵ[2] * β) *
           (1 + erf((sqrt(Pe) / ϵ[3]) * ((1 - ta^ϵ[4]) / (ta^ϵ[4]))))

    # 5. Combine initial and scaling function
    g = convolution(diff([0; g₀]), h)
    return g
end