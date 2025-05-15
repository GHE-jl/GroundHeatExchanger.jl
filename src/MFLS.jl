using SpecialFunctions
using QuadGK


function mfls_single_borehole(
        t::Union{T,AbstractVector{T}}, ks::T, Cs::T, rb::T, H::T, D::T, vD::T) where {T <: Real}
    """
        mfl_single_borehole(t, ks, Cs, rb, H, D, vD, xy)

    Compute the moving finite line source (MFLS) model of Guo et al. 2021, which integrates the 
    buried depth and groundwater flow (with direction).The output is a g-function that requires a
    heat load per unit of borehole length [W/m] to provide the borehole wall temperature.
    Inputs:
        - t: Time vector [s]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - rb: Borehole radius [m]
        - H: Borehole depth [m]
        - D: Buried depth [m]
        - vD: Uniform Darcy velocity [m/s]
        - xy: Matrix of borehole coordinates [nb x 2], where nb is the number of borehole.
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the borehole [-]
    Reference:
    Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite 
    line source model for vertical borehole heat exchangers—A new solution. Energy and Buildings, 
    214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
    """

    # Set initial parameters
    nt = length(t)                          # Number of element in the time vector
    g = zeros(nt)                           # Preallocation of the borehole wall temperature
    Foₛ = similar(g)                        # Preallocation for the Fourier number

    # Dimensionless parameters
    α = ks / Cs                             # Ground thermal diffusivity [m^2/s]
    Rᵦ = rb / H                             # Dimensionless radius
    d = D / H                               # Dimensionless buried dept
    Foₛ = 1 ./ sqrt.(4 * α * t / (H^2))     # 1 over square root of Fourier number
    U = vD * 999.7 * 4190 / Cs              # Water volumetric heat capacity at 10 degC
    Pe = U * H / α                          # Peclet number
    I₀ = besseli(0, Rᵦ * Pe / 2)            # Modified Bessel of the first kind with zero-order

    # Compute the MFLS
    for i in 1:nt
        integral, _ = quadgk(
            s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵦ^2 * s^2)) * integrand_mfls(s, d) / (s^2),
                Foₛ[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g * I₀ / (4 * π * ks)
end

function mfls_borefield_I(
        t::Union{T,AbstractVector{T}}, ks::T, Cs::T, rb::T, H::T, D::T, vD::T,
        xy::Matrix) where {T <: Real}
    """
        mfls_borefield_I(t, ks, Cs, rb, H, D, vD, xy)

    Compute the moving finite line source (MFLS) model of Guo et al. 2021, which integrates the 
    buried depth and groundwater flow (with direction).The output is a g-function that requires a
    heat load per unit of borehole length [W/m] to provide the borehole wall temperature.
    Inputs:
        - t: Time vector [s]
        - ks: Ground thermal conductivity [W/mK]
        - Cs: Ground volumetric specific heat [J/m³K]
        - rb: Borehole radius [m]
        - H: Borehole depth [m]
        - D: Buried depth [m]
        - vD: Uniform Darcy velocity [m/s]
        - xy: Matrix of borehole coordinates [nb x 2], where nb is the number of borehole.
    Output:
        - g: A g-function corresponding to the borehole wall temperature of the borehole [-]
    Reference:
    Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2021). Considering buried depth for vertical borehole 
    heat exchangers in a borehole field with groundwater flow—An extended solution. Energy and 
    Buildings, 235, 110722. https://doi.org/10.1016/j.enbuild.2021.110722
    """

    # Set initial parameters
    nt = length(t)                          # Number of element in the time vector
    nb = size(xy, 1)                        # Number of boreholes
    integral_s = zeros(nt)                  # Preallocation of integral computation
    integral_m = similar(integral_s)        # Preallocation of integral computation
    gᵢⱼ = zeros(nt)                         # Preallocation of MFLS at r>rb
    Foₛ = similar(integral_s)               # Preallocation for the Fourier number

    # Dimensionless parameters
    α = ks / Cs                             # Ground thermal diffusivity [m^2/s]
    Rᵦ = rb / H                             # Dimensionless radius
    d = D / H                               # Dimensionless buried dept
    Foₛ = 1 ./ sqrt.(4 * α * t / (H^2))     # 1 over square root of Fourier number
    U = vD * 999.7 * 4190 / Cs              # Water volumetric heat capacity at 10 degC
    Pe = U * H / α                          # Peclet number
    I₀ = besseli(0, Rᵦ * Pe / 2)            # Modified Bessel of the first kind with zero-order

    # Compute ΔX and R and ϕ for the spatial superposition
    
    ΔX = zeros(nb, nb)
    R = similar(ΔX)
    ϕ = similar(ΔX)

    for i in 1:nb
        ΔX[:, i] = abs.(xy[:, 1] .- xy[i, 1]) ./ H
        R[:, i] = sqrt.((xy[:, 1] .- xy[i, 1]) .^ 2 .+ (xy[:, 2] .- xy[i, 2]) .^ 2) ./ H
        ϕ[:, i] = atan.(xy[:, 1] .- xy[i, 1], xy[:, 2] .- xy[i, 2])
    end

    ΔXᵥ = reshape(ΔX, nb * nb)
    #ΔXᵤ = unique(ΔXᵥ)
    #ΔXᵢ = indexin(ΔXᵥ, ΔXᵤ)

    Rᵥ = reshape(R, nb * nb)
    #Rᵤ = unique(Rᵥ)
    #Rᵢ = indexin(Rᵥ, Rᵤ)

    ϕᵥ = reshape(ϕ, nb * nb)
    #ϕᵤ = unique(ϕᵥ)
    #ϕᵢ = indexin(ϕᵥ, ϕᵤ)

    # Compute the MFLS at the borehole i
    for i in 1:nt
        temp, _ = quadgk(
            s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵦ^2 * s^2)) *
                 integrand_mfls(s, d) / (s^2),
            Foₛ[i], Inf, rtol = 1e-6)
        integral_s[i] = temp
    end
    gᵢ = integral_s * I₀ / (4 * π * ks)

    # Compute the MFLS at all distances of the borefield
    for j in 1:(nb * nb)
        if Rᵥ[j] != 0
        for i in 1:nt
            temp, _ = quadgk(
                s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵥ[j]^2 * s^2)) * integrand_mfls(s, d) /
                     (s^2),
                Foₛ[i], Inf, rtol = 1e-6)
            integral_m[i] = exp(ΔXᵥ[j] * Pe^2 * Rᵦ * cos(ϕ[j]) / 4) .* temp / (4 * π * ks)
        end
        gᵢⱼ = gᵢⱼ .+ integral_m[:]
        end
    end

    return (gᵢ .+ gᵢⱼ) ./ (nb)
end

function ierf(x::T) where {T <: AbstractFloat}
    """
        ierf(x)
       
    Inverse "erf" function used in the FLS model
    """
    return x * erf(x) - 1 / sqrt(π) * (1 - exp(-x^2))
end

function integrand_mfls(s::T, d::T) where {T <: AbstractFloat}
    """
        integrand_fls(s, d)

    Integrand of the FLS model. Assumes constant heat flux boundary condition.
    """
    return (2 * ierf(s) + 2 * ierf(s + 2 * d * s) - ierf(2 * s + 2 * d * s) -
            ierf(2 * d * s))
end