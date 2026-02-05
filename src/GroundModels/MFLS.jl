using SpecialFunctions: erf, besseli
using QuadGK: quadgk

"""
    _ierf(x)

Inverse "erf" function used in the moving finite line source model.
"""
function _ierf(x::T) where {T<:AbstractFloat}
    return x * erf(x) - inv(sqrt(T(π))) * (one(T) - exp(-x^2))
end

"""
    _mfls_integrand(s, rb, H, D)

Computes the integrand of the finite line source model.
"""
function _mfls_integrand(s::T, d::T, Rᵦ::T, Pe::T) where {T<:AbstractFloat}
    # Calculate terms
    term1 = 2 * _ierf(s)
    term2 = 2 * _ierf(s + 2 * d * s)
    term3 = _ierf(2 * s + 2 * d * s)
    term4 = _ierf(2 * d * s)
    
    return exp(-Pe^2 / (16 * s^2) - (Rᵦ^2 * s^2)) * (term1 + term2 - term3 - term4) / s^2
end

"""
    _mfls(t, ks, Cs, Cf, r, rb, H, D, vD)

Kernel function for the moving finite line source model based on Guo et al. (2020). The response 
function is based on an impulse of 1 W/m.
"""
function _mfls(t::T, ks::T, Cs::T, Cf::T, r::T, rb::T, H::T, D::T, vD::T) where {T<:AbstractFloat}
    # Initial parameters
    α = ks / Cs
    d = D / H
    U = vD * Cf / Cs
    Pe = U * H / α

    # Determine how to compute the integrand based on radius
    if r < rb
        Rᵦ = rb / H
        I = besseli(zero(T), Rᵦ * Pe / 2)
    else
        Rᵦ = r / H
        I = exp(Rᵦ * Pe / 2)
        #TODO Check, there should be an angle to compute here (cos(θ))
    end
    
    # Numerical integration
    lower_lim = inv(sqrt(4 * α * t / H^2)) # Integration limit (1/√Fo)
    integral, _ = quadgk(s -> _mfls_integrand(s, d, Rᵦ, Pe), lower_lim, T(Inf), rtol = T(1e-6))
    
    return (integral * I) / (4 * T(π) * ks)
end

"""
    mfls(t, ks, Cs, Cf, r, rb, H, D, vD)

Compute the moving finite line source (MFLS) model of Guo et al. (2020), which integrates the 
buried depth and groundwater flow (with direction). The output is a g-function that requires a
heat load per unit of borehole length [W/m] to provide the borehole wall temperature.
# Arguments
    - `t`: Time vector [s]
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cf`: Groundwater volumetric specific heat [J/m³K]
    - `r`: Radius at which to computed [m]
    - `rb`: Borehole radius [m]
    - `H`: Borehole depth [m]
    - `D`: Buried depth [m]
    - `vD`: Uniform Darcy velocity [m/s]
        - Must not be zero, set to low value for impervious (1e-12)
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the borehole [°Cm/W]
# Reference
    - Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving
        finite line source model for vertical borehole heat exchangers—A new solution. Energy and 
        Buildings, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
# Example
    mfls(60:60:3600, 3.0, 2e6, 4.2e6, 0.076, 0.076, 150, 4, 1e-6)
    mfls(60:60:3600, 3.0, 2e6, 4.2e6, 5, 0.076, 150, 4, 1e-6)
"""
function mfls(t::Real, ks::Real, Cs::Real, Cf::Real, r::Real, H::Real, D::Real, vD::Real)
    T = float(promote_type(typeof(t), typeof(ks), typeof(Cs), typeof(Cf), typeof(r), typeof(rb), 
        typeof(H), typeof(D), typeof(vD)))
    return _mfls(T(t), T(ks), T(Cs), T(Cf), T(r), T(rb), T(H), T(D), T(vD))
end

function mfls(t::AbstractVector{<:Real}, ks::Real, Cs::Real, Cf::Real, r::Real, H::Real, D::Real,
    vD::Real)
    # Check type
    T = float(promote_type(eltype(t), typeof(ks), typeof(Cs), typeof(Cf), typeof(r), typeof(rb),
        typeof(H), typeof(D), typeof(vD)))
    t_T  = convert(Vector{T}, t)
    
    # Preallocate and MILS
    g = similar(t_T)
    @inbounds @simd for i in eachindex(t_T)
        g[i] = _mfls(t_T[i], T(ks), T(Cs), T(Cf), T(r), T(rb), T(H), T(D), T(vD))
    end
    return g
end

function mfls!(g::AbstractVector{T}, t::AbstractVector, ks::Real, Cs::Real, Cf::Real, r::Real, 
    H::Real, D::Real, vD::Real) where {T<:AbstractFloat}
    # Check for same vector length
    @assert length(g) == length(t)

    # Convert parameters to T once
    ks_T, Cs_T, Cf_T, r_T, rb_T = T(ks), T(Cs), T(Cf), T(r), T(rb)
    H_T, D_T, vD_T = T(H), T(D), T(vD)
    
    # MILS
    @inbounds @simd for i in eachindex(g, t)
        g[i] = _mfls(T(t[i]), ks_T, Cs_T, Cf_T, r_T, rb_T, H_T, D_T, vD_T)
    end
    return g
end

function mfls_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real,
    H::Real, D::Real, vD::Real)
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

    # Compute the MFLS R
    for i in 1:nt
        integral, _ = quadgk(
            s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵦ^2 * s^2)) * integrand_mfls(s, d) / (s^2),
                Foₛ[i], Inf, rtol = 1e-6)
        g[i] = integral
    end
    return g * I₀ / (4 * π * ks)
end

function mfls_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real,
    D::Real, vD::Real, xy::Matrix)
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
        temp, _ = quadgk(s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵦ^2 * s^2)) * integrand_mfls(s, d) / 
            (s^2), Foₛ[i], Inf, rtol = 1e-6)
        integral_s[i] = temp
    end
    gᵢ = integral_s * I₀ / (4 * π * ks)

    # Compute the MFLS at all distances of the borefield
    for j in 1:(nb * nb)
        if Rᵥ[j] != 0
        for i in 1:nt
            temp, _ = quadgk(s -> exp(-Pe^2 / (16 * (s^2)) - (Rᵥ[j]^2 * s^2)) * integrand_mfls(s, d)
             / (s^2), Foₛ[i], Inf, rtol = 1e-6)
            integral_m[i] = exp(ΔXᵥ[j] * Pe^2 * Rᵦ * cos(ϕ[j]) / 4) .* temp / (4 * π * ks)
        end
        gᵢⱼ = gᵢⱼ .+ integral_m[:]
        end
    end

    return (gᵢ .+ gᵢⱼ) ./ (nb)
end