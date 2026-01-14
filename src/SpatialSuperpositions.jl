"""
Collection of functions that allows to perform spatial superposition of multiple borehole heat
exchanger (BHE) to form a ground heat exchanger (GHE). Current methods are:
    - Bloc matrix (Dusseault et al., 2018) (constant heat flux)
    - Enhanced successive flux (Nguyen and pasquier, 2021) (constant heat flux)
"""

using LinearAlgebra
using FFTW

"""
    borefield_radius(xy)

Function that computes a radius matrix, vector, unique values and indices of a borefield given
the coordinates of each borehole.
# Arguments
    - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        E.g.: [0 0] (to have a matrix input).
    - rb: Borehole radius (1x1) [m]
# Output
    - r: Radius of the borefield (1x1) [m]
"""
function borefield_radius(xy::AbstractArray{<:Real}, rb::Real)
    nb = size(xy, 1)    # Number of boreholes
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))
    r = r + Diagonal(rb * ones(nb))
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)
    return r, rᵥ, rᵤ, rᵢ, nb
end

"""
    g_matrix(t, ks, Cs, rb, H, D, xy)

Function creating a matrix named "gg" that has dimensions [nt x nr], where nt is the length of 
time steps and nr is the number of radius for the borefield. The matrix can subsequently be used
in spatial superposition, such as with the functions `bloc_matrix()` or `successive_flux()` from
this file.
# Arguments
    - t: Time vector (nt x 1) [s]
    - ks: Soil thermal conductivity (1x1) [W/mK]
    - Cs: Soil volumetric specific heat (1x1) [J/m³K]
    - rb: Borehole radius (1x1) [m]
    - H: Borehole depth (1x1) [m]
    - D: Borehole burried depth (1x1) [m]
    - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
# Output
    - gₘ: A matrix of g-function computed at different time and radius (nt x nr) [-]
"""
function g_matrix(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real,
    D::Real, xy::AbstractArray{<:Real})

    # Evaluate a matrix of radius
    ~, ~, rᵤ, ~, ~ = borefield_radius(xy, rb)

    # Compute temperature array depending on the model specified
    gₘ = Matrix{Float64}(undef, length(t), length(rᵤ))

    for i in eachindex(rᵤ)
        gₘ[:, i] = fls(t, ks, Cs, rᵤ[i], H, D)
    end    
    return gₘ
end

"""
    bloc_matrix(t, ks, Cs, rb, H, D, xy)

Function that computes the spatial superposition of a borefield using the bloc matrix approach
of Dusseault et al. (2018) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - t: Time vector (nt x 1) [s]
    - ks: Soil thermal conductivity (1x1) [W/mK]
    - Cs: Soil volumetric specific heat (1x1) [J/m³K]
    - rb: Borehole radius (1x1) [m]
    - H: Borehole depth (1x1) [m]
    - D: Borehole burried depth (1x1) [m]
    - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        E.g.: [0 0] (to have a matrix input).
# Output
    - g: g-function of the borefield spatial superposition [-]
# Reference
    Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient 
    g-function construction. Renewable Energy, 121, 249–260. 
    https://doi.org/10.1016/j.renene.2017.12.092
"""
function bloc_matrix(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real,
    D::Real, xy::AbstractArray{<:Real})
    # Basic parameters
    nt = length(t)
    ~, rᵥ, rᵤ, rᵢ, nb = borefield_radius(xy, rb)

    # Compute the ground model for all different radius of the borefield
    gₘ = Matrix{Float64}(undef, length(t), length(rᵤ))
    for i in eachindex(rᵤ)
        gₘ[:, i] = fls(t, ks, Cs, rᵤ[i], H, D)
    end
    # gₘ2 = g_matrix(t, ks, Cs, rb ,H, D, xy)

    # Compute temperature array depending on the model specified
    gg = Matrix{Float64}(undef, length(rᵥ), nt)
    gg = gₘ[:, rᵢ]'

    # Building the convolution matrix
    gᵢ = zeros(nt, nb, nb)
    G = zeros(nb * nt, nb * nt)

    for i in 1:nt
        for j in 1:nt
            if i <= j
                gᵢ[j, :, :] = reshape(gg[:, j], (1, nb, nb))
                # gᵢ[j, :, :] = reshape(gₘ[:, j], (1, nb, nb))
            end
        end

        for ii in 1:nb
            for jj in 1:nb
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gᵢ[i:end, ii, jj]
            end
        end
    end

    # Create inputs to solve the linear system
    Gₕ = [[G; repeat(I(nt), 1, nb)] [repeat(I(nt), nb, 1); zeros(nt, nt)]]
    b = zeros((nb + 1) * nt)
    b[nb*nt+1] = 1
    g = similar(b)

    # Solve the linear system
    g = Gₕ \ b

    # Output the transfer function
    return -g[nb*nt+1:end]
end

"""
    successive_flux(t, ks, Cs, rb, H, D, xy)

Iteratively solve spatial superposition for a borefield using the successive flux approach of
Nguyen and Pasquier (2021) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - t: Time vector (nt x 1) [s]
    - ks: Soil thermal conductivity (1x1) [W/mK]
    - Cs: Soil volumetric specific heat (1x1) [J/m³K]
    - rb: Borehole radius (1x1) [m]
    - H: Borehole depth (1x1) [m]
    - D: Borehole burried depth (1x1) [m]
    - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
# Output
    - gₛ: g-function of the borefield spatial superposition [-]
# Reference
    Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function 
    construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. 
    https://doi.org/10.1016/j.renene.2020.10.074
"""
function successive_flux(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real,
    H::Real, D::Real, xy::AbstractArray{<:Real})
    # Basic parameters
    nt = length(t)      # Number of time steps
    pad = 2 * nt - 1    # Padding for FFT
    r, rᵥ, rᵤ, rᵢ, nb = borefield_radius(xy, rb)

    # Creating an interpolator if there is either too many time steps or boreholes
    # nt > 100 ? t_ = exp10.(range(log10(60.), log10(3600. * 24*365*100), length=100)) : t_ = t
    # nb > 100 ? r_ = exp10.(range(log10(rb), log10(100), 100)) : r_ = rᵤ
    
    # Compute the ground model for all different radius of the borefield
    gₘ = Matrix{Float64}(undef, length(t), length(rᵤ))
    for i in eachindex(rᵤ)
        gₘ[:, i] = fls(t, ks, Cs, rᵤ[i], H, D)
    end
    #g_int = Interpolator(t, gₘ)
    #gₘ = g_matrix(t, ks, Cs, rb ,H, D, xy)

    # Setup matrices for successive flux
    #Temp = gFLS_int(repeat(rr, inner=(1, nt_hour)), repeat(tspan_hour, outer=(length(rr), 1)))
    Gg = reshape(gₘ, (nt, nb, nb))
    G = permutedims(Gg, (2, 3, 1)) # vector k
    g = permutedims(Gg, (3, 1, 2)) # vector j

    # Call gfunc for different time spans
    n = size(G, 1)

    # Initialize the matrices
    GG = zeros(n + 1, n + 1, nt)
    GG[1:n, 1:n, :] .= G
    GG[n+1, 1:n, :] .= 1
    GG[1:n, n+1, :] .= 1

    # Solve the linear system
    b = [zeros(n); 1]
    x = zeros(nt, n + 1)
    for i in 1:nt
        x[i, :] = GG[:, :, i] \ b
    end

    # Setup convolution
    # g_pad = zeros(pad)
    # g_pad[1:nt] .= gₘ
    g_fft = rfft([g, pad], 2) # Check if I have to use g or gₘ!!!

    # Set convergence criteria
    err_1 = 10.0
    err_2 = Inf
    rel_err = Inf
    k = 0

    # Iteratively solve the system
    while rel_err > 0.15 && err_1 > 1e-3 && err_1 < err_2
        f = [x[1, 1:n]'; diff(x[:, 1:n], dims=1)']
        h = irfft(rfft([f, pad], 2) .* g_fft, 2)
        hh = sum(h, dims=1)

        Tf = hh[1:nt, :]
        mTf = sum(x[:, 1:n] .* Tf, dims=2)
        CC = (Tf .- mTf) ./ mTf
        x[:, 1:n] .-= CC .* x[:, 1:n]

        k += 1
        err = maximum(abs.(CC))

        rel_err = (err_1 - err) / err_1
        err_2 = err_1
        err_1 = err
    end

    x[:, n + 1] .= mTf
    return x
end