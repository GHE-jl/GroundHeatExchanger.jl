using LinearAlgebra
using DSP: conv
includet("ground_models/finite_line_source.jl")
includet("utils.jl")

"""
    bloc_matrix(g)
    bloc_matrix(t, H, rb, D, ks, Cs, xy)

Function that computes the spatial superposition of a borefield using the bloc matrix approach
of Dusseault et al. (2018) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
    - `t`: Time vector (nt x 1) [s]
    - `H`: Borehole depth (1x1) [m]
    - `rb`: Borehole radius (1x1) [m]
    - `D`: Borehole burried depth (1x1) [m]
    - `ks`: Soil thermal conductivity (1x1) [W/mK]
    - `Cs`: Soil volumetric specific heat (1x1) [J/m³K]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        - Can be computed with `borefield_xy()` from utils.jl.
# Output
    - `g`: g-function of the borefield spatial superposition [-]
# Reference
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient 
        g-function construction. Renewable Energy, 121, 249–260. 
        https://doi.org/10.1016/j.renene.2017.12.092
"""
function bloc_matrix(gm::AbstractArray{<:Real})
    # Basic parameters
    nt, nb = size(gm)

    # Building the convolution matrix
    G = zeros(nt * nb, nt * nb)
    for i in 1:nt
        for ii in 1:nb
            for jj in 1:nb
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gm[i:end, ii, jj]
            end
        end
    end

    # Create inputs to solve the linear system
    Gₕ = [[G; repeat(I(nt), 1, nb)] [repeat(I(nt), nb, 1); zeros(nt, nt)]]
    b = zeros((nb + 1) * nt)
    b[nb*nt+1] = 1
    g = similar(b)

    # Solve the linear system
    sol = Gₕ \ b
    # Output the transfer function
    return -sol[nb*nt+1:end]
end
function bloc_matrix(t, H, rb, D, ks, Cs, xy)
    # Compute the radius matrix of the borefield
    r, _, _, _, _, _ = borefield_radius(xy, rb)

    # Compute the ground model for all different radius of the borefield
    g = fls(t, H, r, D, ks, Cs) # TODO: at some point add more models.

    # Call the bloc matrix function
    return bloc_matrix(g)
end

"""
    successive flux(g)
    successive_flux(t, H, rb, D, ks, Cs, xy)

Iteratively solve spatial superposition for a borefield using the successive flux approach of
Nguyen and Pasquier (2021) to obtain the g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
    - `t`: Time vector (nt x 1) [s]
    - `H`: Borehole depth (1x1) [m]
    - `rb`: Borehole radius (1x1) [m]
    - `D`: Borehole burried depth (1x1) [m]
    - `ks`: Soil thermal conductivity (1x1) [W/mK]
    - `Cs`: Soil volumetric specific heat (1x1) [J/m³K]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        - Can be computed with `borefield_xy()` from utils.jl.
# Output
    - `g`: g-function of the borefield spatial superposition [-]
# Reference
    - Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function 
        construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. 
        https://doi.org/10.1016/j.renene.2020.10.074
"""
function successive_flux(g::AbstractArray{<:Real,3})
    # Basic parameters
    nt, nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    nb = nb1

    # First estimation of g-function using block matrix (Eq. 20)
    GG = zeros(eltype(g), nt, nb + 1, nb + 1)
    @views GG[:, 1:nb, 1:nb] .= g
    @views GG[:, nb + 1, 1:nb] .= 1
    @views GG[:, 1:nb, nb + 1] .= 1

    b = zeros(eltype(g), nb + 1)
    b[end] = 1

    x = zeros(eltype(g), nt, nb)
    for it in 1:nt
        sol = GG[it, :, :] \ b
        x[it, :] .= sol[1:nb]
    end

    # Successive flux estimation
    e1 = 10.0
    e2 = Inf
    e3 = Inf
    k  = 0
    kmax = 100
    gi = zeros(eltype(g), nt)

    while e3 > 0.15 && e1 > 1e-3 && e1 < e2 && k < kmax
        k += 1                                  # Iteration counter
        f = vcat(x[1:1, :], diff(x, dims=1))    # Step fluxes (Eq. 4)
        # Temperature response via pairwise convolutions (Eq. 6, 8)
        hh = zeros(eltype(g), nt, nb)
        for j in 1:nb, i in 1:nb
            hh[:, i] .+= conv(f[:, j], g[:, i, j])[1:nt]
        end
        gi = vec(sum(x .* hh, dims=2))          # Eq, 11 ĥ
        c = hh ./ gi .- 1                       # Eq. 12
        x .*= (1 .- c)                          # Eq. 16 (or 15?)
        # Convergence check
        err = maximum(abs, c)
        e3 = abs((e1 - err) / e1)
        e2 = e1
        e1 = err
    end
    return gi
end
function successive_flux(t, H, rb, D, ks, Cs, xy)
    # Compute the radius matrix of the borefield
    r, _, _, _, _, _ = borefield_radius(xy, rb)

    # Compute the ground model for all different radius of the borefield
    g = fls(t, H, r, D, ks, Cs) # TODO: at some point add more models.

    # Call the bloc matrix function
    return successive_flux(g)
end

function successive_flux_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real,
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

function bloc_matrix_old(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real,
    H::Real, D::Real, xy::AbstractArray{<:Real})
    # Basic parameters
    nt = length(t)
    _, rᵥ, rᵤ, rᵢ, _, nb = borefield_radius(xy, rb)

    # Compute the ground model for unique radius of the borefield
    gₘ = Matrix{Float64}(undef, nt, length(rᵤ))
    for i in eachindex(rᵤ)
        gₘ[:, i] = fls(t, ks, Cs, rᵤ[i], H, D)
    end
    # gₘ2 = g_matrix(t, ks, Cs, rb ,H, D, xy)

    # Create matrix of g-function for all radius of the borefield
    gg = Matrix{Float64}(undef, nt, length(rᵥ))     # Create matrix for all radius of the borefield
    gg = gₘ[:, rᵢ]                                  # Fill depending on the indices of unique radius

    # Building the convolution matrix
    gᵢ = zeros(nt, nb, nb)
    G = zeros(nb * nt, nb * nt)

    for i in 1:nt
        for j in 1:nt
            if i <= j
                gᵢ[j, :, :] = reshape(gg[j, :], (1, nb, nb)) # Fill a 3D matrix of g-functions
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