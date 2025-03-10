"""
Collection of functions that allows to perform spatial superposition of multiple borehole heat
exchanger (BHE) to form a ground heat exchanger (GHE). Current methods are:
    - Bloc matrix (Dusseault et al., 2018)
    - Enhanced successive flux (Nguyen and pasquier, 2021)
"""

using LinearAlgebra
using PCHIPInterpolation
using FFTW

mutable struct GHE_param
    """
    Structure including all parameters required to simulate a GHE response.
    """
    t::Vector{Float64}      # Time range (log)
    H::Float64              # Borehole depth
    D::Float64              # Borehole buried depth
    s::Float64              # Shank spacing (s/2 is the half-shank spacing)
    rb::Float64             # Borehole radius
    ri::Float64             # Pipe inlet radius
    ro::Float64             # Pipe outlet radius
    ks::Float64             # Ground thermal conductivity
    kg::Float64             # Grout thermal conductivity
    kp::Float64             # Pipe thermal conductivity
    kf::Float64             # Fluid thermal conductivity
    Cs::Float64             # Ground volumetric specific heat
    Cg::Float64             # Grout volumetric specific heat
    Cp::Float64             # Pipe volumetric specific heat
    Cf::Float64             # Fluid volumetric specific heat
    ρs::Float64             # Groud density
    ρg::Float64             # Grout density
    ρp::Float64             # Pipe density
    ρf::Float64             # Fluid density
    V::Float64              # Circulating flow rate
    vD::Float64             # Groundwater flow
end

function gfunc_matrix(params::GHE_param, xy::Matrix{T}, model::String) where T<:Real
    """
    Function creating a matrix named "gg" that has dimensions [nt x nr], where nt is the length of 
    time steps and nr is the number of radius for the borefield. The matrix can subsequently be used
    in spatial superposition.
    Inputs:
        - t: Time vector (nt x 1) [s]
        - rb: Borehole radius [m]
        - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        - model: Name of the numerical model used to generate the interpolant [string]
    Output:
        - gg: A matrix of g-function computed at different time and radius (nt x nr) [-]
    """
    # Basic parameters
    nt = length(params.t)
    nb = size(xy, 1)    # Number of boreholes

    # Evaluate a matrix of radius
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))
    r = r + Diagonal(params.rb * ones(nb))
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    #rᵢ = indexin(rᵥ, rᵤ)

    # Compute temperature array depending on the model specified
    g = Matrix{Float64}(undef, nt, length(rᵤ))

    if model == "fls" || "FLS"
        # Compute the model for each radius
        for i in eachindex(rᵤ)
            g[:, i] = fls(params.t, params.ks, params.Cs, rᵤ[i], params.H, params.D)
        end
    end
    return g
end

function bloc_matrix(params::GHE_param, g::Matrix{T}, xy::Matrix{T}) where {T<:Real}
    """
    Function that computes the spatial superposition of a borefield using the bloc matrix approach
    of Dusseault et al. (2018). The g-functions at every radius must already be computed with
    for example, the function `gfunc_matrix` from this file. This function only produce the spatial
    superposition from the matrix of g-function.
    Inputs:
        - t: Time vector (nt x 1) [s]
        - g: Matrix of g-function for all time steps and all radius (nt x nr) [-]
        - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
    Output:
        - gₛ: g-function of the borefield spatial superposition [-]
    Reference:
    Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient 
    g-function construction. Renewable Energy, 121, 249–260. 
    https://doi.org/10.1016/j.renene.2017.12.092
    """

    # Basic parameters
    nt = length(params.t)
    nb = size(xy, 1)    # Number of boreholes

    # Evaluate a matrix of radius
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))
    r = r + Diagonal(params.rb * ones(nb))
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)

    # Compute temperature array depending on the model specified
    gg = Matrix{Float64}(undef, nt, length(rᵥ))
    gg = g[:, rᵢ]'
    println(size(gg))

    # Building the convolution matrix
    gᵢ = fill(0.0, nt, nb, nb)
    G = fill(0.0, length(rᵤ) * nt, length(rᵤ) * nt)

    for i in 1:nt
        for j in 1:nt
            if i <= j
                gᵢ[j, :, :] = reshape(gg[:, j], (1, nb, nb))
            end
        end

        for ii in 1:nb
            for jj in 1:nb
                println(ii)
                println(jj)
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gᵢ[i:end, ii, jj]
            end
        end
    end

    ## Create inputs to solve the linear system
    Gₕ = [[G; repeat(I(nt), 1, nb)] [repeat(I(nt), nb, 1); zeros(nt, nt)]]
    b = zeros((nb + 1) * nt)
    b[nb*nt+1] = 1
    gₛ = similar(b)

    # Solve the linear system
    gₛ = Gₕ \ b
    gₛ = -nb * gₛ[nb*nt+1:end]
    return gₛ
end

function successive_flux(t::Vector{T}, g::Matrix{T}, rb::T, xy::Matrix{T}) where {T<:Real}
    """
    Iteratively solve spatial superposition for a borefield using the successive flux approach of
    Nguyen and Pasquier (2021). The required g-functions are computed by, for example, the function
    `gfunc_matrix` from this file. This function only produce the spatial
    superposition from the matrix of g-function.
    Inputs:
        - t: Time vector (nt x 1) [s]
        - g: Matrix of g-function for all time steps and all radius (nt x nr) [-]
        - rb: Borehole radius [m]
        - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
    Output:
        - gₛ: g-function of the borefield spatial superposition [-]
    Reference:
    Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function 
    construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. 
    https://doi.org/10.1016/j.renene.2020.10.074
    """

    # Basic parameters
    nt = length(t)
    nb = size(xy, 1)    # Number of boreholes

    # Evaluate a matrix of radius
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))
    r = r + Diagonal(rb * ones(nb))
    rᵥ = reshape(r, nb * nb)
    rᵤ = unique(rᵥ)
    rᵢ = indexin(rᵥ, rᵤ)

    g_int = Interpolator(t, g)

    # Define constants
    tfinal = tspan[end]
    nt_hour = 24
    nt_day = 30
    nt_month = 12
    nt_year = 10

    # Determine nt_decade
    if tfinal / 3600 / 8760 / 10 > 2
        nt_decade = ceil(Int, tfinal / 3600 / 8760 / 10)
    else
        nt_decade = 2
    end

    # Time spans
    tspan_hour = range(3600, step=3600, length=nt_hour)
    tspan_day = tspan_hour[end] .* range(1, step=1, length=nt_day)
    tspan_month = tspan_day[end] .* range(1, step=1, length=nt_month)
    tspan_year = tspan_month[end] .* range(1, step=1, length=nt_year)
    tspan_decade = tspan_year[end] .* range(1, step=1, length=nt_decade)

    # Combine time spans into a single array
    tspan_temp = vcat(
        tspan_hour,
        tspan_day[2:end],
        tspan_month[2:end],
        tspan_year[2:end],
        tspan_decade[2:end]
    )

    # Hour calculations
    Temp = gFLS_int(repeat(rr, inner=(1, nt_hour)), repeat(tspan_hour, outer=(length(rr), 1)))
    Gg_hour = reshape(Temp[Ind_c, :]', (nt_hour, n, n))
    G_hour = permutedims(Gg_hour, (2, 3, 1)) # vector k
    g_hour = permutedims(Gg_hour, (3, 1, 2)) # vector j

    # Day calculations
    Temp = gFLS_int(repeat(rr, inner=(1, nt_day)), repeat(tspan_day, outer=(length(rr), 1)))
    Gg_day = reshape(Temp[Ind_c, :]', (nt_day, n, n))
    G_day = permutedims(Gg_day, (2, 3, 1)) # vector k
    g_day = permutedims(Gg_day, (3, 1, 2)) # vector j

    # Month calculations
    Temp = gFLS_int(repeat(rr, inner=(1, nt_month)), repeat(tspan_month, outer=(length(rr), 1)))
    Gg_month = reshape(Temp[Ind_c, :]', (nt_month, n, n))
    G_month = permutedims(Gg_month, (2, 3, 1)) # vector k
    g_month = permutedims(Gg_month, (3, 1, 2)) # vector j

    # Year calculations
    Temp = gFLS_int(repeat(rr, inner=(1, nt_year)), repeat(tspan_year, outer=(length(rr), 1)))
    Gg_year = reshape(Temp[Ind_c, :]', (nt_year, n, n))
    G_year = permutedims(Gg_year, (2, 3, 1)) # vector k
    g_year = permutedims(Gg_year, (3, 1, 2)) # vector j

    # Decade calculations
    Temp = gFLS_int(repeat(rr, inner=(1, nt_decade)), repeat(tspan_decade, outer=(length(rr), 1)))
    Gg_decade = reshape(Temp[Ind_c, :]', (nt_year, n, n))
    G_decade = permutedims(Gg_decade, (2, 3, 1)) # vector k
    g_decade = permutedims(Gg_decade, (3, 1, 2)) # vector j

    # Call gfunc for different time spans
    x_hour, _ = gfunc(tspan_hour, G_hour, g_hour)
    x_day, _ = gfunc(tspan_day, G_day, g_day)
    x_month, _ = gfunc(tspan_month, G_month, g_month)
    x_year, _ = gfunc(tspan_year, G_year, g_year)
    x_decade, _ = gfunc(tspan_decade, G_decade, g_decade)

    # Combine results
    g_temp = vcat(
        x_hour[:, end],
        x_day[2:end, end],
        x_month[2:end, end],
        x_year[2:end, end],
        x_decade[2:end, end]
    )

    gₛ = pchip_interpolation(tspan_temp, -g_temp, tspan)
    return gₛ
end

function successibe_assembly(t, G, g)
    """
    Function that assemble the successive flux of Nguyen and Pasquier (2021). Used in the 
    `successive_flux` function from this file.
    """
    nt = length(t)
    n = size(G, 1)

    GG = zeros(n + 1, n + 1, nt)
    GG[1:n, 1:n, :] .= G
    GG[n+1, 1:n, :] .= 1
    GG[1:n, n+1, :] .= 1

    b = [zeros(n); 1]
    x = zeros(nt, n + 1)
    for i in 1:nt
        x[i, :] = GG[:, :, i] \ b
    end

    pad0 = 2 * nt
    if log2(nt) < 16
        p = ceil(Int, log2(pad0))
        pad0 = 2^p
    end
    g_fft = fft([g; pad0], 2)

    err_1 = 10.0
    err_2 = Inf
    rel_err = Inf
    k = 0

    while rel_err > 0.15 && err_1 > 1e-3 && err_1 < err_2
        f = [x[1, 1:n]'; diff(x[:, 1:n], dims=1)']
        h = ifft(fft([f; pad0], 2) .* g_fft, 2)
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

    x[:, n+1] .= -mTf

    return x, nothing  # 'l' is not used in the function, so we return nothing
end