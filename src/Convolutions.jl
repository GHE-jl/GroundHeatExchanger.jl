using FFTW
FFTW.set_num_threads(Sys.CPU_THREADS)
using Clustering

function convolution(f::Union{T,AbstractVector{T}}, g::Union{T,AbstractVector{T}}) where {T<:Real}
    """
        convolution(f, g)
    
    Function that solves a convolution product in the spectral domain. Padding is used to
    avoid circular convolution.
    Inputs:
        - f: Incremental load functions (impulse function) (nₜ x 1)
        - g: Transfer function (or g-function) (nₜ x 1)
    Output:
        - Convolved signal
    """
    n = length(f)
    pad = 2 * n - 1

    # Preallocate arrays
    f_pad = zeros(T, pad)
    g_pad = zeros(T, pad)
    #copyto!(f_pad, f)
    #copyto!(g_pad, g)
    f_pad[1:n, :] .= f
    g_pad[1:n, :] .= g

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_fg = F_f .* F_g

    # Inverse FFT to get convolution results
    y = irfft(F_fg, pad)[1:n]
    return y
end

function step_signal(x::AbstractVector{AbstractFloat}, steps::Integer)
    """
        step_signal(x, steps)
    
    Function that separated a vector signal "x" in n number of constant steps. This applies mainly
    to help interprete noisy signal into constant values based on average abrupt changes.
    Note: The k-means algorithm used to idenfity the changes can be unstable when too few steps are
    used.
    Inputs:
        - x: A vector
        - steps: The number of constant steps wanted in the output step-constant signal
    Output:
        - A step-constant signal
    """
    # Find group of same data using a k-mean cluster algorithm from Clustering.jl
    x_mat = reshape(x, 1, :)
    result = kmeans(x_mat, steps)
    #result = kmedoids(x_mat, steps)

    # Compute mean for each cluster
    means = [mean(x[result.assignments .== k]) for k in 1:steps]

    # Assign each value its cluster's mean
    return means[result.assignments]
end

function state_transitions(vectors::AbstractArray...)
    """
        state_transitions(vectors::AbstractArray...)
    
    Function that allows to find indices and states changes for vectors of values. Specifically, 
    this function is used in non-stationary operation to identify the operating conditions for all 
    values between either one or a serie of input vectors. The function will identify position of 
    state changes and state index for a set of operating conditions.
    Inputs:
        - vectors: any number of input vector for all parameter affecting the operating conditions
    Outputs:
        - ind: Indices of state change on the vectors
        - s: State index for each segment delimited by the indices
    """
    # Check to have at least one vector
    n_vector = length(vectors)
    n_vector > 0 || throw(ArgumentError("At least one vector must be provided"))

    # Ensure all vectors have the same length
    n = length(vectors[1])
    for v in vectors
        length(v) == n || throw(ArgumentError("All vectors must have the same length"))
    end

    # Stack vectors into an n × k matrix
    mat = hcat(vectors...)

    # Find index for each row
    rows = unique(eachrow(mat))
    # index = map(r -> findfirst(isequal(r), unique(rows)), rows)
    dict = Dict{typeof(rows[1]), Int}()
    for (i, r) in enumerate(rows)
        dict[r] = i
    end
    index = [dict[r] for r in eachrow(mat)]

    # Find indices of changes and state value and indice of unique state changes for g-function
    change = diff(index) .!= 0
    #ind = findall(change) .- 1
    ind = findall(change)
    state = index[[true; change]]
    ind_unique = unique(i -> index[i], eachindex(index))
    
    return [1; ind], state, ind_unique
end

function convolution_ns(Q::AbstractVector{T}, g::AbstractArray{T}, ind::AbstractVector{<: Integer}, 
    s::AbstractVector{<: Integer}) where {T <: AbstractFloat}
    """
        convolution_ns(Q, g, ind, s)

    Function that performs a non-stationary convolution based on varying operating conditions in the
    simulation of GHE. This script follows the implementation of Beaudry et al (2024).
    Inputs:
        - Q: Thermal load vector (nt x 1) [W, W/m, °C]
        - g: Set of transfer functions (nt x ns) [°C/W, °Cm/W, -]
        - ind: Time vector of state transition (ns-1 x 1) [-]
        - s: State index (ns x 1) [-]
    Outputs:
        - T: Temperature vector (nt x 1)
    Reference:
        Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental
        validation of non-stationary convolutions for the fast simulation of time-variant flowrates
        in ground heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
        https://doi.org/10.1080/23744731.2023.2279468
    """
    # Basic inputs
    n = length(Q)
    index_count = diff([ind; n + 1])
    
    # Repeat each state index by corresponding count
    state_vec = vcat([fill(s[i], index_count[i]) for i in eachindex(s)]...)

    # Initialize Q_mu, size n × size(g,2)
    Q_s = zeros(T, n, size(g, 2))
    
    # Assign Q elements to Q_s according to state_vec
    for i in 1:size(g, 2)
        indt = findall(state_vec .== i)
        Q_s[indt, i] = Q[indt]
    end

    # Compute f_s as difference along rows, with first row the same as Q_s(1,:)
    f_s = vcat(Q_s[1:1, :], diff(Q_s, dims=1))

    # Compute zero-padding length
    pad = 2 * n - 1

    # Prepare padded arrays for fft
    f_pad = zeros(T, pad, size(g, 2))
    g_pad = zeros(T, pad, size(g, 2))

    # Copy data into padded arrays
    f_pad[1:n, :] .= f_s
    g_pad[1:size(g, 1), :] .= g
    #copyto!(f_pad, f_s)
    #copyto!(g_pad, g)
    
    # Perform FFTs along columns
    #F_f = rfft.(eachcol(f_pad))
    #F_g = rfft.(eachcol(g_pad))
    fft_plan = plan_rfft.(eachcol(f_pad))
    F_f = fft_plan .* eachcol(f_pad)
    F_g = fft_plan .* eachcol(g_pad)
    
    # Element-wise multiply FFTs and inverse FFT the results
    H_cols = map((Ff, Fg) -> irfft(Ff .* Fg, pad)[1:n], F_f, F_g)
    
    # Convert the vector of vectors back to matrix, each column is the result
    H = hcat(H_cols...)

    # Sum over the first n rows, along columns (dim=2), results in vector length n
    dT = sum(H[1:n, :], dims=2)

    # Return the result as a vector
    return vec(dT)
end