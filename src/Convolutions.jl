using FFTW
FFTW.set_num_threads(Sys.CPU_THREADS)

# Dictionary to cache FFT plans by length
const FFT_PLANS = Dict{Int, FFTW.rFFTWPlan{Float64, true, 1}}()

# Helper to get or create a plan for a given length
function get_fft_plan(len::Int)
    get!(FFT_PLANS, len) do
        plan_rfft(zeros(len))
    end
end

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
    total_length = 2 * n - 1
    optimal_length = nextprod([2, 3, 5, 7], total_length)

    # Preallocate arrays
    f_pad = zeros(T, optimal_length)
    g_pad = zeros(T, optimal_length)
    copyto!(f_pad, f)
    copyto!(g_pad, g)

    # Precompute FFT plans
    #fft_plan = get_fft_plan(optimal_length)
    fft_plan = plan_rfft(f_pad)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_fg = F_f .* F_g

    # Inverse FFT to get convolution results
    y = irfft(F_fg, optimal_length)[1:n]
    return y
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
    #index = indexin(eachrow(mat), unique(eachrow(mat)))
    rows = collect(eachrow(mat))
    index = map(r -> findfirst(isequal(r), unique(rows)), rows)

    # Find indices of changes and state value
    change = diff(index) .!= 0
    ind = findall(change) .- 1
    state = index[[true; change]]
    
    return ind, state
end

function convolution_ns(Q::AbstractVector{T}, g::AbstractArray{T}, ind::AbstractVector{T}, 
    s::AbstractVector{T}) where {T <: Real}
    """
        convolution_ns(Q, g, ind, s)

    Function that performs a non-stationary convolution based on varying operating conditions in the
    simulation of GHE. This script follows the implementation of Beaudry et al (2024).
    Inputs:
        - Q: Thermal load vector (nₜ x 1) [W or W/m or °C]
        - g: Set of transfer functions (nₜ x nₛ) [°C/W °Cm/W, -]
        - ind: Time vector of state transition (nₜ x 1) [-]
        - s: State index (nₜ x 1) [-]
    Outputs:
        - T: Temperature vector (nₜ x 1)
    Reference:
        Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental
        validation of non-stationary convolutions for the fast simulation of time-variant flowrates
        in ground heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
        https://doi.org/10.1080/23744731.2023.2279468
    """
    # Basic inputs
    n = length(Q)
    index_count = diff([0; ind; n])
    
    # Repeat each state index by corresponding count
    state_vec = repeat(s, inner=index_count)
    
    
    # Initialize Q_mu, size n × size(g,2)
    Q_mu = zeros(eltype(Q), n, size(g, 2))
    
    # Assign Q elements to Q_mu according to state_vec
    for i in 1:size(g, 2)
        inds = findall(state_vec .== i)
        Q_mu[inds, i] = Q[inds]
    end

    # Compute f_mu as difference along rows, with first row the same as Q_mu(1,:)
    f_mu = vcat(Q_mu[1:1, :], diff(Q_mu, dims=1))

    # Compute zero-padding length (twice n, rounded up to nextprod of [2,3,5,7])
    pad0 = 2 * n
    pad0_opt = nextprod([2, 3, 5, 7], pad0)

    # Prepare padded arrays for fft
    f_pad = zeros(eltype(Q_mu), pad0_opt, size(g, 2))
    g_pad = zeros(eltype(g), pad0_opt, size(g, 2))

    # Copy data into padded arrays
    f_pad[1:n, :] .= f_mu
    g_pad[1:size(g, 1), :] .= g
    
    # Perform FFTs along columns (dim 1)
    F_f = fft.(eachcol(f_pad))
    F_g = fft.(eachcol(g_pad))
    
    # Element-wise multiply FFTs and inverse FFT the results
    H_cols = map((Ff, Fg) -> ifft(Ff .* Fg), F_f, F_g)
    
    # Convert the vector of vectors back to matrix, each column is the result
    H = hcat(H_cols...)

    # Sum over the first n rows, along columns (dim=2), results in vector length n
    dT_NSC_COR = sum(H[1:n, :], dims=2)

    # Return the result as a vector
    return vec(dT_NSC_COR)
end