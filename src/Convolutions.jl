using FFTW
FFTW.set_num_threads(Sys.CPU_THREADS)

"""
    convolution(f, g)
    convolution(Q, g)

Performs a temporal superposition through a convolution product solved in the spectral domain using
fast fourier transform. Padding is used to avoid circular convolution.
# Arguments
    - `f`: Incremental thermal load function (nₜ x 1) [W, W/m, °C] (incremental_func())
    - `Q`: Load functions (nₜ x 1) [W, W/m, °C]
    - `g`: Transfer function (or g-function) (nₜ x 1) [°C/W, °Cm/W, -]
# Output
    - Convolved signal. Typically temperature response (nₜ x 1) [°C].
# Reference
    - Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for 
        geothermal ground-loop heat exchanger systems. Geothermics, 37(6), 651–665.
        https://doi.org/10.1016/j.geothermics.2008.08.003
    - Pasquier, P., & Marcotte, D. (2013). Efficient computation of heat flux signals to ensure the
        reproduction of prescribed temperatures at several interacting heat sources. Applied Thermal
        Engineering, 59(1–2), 515–526. https://doi.org/10.1016/j.applthermaleng.2013.06.018
"""
function convolution(f::AbstractVector{T}, g::AbstractVector{T}) where {T<:AbstractFloat}
    # Ensure f and g are vectors of the same length
    n = length(f)
    n == length(g) || throw(ArgumentError("Vectors f and g must have the same length"))

    # Preallocate and replace arrays
    pad = 2 * n - 1
    f_pad = zeros(T, pad)
    g_pad = zeros(T, pad)
    f_pad[1:n] .= f     # If ever, I changed from f_pad[1:n, :]
    g_pad[1:n] .= g     # If ever, I changed from f_pad[1:n, :]

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Inverse FFT to get convolution results
    y = irfft(F_f .* F_g, pad)[1:n]
    return y
end
function convolution(Q::AbstractVector{T}, g::AbstractVector{T}) where {T<:AbstractFloat}
    f = incremental_func(Q)
    return convolution(f, g)
end

"""
    incremental_func(Q)

Function to compute the incremental function for convolution.
# Arguments
    - `Q`: Load functions (nₜ x 1) [W, W/m, °C]
# Output
    - `f`: Incremental thermal load function (nₜ x 1) [W, W/m, °C]
"""
function incremental_func(Q::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(Q)
    f = zeros(T, n)
    f[1] = Q[1]
    @inbounds for t in 2:n
        f[t] = Q[t] - Q[t-1]
    end
    return f
end

"""
    convolution_ns(f, g, state_vec)
    convolution_ns(f, g, ind, state)
    convolution_ns(Q, g, state_vec)
    convolution_ns(Q, g, ind, state)

Function that performs temporal superposition for varying operating conditions through a
non-stationary convolution solved in the spectral domain. The varying operating conditions are
defined by the `state_vec`, which assigns each time step to a specific state. Each state corresponds
to a specific column in the transfer function matrix `g`. Padding is used to avoid circular
convolution.
# Arguments
    - `f`: Incremental thermal load function (nₜ x nₛ) [W, W/m, °C] (incremental_func_ns())
    - `Q`: Load functions (nₜ x 1) [W, W/m, °C]
    - `g`: Matrix of transfer functions (or g-function) (nₜ x nₛ) [°C/W, °Cm/W, -]
        - `nₛ`: Number of states (operating conditions).
    - `state_vec`: State index for each time step (nₜ x 1) [-], with values from 1 to nₛ.
    - `ind`: Indices of state change on the vectors [-]
    - `state`: State index for each segment delimited by the indices [-]
# Output
    - Convolved signal. Typically temperature response (nₜ x 1) [°C].
# Reference
    - Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental
        validation of non-stationary convolutions for the fast simulation of time-variant flowrates
        in ground heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
        https://doi.org/10.1080/23744731.2023.2279468
"""
function convolution_ns(f::AbstractMatrix{T}, g::AbstractMatrix{T},
     state_vec::AbstractVector{<:Integer}) where {T<:AbstractFloat}
    # Check if all vectors have the same length
    n, ns = size(f)                      # Number of time steps and states
    n == length(state_vec) || throw(ArgumentError("state_vec length mismatch"))
    n == size(g, 1) || throw(ArgumentError("g row mismatch"))
    ns == size(g, 2) || throw(ArgumentError("g column mismatch"))

    # Setup the loop for non-stationary convolution in the spectral domain
    pad = 2 * n - 1                     # Padding length for convolution
    dT = zeros(T, n)                    # Initialize output vector with zeros
    f_mask = zeros(T, pad)              # Padded input vector
    g_pad  = zeros(T, pad)              # Padded transfer function vector
    fft_plan = plan_rfft(f_mask)        # FFT plan for input vector

    # Loop over states and perform convolution for each state
    for i in 1:ns
        # Build f
        fill!(f_mask, zero(T))
        @views f_mask[1:n] .= f[:, i]

        # Build g
        fill!(g_pad, zero(T))
        @views g_pad[1:n] .= g[:, i]

        # Convolution in the spectral domain (multiplication)
        Ff = fft_plan * f_mask
        Fg = fft_plan * g_pad

        # Inverse FFT to get convolution results for this state
        H = irfft(Ff .* Fg, pad)

        @inbounds for t in 1:n
            dT[t] += H[t]
        end
    end
    return dT
end
function convolution_ns(f::AbstractMatrix{T}, g::AbstractArray{T}, ind::AbstractVector{<:Integer},
    state::AbstractVector{<:Integer}) where {T<:AbstractFloat}
    # Convert indices and states to state vector
    state_vec = state_vector(ind, state, length(Q))

    # Call internal convolution function
    return convolution_ns(f, g, state_vec)
end
function convolution_ns(Q::AbstractVector{T}, g::AbstractArray{T},
    state_vec::AbstractVector{<:Integer}) where {T<:AbstractFloat}
    # Compute the incremental function from the load vector Q and state vector
    f = impulse_func_ns(Q, state_vec)

    # Call internal convolution function
    return convolution_ns(f, g, state_vec)
end
function convolution_ns(Q::AbstractVector{T}, g::AbstractArray{T}, ind::AbstractVector{<:Integer},
    state::AbstractVector{<:Integer}) where {T<:AbstractFloat}
    # Convert indices and states to state vector
    state_vec = state_vector(ind, state, length(Q))

    # Compute the incremental function from the load vector Q and state vector
    f = impulse_func_ns(Q, state_vec)

    # Call internal convolution function
    return convolution_ns(f, g, state_vec)
end

"""
    impulse_func_ns(Q, state_vec)

Function to compute the incremental function for non-stationary convolution. The incremental 
function is computed for each state separately, by masking the load vector Q with the state vector.
# Arguments
    - `Q`: Load functions (nₜ x 1) [W, W/m, °C]
    - `state_vec`: State index for each time step (nₜ x 1) [-], with values from 1 to nₛ.
# Output
    - `f`: impulse thermal load response (nₜ x nₛ) [W, W/m, °C]
"""
function impulse_func_ns(Q::AbstractVector{T}, state_vec::AbstractVector{<:Integer}
    ) where {T<:AbstractFloat}
    n = length(Q)
    ns = maximum(state_vec)
    fs = zeros(T, n, ns)
    Q_mask = zeros(T, n)

    for i in 1:ns
        @inbounds for t in 1:n
            Q_mask[t] = (state_vec[t] == i) ? Q[t] : zero(T)
        end

        fs[1,i] = Q_mask[1]
        @inbounds for t in 2:n
            fs[t,i] = Q_mask[t] - Q_mask[t-1]
        end
    end
    return fs
end

"""
    state_vector(ind, state, n)
    state_vector(vectors)

Helper function to create a state vector of length `n` from the indices of state changes (ind) and 
state values (state).
# Arguments
    - `ind`: Indices of state change on the vectors
    - `state`: State index for each segment delimited by the indices
    - `n`: Length of the output state vector
    - `vectors`: any number of input vector for all parameter affecting the operating conditions
        - E.g., `state_vector(x, y, z)`, were the vectors `x`, `y`, `z` are the same legnth.
# Output
    - `state_vec`: State index for each time step (n x 1) [-]
# Example
    state_vector([1, 4, 6], [1, 2, 3], 6) # Output: [1, 1, 1, 2, 2, 3]
"""
function state_vector(ind, state, n)
    index_count = diff([ind; n + 1])
    state_vec = vcat([fill(state[i], index_count[i]) for i in eachindex(state)]...)
    return state_vec
end
function state_vector(vectors::AbstractArray...)
    ind, state, _ = state_indices(vectors...)
    return state_vector(ind, state, length(vectors[1]))
end

"""
    state_indices(vectors...)
    state_indices(matrix)

Function that allows to find indices and states changes for vectors of values. 
This function is used as a helper function in non-stationary convolution to identify the operating
conditions at all time steps for values between either one or a serie of input vectors. The function
 will identify the location of state changes and state index for a set of operating conditions.
The indices of state changes correspond to the last time step of a given state. To have the value
of the state change at the first time step of the new state, use the `ind_value` output.
# Arguments
    - `vectors...`: any number of input vector for all parameter affecting the operating conditions
        - E.g., `state_transitions(x, y, z)`, were the vectors `x`, `y`, `z` are the same legnth.
    - `matrix`: A matrix where each column is a parameter affecting the operating conditions.
# Outputs
    - `ind`: Indices of state change on the vectors
    - `state`: State index for each segment delimited by the indices
    - `ind_unique`: Indices of unique state values (if required)
# Example
    state_indices([1, 1, 1, 2, 2, 3]) # Output: ([1, 4, 6], [1, 2, 3], [1, 4, 6])
    state_indices([1, 1, 1, 2, 2, 2, 1, 1],[3, 3, 4, 4, 5, 5, 3, 3])
    # Output: ([1, 3, 4, 5, 7], [1, 2, 3, 4, 1], [1, 3, 4, 5])
    # Last state is repeated, but only the first occurrence is kept in `ind_unique`
"""
function state_indices(vectors::AbstractArray...)
    # Ensure all vectors have the same length
    n = length(vectors[1])
    for v in vectors
        length(v) == n || throw(ArgumentError("All vectors must have the same length"))
    end

    # Stack vectors into an n × k matrix
    mat = hcat(vectors...)

    # Find index for each row
    rows = unique(eachrow(mat))
    dict = Dict{typeof(rows[1]), Int}()
    for (i, r) in enumerate(rows)
        dict[r] = i
    end
    index = [dict[r] for r in eachrow(mat)]

    # Find indices of changes and state value and indice of unique state changes for g-function
    change = diff(index) .!= 0
    # ind = findall(change)       # Indices of state change
    ind = findall(change) .+ 1  # Indices of state change, value at first time step of new state
    state = index[[true; change]]
    ind_value = unique(i -> index[i], eachindex(index))
    
    return [1; ind], state, ind_value
end
function state_indices(matrix::AbstractMatrix)
    return state_indices(eachcol(matrix)...)
end