using FFTW

function convolution(f::AbstractVector{T}, g::AbstractVector{T}) where {T<:Real}
    """
        convolution(f, g)
    
    Function that solves a convolution product in the spectral domain. Padding is used to
    avoid circular convolution.
    Inputs:
        - f: incremental load functions (impulse function)
        - g: transfer function (or g-function)
    Output:
        - Convolved signal
    """
    n = length(f)
    total_length = 2 * n - 1

    # Preallocate arrays
    f_pad = zeros(T, total_length)
    g_pad = zeros(T, total_length)
    y = zeros(n)

    # Fill preallocated arrays
    copyto!(f_pad, f)
    copyto!(g_pad, g)

    # Precompute FFT plans
    fft_plan = plan_rfft(f_pad)

    # Compute FFTs
    F_f = fft_plan * f_pad
    F_g = fft_plan * g_pad

    # Element-wise multiplication (in-place)
    F_f .*= F_g

    # Inverse FFT to get convolution resuls
    y = irfft(F_f, total_length)
    return y[1:n]
end