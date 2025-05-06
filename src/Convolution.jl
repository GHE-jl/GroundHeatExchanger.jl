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
        - f: incremental load functions (impulse function)
        - g: transfer function (or g-function)
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

    # Inverse FFT to get convolution resuls
    y = irfft(F_fg, optimal_length)[1:n]
    return y
end

