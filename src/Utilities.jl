using PCHIPInterpolation, Colors


function set_nodes(nt::Int, n₀::Int)
    """
        set_nodes(nt, n₀)
    
    Function that sets a logarithmic progression of node positions on a transfer function.
    Inputs:
        - nt: Total number of data in the input vectors [-]
        - n₀: User defined number of nodes on the transfer function [-]
    Output:
        - id: A vector of length "n₀" of node positions on the transfer function [-]
    """

    n_tmp = n₀ - 1
    id = Vector{Integer}(undef, n_tmp)

    while length(id) < n₀
        empty!(id)
        for x in range(0, stop=log10(nt), length=n_tmp)
            push!(id, round(Int, exp10(x)))
        end
        unique!(id)
        n_tmp += 1
    end
    return id
end

function pchip_interpolation(tᵢ::AbstractVector{T}, gᵢ::Vector{T},
    t::AbstractVector{T}) where T<:Real
    """
        pchip_interpolation(tᵢ, gᵢ, t)
    """
    interp = Interpolator(tᵢ, gᵢ)
    g = interp.(t)
    return g
end

function fig_color()
    """
        fig_color()
    
    Outputs a palette of predefined colors for figures.
    """
    col = [
        RGB(30/255, 144/255, 205/255),      # Navy Blue
        RGB(255/255, 140/255, 0/255),       # Dark Orange
        RGB(50/255, 157/255, 13/255),       # Kinda Forest Green
        RGB(205/255, 175/255, 0/255),       # Gold yellow
        RGB(128/255, 0/255, 128/255),       # Purple
        RGB(205/255, 0/255, 0/255),         # Medium Red
        RGB(0/255, 0/255, 205/255),         # Medium Blue
        RGB(0/255, 205/255, 0/255),         # Medium green
        RGB(105/255, 105/255, 105/255),     # DimGrey1
        RGB(60/255, 60/255, 60/255)         # DimGrey2
    ]
    return col
end