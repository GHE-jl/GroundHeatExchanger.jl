using PCHIPInterpolation


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

function pchip_interpolation(tᵢ::Vector{T}, gᵢ::Vector{T}, t::Vector{T}) where T<:Real
    """
        pchip_interpolation(tᵢ, gᵢ, t)
    """
    interp = Interpolator(tᵢ, gᵢ)
    g = interp.(t)
    return g
end