"""
    AbstractANN

Model tag selecting which trained short-term ANN [`short_term_nodes`](@ref),
[`short_term_response`](@ref) and [`outlet_transfer_function`](@ref) evaluate.
"""
abstract type AbstractANN end

"""
    PublishedANN()

Short-term ANN of Pasquier, Zarrella & Labib (2018): 8 physical inputs (pipe properties fixed),
85 time nodes, 7-day validity horizon. See `short_term_response(::PublishedANN, ...)` for the
validity ranges. Implemented in `ann_published.jl`.
"""
struct PublishedANN <: AbstractANN end

"""
    DeepANN()

Short-term ANN of Pasquier & Marcotte (2020): pipe properties (`kp`, `Cf`, `ri`, `ro`) are free
parameters with their own (wider) validity ranges, 95 time nodes, 21-day validity horizon. This is
the **default** model tag for [`outlet_transfer_function`](@ref). See
`short_term_nodes(::DeepANN, ...)` for the validity ranges. Only the borehole-outlet response is
exposed (the underlying network can in principle report the temperature at any probe location
along the U-tube). Implemented in `ann_deep.jl`.
"""
struct DeepANN <: AbstractANN end

# Min–max scaler holding the offset, gain and target minimum used to normalise a network's input
# and de-normalise its output. Shared by every `AbstractANN` model.
struct _MinMaxScaler
    xoffset::Vector{Float64}
    gain::Vector{Float64}
    ymin::Int8
end

"""
    _mapminmax_apply(x, scaler)

Forward min–max scaling (MATLAB `mapminmax_apply`): `y = (x - xoffset) * gain + ymin`. First
step of a network's simulation.
"""
function _mapminmax_apply(x, scaler)
    y = x .- scaler.xoffset
    y .*= scaler.gain
    y .+= scaler.ymin
    return y
end

"""
    _mapminmax_reverse(y, scaler)

Inverse min–max scaling (MATLAB `mapminmax_reverse`): `x = (y - ymin) / gain + xoffset`. Last
step of a network's simulation.
"""
function _mapminmax_reverse(y, scaler)
    x = y .- scaler.ymin
    x ./= scaler.gain
    x .+= scaler.xoffset
    return x
end
