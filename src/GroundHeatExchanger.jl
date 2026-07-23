module GroundHeatExchanger

using DSP
using PCHIPInterpolation
using BoreholeResistance
using GroundResponse: AbstractGroundModel, ILSModel, ICSModel, FLSModel, MILSModel,
    MFLSModel, ground_response, ils, ics, fls, mils, mfls,
    successive_flux, bloc_matrix, borefield_geometry, borefield, borefield_rectangle,
    borefield_line, borefield_circle, borefield_L, borefield_U, borefield_open_rectangle

# Utilities (pchip_interpolation, GHE, ground_load_profile)
include("utils.jl")

# Artificial neural network (abstraction, then one file per concrete model)
include("ann_model.jl")
include("ann_published.jl")
include("ann_deep.jl")

# Temporal superposition (FFT-based convolution)
include("temporal_superposition.jl")

# Temperature simulation (uses convolution; g-functions come from GroundResponse.ground_response)
include("temperature_simulation.jl")

# Temporal superposition
export convolution, convolutionf, impulse_func, convolution_ns!,
    convolution_ns, impulse_func_ns, state_vector, state_indices

# Temperature simulation
export fluid_temperature, outlet_temperature, inlet_temperature, outlet_transfer_function

# Short-term ANN model dispatch and models
export AbstractANN, PublishedANN, DeepANN

# Short-term ANN (Pasquier, Zarrella & Labib, 2018 — `PublishedANN`; Pasquier & Marcotte, 2020 —
# `DeepANN`, the default)
export short_term_response, short_term_nodes

# Interpolation utilities
export pchip_interpolation

# Package utilities
export GHE, ground_load_profile

# Re-exports from BoreholeResistance
export water_k, water_cp, water_ρ, water_μ
export Reynolds, Prandtl, Nusselt, Nusselt_annulus, convection_coefficient
export friction_factor_Colebrook_White, friction_factor_Tkachenko_Mileikovskyi
export resistance_fluid, resistance_pipe
export resistance_ULoop_borehole, resistance_ULoop_total_internal, resistance_ULoop_effective
export resistance_coaxial, resistance_coaxial_effective

# Re-exports from GroundResponse
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel
export ground_response
export ils, ics, fls, mils, mfls
export successive_flux, bloc_matrix
export borefield_geometry
export borefield, borefield_rectangle, borefield_line, borefield_circle,
    borefield_L, borefield_U, borefield_open_rectangle

end
