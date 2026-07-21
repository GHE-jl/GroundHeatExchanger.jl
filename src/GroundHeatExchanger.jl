module GroundHeatExchanger

using DSP
using PCHIPInterpolation
using BoreholeResistance
# ground_response is intentionally excluded from selective import — GHE defines its own wrapper
# that adds optional PCHIP compression. GroundResponse.ground_response is called directly inside.
import GroundResponse
using GroundResponse: AbstractGroundModel, ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel,
    ils, ics, fls, mils, mfls, successive_flux, bloc_matrix,
    borefield_radius, borefield, borefield_rectangle, borefield_line, borefield_circle,
    borefield_L, borefield_U, borefield_open_rectangle

# Temporal superposition (FFT-based convolution)
include("temporal_superposition.jl")

# Utilities — must precede temperature_simulation.jl (pchip_interpolation, set_nodes)
include("utils.jl")

# Temperature simulation and g-function wrappers
include("temperature_simulation.jl")

# Temporal superposition
export convolution, convolutionf, impulse_func, convolution_ns!,
    convolution_ns, impulse_func_ns, state_vector, state_indices

# Temperature simulation
export fluid_temperature, outlet_temperature, inlet_temperature
export ground_response  # GHE wrapper: GroundResponse.ground_response + optional PCHIP compression

# Interpolation utilities
export pchip_interpolation, set_nodes

# Package utilities
export GHE, heat_load_profile, head_loss_Darcy_Weisbach

# Re-exports from BoreholeResistance
export water_k, water_cp, water_ρ, water_μ
export Reynolds, Prandtl, Nusselt, Nusselt_annulus
export friction_factor_Colebrook_White, friction_factor_Tkachenko_Mileikovskyi
export resistance_fluid, resistance_pipe
export resistance_ULoop_borehole, resistance_ULoop_total_internal,
    resistance_ULoop_effective

# Re-exports from GroundResponse
export AbstractGroundModel
export ILSModel, ICSModel, FLSModel, MILSModel, MFLSModel
export ils, ics, fls, mils, mfls
export successive_flux, bloc_matrix
export borefield_radius
export borefield, borefield_rectangle, borefield_line, borefield_circle,
    borefield_L, borefield_U, borefield_open_rectangle

end
