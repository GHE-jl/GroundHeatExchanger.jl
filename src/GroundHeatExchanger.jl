module GroundHeatExchanger

using DSP
using PCHIPInterpolation
using LinearAlgebra

# Temporal superposition
include("temporal_superposition.jl")

# Temperature simulations
include("temperature_simulation.jl")

# Utilities (GHE default parameters, heat load profile)
include("utils.jl")

# Temporal superposition exports
export convolution, convolutionf, impulse_func,
    convolution_ns, impulse_func_ns, state_vector, state_indices

# Temperature simulation exports
export fluid_temperature, outlet_temperature, inlet_temperature

# Utility exports
export GHE, heat_load_profile

end # module GroundHeatExchanger
