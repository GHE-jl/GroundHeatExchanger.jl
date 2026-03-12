module GroundHeatExchanger

using Revise

# Infinite line source of Ingersol (1948)
includet("ground_models/infinite_line_source.jl")
# Infinite cylindrical source of Ingersol (1959)
includet("ground_models/infinite_cylindrical_source.jl")
# Finite line source of Claesson and Javed (2011)
includet("ground_models/finite_line_source.jl")
# Moving finite line source of Pasquier et Lamarche (2022)
includet("ground_models/moving_infinite_line_source.jl")
# Moving finite line source of Guo et al. (2021)
includet("ground_models/moving_finite_line_source.jl")
# Standing column well model of Jacques et al. (2025)
includet("ground_models/beta_infinite_line_source.jl")

# Temporal superposition
includet("temporal_superposition.jl")

# Spatial superpositions techniques
includet("spatial_superposition.jl")

# Include thermal resistance evaluation
includet("borehole_thermal_resistance/resistance_fluid.jl")
includet("borehole_thermal_resistance/resistance_pipe.jl")
includet("borehole_thermal_resistance/resistance_borehole.jl")

# Include other varied functions
includet("utils.jl")

# Ground model export for closed-loop GHEs
export ils, ics, fls, mils, mfls

# Ground model export for standing column wells (probably to be moved to a separate module)
export βils_outlet, βils, Rb_SCW, convergence_flow

# Temporal superposition export
export convolution, convolutionf, impulse_func, 
    convolution_ns, impulse_func_ns, state_vector, state_indices

# Spatial superposition export
export bloc_matrix, successive_flux

# Thermal resistance export
export Reynold, Prandtl, Nusselt, friction_factor_Colebrook_White, resistance_fluid,
    resistance_pipe,
    resistance_borehole_multipole, resistance_total_internal_multipole,
    resistance_borehole_effective

# Utils export
export pchip_interpolation, set_nodes, borefield_xy, borefield_radius,
    water_ρ, water_cp, water_k, water_μ,
    head_loss_Darcy_Weisbach,
    GHE

# Temperature simulations
export g_model, T_f

"""
    g_model(t, ks, Cs, rb, H, D, xy)

Function that allows to generate a transfer function efficiently depending on the model required,
the number of time steps and arrangement of the borefield. Interpolations are performed if there are
more than 150 time steps, and more than 50 radius.
# Arguments
    - t: Time vector (nt x 1) [s]
    - ks: Soil thermal conductivity (1x1) [W/mK]
    - Cs: Soil volumetric specific heat (1x1) [J/m³K]
    - rb: Borehole radius (1x1) [m]
    - H: Borehole depth (1x1) [m]
    - D: Borehole burried depth (1x1) [m]
    - xy: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
# Output
    - g: The transfer function of the GHE (nt x 1) [°Cm/W]
"""
function g_model(t::Union{Real, AbstractVector{<:Real}}, ks::Real, Cs::Real, rb::Real, H::Real,
    D::Real, xy::AbstractArray{<:Real}; model::String="fls")
    # Set nodes if time vector is too long
    nt = length(t)
    if nt >= 150
        s = set_nodes(nt, 150)
    else
        s = range(1, nt)
    end

    # Select if spatial superposition is required or not
    if size(xy, 1) > 1
        gₛ = bloc_matrix(t[s], ks, Cs, rb, H, D, xy)
    else
        if model == "fls" || "FLS"
            gₛ = fls(t[s], ks, Cs, rb, H, D)
        end
    end

    # Interpolate to the right length
    if nt > 150
        return pchip_interpolation(t[s], gₛ, t)
    else
        return gₛ
    end
end

"""
    T_f(g, Q, Rbₑ, T₀)

Function that computes the average temperature between the inlet and outlet pipes using a borehole 
wall g-functions.
# Arguments
    - g: g-function for the borefield (nt x 1) [°Cm/W]
    - q: Ground loads (nt x 1) [W/m]
    - Rbₑ: Effective borehole thermal resistance (1 x 1) [mK/W]
    - T₀: Undisturbed ground temperature (1 x 1) [°C]
# Output
    - Tf: Average temperatut between the inlet and outlet pipes (nt x 1) [°C]
"""
function T_f(g::Union{Real, AbstractVector{<:Real}}, q::Union{Real, AbstractVector{<:Real}},
    Rbₑ::Real, T₀::Real)
    # Mean fluid temperature
    @. Tf = T₀ + (q * Rbₑ) + convolution(diff([0; q]), g)    
    return Tf
end
end