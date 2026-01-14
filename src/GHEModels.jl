"""
Module that is used to define various analytical models used to simulate the temperature of 
ground heat exchanger (GHE). The approach used is to compute a g-function (temperature at 
the borehole wall) using the analytical model, and then using a convolution approach to 
simulate the borehole outlet temperature.

Features that are used include spatial superposition, computing the borehole thermal 
resistance Rb through the first-order multipole method, ...

Author: Gabriel Dion (dion.gabriel100@gmail.com)
Date: 2026-01
Julia version: 1.11.3
"""

module GHEModels

using Revise

# Include files to make the package
includet("GroundModels.jl")
includet("SpatialSuperpositions.jl")
includet("ThermalResistances.jl")
includet("Convolutions.jl")
includet("Utilities.jl")

# Ground model export for closed-loop GHEs
export ils, ics, fls, mfls

# Ground model export for standing column wells
export βils

# Spatial superposition export
export borefield_radius, g_matrix, bloc_matrix, successive_flux

# Thermal resistance export
export R_f, R_p, R_b_zeroth_order_multipole, R_b_first_order_multipole, R_b,
    R_a_first_order_multipole, R_bₑ

# Convolution export
export convolution, step_signal, state_transitions, convolution_ns

# Other export
export set_nodes, pchip_interpolation

# Temperature simulations
export g_model,
    T_f

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