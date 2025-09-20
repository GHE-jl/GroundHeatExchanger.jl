"""
Module that is used to define various analytical models used to simulate the temperature of 
ground heat exchanger (GHE). The approach used is to compute a g-function (temperature at 
the borehole wall) using the analytical model, and then using a convolution approach to 
simulate the borehole outlet temperature.

Features that are used include spatial superposition, computing the borehole thermal 
resistance Rb through the first-order multipole method, ...

Author: Gabriel Dion (dion.gabriel100@gmail.com)
Date: 2025-03
Julia version: 1.11.3
"""

module GHEModels

# Include files to make the package
include("GroundModels.jl")
include("MFLS.jl")
include("SpatialSuperposition.jl")
include("Convolutions.jl")
include("ThermalResistances.jl")
include("Utilities.jl")

# Ground model export
export ils, ics, fls, ilsβ
export mfls_single_borehole, mfls_borefield_I

# Spatial superposition export
export borefield_radius, g_matrix, bloc_matrix, successive_flux

# Thermal resistance export
export R_f,
    R_p,
    R_b_zeroth_order_multipole,
    R_b_first_order_multipole,
    R_b,
    R_a_first_order_multipole,
    R_bₑ

# Other export from included files
export convolution,
    set_nodes,
    pchip_interpolation

# Temperature simulations
export ghe_model,
    building_to_ground_loads,
    outlet_temperature

function ghe_model(t::AbstractVector{T}, ks::T, Cs::T, rb::T, H::T, D::T, xy::Array{T}; 
    model::String="fls") where T<:Real
    """
    Allow to run different model with various inputs. (Very poor definition right here)
    """

    # Set nodes if time vector is too long
    nt = length(t)
    if nt >= 150
        s = set_nodes(nt, 150)
    else
        s = range(1,nt)
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

function outlet_temperature(g::Vector{T}, Q::Vector{T}, V::Vector{T}, H::T, Rb::T, T₀::T
    ) where T<:Real
    """
    Function that computes the borehole inlet and outlet temperature using g-functions.
    Inputs:
        - g: g-function for the borefield [-]
        - Q: Ground loads [W]
        - V: Circulating flow rate [m^3/s]
        - H: Borehole depth [m]
        - Rb: Borehole thermal resistance (first-order multipole) [mK/W]
        - T₀: Undisturbed ground temperature [degC]
    """

    # Mean fluid temperature
    Tf = T₀ .+ Q*Rb/H .+ convolution(diff([0;Q/H]), g)

    # Temperature variation between inlet and outlet
    dT = Q/(V*999.7*4190)

    # Outlet and inlet temperature
    Tout = Tf .- dT_2
    Tin = Tf .+ dT/2
    
    return Tin, Tout
end
end