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
includet("Models.jl")
includet("SpatialSuperposition.jl")
includet("Convolution.jl")
includet("ThermalResistance.jl")
includet("Utilities.jl")

# Ground model export
export ils, fls, mfls, scwm

# Spatial superposition export
export GHE_param, gfunc_matrix, bloc_matrix, successive_flux

# Other export from included files
export convolution,
    Rb_first_order_multipole,
    Ra_first_order_multipole,
    set_nodes,
    pchip_interpolation

# Temperature simulations
export ghe_model,
    building_to_ground_loads,
    outlet_temperature

function ghe_model(t::Vector{T}, kg::T, Cg::T, rb::T, H::T, D::T, xy::Array{T}; 
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
    if size(xy,1) > 1
        gₛ = spatial_sup(t[s], kg, Cg, rb, H, D, xy, model)
    else
        if model == "fls" || "FLS"
            gₛ = fls(t[s], kg, Cg, rb, H, D)
        end
    end

    # Interpolate to the right length
    if nt > 150
        return pchip_interpolation(t[s], gₛ, t)
    else
        return gₛ
    end
end

function building_to_ground_loads(Qb::Vector{T1}, COPh::T2, COPc::T2,
    pc_peakh=1, pc_peakc=1) where {T1<:Real, T2<:Real}
    """
    Converts building loads to ground loads, as a function of COP for both heating (COPh) and 
    cooling (COPc). Option inputs can be used to specify the percentage (pc) of the peak loads (for
    both heating and cooling) that has to be covered by the geothermal system. The default is 100%
    coverage.
    """

    # Cut the loads to the percentage of peak coverage.
    Qb[Qb .< pc_peakh*minimum(Qb)] .= pc_peakh*minimum(Qb)
    Qb[Qb .> pc_peakc*maximum(Qb)] .= pc_peakc*maximum(Qb)

    # Convert building loads (Qb) to ground loads (Qg)
    Qgh = Qb[Qb .<= 0] * (1 - (1/COPh))
    Qgc = Qb[Qb .> 0] * (1 + (1/COPc))
    Qg = Qgh + Qgc
    return Qg
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