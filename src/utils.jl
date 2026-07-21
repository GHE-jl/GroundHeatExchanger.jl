"""
    pchip_interpolation(tᵢ, vᵢ, t)

Function that performs the complete interpolation of a vector using the PCHIP interpolation method.
# Arguments
    - `tᵢ`: The id on which to interpolate
    - `vᵢ`: The vector to interpolate, sampled at `tᵢ`
    - `t`: The new interpolated vector sample
# Output
    - `v`: The new interpolated vector
"""
function pchip_interpolation(tᵢ::AbstractVector{<:Real}, vᵢ::AbstractVector{<:Real},
    t::AbstractVector{<:Real})
    interp = Interpolator(tᵢ, vᵢ)
    v = interp.(t)
    return v
end

"""
    GHE()

Function that defines a set of parameters to be used in the testing scripts for the
GroundHeatExchanger.jl package. The parameters are defined with typical values for ground heat
exchanger applications.
# Output
    - `t`: Time vector [s]
    - `H`: Borehole depth [m]
    - `D`: Borehole buried depth [m]
    - `s`: Shank spacing (s/2 is the half-shank spacing) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inlet radius [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK]
    - `Cs`: Ground volumetric specific heat [J/m³K]
    - `Cg`: Grout volumetric specific heat [J/m³K]
    - `Cp`: Pipe volumetric specific heat [J/m³K]
    - `Cf`: Fluid volumetric specific heat [J/m³K]
    - `ρs`: Ground density [kg/m³]
    - `ρg`: Grout density [kg/m³]
    - `ρp`: Pipe density [kg/m³]
    - `ρf`: Fluid density [kg/m³]
    - `μf`: Fluid dynamic viscosity [Pa·s or kg/(m·s)]
    - `vD`: Darcy velocity in the ground (groundwater flow) in m/s
    - `V`: Fluid flow rate [m³/s]
"""
function GHE()
    # Fluid properties (water_k/cp/ρ/μ) come from BoreholeResistance, re-exported by GHE.
    t = exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length = 500)) # Time (log) [s]
    H = 150.0                       # Borehole depth [m]
    D = 2.0                         # Borehole buried depth [m]
    s = 0.05                        # Shank spacing (s/2 is the half-shank spacing) [m]
    rb = 0.08                       # Borehole radius [m]
    ro = 0.022                      # Pipe outlet radius [m]
    ri = 0.017                      # Pipe inlet radius [m]
    T0 = 10.0                       # Undisturbed ground temperature [°C]
    ks = 3.0                        # Ground thermal conductivity [W/mK]
    kg = 1.6                        # Grout thermal conductivity [W/mK]
    kp = 0.4                        # Pipe thermal conductivity [W/mK]
    kf = water_k(T0)                # Fluid thermal conductivity [W/mK]
    Cs = 2.11e6                     # Ground volumetric specific heat [J/m³K]
    Cg = 2.25e6                     # Grout volumetric specific heat [J/m³K]
    Cp = 1.9e6                      # Pipe volumetric specific heat [J/m³K]
    Cf = water_cp(T0) * water_ρ(T0) # Fluid volumetric specific heat [J/m³K]
    ρs = 1000.0                     # Ground density [kg/m³]
    ρg = 1000.0                     # Grout density [kg/m³]
    ρp = 1000.0                     # Pipe density [kg/m³]
    ρf = water_ρ(T0)                # Fluid density [kg/m³]
    μf = water_μ(T0)                # Fluid dynamic viscosity [Pa.s or kg/m/s]
    vD = 1e-7                       # Darcy velocity in [m/s]
    V = 30/6e4                      # Fluid flow rate in m³/s (30 L/s)
    return t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, vD, V
end

"""
    ground_load_profile(t, A, B, C, D, E, F, G)

Function that defines a ground heat load profile for testing the convolution of the g-function with
a heat load profile. See Eqs. 7 and 8 of Bernier 2004.
# Arguments
    - `t`: Time vector or scalar *in hours* [h]
    - `A`: Amplitude parameter [W] (default: 2000)
    - `B`: Phase shift parameter [hours] (default: 2190)
    - `C`: Harmonic coefficient [hours] (default: 80)
    - `D`: Period parameter [-] (default: 2)
    - `E`: Modulation amplitude [-] (default: 0.01)
    - `F`: Phase offset [hours] (default: 0)
    - `G`: Phase parameter [-] (default: 0.95)
# Output
    - `Q`: Ground heat load profile [W]
# Reference
    - Bernier, M. A., Pinel, P., Labib, R., & Paillot, R. (2004). A Multiple Load Aggregation
        Algorithm for Annual Hourly Simulations of GCHP Systems. HVAC&R Research, 10(4), 471–487.
        https://doi.org/10.1080/10789669.2004.10391115
"""
function ground_load_profile(t, A=2000, B=2190, C=80, D=2, E=0.01, F=0, G=0.95)
    function f(t)
        harmonics = zeros(size(t))
        for n in 1:3
            coeff = 1 / (π * n) * (cos(C * π * n / 84) - 1)
            harmonics .+= coeff .* sin.(π * n .* (t .- B) / 84)
        end
        return A .* sin.(π/12 .* (t .- B)) .* sin.(π / 4380 .* (t .- B)) .*
               ((168 - C) / 168 .+ harmonics)
    end
    Q = f(t) .+ (-1) .^ floor.(D / 8760 .* (t .- B)) .* abs.(f(t)) .+
        E .* (-1) .^ floor.(D / 8760 .* (t .- B)) .* sign.(cos.(D .* π / 4380 .* (t .- F)) .+ G)
    return Q
end