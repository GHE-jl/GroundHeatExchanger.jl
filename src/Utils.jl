using PCHIPInterpolation
using LinearAlgebra

"""
    pchip_interpolation(tᵢ, gᵢ, t)

Function that performs the complete interpolation of a vector using the PCHIP interpolation method.
# Arguments
    - tᵢ: The id on which to interpolate
    - vᵢ: The vector to interpolate, sampled at tᵢ
    - t: The new interpolated vector sample
# Output
    - v: The new interpolated vector
"""
function pchip_interpolation(tᵢ::AbstractVector{<:Real}, vᵢ::AbstractVector{<:Real},
    t::AbstractVector{<:Real})
    interp = Interpolator(tᵢ, vᵢ)
    v = interp.(t)
    return v
end

"""
    set_nodes(nt, n₀)

Function that sets a logarithmic progression of node positions on a transfer function.
# Arguments
    - nt: Total number of data in the input vectors [-]
    - n₀: User defined number of nodes on the transfer function [-]
# Output
    - id: A vector of length "n₀" of node positions on the transfer function [-]
"""
function set_nodes(nt::Real, n₀::Integer)
    # Basic inputs
    n_tmp = n₀ - 1
    id = Vector{Integer}(undef, n_tmp)
    # Fill the vector with node positions
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

"""
    borefield_xy(nx, ny, B)
    borefield_xy(nx, ny, Bx, By)

Function that generates the coordinates of a rectangular borefield given the number of boreholes in
the `x` and `y` directions and the spacing between them. For a square borefield, `Bx` and `By` are
equal.
# Arguments
    - `nx`: Number of boreholes in the x direction [-]
    - `ny`: Number of boreholes in the y direction [-]
    - `B`: Spacing between boreholes in both directions [m]
    - `Bx`: Spacing between boreholes in the x direction [m]
    - `By`: Spacing between boreholes in the y direction [m]
# Output
    - `xy`: Matrix of borehole coordinates (nb x 2) [m]
"""
function borefield_xy(nx::Integer, ny::Integer, B::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* B .- B
    return xy
end
function borefield_xy(nx::Integer, ny::Integer, Bx::Real, By::Real)
    xy = hcat([[i, j] for i in 1:nx for j in 1:ny]...)' .* [Bx By] .- [Bx By]
    return xy
end

"""
    borefield_radius(xy, rb)

Function that computes a radius matrix, vector, unique values and indices of a borefield given
the coordinates of each borehole. This helps to compute the g-function of a borefield using spatial
superposition, as it allows to compute the ground thermal response at each unique radius of the
borefield.
# Arguments
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nb x 2) [m]
        - E.g.: [0 0] (to have a matrix input).
    - `rb`: Borehole radius [m]
# Outputs
    - `r`: Matrix of the borefield radius (nb x nb) [m]
    - `rᵥ`: Vector of the borefield radius ((nb x nb) x 1) [m]
    - `rᵤ`: Unique values of the borefield radius (nbᵤ x 1) [m]
    - `rᵢ`: Indices of the unique radius values (nb*nb x 1) [m]
    - `θ`: Angle of the boreholes in the borefield from the origin (0,0) (nb x nb) [°]
    - `nb`: Number of boreholes [-]
"""
function borefield_radius(xy::AbstractArray{<:Real}, rb::Real)
    # Analyse the radius of the borefield
    r = sqrt.(sum(abs2, xy, dims=2) .+ sum(abs2, xy, dims=2)' .- 2 * (xy * xy'))

    # Outputs
    nb = size(xy, 1)                # Number of boreholes
    r = r + Diagonal(rb * ones(nb)) # Matrix of the borefield radius (nb x nb) [m]
    rᵥ = reshape(r, nb * nb)        # Vector of the borefield radius (nb x 1) [m]
    rᵤ = unique(rᵥ)                 # Unique values of the borefield radius (nbᵤ x 1) [m]
    rᵢ = indexin(rᵥ, rᵤ)            # Indices of the unique radius values (nb*nb x 1) [m]
    θ = atan.(xy[:, 2], xy[:, 1])   # Angle of the boreholes from the origin (0,0) (nb x nb) [°]
    return r, rᵥ, rᵤ, rᵢ, θ, nb
end

"""
    water_k(T::Real)

Water thermal conductivity k(T), 0 ≤ T ≤ 99.6°C at 1 bar. The polynomial equation is a fit to the 
data from the Engineering Toolbox. The values can be validates with a temperature vector 
`T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - k: Thermal conductivity [W/mK]
# Reference
    - The Engineering ToolBox (2018). Thermal Conductivity of Water: Temperature and Pressure Data. 
        [online] Available at: https://www.engineeringtoolbox.com/water-liquid-gas-thermal-
        conductivity-temperature-pressure-d_2012.html [Accessed 2026-01-14].
"""
function water_k(T::Real)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_k` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 0.5557250521318174 + 0.002490814640452007*T - 2.1170044416971473e-5*T^2 + 
        1.285515973680875e-7*T^3 - 4.546428806458628e-10*T^4 + 9.750314739837196e-14*T^5
end

"""
    water_cp(T::Real)

Water specific heat capacity cp(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the isobaric 
specific heat capacity data from the Engineering Toolbox. The values can be validates with a 
temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - cp: Specific heat capacity [J/(kg·K)]
# Reference
    - The Engineering ToolBox (2004). Specific Heat Capacity of Water: Temperature-Dependent Data 
        and Calculator. [online] Available at: 
        https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html 
        [Accessed 2026-01-14].
"""
function water_cp(T::Real)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_cp` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 4219.849078078278 - 3.266686616602623*T + 0.09969277880041719*T^2 - 
        0.0014860911377001344*T^3 + 1.161963811563561e-5*T^4 - 3.5034316470844105e-8*T^5
end

"""
    water_ρ(T::Real)

Water density ρ(T), 0 ≤ T ≤ 100°C at 1 atm. The polynomial equation is a fit to the data from the 
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - ρ: Density [kg/m³]
# Reference
    - The Engineering ToolBox (2003). Water Density, Specific Weight and Thermal Expansion 
        Coefficients - Temperature and Pressure Dependence. [online] Available at: 
        https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html 
        [Accessed 2026-01-14].
"""
function water_ρ(T::Real)::Float64
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_ρ` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 999.8475436930158 + 0.06180756931966996*T - 0.008309049138917115*T^2 + 
        6.35713412865478e-5*T^3 - 3.8404497053894326e-7*T^4 + 1.0249871031879443e-9*T^5
end

"""
    water_μ(T::Real)

Dynamic viscosity μ(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the data from the 
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - μ: Dynamic viscosity [Pa·s or kg/(m·s)]
# Reference
    - The Engineering ToolBox (2004). Water - Dynamic and Kinematic Viscosity at Various 
        Temperatures and Pressures. [online] Available at: https://www.engineeringtoolbox.com/water-
        dynamic-kinematic-viscosity-d_596.html [Accessed 2026-01-14].
"""
function water_μ(T::Real)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_μ` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 0.001790966556989398 - 5.965082369793418e-5*T + 1.3185191782991122e-6*T^2 - 
        1.8236868027209892e-8*T^3 + 1.3644271817518522e-10*T^4 - 4.137645533574321e-13*T^5
end

"""
    head_loss_Darcy_Weisbach(L, D, V, f)

Function that calculates the head loss in a pipe using the Darcy-Weisbach equation.
# Arguments
    - `L`: Length of the pipe [m]
    - `r`: Pipe inside or annulus (r = rb - ro) radius [m]
    - `V`: Velocity of the fluid in the pipe [m/s]
    - `f`: Darcy friction factor [-] (from `friction_factor_Colebrook_White()`)
# Output
    - `h`: Head loss [m]
"""
function head_loss_Darcy_Weisbach(L::Real, r::Real, V::Real, f::Real)
    # Friction factor
    return f * (L / (2 * r)) * (V^2 / 2)
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
    - `ϵ`: Roughness of pipes [m]
    - `vD`: Darcy velocity in the ground (groundwater flow) in m/s
    - `V`: Fluid flow rate [m³/s]
# Example
    t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
"""
function GHE()
    # Define paremeters
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
    ϵ = 1e-5                        # roughness of pipes [m]
    vD = 1e-7                       # Darcy velocity in [m/s]
    V = 30/6e4                      # Fluid flow rate in m³/s (30 L/s)
    return t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V
end

"""
    heat_load_profile(t, A, B, C, D, E, F, G)

Function that defines a heat load profile for testing the convolution of the g-function with a heat
load profile. See Eqs. 7 and 8 of Bernier 2004.
# Arguments
    - `t`: Time vector or scalar [s]
    - `A`: Amplitude parameter [W] (default: 10246)
    - `B`: Phase shift parameter [hours] (default: 2409)
    - `C`: Harmonic coefficient [hours] (default: 64)
    - `D`: Period parameter [-] (default: 2)
    - `E`: Modulation amplitude [-] (default: 0.01)
    - `F`: Phase offset [hours] (default: 0)
    - `G`: Phase parameter [-] (default: 0.95)
# Output
    - `Q`: Heat load profile [W]
# Reference
    - Bernier, M. A., Pinel, P., Labib, R., & Paillot, R. (2004). A Multiple Load Aggregation
        Algorithm for Annual Hourly Simulations of GCHP Systems. HVAC&R Research, 10(4), 471–487.
        https://doi.org/10.1080/10789669.2004.10391115
"""
function heat_load_profile(t, A=2000, B=2190, C=80, D=2, E=0.01, F=0, G=0.95)
    # Define the base heat load function that handles both scalars and vectors
    function f(t)
        # Initialize harmonics sum
        harmonics = zeros(size(t))
        for n in 1:3
            coeff = 1 / (π * n) * (cos(C * π * n / 84) - 1)
            harmonics .+= coeff .* sin.(π * n .* (t .- B) / 84)
        end
        
        return A .* sin.(π/12 .* (t .- B)) .* sin.(π / 4380 .* (t .- B)) .* 
               ((168 - C) / 168 .+ harmonics)
    end
    
    # Compute the heat load profile with modulation
    Q = f(t) .+ (-1) .^ floor.(D / 8760 .* (t .- B)) .* abs.(f(t)) .+ 
        E .* (-1) .^ floor.(D / 8760 .* (t .- B)) .* sign.(cos.(D .* π / 4380 .* (t .- F)) .+ G)
    return Q
end