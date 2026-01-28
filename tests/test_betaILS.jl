"""
Script testing the enhanced \beta-ILS model from GHEModels.jl using a single borehole.
"""

using BenchmarkTools
using CairoMakie
includet("../src/GHEModels.jl")
using .GHEModels

# Define physical parameters for the SCW
H = 500.0                       # Borehole depth [m]
Hp = 495.0                      # Pipe length [m]
ri = 0.02434                    # Inner radius of the pipe [m]
ro = 0.03016                    # Outer radius of the pipe [m]
rb = 0.083                      # Borehole radius [m]
rm = 100.0                      # Influence radius of the borehole [m]

# Define ground properties
T0 = 8.0                        # Initial temperature of the ground [°C]
ks = 2.74                       # Thermal conductivity of the ground [W/m·K]
kp = 0.42                       # Thermal conductivity of the pipe [W/m·K]
Cs = 2e6                        # Volumetric heat capacity of the ground [J/m³·K]
Cf = 4.2e6                      # Volumetric heat capacity of the fluid [J/m³·K] Should not be used.
Ku = 1e-7                       # Horizontal hydraulic conductivity of the ground [m/s]
Kf = 1e-5                       # Fracture hydraulic conductivity of the ground [m/s]

K = [495.0 Ku; 5.0 Kf]          # Hydraulic cond. with layer thicknesses and K values [m, m/s]
k = [500.0 ks Cs]               # Thermal properties with layer thicknesses [m, W/m·K, J/m³·K]

# Define operating conditions for the SCW
t = 60.0:60.0:3600.0 * 24 * 30  # Time vector (1 month, 1-minute steps, in seconds)
n = length(t)
Qp = ones(n) * 1000.0           # Constant heat emitted by the pump [W]
Qg = ones(n) * 30000.0 + Qp     # Constant heat injected from the ground [W]
Qg[1:60 * 12] = zeros(60 * 12)  # No injection for the first 12 hours

# Pumping flow rate [m³/s]
Vp = ones(n) * 175 / 60000
Vp[60 * 24 * 15:end] .= 125 / 60000

# Bleed flow rate [m³/s]
Vb = ones(n) * 17.5 / 60000
Vb[1:60 * 24 * 10] .= 1e-9      # Avoid value of zero for the β-ILS model
Vb[60 * 24 * 10 + 1:60 * 24 * 20] .= 17.5 / 60000 / 2

# Other parameters used for calculation
h_recircu = 0.22                # Observed drawdown due to pumping [m]
r_influence = 100               # Influence radius of pumping [m]

# Pre-processing of the input data
ind, state, ind_unique = state_transitions(Vp, Vb)
V = Vp[ind_unique]              # Unique pumping flow rates [m³/s]
B = Vb[ind_unique]              # Unique bleed flow rates [m³/s]
ns = length(ind)                # Number of states

# Evaluate Peclet number for each state
# At = B / 2 / π ./ H             # Fluid planar rate in the pipe [m^2/s]
# α = ks ./ (ones(ns) * Cs)       # Diffusivity (m^2/s)
# Pe = At / α                     # Peclet number (-)

# Evaluate the thermal resistance
Rb, Rv, f = zeros(length(ind)), zeros(length(ind)), zeros(length(ind)) # Preallocation
for i in 1:length(ind)
    Rb[i], Rv[i], f[i] = Rb_SCW(V[i], kp, rb, ro, ri, H, Hp, 8.0)
end

KK = effective_K(K, H)

### Testing ###
g1 = _βils(t, k, rb, H, V[1], B[1] / V[1], f[1], K = K)

# Evaluate the transfer functions for each state
g = Matrix{Float64}(undef, n, ns)
for i in 1:length(ind)
    g[:, i] = βils_outlet(t, k, kp, rb, ro, ri, H, Hp, V[i], B[i] ./ V[i], T0; K = K)
end

# Plot the transfer functions
f = Figure(; size = (17 * 96 / 2.54, 12 * 96 / 2.54))
ax = Axis(f[1, 1], xlabel = L"$t$ (s)", ylabel = "Transfer function (-)", xscale = log10)
for i in 1:ns
lines!(ax, t, g[:, i], color = :blue, linewidth = 1.5,
    label = "State $i: V=$(round(V[i]*6e4)) L/min, B=$(round(B[i]*6e4)) L/min")
end
axislegend(ax; position = :rt)
