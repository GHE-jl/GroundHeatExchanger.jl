# Demonstration of FFT-based stationary temporal superposition.
# Ground response is computed with the FLS model from GroundResponse.jl.

using CairoMakie
using GroundHeatExchanger

# Parameters
_, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
t = 60.0:60:3600*24*6   # 6 days, 1-minute steps [s]
th = t ./ 3600            # hours (for plotting)

# g-function (FLS)
g = fls(t, rb, H, D, ks, Cs)

# Heat load with random noise
q = [12000.0 * ones(60*24*2); 8000.0 * ones(60*24*2); 10000.0 * ones(60*24*2)] / H
q_rnd = q .+ rand(length(t))

# Temporal superposition
ΔT = convolution(q_rnd, g)
ΔT2 = convolutionf(impulse_func(q_rnd), g)
@assert ΔT ≈ ΔT2  "FFT-based convolutions match"

# Plot
fig = Figure(size=(900, 400))
ax1 = Axis(fig[1, 1], xlabel="Time (h)", ylabel="q (W/m)", title="Heat load with noise")
lines!(ax1, th, q_rnd)
ax2 = Axis(fig[1, 2], xlabel="Time (h)", ylabel="ΔT (°C)", title="Temperature response")
lines!(ax2, th, ΔT)
display(fig)
