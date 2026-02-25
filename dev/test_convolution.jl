"""
Script that tests the convolution operation on a simple example.
Notes:
    - For a heating power profile with constant value (such as 1 W/m or 1000 W/m) across all time
        steps, each method of convolution yields the same thermal response. In this case, the
        fastest approach remains the fft-based convolution, no matter the vector size. The matrix
        approach is the slowest across all vector sizes due to its O(n²) complexity.
    - For a heating power profile with varying values across time steps, the convolution results...
"""

using BenchmarkTools
using CairoMakie

includet("../src/ground_models/infinite_line_source.jl")
includet("../src/temporal_superposition.jl")
includet("../src/utils.jl")

# Define varying step heat flux vectors with even and uneven time steps
t_var = [0, 10, 15, 40, 70, 85] * 3600          # Time steps in seconds (uneven)
# t_var = [10:10:60] * 3600
q_var = [50.0, 30.0, 40.0, 20.0, 60.0, 30.0]    # Heat flux in W/m (uneven)

t_const = 1:60:t_var[end]                       # Interpolate for constant time steps
dt = diff([t_var; t_var[end]] .÷ 60)                     # Time step sizes
q_const = vcat([fill(q_var[i], dt[i]) for i in eachindex(q_var)]...)

# Define response functions using the infinite line source model
g_var = ils(t_var, 3.0, 2e6, 0.076)
g_const = ils(t_const, 3.0, 2e6, 0.076)

# Perform convolution using the three methods
T_var_time = convolution_time(q_var, g_var)
T_var_matrix = convolution_matrix(q_var, g_var)
T_var_fft = convolution(q_var, g_var)

T_const_time = convolution_time(q_const, g_const)
T_const_matrix = convolution_matrix(q_const, g_const)
T_const_fft = convolution(q_const, g_const)

# Plot results
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Thermal load (W/m)")
scatter!(ax, t_var, q_var, markersize=10)
lines!(ax, t_const, q_const)
ax = Axis(fig[1, 2], xlabel="Time (s)", ylabel="g-function (°Cm/W)")
scatter!(ax, t_var, g_var, markersize=10)
lines!(ax, t_const, g_const)
ax = Axis(fig[2, 1:2], xlabel="Time (s)", ylabel="Temperature (°C)")
lines!(ax, t_var, T_var_time, color=:blue, linewidth=3, linestyle=:solid, label="var_t")
lines!(ax, t_const, T_const_time, color=:red, linewidth=3, linestyle=:solid, label="const_t")
lines!(ax, t_var, T_var_matrix, color=:orange, linewidth=3, linestyle=:dash, label="const_m")
lines!(ax, t_const, T_const_matrix, color=:purple, linewidth=3, linestyle=:dash, label="var_m")
lines!(ax, t_var, T_var_fft, linewidth=5, color=:green, linestyle=:dot, label="var_f")
lines!(ax, t_const, T_const_fft, color=:pink, linewidth=5, linestyle=:dot, label="const_f")
axislegend(ax, position=:rb, orientation=:horizontal, nbanks=2)
display(fig)
