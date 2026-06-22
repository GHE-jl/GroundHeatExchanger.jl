# Pipe head loss calculation for a single U-loop GHE using the Darcy-Weisbach equation.
# Demonstrates how to compute pressure drop across borehole pipes for pump sizing.

using GroundHeatExchanger

# Parameters from GHE()
t, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()

kf = water_k(T0)

V̇ = V / (π * ri^2)            # Mean fluid speed in pipe [m/s]
Re = Reynolds(V̇, ri, ρf, μf)
f  = friction_factor_Colebrook_White(Re, ri, ϵ)

# Borehole circuit: fluid travels down + up = 2H
L_circuit = 2 * H
Δh = head_loss_Darcy_Weisbach(L_circuit, ri, V̇, f)

println("=== Head loss — single U-loop borehole ===")
println("  Flow rate : $(round(V * 1e6 * 60, digits=1)) L/min")
println("  Re        = $(round(Int, Re))")
println("  f (CW)    = $(round(f, digits=5))")
println("  Δh        = $(round(Δh, digits=3)) m  (L = $(round(Int, L_circuit)) m)")

# Friction factors: Colebrook-White vs Tkachenko-Mileikovskyi must agree within 1%
f_TM = friction_factor_Tkachenko_Mileikovskyi(Re, ri, ϵ)
@assert abs(f - f_TM) / f < 0.01  "CW and TM friction factors must agree within 1%"

# Flow-rate sweep
println("=== Head loss sweep ===")
println(rpad("V [L/min]", 12), " ", rpad("Re", 8), " ", rpad("f", 10), " ", "Δh [m]")
for Q_Lmin in [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]
    V_i  = Q_Lmin / 1000 / 60
    V̇_i = V_i / (π * ri^2)
    Re_i = Reynolds(V̇_i, ri, ρf, μf)
    f_i  = friction_factor_Colebrook_White(Re_i, ri, ϵ)
    Δh_i = head_loss_Darcy_Weisbach(L_circuit, ri, V̇_i, f_i)
    println(rpad(string(round(Q_Lmin, digits=1)), 12), " ",
            rpad(string(round(Int, Re_i)), 8), " ",
            rpad(string(round(f_i, digits=5)), 10), " ",
            round(Δh_i, digits=3))
end

@assert Δh > 0  "Head loss must be positive"
println("\nAll assertions passed.")
