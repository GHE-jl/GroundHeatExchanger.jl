# Borehole-outlet (T_out / source-side entering water) transfer function: the short-term ANN
# response (PublishedANN or DeepANN) joined to a long-term ground model (FLS) through the
# effective borehole resistance Rb*.
# The script reproduces two figures of the 2018 paper with PublishedANN, then compares both
# short-term ANNs on a multi-borehole field:
#   Figure 1 → Fig. 2a: five short-term transfer functions (Table 1 cases) vs t/tₛ
#   Figure 2 → Fig. 3 : short-term joined to long-term for a single borehole (outlet + inlet)
#   Figure 3          : borefield outlet transfer function — PublishedANN vs DeepANN, same input
# Reference
#   Pasquier, P., Zarrella, A., & Labib, R. (2018). Application of artificial neural networks
#   to near-instant construction of short-term g-functions. Applied Thermal Engineering, 143,
#   910–921.
#   Pasquier, P., & Marcotte, D. (2020). Robust identification of volumetric heat capacity and
#   analysis of thermal response tests by Bayesian inference with correlated residuals. Applied
#   Energy, 261, 114394.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundHeatExchanger
using BoreholeResistance

# Pipe and fluid properties shared by every case (Table 1 footnote)
kp = 0.4; Cp = 1.9e6; Cf = 4.2e6; ri = 0.017; ro = 0.022

# Figure 1 — Fig. 2a: short-term transfer functions for the five Table 1 cases
#            ks    Cs     kg    Cg     rb     H      V̇ (L/min)  s
cases = [
    (0.75, 1.8e6, 0.75, 1.8e6, 0.06, 125.0, 21.0, 0.064),   # Case 1
    (0.75, 2.2e6, 2.75, 2.4e6, 0.07, 175.0, 29.0, 0.064),   # Case 2
    (3.75, 2.2e6, 2.75, 1.8e6, 0.07, 175.0, 29.0, 0.064),   # Case 3
    (3.75, 2.4e6, 0.75, 2.4e6, 0.06, 125.0, 29.0, 0.064),   # Case 4
    (3.75, 2.4e6, 2.75, 2.4e6, 0.07, 195.0, 21.0, 0.084),   # Case 5
]

dt = 15.0
tf = 7 * 24 * 3600.0

fig1 = Figure()
ax1 = Axis(fig1[1, 1];
    title  = "Fig. 2a — short-term transfer functions (Table 1)",
    xlabel = "t / tₛ  (-)", ylabel = "g  (-)", xscale = log10)
for (i, (ks, Cs, kg, Cg, rb, H, Vlpm, s)) in enumerate(cases)
    local V̇ = Vlpm / 1000 / 60
    local t, g = short_term_response(PublishedANN(), collect(dt:dt:tf), rb, ri, ro, H, s, V̇, ks, Cs,
        kg, Cg, kp, Cp, Cf)
    local ts_char = H^2 / (9 * ks / Cs)
    lines!(ax1, t ./ ts_char, g; label = "Case $i", linewidth = 1.3)
end
axislegend(ax1; position = :lt)
display(fig1)

# Figure 2 — Fig. 3: short- and long-term combination for a single borehole
# TRT case of Table 1
ks = 2.13
Cs = 2.0e6
kg = 1.65
Cg = 2.0e6
rb = 0.08
H  = 152.4
V̇  = 23.7 / 1000 / 60
s  = 0.058
D = 4.0

# Effective borehole resistance Rb* (the paper's "equivalent" resistance), first-order
# multipole, with water properties at a representative 10 °C.
Tmean = 10.0
kf = water_k(Tmean); cf = water_cp(Tmean); ρf = water_ρ(Tmean); μf = water_μ(Tmean)
Rbₑ = resistance_ULoop_effective(V̇, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
println("Effective borehole resistance Rb* = ", round(Rbₑ; sigdigits = 4), " m·K/W")

# Long-term ground model and a log-spaced time vector (1 min → 50 yr)
m = FLSModel(H, D, ks, Cs)
t = 10 .^ range(log10(60), log10(50 * 365 * 24 * 3600); length = 400)
xy = borefield(:rectangle, 4, 3, 5.0)

# Combined outlet transfer function, and its short-/long-term components for context.
# `interp = true` is the documented correctness requirement for the borefield temporal solver
# (successive_flux) on a non-uniform, log-spaced `t`.
nb = size(xy, 1)
g_out = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s, Rbₑ, m;
    xy = xy, interp = true, model = PublishedANN())
t_st, g_st = short_term_response(PublishedANN(), collect(dt:dt:tf), rb, ri, ro, H, s, V̇, ks, Cs,
    kg, Cg, kp, Cp, Cf)
# Raw long-term branch (no splice/shift), on the same per-borehole basis as the combined curve:
# the field ground response is normalised to 1 W/m of *total* field heat rate, so scale by nb.
g_lt = (nb .* ground_response(t, rb, xy, m; interp = true) .+ Rbₑ) .* (V̇ * Cf / H)

tc = 7 * 24 * 3600.0

fig2 = Figure()
ax2 = Axis(fig2[1, 1];
    title  = "Fig. 3 — short-term joined to long-term (4×3 borefield)",
    xlabel = "t (y)", ylabel = "g  (-)", xscale = log10)
lines!(ax2, t_st / (3600*24*365*50), g_st; label = "short-term (ANN)", color = :blue)
lines!(ax2, t / (3600*24*365*50), g_lt; label = "FLS long-term (per borehole)", color = :orange,
    linestyle = :dash)
lines!(ax2, t / (3600*24*365*50), g_out; label = "combined — outlet", color = :black,
    linestyle = :dash, linewidth = 2)
lines!(ax2, t / (3600*24*365*50), g_out .+ 1; label = "combined — inlet", color = :red,
    linewidth = 2)
vlines!(ax2, tc / (3600*24*365*50); color = :gray, linestyle = :dot)
axislegend(ax2; position = :lt)
display(fig2)

# Figure 3 — Borefield outlet transfer function: PublishedANN vs DeepANN, same physical input
# Same TRT case and 4×3 borefield as Figure 2, evaluated with both short-term ANNs (pipe
# properties left at the PublishedANN's fixed values, which also lie within DeepANN's wider
# ranges). Each ANN only covers its own short-term horizon (7 days for PublishedANN, 21 days for
# DeepANN) before joining the FLS long-term branch, so the two combined curves are expected to
# diverge between the two contact times rather than match exactly.
tc_deep = 21 * 24 * 3600.0

g_out_published = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s, Rbₑ, m;
    xy = xy, interp = true, model = PublishedANN())
g_out_deep = outlet_transfer_function(t, ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s, Rbₑ, m;
    xy = xy, interp = true, model = DeepANN())

fig3 = Figure()
ax3 = Axis(fig3[1, 1];
    title  = "Borefield outlet transfer function — PublishedANN vs DeepANN (4×3 borefield)",
    xlabel = "t (y)", ylabel = "g  (-)", xscale = log10)
lines!(ax3, t / (3600*24*365*50), g_out_published; label = "PublishedANN", color = :blue,
    linewidth = 2)
lines!(ax3, t / (3600*24*365*50), g_out_deep; label = "DeepANN", color = :red,
    linestyle = :dash, linewidth = 2)
vlines!(ax3, tc / (3600*24*365*50); color = :blue, linestyle = :dot)
vlines!(ax3, tc_deep / (3600*24*365*50); color = :red, linestyle = :dot)
axislegend(ax3; position = :lt)
display(fig3)
