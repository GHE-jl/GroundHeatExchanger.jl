"""
Script testing the analytical model for standing column well (SCW) from GHEModels.jl.
"""

includet("../src/GHEModels.jl")

using Plots
using .GHEModels

function test_scwm()
    # Parameters
    t = collect(exp10.(range(log10(60.0), log10(3600 * 24 * 365 * 100), length=500)))
    ks = 3.
    Cs = 2.11e6
    rb = 0.08
    H = 150.0
    V = 30. / 60000
    β = 0.05

    # Run model
    g = scwm(t, ks, Cs, rb, H, V, β)

    # Convolution to see temperatures
    q = 10000 * ones(length(t)) / H
    T = convolution(diff([0; q]), g) .+ 10

    return t, g, T
end

t, g, T = test_scwm()

plot(t, g; xscale=:log10)