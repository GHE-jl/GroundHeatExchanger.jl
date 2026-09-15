# Validation of GroundHeatExchanger.jl (and the BoreholeResistance / GroundResponse functions it
# re-exports) against the worked examples of Lamarche (2023).
#
# Reference
#   - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application.
#       Springer Nature Switzerland. https://books.google.ca/books?id=3erMzwEACAAJ

using Test
using GroundHeatExchanger

@testset "Validation — Lamarche (2023)" begin

    # Reynolds number — Chap. 2, p.20, ex. 2.1
    @testset "Reynolds number — Chap. 2, ex. 2.1" begin
        ρf = 997                      # density [kg/m³]
        μf = 8.54e-4                  # viscosity [kg/m⋅s]
        r  = 0.01                     # pipe radius [m]
        Ac = π/4 * (0.04^2 - 0.02^2)  # annulus area [m²]
        V̇  = 0.5 / (ρf * Ac)          # fluid speed [m/s]

        Re_ref = 12428
        Re = Reynolds(V̇, r, ρf, μf)
        @test Re_ref ≈ Re rtol = 1e-3
    end

    # ILS and ICS transfer functions — Chap. 4, p.60, ex. 4.1
    @testset "ILS and ICS — Chap. 4, ex. 4.1" begin
        q0 = -40           # heat load [W/m]
        α  = 0.1/(24*3600) # thermal diffusivity [m²/s]
        t  = 250*24*3600   # time [s]
        r  = 2             # radial position [m]
        rc = 0.1           # cylinder radius [m]
        ks = 2.5           # thermal conductivity [W/m-K]
        Cs = ks/α          # volumetric heat capacity [J/m³-K]

        ΔT_ils = -q0 * ils(t, r, ks, Cs)
        ΔT_ics = -q0 * ics(t, r, rc, ks, Cs)

        @test 3.414 ≈ ΔT_ils rtol = 1e-3
        @test 3.412 ≈ ΔT_ics rtol = 1e-3
    end

    # ILS, ICS and FLS g-functions, and the resulting step-load transfer function — Chap. 4, p.66,
    # ex. 4.2
    @testset "ILS, ICS and FLS — Chap. 4, ex. 4.2" begin
        α  = 0.1 / (24*3600)          # thermal diffusivity [m²/s]
        ks = 2                        # thermal conductivity [W/m-K]
        Cs = ks/α                     # volumetric heat capacity [J/m³-K]
        rc = 0.075                    # cylinder radius [m]
        H  = 100                      # borehole depth [m]
        D  = 0                        # buried depth [m]
        q  = [1000, 1800, -200] / 100 # heat loads [W/m]
        t  = 3600 * [1, 3, 5]         # time steps [s]

        g_ics = [
            ks * ics(t[3], rc, rc, ks, Cs),
            ks * ics(t[3] - t[1], rc, rc, ks, Cs),
            ks * ics(t[3] - t[2], rc, rc, ks, Cs),
        ]
        g_ils = [
            ks * ils(t[3], rc, ks, Cs),
            ks * ils(t[3] - t[1], rc, ks, Cs),
            ks * ils(t[3] - t[2], rc, ks, Cs),
        ]
        g_fls = [
            ks * fls(t[3], rc, H, D, ks, Cs),
            ks * fls(t[3] - t[1], rc, H, D, ks, Cs),
            ks * fls(t[3] - t[2], rc, H, D, ks, Cs),
        ]

        @test [0.198, 0.1849, 0.1469]  ≈ g_ics rtol = 1e-3
        @test [0.1739, 0.1574, 0.1085] ≈ g_ils rtol = 1e-3
        @test [0.1737, 0.1573, 0.1085] ≈ g_fls rtol = 1e-3

        # Superposed step-load transfer function, ΔT = -1/ks * Σ Δq·g
        transfer_function(g) = -1/ks * (q[1]*g[1] + (q[2] - q[1])*g[2] + (q[3] - q[2])*g[3])

        @test -0.26 ≈ transfer_function(g_ics) rtol = 1e-2
        @test -0.41 ≈ transfer_function(g_ils) rtol = 1e-2
        @test -0.41 ≈ transfer_function(g_fls) rtol = 1e-2
    end

    # Water properties and single U-loop borehole resistance — Chap. 5, p.75, ex. 5.1
    @testset "Water properties and Rb (single U-loop) — Chap. 5, ex. 5.1" begin
        kf, cf, ρf, μf = fluid_property(30, :water) # k [W/m-K], cp [J/kg-K], ρ [kg/m³], μ [kg/m-s]
        @test 0.614  ≈ kf rtol = 1e-2
        @test 8.0e-4 ≈ μf rtol = 1e-2
        @test 996    ≈ ρf rtol = 1e-2

        s  = 0.05               # distance between the 2 legs of the U-tube [m]
        rb = 0.075               # borehole radius [m]
        ro = 0.021                # outer radius of the pipe [m]
        ri = 0.0172               # inner radius of the pipe [m]
        V̇  = 0.0004 / (π*ri^2)    # fluid speed in the pipe [m/s]
        ks = 2.5                  # soil conductivity [W/m-K]
        kg = 1.7                  # grout conductivity [W/m-K]
        kp = 0.4                  # pipe conductivity [W/m-K]

        Rf  = resistance_fluid(V̇, ri, kf, cf, ρf, μf, 5e-6)                              # [m-K/W]
        Rp  = resistance_pipe(ro, ri, kp)                                                # [m-K/W]
        Rg  = resistance_ULoop_borehole(s, rb, ro, ks, kg, 0, 0, nLoop=1, order=0)       # [m-K/W]
        Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0)     # line source
        Rb₁ = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1)     # multipole

        @test 0.00415 ≈ Rf  rtol = 1e-2
        @test 0.0794  ≈ Rp  rtol = 1e-2
        @test 0.07845 ≈ Rg  rtol = 1e-2
        @test 0.120   ≈ Rb  rtol = 1e-2
        @test 0.120   ≈ Rb₁ rtol = 1e-2
    end

    # Grout and borehole resistance with the line-source method, for U-tube configurations A
    # (touching pipes), B (evenly spaced) and C (pipes against the borehole wall) — Chap. 5, p.76,
    # ex. 5.2
    @testset "Rg and Rb (configs A/B/C) — Chap. 5, ex. 5.2" begin
        rb = 0.075                      # borehole radius [m]
        ro = 0.033/2                    # outer radius of the pipe [m]
        ri = 0.027/2                    # inner radius of the pipe [m]
        ks = 2.5                        # soil conductivity [W/m-K]
        kg = 1.0                        # grout conductivity [W/m-K]
        kp = 0.4                        # pipe conductivity [W/m-K]
        Rp = resistance_pipe(ro, ri, kp)              # pipe resistance [m-K/W]
        Rf = 0.08 - Rp                                # fluid resistance [m-K/W]
        s  = [2*ro, 2*rb / 3, 2*(rb-ro)]              # leg spacing for cases A, B, C [m]

        Rg = [resistance_ULoop_borehole(si, rb, ro, ks, kg, 0, 0, nLoop=1, order=0) for si in s]
        Rb = [resistance_ULoop_borehole(si, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0) for si in s]

        @test [0.186, 0.152, 0.069] ≈ Rg rtol = 1e-2
        @test [0.226, 0.192, 0.109] ≈ Rb rtol = 1e-2
    end

    # Internal resistance with the line-source and multipole methods, same U-tube configurations
    # A/B/C — Chap. 5, p.77, ex. 5.3
    @testset "Internal resistance Ra (configs A/B/C) — Chap. 5, ex. 5.3" begin
        rb = 0.075                              # borehole radius [m]
        ro = 0.033/2                            # outer radius of the pipe [m]
        ri = 0.027/2                            # inner radius of the pipe [m]
        ks = 2.5                                # soil conductivity [W/m-K]
        kg = 1.0                                # grout conductivity [W/m-K]
        kp = 0.4                                # pipe conductivity [W/m-K]
        Rp = resistance_pipe(ro, ri, kp)                       # pipe resistance [m-K/W]
        Rf = 0.08 - Rp                                         # fluid resistance [m-K/W]
        s  = [2*ro + 0.001, 2*rb / 3, 2*(rb-ro)]               # leg spacing for cases A, B, C [m]

        Rals   = [resistance_ULoop_total_internal(si, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0,
                      network="diagonal") for si in s]
        Ramult = [resistance_ULoop_total_internal(si, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1,
                      network="diagonal") for si in s]

        @test [0.3761, 0.482, 0.592] ≈ Rals   rtol = 1e-2  # line-source method
        @test [0.353, 0.475, 0.592]  ≈ Ramult rtol = 1e-2  # multipole method
    end

    # Double U-loop borehole and internal resistance with the line-source method — Chap. 5, p.80,
    # ex. 5.4
    @testset "Double U-loop Rb and Ra — Chap. 5, ex. 5.4" begin
        rb = 0.075                       # borehole radius [m]
        ro = 0.033/2                     # outer radius of the pipe [m]
        ri = 0.027/2                     # inner radius of the pipe [m]
        ks = 2.5                         # soil conductivity [W/m-K]
        kg = 1.0                         # grout conductivity [W/m-K]
        kp = 0.4                         # pipe conductivity [W/m-K]
        Rp = resistance_pipe(ro, ri, kp) # pipe resistance [m-K/W]
        Rf = 0.08 - Rp                   # fluid resistance [m-K/W]
        s  = 2*rb / 3                    # leg spacing, config B [m]

        Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=2, order=0)
        Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf, nLoop=2, order=0,
            network="diagonal")

        @test 0.156 ≈ Rb rtol = 1e-2
        # The diagonal-network line-source formula runs ~20% above the book's value for Ra.
        @test 0.241 ≈ Ra rtol = 0.22
    end

    # Borehole temperature from a linear (ILS) approximation, plus the associated grout, borehole
    # and internal resistances — Chap. 5, p.84, ex. 5.5
    @testset "Linear approximation (Tb, Rg, Rb, Ra) — Chap. 5, ex. 5.5" begin
        t  = 100*3600  # time [s]
        T₀ = 10        # initial temperature [°C]
        rb = 0.075     # borehole radius [m]
        ro = 0.0335/2  # outer radius of the pipe [m]
        ks = 2.1       # soil conductivity [W/m-K]
        kg = 1.5       # grout conductivity [W/m-K]
        kp = 0.4       # pipe conductivity [W/m-K]
        Cg = 2.43e6    # ground volumetric heat capacity [J/m³-K]
        s  = 2*0.055   # distance between the 2 legs of the U-tube [m]
        q  = -50       # heat load [W/m]

        Tb = T₀ + q * ils(t, rb, ks, Cg)
        @test 0.856 ≈ Tb rtol = 1e-2

        V̇ = 0.00019 / (π*ro^2)                            # fluid speed [m/s]
        kf, _, _, μf = fluid_property(T₀, :water)          # k [W/m-K], μ [kg/m⋅s]
        cf = 4200                                          # fluid volumetric heat capacity [J/m³-K]
        ρf = 1000                                          # fluid density [kg/m³]
        Rf = resistance_fluid(V̇, ro, kf, cf, ρf, μf, 5e-6) # fluid resistance [m-K/W]
        Rp = 0.076 - Rf                                    # pipe resistance [m-K/W]

        Rg   = resistance_ULoop_borehole(s, rb, ro, ks, kg, 0, 0, nLoop=1, order=0)
        Rb   = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0)
        Rals = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=0,
            network="diagonal")

        @test 0.0562 ≈ Rg   rtol = 1e-2
        @test 0.094  ≈ Rb   rtol = 1e-2
        @test 0.509  ≈ Rals rtol = 1e-2
    end

end
