using Test
using GroundHeatExchanger

@testset "GroundHeatExchanger.jl" begin

    # Shared fixture — GHE() provides realistic parameters; Rb is used as a scalar fixture.
    _, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, ϵ, vD, V = GHE()
    cf = water_cp(T0)
    Rb = resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

    @testset "set_nodes" begin
        id = set_nodes(1000, 50)
        @test length(id) == 50
        @test id[1]   == 1
        @test id[end] == 1000
        @test issorted(id)
        @test allunique(id)
    end

    @testset "pchip_interpolation" begin
        tᵢ  = [1.0, 2.0, 4.0, 8.0, 10.0]
        vᵢ  = 3.0 .* tᵢ .- 1.0
        t_q = collect(1.0:0.25:10.0)
        @test isapprox(pchip_interpolation(tᵢ, vᵢ, t_q), 3.0 .* t_q .- 1.0, atol=1)
        @test isapprox(pchip_interpolation(tᵢ, vᵢ, tᵢ), vᵢ, atol=1)
    end

    @testset "GHE" begin
        @test H > 0 && D > 0 && rb > 0 && ks > 0
        @test rb > ro > ri > 0
        @test V > 0 && Cf > 0
    end

    @testset "heat_load_profile" begin
        t_h = collect(0.0:1.0:8760.0)
        Q_h = heat_load_profile(t_h)
        @test length(Q_h) == length(t_h)
        @test !all(iszero, Q_h)
        @test heat_load_profile(0.0) isa Real    # scalar input
    end

    @testset "head_loss_Darcy_Weisbach" begin
        L_p, r_p, u_p, f_p = 100.0, 0.02, 1.0, 0.02
        Δh = head_loss_Darcy_Weisbach(L_p, r_p, u_p, f_p)
        @test Δh ≈ f_p * (L_p / (2r_p)) * (u_p^2 / (2 * 9.81))
        @test head_loss_Darcy_Weisbach(2L_p, r_p, u_p, f_p) ≈ 2Δh   # linear in L
        @test head_loss_Darcy_Weisbach(L_p, r_p, 2u_p, f_p) ≈ 4Δh   # quadratic in u
    end

    @testset "impulse_func" begin
        q_ramp  = cumsum(rand(Float64, 30))
        f_ramp  = impulse_func(q_ramp)
        @test length(f_ramp) == 30
        @test cumsum(f_ramp) ≈ q_ramp       # cumulative sum of increments = load

        q_const = fill(5.0, 20)
        f_const = impulse_func(q_const)
        @test f_const[1]        == 5.0
        @test all(iszero, f_const[2:end])    # constant load: only first increment is non-zero
    end

    @testset "convolution" begin
        n = 30
        q = cumsum(rand(Float64, n))
        g = collect(Float64, 1:n) ./ n

        # Convolution with ones returns the original load
        @test convolution(q, ones(Float64, n)) ≈ q

        # Constant load: f = [q₀, 0, …], so convolution reduces to q₀·g
        q_c = fill(3.0, n)
        @test convolution(q_c, g) ≈ q_c[1] .* g

        # Agreement with the time-domain definition for small n
        function time_domain_conv(q, g)
            f = impulse_func(q)
            n = length(q); y = zeros(n)
            for i in 1:n, j in 1:i
                y[i] += f[j] * g[i - j + 1]
            end
            return y
        end
        n_s = 12
        q_s = cumsum(rand(Float64, n_s))
        g_s = collect(Float64, 1:n_s)
        @test convolution(q_s, g_s) ≈ time_domain_conv(q_s, g_s)
    end

    @testset "convolutionf" begin
        n = 25
        q = cumsum(rand(Float64, n))
        g = collect(Float64, 1:n) ./ n
        @test convolutionf(impulse_func(q), g) ≈ convolution(q, g)
        @test length(convolutionf(impulse_func(q), g)) == n
    end

    @testset "impulse_func_ns" begin
        n_ns = 20; n_h = n_ns ÷ 2
        q_ns = cumsum(rand(Float64, n_ns))

        # Single state: matches stationary impulse_func
        f_ns1 = impulse_func_ns(q_ns, ones(Int, n_ns))
        @test size(f_ns1) == (n_ns, 1)
        @test vec(f_ns1) ≈ impulse_func(q_ns)

        # Two states: masking and transition impulses
        sv2 = vcat(fill(1, n_h), fill(2, n_ns - n_h))
        f2  = impulse_func_ns(q_ns, sv2)
        @test size(f2) == (n_ns, 2)
        @test all(iszero, f2[1:n_h, 2])            # state 2 inactive in first half
        @test f2[n_h + 1, 2] ≈  q_ns[n_h + 1]     # startup impulse for state 2
        @test f2[n_h + 1, 1] ≈ -q_ns[n_h]          # shutoff impulse for state 1
    end

    @testset "state_vector and state_indices" begin
        # state_vector: docstring example
        @test state_vector([1, 4, 6], [1, 2, 3], 6) == [1, 1, 1, 2, 2, 3]

        # state_indices → state_vector round-trip
        v = [1.0, 1.0, 1.0, 2.0, 2.0, 3.0]
        ind_v, state_v, _ = state_indices(v)
        @test state_vector(ind_v, state_v, length(v)) == [1, 1, 1, 2, 2, 3]

        # Re-entry of a repeated state is consistent (same value → same state index)
        v2 = [1.0, 1.0, 2.0, 2.0, 1.0, 1.0]
        ind2, state2, _ = state_indices(v2)
        sv2 = state_vector(ind2, state2, length(v2))
        @test sv2[1] == sv2[5]    # same value in v2 → same state
        @test sv2[1] != sv2[3]    # different value → different state
    end

    @testset "convolution_ns — overloads agree" begin
        n = 40; n_h = n ÷ 2
        q   = cumsum(rand(Float64, n))
        sv  = vcat(fill(1, n_h), fill(2, n - n_h))
        g1  = collect(Float64, 1:n) ./ n
        g2m = hcat(g1, reverse(g1))   # two distinct transfer functions
        f2  = impulse_func_ns(q, sv)

        # Single state: NS convolution equals stationary convolution
        f_single = impulse_func_ns(q, ones(Int, n))
        @test vec(convolution_ns(f_single, reshape(g1, n, 1))) ≈ convolution(q, g1)

        # All four overloads agree
        dT_ref = similar(q); convolution_ns!(dT_ref, f2, g2m)
        @test convolution_ns(f2, g2m)  ≈ dT_ref
        @test convolution_ns(q, g2m, sv) ≈ dT_ref
        ind_ns, st_ns, _ = state_indices(Float64.(sv))
        @test convolution_ns(q, g2m, ind_ns, st_ns) ≈ dT_ref
    end

    t_sim  = collect(3600.0:3600:3600 * 24 * 30)   # 30 days, hourly
    n_sim  = length(t_sim)
    model  = FLSModel(H, D, ks, Cs)
    g_sim  = ground_response(t_sim, rb, [0.0 0.0], model)
    Q_sim  = heat_load_profile(t_sim ./ 3600)
    q_sim  = Q_sim ./ H

    @testset "ground_response wrapper" begin
        g_full  = ground_response(t_sim, rb, [0.0 0.0], model; n_nodes=0)
        g_pchip = ground_response(t_sim, rb, [0.0 0.0], model; n_nodes=30)
        @test length(g_full)  == n_sim
        @test length(g_pchip) == n_sim
        @test issorted(g_full)   # FLS g-function is monotone increasing
        rmse = sqrt(sum((g_full .- g_pchip).^2) / n_sim)
        @test rmse < 1e-2        # PCHIP interpolation error is small
    end

    @testset "fluid_temperature" begin
        # Zero load → Tf = T0 everywhere
        @test fluid_temperature(t_sim, zeros(n_sim), g_sim, T0, ks, Rb) ≈ fill(T0, n_sim)

        # Zero g → Tf = T0 + q·Rb (no ground conduction)
        @test fluid_temperature(t_sim, q_sim, zeros(n_sim), T0, ks, Rb) ≈ T0 .+ q_sim .* Rb

        # h-formulation equivalence: T0 + conv(q, Rb + g/(2π·ks)) == fluid_temperature
        h_nom = Rb .+ g_sim ./ (2π * ks)
        @test T0 .+ convolution(q_sim, h_nom) ≈ fluid_temperature(t_sim, q_sim, g_sim, T0, ks, Rb)

        # All overloads agree
        Tf_g = fluid_temperature(t_sim, q_sim, g_sim, T0, ks, Rb)
        Tf_m = fluid_temperature(t_sim, q_sim, model, rb, T0, ks, Rb)
        Tf_Q = fluid_temperature(t_sim, Q_sim, g_sim, H, T0, ks, Rb)
        @test Tf_g ≈ Tf_m
        @test Tf_g ≈ Tf_Q
    end

    @testset "outlet and inlet temperature" begin
        Tf_sim = fluid_temperature(t_sim, q_sim, g_sim, T0, ks, Rb)
        Tout   = outlet_temperature(Tf_sim, Q_sim, V, Cf)
        Tin    = inlet_temperature(Tf_sim, Q_sim, V, Cf)

        # Mean of Tin and Tout equals Tf
        @test (Tin .+ Tout) ./ 2 ≈ Tf_sim

        # Energy balance: Q = V·C·(Tin − Tout)
        @test Tin .- Tout ≈ Q_sim ./ (V * Cf)

        # q/H overloads agree with Q overloads
        @test outlet_temperature(Tf_sim, q_sim, H, V, Cf) ≈ Tout
        @test inlet_temperature(Tf_sim, q_sim, H, V, Cf)  ≈ Tin
    end

end
