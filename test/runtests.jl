using Test
using GroundHeatExchanger

@testset "GroundHeatExchanger.jl" begin

    # Shared fixture — GHE() provides realistic parameters; Rb is used as a scalar fixture.
    _, H, D, s, rb, ro, ri, T0, ks, kg, kp, kf, Cs, Cg, Cp, Cf, ρs, ρg, ρp, ρf, μf, vD, V = GHE()
    cf = water_cp(T0)
    Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

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

    @testset "ground_load_profile" begin
        t_h = collect(0.0:1.0:8760.0)
        Q_h = ground_load_profile(t_h)
        @test length(Q_h) == length(t_h)
        @test !all(iszero, Q_h)
        @test ground_load_profile(0.0) isa Real    # scalar input
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
    Q_sim  = ground_load_profile(t_sim ./ 3600)
    q_sim  = Q_sim ./ H

    @testset "ground_response" begin
        g_full  = ground_response(t_sim, rb, [0.0 0.0], model; interp=false)
        g_pchip = ground_response(t_sim, rb, [0.0 0.0], model; interp=true)
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
        @test Tf_g ≈ Tf_m
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


    # ---------------------------------------------------------------------------------------
    # short_term_ann — Pasquier, Zarrella & Labib (2018)
    # Validated against the reference MATLAB implementation (ANN_gfunction.m). Five parameter
    # sets spanning the ANN training range (Table 2): min, three interpolated, max. Reference
    # g-functions are the 85-point raw ANN output (before PCHIP interpolation) recomputed with
    # ANN_gfunction.m for the same inputs.
    # ---------------------------------------------------------------------------------------
    @testset "short_term_ann (Pasquier, Zarrella & Labib, 2018)" begin

        # 1 - First data set (min values)
        ks_1=0.5
        Cs_1=1.7e6
        kg_1=0.5
        Cg_1=1.7e6
        kp_1=0.4
        Cp_1=1.9e6
        Cf_1=4.2e6
        ri_1=0.017
        ro_1=0.022
        rb_1= 2ro_1+2e-3 
        H_1=110
        V̇_1=3.34e-4
        s_1=2*(0.02(rb_1-2ro_1)+ro_1)
        dt_1=15 
        tf_1=15
        
        # 2 - Second data set 
        ks_2=1.375
        Cs_2=1.925e6
        kg_2=1.125
        Cg_2=1.925e6
        kp_2=0.4
        Cp_2=1.9e6
        Cf_2=4.2e6
        ri_2=0.017
        ro_2=0.022
        rb_2= 0.0595
        H_2=132.5
        V̇_2=3.75e-4
        s_2=0.06272
        dt_2=15
        tf_2=151211.25
        
        # 3 - Third data set
        ks_3=2.25
        Cs_3=2.15e6
        kg_3=1.75
        Cg_3=2.15e6
        kp_3=0.4
        Cp_3=1.9e6
        Cf_3=4.2e6
        ri_3=0.017
        ro_3=0.022
        rb_3= 0.073
        H_3=155.0
        V̇_3=4.17e-4
        s_3=0.08
        dt_3=43207.5
        tf_3=302407.5
        
        # 4 - Fourth data set
        ks_4=3.125
        Cs_4=2.375e6
        kg_4=2.375
        Cg_4=2.375e6
        kp_4=0.4
        Cp_4=1.9e6
        Cf_4=4.2e6
        ri_4=0.017
        ro_4=0.022
        rb_4= 0.087
        H_4=177.5
        V̇_4=4.59e-4
        s_4=0.09728
        dt_4=64803.75
        tf_4=453603.75
        
        # 5 - Fifth data set (max values)
        ks_5=4
        Cs_5=2.0e6
        kg_5=3
        Cg_5=2.0e6
        kp_5=0.4
        Cp_5=1.9e6
        Cf_5=4.2e6
        ri_5=0.017
        ro_5=0.022
        rb_5=0.1
        H_5=200
        V̇_5=5.0e-4
        s_5=2*(0.98(rb_5-2ro_5)+ro_5)
        dt_5= 86400
        tf_5= 604800

        # MATLAB reference outputs (ANN_gfunction.m, 85-point raw g per case)
        g_raw_1m = [     
           0.000001116570953,
           0.000020782328099,
           0.000085105450193,
           0.000242350790129,
           0.000536214331583,
           0.000857447852279,
           0.001239359543187,
           0.001663643708388,
           0.002111063003516,
           0.002621520355850,
           0.003726658767974,
           0.004277631271503,
           0.005483372979927,
           0.006173439282070,
           0.007622886911145,
           0.009284143698116,
           0.012039272174396,
           0.014881849202524,
           0.016979961243905,
           0.023074374955951,
           0.040932724216621,
           0.092645777244742,
           0.185549848016775,
           0.328960380098967,
           0.421953278012390,
           0.483398862090885,
           0.537034249753091,
           0.596014982981123,
           0.675522770957572,
           0.795638629890888,
           0.889842419951072,
           0.999402441377715,
           1.108016496071149,
           1.231863055897553,
           1.351994315774868,
           1.499123267753769,
           1.631318716638149,
           1.784927203973671,
           1.934195780940208,
           2.104169937672967,
           2.262450581462244,
           2.444842049265199,
           2.617746475414431,
           2.811796340826628,
           2.993556149864516,
           3.202613629962255,
           3.391188077253676,
           3.609584739748685,
           3.808840963597227,
           4.034992469481594,
           4.240775826954676,
           4.473906626144781,
           4.685452884694082,
           4.926202450683492,
           5.143569844596157,
           5.388111967381686,
           5.608767523120203,
           5.858165826552790,
           6.081252239551037,
           6.332297306178447,
           6.557641451138187,
           6.810860244133304,
           7.036502079971269,
           7.288542658779940,
           7.514102705081680,
           7.765411410354652,
           7.989462037553052,
           8.238364119215735,
           8.460131934339381,
           8.707622110564053,
           8.927600707786937,
           9.172478947622594,
           9.390514319447755,
           9.633150590645323,
           9.849400502241481,
          10.090224079059968,
          10.305244624504789,
          10.543806538146432,
          10.756917201650325,
          10.994808990065501,
          11.206790948239503,
          11.443140604785391,
          11.654277531809850,
          11.774675019752896,
          11.887947716165250
        ]
        g_raw_2m = [
           0.000001115977860,
           0.000008458791322,
           0.000001128072118,
                           0,
                           0,
                           0,
           0.000006338663738,
           0.000081884175225,
           0.000196228128120,
           0.000344550703031,
           0.000648855436974,
           0.000799919537598,
           0.001177156539813,
           0.001380823907089,
           0.001836287902005,
           0.002339095822345,
           0.003197267423380,
           0.004203319711232,
           0.004848304013707,
           0.006601616040646,
           0.011263713648270,
           0.030573976753514,
           0.083203108623058,
           0.198028588636056,
           0.284339372048606,
           0.334109310912200,
           0.373092684679785,
           0.410544425364139,
           0.450116290512749,
           0.520694132732791,
           0.584189606535683,
           0.651888925475764,
           0.710847105606868,
           0.780341303827392,
           0.847455178494192,
           0.925512869888060,
           0.994511885369821,
           1.072337966044934,
           1.145891824426999,
           1.227458345763139,
           1.301287292898951,
           1.384150771106126,
           1.460652740329491,
           1.544189837135436,
           1.620786024446998,
           1.706657250861657,
           1.782456065846119,
           1.868504232214574,
           1.945599511483201,
           2.031702027936047,
           2.108801064557226,
           2.194993428847679,
           2.272181029193733,
           2.358772301721867,
           2.436257464945800,
           2.522706163476572,
           2.599955270510014,
           2.686530946426337,
           2.763526238090192,
           2.849590662992047,
           2.926361206359488,
           3.012128240066543,
           3.088415825410388,
           3.173350543855911,
           3.249147765581541,
           3.333467420369248,
           3.408558155055272,
           3.492145095161220,
           3.566613739895278,
           3.649613112079702,
           3.723507298177323,
           3.805821152824054,
           3.879269653161861,
           3.961128246108818,
           4.034197192969684,
           4.115604595028980,
           4.188317477070671,
           4.269429704544772,
           4.341830466992732,
           4.422671646128210,
           4.494873688151126,
           4.575440902177598,
           4.647456214200107,
           4.688546498153606,
           4.727242546359512
        ]
        g_raw_3m = [
            0.000001115793989,
            0.000004920134532,
                           0,
                           0,
                           0,
                           0,
                           0,
           0.000009454057975,
           0.000038497761107,
           0.000079621978218,
           0.000167771591269,
           0.000206695130889,
           0.000337517330549,
           0.000409026555273,
           0.000585720306925,
           0.000787627291345,
           0.001123700473710,
           0.001509225605831,
           0.001773340114998,
           0.002521511883032,
           0.004135775326935,
           0.011804553376055,
           0.040645592775917,
           0.128784185339416,
           0.211989452646398,
           0.259051323219539,
           0.291145377455341,
           0.320330295751946,
           0.345899983221785,
           0.389995275950932,
           0.436826805438582,
           0.487909013432777,
           0.527742469078487,
           0.572703480610697,
           0.618470290202983,
           0.669678676325729,
           0.714374811667430,
           0.764560526752670,
           0.811362634029197,
           0.862835988146816,
           0.909006134709423,
           0.960370358088265,
           1.007394222158855,
           1.058335685134857,
           1.104664032484845,
           1.156234908160187,
           1.201479838082709,
           1.252512439657222,
           1.297981999577269,
           1.348536820535001,
           1.393640821224185,
           1.443870889407203,
           1.488751011498426,
           1.539027938502423,
           1.583947698791546,
           1.634004288373935,
           1.678737640279749,
           1.728849756248192,
           1.773430600369059,
           1.823252468989272,
           1.867729935151286,
           1.917415738863161,
           1.961625744534698,
           2.010861936163202,
           2.054835574136186,
           2.103750314565380,
           2.147349669988853,
           2.195869324029892,
           2.239103113588846,
           2.287315808269161,
           2.330256977543068,
           2.378113084163656,
           2.420823424979291,
           2.468412820956364,
           2.510925181778493,
           2.558272587627312,
           2.600596904462089,
           2.647825114962506,
           2.689986543243635,
           2.737079430686568,
           2.779149757025833,
           2.826119382469891,
           2.868106267874902,
           2.892063411891538,
           2.914641430166469
        ]
        g_raw_4m = [
           0.000001115739243,
           0.000003976843701,
                           0,
                           0,
           0.000025918240050,
           0.000033425142371,
           0.000038475151960,
           0.000024879175408,
           0.000015880862659,
           0.000012601057462,
           0.000026021112071,
           0.000023575458457,
           0.000067360852723,
           0.000097516145504,
           0.000167883063943,
           0.000260390922080,
           0.000409268619823,
           0.000575703336745,
           0.000740186497769,
           0.001054670817416,
           0.001598932582411,
           0.005075552754506,
           0.020763753783113,
           0.086817858310772,
           0.165372470999545,
           0.213276407804316,
           0.241780305133472,
           0.265247935429876,
           0.284152364001830,
           0.313439717895375,
           0.349265366045489,
           0.390476376480662,
           0.420071346908411,
           0.451599888179016,
           0.485357832232211,
           0.522483631144778,
           0.554375878186059,
           0.590338917362047,
           0.623733824531577,
           0.660388943720859,
           0.693217871344164,
           0.729653537540839,
           0.762951473890161,
           0.798945324416238,
           0.831589738777576,
           0.867790511677107,
           0.899476396403140,
           0.935083093224930,
           0.966722718355832,
           1.001806331148140,
           1.033054474109667,
           1.067770228757087,
           1.098752887171345,
           1.133433165734059,
           1.164430700201379,
           1.198953396123335,
           1.229833686596721,
           1.264414046266514,
           1.295218483760938,
           1.329667519758504,
           1.360429242763522,
           1.394843877378452,
           1.425476653509231,
           1.459588168835318,
           1.490091191459177,
           1.524032804537892,
           1.554299790323954,
           1.587973352610745,
           1.617985563086897,
           1.651471155863977,
           1.681294102726100,
           1.714526362320836,
           1.744180982232893,
           1.777210091395886,
           1.806710762027635,
           1.839579434377721,
           1.868933617982615,
           1.901677497565531,
           1.930921189583336,
           1.963578461764204,
           1.992734692532257,
           2.025278734223289,
           2.054372995729974,
           2.070978900129203,
           2.086618017297788
        ]
        g_raw_5m = [
            0.000001115775265,
            0.000004879898726,
            0.000021521501362,
            0.000045652821692,
            0.000012468017308,
                           0,
                           0,
                           0,
                           0,
                           0,
                           0,
           0.000025774987848,
           0.000088183922578,
           0.000124815044826,
           0.000139123634165,
           0.000148109086153,
           0.000135399781494,
           0.000025932249457,
                           0,
           0.000215466180049,
           0.000338057011357,
           0.001620636594693,
           0.010660225221072,
           0.059019613290934,
           0.132083467239161,
           0.184069016549541,
           0.214854286131434,
           0.236813047292793,
           0.252323590270083,
           0.272832581073185,
           0.300021956626774,
           0.334614863078190,
           0.357844985503260,
           0.380022328599366,
           0.403270716730318,
           0.428923458511423,
           0.449745065200867,
           0.473088376214690,
           0.494496585306366,
           0.517935745494408,
           0.539042482639495,
           0.562693537967501,
           0.584720043995369,
           0.609083400101245,
           0.631615114279325,
           0.657155548327890,
           0.679999849905383,
           0.706064499185493,
           0.729621033148482,
           0.755933867945722,
           0.779526106301970,
           0.805771035292445,
           0.829219131462613,
           0.855383411866505,
           0.878719371150520,
           0.904623076550657,
           0.927732625988247,
           0.953519331192912,
           0.976470576128820,
           1.002199427965497,
           1.025131757632217,
           1.050946757981725,
           1.073903864619438,
           1.099523499446625,
           1.122419030176211,
           1.147943198348941,
           1.170706545834038,
           1.196067695798786,
           1.218707468001210,
           1.243875617227930,
           1.266331201129403,
           1.291372327983105,
           1.313697466674097,
           1.338594572075022,
           1.360788557715543,
           1.385625358884775,
           1.407848239013908,
           1.432584303559083,
           1.454745376165156,
           1.479426696367302,
           1.501531775340319,
           1.526233059679295,
           1.548312414821026,
           1.560933948284209,
           1.572839212396244
        ]

        datasets = [
            (ks_1, Cs_1, kg_1, Cg_1, kp_1, Cp_1, Cf_1, ri_1, ro_1, rb_1, H_1, V̇_1, s_1),
            (ks_2, Cs_2, kg_2, Cg_2, kp_2, Cp_2, Cf_2, ri_2, ro_2, rb_2, H_2, V̇_2, s_2),
            (ks_3, Cs_3, kg_3, Cg_3, kp_3, Cp_3, Cf_3, ri_3, ro_3, rb_3, H_3, V̇_3, s_3),
            (ks_4, Cs_4, kg_4, Cg_4, kp_4, Cp_4, Cf_4, ri_4, ro_4, rb_4, H_4, V̇_4, s_4),
            (ks_5, Cs_5, kg_5, Cg_5, kp_5, Cp_5, Cf_5, ri_5, ro_5, rb_5, H_5, V̇_5, s_5),
        ]
        references = [g_raw_1m, g_raw_2m, g_raw_3m, g_raw_4m, g_raw_5m]

        @testset "Case $i vs MATLAB reference" for (i, (p, ref)) in enumerate(zip(datasets, references))
            _, g_raw = short_term_nodes(p...)
            @test length(g_raw) == 85
            @test isapprox(g_raw, ref)
        end

        @testset "short_term_response matches short_term_nodes at the native ANN times" begin
            ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s = datasets[3]
            ts_nodes, g_nodes = short_term_nodes(ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s)
            _, g_req = short_term_response(collect(Float64, ts_nodes), rb, ri, ro, H, s, V̇,
                ks, Cs, kg, Cg, kp, Cp, Cf)
            @test isapprox(g_req, g_nodes)
        end

        @testset "clamp policy on out-of-range inputs" begin
            ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, V̇, s = datasets[3]
            @test_throws ArgumentError short_term_response(3600.0, rb, ri, ro, H, s, V̇, 100.0, Cs,
                kg, Cg, kp, Cp, Cf; clamp = false)
            _, g_clamped = short_term_response(3600.0, rb, ri, ro, H, s, V̇, 100.0, Cs, kg, Cg, kp,
                Cp, Cf; clamp = true)
            @test isfinite(g_clamped)
        end
    end


end
