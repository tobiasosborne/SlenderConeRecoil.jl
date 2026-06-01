using Test
using SlenderConeRecoil

@testset "Capillary-wave phase and envelope diagnostics" begin
    epsilon = 0.1
    lambda = 2.0
    decay = 0.05
    phase = 0.3
    xi = collect(range(0.0, 20.0, length=801))
    excess = exp.(-decay .* xi) .* sin.(2pi .* xi ./ lambda .+ phase)
    S = epsilon .* xi .+ excess

    @testset "Resolved local implementation wave train" begin
        diagnostics = wave_diagnostics(xi, S; epsilon=epsilon)

        @test diagnostics.successful
        @test diagnostics.status == :ok
        @test diagnostics.profile_kind == :sampled_radius
        @test diagnostics.source_status == "IMPL-inferred"
        @test diagnostics.quantitative_source_status == "IMPL-inferred"
        @test diagnostics.source_basis.qualitative.source_id ==
              "DK2001-figure1-oscillation-qualitative"
        @test diagnostics.source_basis.qualitative.source_status == "C2001"
        @test diagnostics.source_basis.qualitative.fact_type == :qualitative
        @test diagnostics.zero_crossings.count >= 18
        @test diagnostics.extrema.crest_count >= 9
        @test diagnostics.extrema.trough_count >= 9
        @test diagnostics.wavelength.zero_crossing.median ≈ lambda rtol=0.03
        @test diagnostics.phase.status == :ok
        @test diagnostics.phase.zero_crossing_phase[2] -
              diagnostics.phase.zero_crossing_phase[1] ≈ pi
        @test diagnostics.envelope.status == :ok
        @test diagnostics.envelope.decay_rate ≈ decay rtol=0.12
        @test diagnostics.envelope.decay_length ≈ 1 / decay rtol=0.12
        @test diagnostics.resolution.status == :sufficient_resolution
        @test diagnostics.resolution.min_samples_per_wavelength > 50
        @test diagnostics.qualitative.oscillatory
        @test diagnostics.qualitative.modulated_on_longer_scale
        @test diagnostics.qualitative.source_consistency ==
              :consistent_with_C2001_qualitative
    end

    @testset "Solution overloads use the documented excess convention" begin
        Sxi = fill(epsilon, length(xi))
        Sxixi = zeros(length(xi))
        U = zeros(length(xi))
        inner = InnerSolution(xi, S, Sxi, Sxixi, U, first(xi), first(S), 0.0)
        outer = OuterSolution(xi, excess, zeros(length(xi)), zeros(length(xi)),
                              zeros(length(xi)), epsilon)
        composite = CompositeSolution(xi, S, U, epsilon)

        inner_diag = wave_diagnostics(inner; epsilon=epsilon)
        outer_diag = wave_diagnostics(outer)
        composite_diag = wave_diagnostics(composite)

        @test inner_diag.profile_kind == :inner_solution
        @test inner_diag.epsilon_source == :provided
        @test outer_diag.profile_kind == :outer_solution
        @test outer_diag.epsilon_source == :outer_solution
        @test composite_diag.profile_kind == :composite_solution
        @test composite_diag.epsilon_source == :composite_solution
        @test inner_diag.wavelength.zero_crossing.median ≈
              outer_diag.wavelength.zero_crossing.median rtol=1e-12
        @test composite_diag.envelope.decay_rate ≈ decay rtol=0.12
    end

    @testset "Invalid inputs fail loudly" begin
        @test_throws ArgumentError wave_diagnostics(xi[1:4], S[1:4];
                                                    epsilon=epsilon)
        @test_throws ArgumentError wave_diagnostics([0.0, 1.0, 1.0, 2.0, 3.0],
                                                    S[1:5]; epsilon=epsilon)
        @test_throws ArgumentError wave_diagnostics([0.0, 2.0, 1.0, 3.0, 4.0],
                                                    S[1:5]; epsilon=epsilon)
        bad_S = copy(S[1:5])
        bad_S[3] = Inf
        @test_throws ArgumentError wave_diagnostics(xi[1:5], bad_S;
                                                    epsilon=epsilon)
        @test_throws ArgumentError wave_diagnostics(xi[1:5], S[1:5];
                                                    epsilon=0.0)
        @test_throws ArgumentError wave_diagnostics(xi[1:5], S[1:5];
                                                    epsilon=-0.1)
        @test_throws ArgumentError wave_diagnostics(xi[1:6], S[1:5];
                                                    epsilon=epsilon)
    end

    @testset "Undersampled or nonoscillatory profiles report diagnostic status" begin
        coarse_xi = collect(0.0:1.0:40.0)
        coarse_excess = sin.(2pi .* coarse_xi ./ 4.0)
        coarse_S = epsilon .* coarse_xi .+ coarse_excess
        coarse = wave_diagnostics(coarse_xi, coarse_S; epsilon=epsilon,
                                  min_samples_per_wavelength=8.0)

        @test !coarse.successful
        @test coarse.status == :insufficient_resolution
        @test coarse.zero_crossings.status == :ok
        @test coarse.resolution.status == :insufficient_resolution
        @test coarse.resolution.min_samples_per_wavelength < 8.0

        smooth_xi = collect(range(0.0, 10.0, length=101))
        smooth_S = epsilon .* smooth_xi .+ exp.(-smooth_xi)
        smooth = wave_diagnostics(smooth_xi, smooth_S; epsilon=epsilon)

        @test !smooth.successful
        @test smooth.status == :insufficient_wave_train
        @test smooth.zero_crossings.status == :insufficient_zero_crossings
        @test !smooth.qualitative.oscillatory
        @test smooth.qualitative.source_consistency ==
              :not_established_by_current_profile
    end
end
