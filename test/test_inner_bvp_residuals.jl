using Test
using SlenderConeRecoil
using LinearAlgebra: norm

@testset "Inner BVP residual API" begin
    @testset "Unknown variables are explicit and validated" begin
        unknowns = InnerBVPUnknowns(2.79, 0.28, 0.57)
        @test unknowns.ξ₀ ≈ 2.79
        @test unknowns.S₀ ≈ 0.28
        @test unknowns.Sξξ₀ ≈ 0.57
        @test InnerBVPUnknowns([2.79, 0.28, 0.57]) == unknowns
        @test InnerBVPUnknowns(; ξ₀=2.79, S₀=0.28, Sξξ₀=0.57) == unknowns

        @test_throws ArgumentError InnerBVPUnknowns(-1.0, 0.28, 0.57)
        @test_throws ArgumentError InnerBVPUnknowns(2.79, 0.0, 0.57)
        @test_throws ArgumentError InnerBVPUnknowns(2.79, 0.28, Inf)
        @test_throws ArgumentError InnerBVPUnknowns([2.79, 0.28])
    end

    @testset "Residual components match the current shooting baseline" begin
        unknowns = InnerBVPUnknowns(2.79, 0.28, 0.57)
        residual = inner_bvp_residual(unknowns; ε=0.1, ξ_max=6.0)
        components = inner_bvp_residual_components(residual)
        vector = inner_bvp_residual_vector(residual)

        ξv, Sv, _, Sppv, Uv = SlenderConeRecoil._shoot(2.79, 0.28, 0.57, 6.0)
        expected = [Sv[end] / ξv[end] - 0.1, Uv[end], Sppv[end]]

        @test residual isa InnerBVPResidual
        @test residual.unknowns == unknowns
        @test residual.diagnostics.component_names ==
              (:far_field_slope, :far_field_velocity, :far_field_curvature)
        @test keys(components) == residual.diagnostics.component_names
        @test vector ≈ expected
        @test inner_bvp_residual_vector(unknowns; ε=0.1, ξ_max=6.0) ≈ vector
        @test components.far_field_slope ≈ expected[1]
        @test components.far_field_velocity ≈ expected[2]
        @test components.far_field_curvature ≈ expected[3]
        @test residual.diagnostics.residual_norm ≈ norm(vector)
        @test residual.diagnostics.final_residual_norm ≈ norm(vector)
    end

    @testset "Residual metadata remains source-fidelity honest" begin
        residual = inner_bvp_residual((2.79, 0.28, 0.57); epsilon=0.1,
                                      ξ_max=6.0)

        @test residual.parameters.ε ≈ 0.1
        @test residual.parameters.far_field_targets ==
              (slope=0.1, velocity=0.0, curvature=0.0)
        @test residual.parameters.tip_conditions.slope == 0.0
        @test residual.domain.independent_variable == :ξ
        @test residual.domain.ξ₀ ≈ 2.79
        @test residual.domain.ξ_max ≈ 6.0
        @test residual.mesh.mesh_variable == :ξ
        @test residual.mesh.mesh_points == length(residual.mesh.ξ)
        @test residual.diagnostics.successful
        @test residual.diagnostics.source_status == "IMPL-inferred"
        @test residual.provenance.source_status == "IMPL-inferred"
        @test residual.provenance.canonical_source_doi == "10.1093/imamat/hxm043"
        @test residual.provenance.component_status.far_field_slope ==
              "IMPL-inferred"
        @test residual.provenance.boundary_status.tip_velocity ==
              "IMPL-inferred"
        @test residual.provenance.solver_settings.problem_kind ==
              :inner_bvp_residual
        @test occursin("endpoint shooting residual",
                       residual.provenance.solver_settings.algorithm)
    end

    @testset "Residual evaluation fails loudly on invalid parameters" begin
        unknowns = InnerBVPUnknowns(2.79, 0.28, 0.57)

        @test_throws ArgumentError inner_bvp_residual(unknowns; ε=0.0,
                                                      ξ_max=6.0)
        @test_throws ArgumentError inner_bvp_residual(unknowns; ε=NaN,
                                                      ξ_max=6.0)
        @test_throws ArgumentError inner_bvp_residual(unknowns; ε=0.1,
                                                      ξ_max=2.79)
        @test_throws ArgumentError inner_bvp_residual(unknowns; ε=0.1,
                                                      ξ_max=Inf)
        @test_throws ArgumentError inner_bvp_residual(unknowns; ε=0.1,
                                                      ξ_max=6.0,
                                                      maxiters=0)
    end

    @testset "Legacy solver reports the centralized residual" begin
        sol = solve_inner_bvp(ε=0.1, ξ_max=6.0, newton_iters=0)
        residual = inner_bvp_residual_vector((2.79, 0.28, 0.57);
                                             ε=0.1, ξ_max=6.0)
        @test sol.final_residual ≈ residual
        @test sol.final_residual_norm ≈ norm(residual)
    end
end
