using Test
using SlenderConeRecoil

@testset "Inner midpoint collocation solver" begin
    @testset "Cheap collocation solve returns common result metadata" begin
        result = solve_inner_bvp_collocation(; ε=0.1, ξ_max=6.0,
                                             nodes=8, maxiters=5,
                                             residual_tol=1e-6,
                                             seed_newton_iters=0)

        @test result isa ConeSimilarityResult
        @test result.solution isa InnerSolution
        @test result.provenance.problem_kind == :inner_bvp_collocation
        @test result.provenance.source_status == "IMPL-inferred"
        @test result.provenance.solver_settings.problem_kind ==
              :inner_bvp_collocation
        @test result.diagnostics.problem_kind == :inner_bvp_collocation
        @test result.diagnostics.solver_backend == :midpoint_collocation
        @test result.diagnostics.seed_source == :shooting
        @test result.diagnostics.collocation_converged
        @test result.diagnostics.collocation_residual_norm < 1e-6
        @test result.diagnostics.collocation_residual_length ==
              4 * (8 - 1) + 7
        @test result.diagnostics.collocation_unknowns == 4 * 8 + 3
        @test result.solution.final_residual_norm < 1e-6
        @test result.solution.ξ[1] ≈ result.solution.ξ₀
        @test result.solution.ξ[end] ≈ 6.0
        @test result.solution.S[end] / result.solution.ξ[end] ≈ 0.1 atol=1e-6
        @test abs(result.solution.U[end]) < 1e-6
    end

    @testset "Collocation solve validates inputs and fails loudly" begin
        @test_throws ArgumentError solve_inner_bvp_collocation(; nodes=2)
        @test_throws ArgumentError solve_inner_bvp_collocation(; nodes=8,
                                                               max_nodes=7)
        @test_throws ArgumentError solve_inner_bvp_collocation(; maxiters=0)
        @test_throws ArgumentError solve_inner_bvp_collocation(; residual_tol=0.0)
        @test_throws ArgumentError solve_inner_bvp_collocation(; seed_newton_iters=-1)

        @test_throws ErrorException solve_inner_bvp_collocation(; ε=0.1,
                                                                ξ_max=6.0,
                                                                nodes=8,
                                                                maxiters=1,
                                                                residual_tol=1e-14,
                                                                seed_newton_iters=0)
    end

    @testset "Problem-taking collocation wrapper keeps override metadata" begin
        problem = ConeSimilarityProblem(; ε=0.1, ξ₀=2.79, S₀=0.28,
                                        Sξξ₀=0.57, ξ_max=6.0)
        result = solve_inner_bvp_collocation(problem; nodes=8, maxiters=5,
                                             residual_tol=1e-6,
                                             seed_newton_iters=0)
        @test result.problem.domain.ξ_max == 6.0
        @test result.diagnostics.collocation_nodes == 8
        @test result.diagnostics.collocation_converged
    end
end
