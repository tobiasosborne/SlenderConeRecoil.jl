using Test
using SlenderConeRecoil

@testset "Inner BVP continuation" begin
    @testset "Domain continuation reuses successful step parameters" begin
        result = continue_inner_bvp_domain(; ξ_max_values=[6.0, 6.1],
                                           ε=0.1,
                                           newton_iters=5,
                                           newton_tol=1e-3)

        @test result isa InnerContinuationResult
        @test length(result.steps) == 2
        @test length(result.solutions) == 2
        @test result.diagnostics.problem_kind == :inner_bvp_continuation
        @test result.diagnostics.axis == :ξ_max
        @test result.diagnostics.successful
        @test isempty(result.diagnostics.failed_steps)
        @test result.parameters.ξ_max_values == [6.0, 6.1]
        @test result.provenance.problem_kind == :inner_bvp_continuation
        @test result.provenance.source_status == "IMPL-inferred"
        @test result.provenance.solver_settings.problem_kind ==
              :inner_bvp_continuation

        first_step, second_step = result.steps
        @test first_step.successful
        @test second_step.successful
        @test first_step.requested_value == 6.0
        @test second_step.requested_value == 6.1
        @test first_step.final_unknowns == second_step.initial_guess
        @test first_step.converged
        @test second_step.converged
        @test first_step.diagnostics.source_status == "IMPL-inferred"
        @test second_step.solution === result.solutions[2]
        @test result.diagnostics.final_residual_norm ==
              second_step.residual_norm
    end

    @testset "Generic continuation entry point supports xi_max alias" begin
        result = continue_inner_bvp(; axis=:xi_max, values=[6.0],
                                    ε=0.1, newton_iters=5,
                                    newton_tol=1e-3)

        @test result isa InnerContinuationResult
        @test length(result.steps) == 1
        @test result.steps[1].requested_value == 6.0
        @test result.diagnostics.successful
    end

    @testset "Controlled failure records failed step metadata" begin
        result = continue_inner_bvp_domain(; ξ_max_values=[2.0, 6.0],
                                           throw_on_failure=false)

        @test result isa InnerContinuationResult
        @test length(result.steps) == 1
        @test isempty(result.solutions)
        @test !result.diagnostics.successful
        @test result.diagnostics.failed_steps == [1]
        @test result.diagnostics.steps_requested == 2
        @test result.diagnostics.steps_attempted == 1
        @test result.steps[1].final_unknowns === nothing
        @test !result.steps[1].successful
        @test occursin("ξ_max must be greater than ξ₀",
                       result.steps[1].termination_reason)

        @test_throws ErrorException continue_inner_bvp_domain(;
            ξ_max_values=[2.0],
            throw_on_failure=true,
        )
    end

    @testset "Continuation validates inputs" begin
        @test_throws ArgumentError continue_inner_bvp_domain()
        @test_throws ArgumentError continue_inner_bvp_domain(;
            ξ_max_values=[6.0],
            xi_max_values=[6.0],
        )
        @test_throws ArgumentError continue_inner_bvp_domain(;
            ξ_max_values=Float64[],
        )
        @test_throws ArgumentError continue_inner_bvp_domain(;
            ξ_max_values=[NaN],
        )
        @test_throws ArgumentError continue_inner_bvp(;
            axis=:ε,
            values=[0.1],
        )
        @test_throws ArgumentError continue_inner_bvp(;
            values=[6.0],
            ξ_max_values=[6.0],
        )
        @test_throws ArgumentError continue_inner_bvp_domain(;
            ξ_max_values=[6.0],
            newton_iters=-1,
        )
    end
end
