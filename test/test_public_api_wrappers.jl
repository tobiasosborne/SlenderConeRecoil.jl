using Test
using SlenderConeRecoil

@testset "Problem constructors validate public API inputs" begin
    cone = ConeSimilarityProblem()
    @test cone isa ConeSimilarityProblem
    @test cone.parameters.ε == 0.1
    @test cone.provenance.source_status == "IMPL-inferred"
    @test cone.provenance.source_ledger ==
          "docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md"

    @test ConeSimilarityProblem(; epsilon=0.2).parameters.ε == 0.2
    @test_throws ArgumentError ConeSimilarityProblem(; ε=-0.1)
    @test_throws ArgumentError ConeSimilarityProblem(; ξ₀=2.0, ξ_max=1.0)
    @test_throws ArgumentError ConeSimilarityProblem(; newton_iters=0)

    @test_throws ArgumentError OuterMatchingProblem(nothing)
    @test_throws ArgumentError OuterMatchingProblem(nothing; ξ_match=2.0, ξ_max=1.0)
    @test_throws ArgumentError CompositeProfileProblem(nothing, nothing)
    @test_throws ArgumentError CompositeProfileProblem(nothing, nothing; n_points=1)
    @test_throws ArgumentError CompositeProfileProblem(nothing, nothing; ξ_grid=[1.0, 1.0])
    @test_throws ArgumentError PDEVerificationProblem(; N=2)
    @test_throws ArgumentError PDEVerificationProblem(; t_end=-1.0)
end

@testset "ProblemResult exposes metadata for legacy solutions" begin
    inner = InnerSolution([1.0, 2.0, 3.0, 4.0, 5.0],
                          [0.1, 0.2, 0.3, 0.4, 0.5],
                          fill(0.1, 5), zeros(5), zeros(5),
                          1.0, 0.1, 0.0)
    cone_problem = ConeSimilarityProblem(; ξ₀=1.0, ξ_max=5.0)
    cone_result = ProblemResult(cone_problem, inner)
    @test cone_result isa ConeSimilarityResult
    @test cone_result.solution === inner
    @test cone_result.parameters.ε == 0.1
    @test cone_result.parameters.initial_ξ₀ == 1.0
    @test cone_result.parameters.S₀ == inner.S₀
    @test cone_result.domain.ξ₀ == 1.0
    @test cone_result.mesh.ξ === inner.ξ
    @test cone_result.diagnostics.mesh_points == 5
    @test cone_result.diagnostics.problem_kind == :cone_similarity
    @test !cone_result.diagnostics.successful
    @test cone_result.provenance.source_status == "IMPL-inferred"

    outer = OuterSolution([2.0, 3.0, 4.0], zeros(3), zeros(3), zeros(3),
                          zeros(3), 0.1)
    outer_problem = OuterMatchingProblem(cone_result; ξ_match=2.0, ξ_max=4.0)
    outer_result = ProblemResult(outer_problem, outer)
    @test outer_result isa OuterMatchingResult
    @test outer_result.parameters.ε == 0.1
    @test outer_result.domain.ξ_match == 2.0
    @test outer_result.mesh.ξ === outer.ξ
    @test outer_result.diagnostics.mesh_points == 3
    @test outer_result.diagnostics.problem_kind == :outer_matching

    composite_problem = CompositeProfileProblem(cone_result, outer_result;
                                                ξ_grid=[2.0, 3.0, 4.0],
                                                ξ_match=2.0)
    composite_result = composite_solution(composite_problem)
    @test composite_result isa CompositeProfileResult
    @test composite_result.solution isa CompositeSolution
    @test composite_result.mesh.ξ == [2.0, 3.0, 4.0]
    @test composite_result.parameters.ε == 0.1
    @test haskey(composite_result.diagnostics, :overlap_slope)
    @test composite_result.diagnostics.problem_kind == :composite_profile
    @test composite_result.diagnostics.successful
    @test composite_result.provenance.source_status == "IMPL-inferred"

    pde_problem = PDEVerificationProblem(; N=3, z_min=0.01, z_max=0.1,
                                         t_end=0.0, n_snapshots=1)
    pde_result = solve_pde(pde_problem)
    @test pde_result isa PDEVerificationResult
    @test pde_result.solution isa PDESolution
    @test pde_result.domain.t_end == 0.0
    @test pde_result.mesh.z === pde_result.solution.z
    @test pde_result.mesh.t === pde_result.solution.t_snapshots
    @test pde_result.diagnostics.mesh_points == 3
    @test pde_result.diagnostics.problem_kind == :pde_verification
    @test pde_result.diagnostics.successful
    @test pde_result.provenance.source_status == "IMPL-inferred"
end

@testset "Problem-taking methods are available without changing legacy calls" begin
    @test hasmethod(solve_inner_bvp, Tuple{ConeSimilarityProblem})
    @test hasmethod(solve_outer_matched, Tuple{OuterMatchingProblem})
    @test hasmethod(composite_solution, Tuple{CompositeProfileProblem})
    @test hasmethod(solve_pde, Tuple{PDEVerificationProblem})

    @test solve_inner_bvp === SlenderConeRecoil.solve_inner_bvp
    @test solve_outer_matched === SlenderConeRecoil.solve_outer_matched
    @test composite_solution === SlenderConeRecoil.composite_solution
    @test solve_pde === SlenderConeRecoil.solve_pde
end

@testset "Problem-taking methods report effective override metadata" begin
    pde_problem = PDEVerificationProblem(; ε=0.1, N=3, z_min=0.01, z_max=0.1,
                                         t_end=0.0, n_snapshots=1)
    pde_result = solve_pde(pde_problem; ε=0.2, z_max=0.2)
    @test pde_result.problem.parameters.ε == 0.2
    @test pde_result.parameters.ε == 0.2
    @test pde_result.problem.domain.z_max == 0.2
    @test pde_result.domain.z_max == 0.2
    @test pde_result.provenance.source_status == "IMPL-inferred"
end
