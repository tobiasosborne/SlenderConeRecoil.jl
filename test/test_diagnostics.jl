using Test
using SlenderConeRecoil

function _diagnostic_inner()
    InnerSolution([1.0, 2.0, 3.0, 4.0, 5.0],
                  [0.1, 0.2, 0.3, 0.4, 0.5],
                  fill(0.1, 5), zeros(5), zeros(5),
                  1.0, 0.1, 0.0,
                  [3e-8, 4e-8, 0.0], 5e-8, 2, true,
                  "residual tolerance reached")
end

function _diagnostic_outer()
    diagnostics = (context="mock outer",
                   retcode=:Success,
                   successful=true,
                   endpoint=4.0,
                   requested_endpoint=4.0,
                   saved_points=3)
    OuterSolution([2.0, 3.0, 4.0], [0.01, 0.0, -0.01],
                  zeros(3), zeros(3), [0.0, 0.02, 0.0],
                  0.1, diagnostics)
end

@testset "Diagnostic summaries normalize solver result fields" begin
    inner = _diagnostic_inner()
    summary = diagnostic_summary(inner)
    nt = as_namedtuple(summary)

    @test summary isa DiagnosticSummary
    @test nt.problem_kind == :cone_similarity
    @test nt.successful
    @test nt.residual_norm ≈ 5e-8
    @test nt.final_residual_norm ≈ 5e-8
    @test nt.mesh_points == 5
    @test nt.mesh_strictly_increasing
    @test nt.domain_start ≈ 1.0
    @test nt.domain_end ≈ 5.0
    @test nt.domain_span ≈ 4.0
    @test nt.converged
    @test diagnostics_succeeded(inner)
end

@testset "Retcode diagnostics are explicit and thresholdable" begin
    ok = (context="ok", retcode=:Success, endpoint=4.0,
          requested_endpoint=4.0, saved_points=2)
    failed = (context="failed", retcode=:MaxIters, endpoint=3.0,
              requested_endpoint=4.0, saved_points=2)
    endpoint_miss = (context="miss", retcode=:Success, endpoint=3.0,
                     requested_endpoint=4.0, saved_points=2)

    @test diagnostics_succeeded(ok)
    @test !diagnostics_succeeded(failed)
    @test !diagnostics_succeeded(endpoint_miss)

    outer = _diagnostic_outer()
    summary = diagnostic_summary(outer)
    @test summary.successful
    @test as_namedtuple(summary).retcode == :Success
    @test as_namedtuple(summary).maximum_abs_perturbation ≈ 0.01
    @test as_namedtuple(summary).maximum_abs_velocity ≈ 0.02
end

@testset "ProblemResult diagnostics expose comparable fields" begin
    inner = _diagnostic_inner()
    cone_result = ProblemResult(ConeSimilarityProblem(; ξ₀=1.0, ξ_max=4.0),
                                inner)

    outer = _diagnostic_outer()
    outer_result = ProblemResult(OuterMatchingProblem(cone_result;
                                                      ξ_match=2.0,
                                                      ξ_max=4.0),
                                 outer)

    composite_result = composite_solution(
        CompositeProfileProblem(cone_result, outer_result;
                                ξ_grid=[2.0, 3.0, 4.0],
                                ξ_match=2.0))

    pde_result = solve_pde(PDEVerificationProblem(; ε=0.1, N=4,
                                                  z_min=0.01, z_max=0.2,
                                                  t_end=0.0,
                                                  n_snapshots=1))

    for result in (cone_result, outer_result, composite_result, pde_result)
        diagnostics = result.diagnostics
        @test haskey(diagnostics, :problem_kind)
        @test haskey(diagnostics, :successful)
        @test haskey(diagnostics, :residual_norm)
        @test haskey(diagnostics, :final_residual_norm)
        @test haskey(diagnostics, :retcode)
        @test haskey(diagnostics, :mesh_points)
        @test haskey(diagnostics, :domain_start)
        @test haskey(diagnostics, :domain_end)
        @test haskey(diagnostics, :domain_span)
        @test diagnostics.mesh_points >= 3
        @test diagnostics.domain_span > 0
    end

    @test cone_result.diagnostics.problem_kind == :cone_similarity
    @test cone_result.diagnostics.final_residual_norm < 1e-6
    @test outer_result.diagnostics.problem_kind == :outer_matching
    @test outer_result.diagnostics.successful
    @test composite_result.diagnostics.problem_kind == :composite_profile
    @test composite_result.diagnostics.successful
    @test composite_result.diagnostics.minimum_S > 0
    @test pde_result.diagnostics.problem_kind == :pde_verification
    @test pde_result.diagnostics.successful
    @test pde_result.diagnostics.minimum_radius > 0
    @test pde_result.diagnostics.radius_positivity_margin > 0
    @test pde_result.diagnostics.pde_data_valid
    @test pde_result.diagnostics.grid_strictly_increasing
    @test pde_result.diagnostics.retcode_string == "Success"
    @test pde_result.diagnostics.retcode_successful
    @test pde_result.diagnostics.endpoint_reached
    @test pde_result.diagnostics.saved_time_points == 1
    @test pde_result.diagnostics.initial_area_mass ≈
          pde_result.diagnostics.final_area_mass
    @test pde_result.diagnostics.max_abs_area_mass_balance_residual == 0.0
    @test haskey(pde_result.solution.diagnostics, :similarity_collapse)
    @test pde_result.diagnostics.similarity_collapse_status ==
          :insufficient_snapshots
    @test !pde_result.diagnostics.similarity_collapse_successful
    @test pde_result.diagnostics.similarity_collapse_grid_points == 0
end

@testset "Mesh and domain helpers reject ambiguous NamedTuples" begin
    @test mesh_summary([1.0, 1.5, 2.5]; variable=:ξ).mesh_points == 3
    @test domain_summary([1.0, 1.5, 2.5]; variable=:ξ).domain_span ≈ 1.5
    @test_throws ArgumentError mesh_summary((q=[1.0, 2.0],))
    @test_throws ArgumentError domain_summary((q_min=1.0, q_max=2.0))
end

@testset "PDE diagnostics gate invalid figure inputs" begin
    z = [1.0, 1.5, 2.0]
    times = [0.0, 0.1]
    R = [[1.0, 1.0, 1.0], [1.1, 1.1, 1.1]]
    u = [[0.0, 0.1, 0.2], [0.0, 0.1, 0.2]]
    diagnostics = (context="manual PDE",
                   retcode=:Success,
                   successful=true,
                   endpoint=0.1,
                   requested_endpoint=0.1,
                   saved_points=2)

    valid = PDESolution(z, times, R, u, 0.1, diagnostics)
    invalid_R = deepcopy(R)
    invalid_R[2][2] = 0.0
    invalid = PDESolution(z, times, invalid_R, u, 0.1, diagnostics)

    function require_figure_diagnostics(x)
        diagnostics_succeeded(x) ||
            error("figure diagnostics failed; refusing to write figure metadata")
        true
    end

    @test require_figure_diagnostics(valid)
    @test_throws ErrorException require_figure_diagnostics(invalid)
end
