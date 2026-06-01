# Public problem and result API layered over the current solver functions.
#
# These wrappers intentionally preserve the existing solver entry points and
# return types for legacy calls. Problem-taking methods return ProblemResult
# objects with stable metadata fields for downstream diagnostics and provenance.

export AbstractRecoilProblem,
       ConeSimilarityProblem, OuterMatchingProblem, CompositeProfileProblem,
       PDEVerificationProblem,
       ProblemResult,
       ConeSimilarityResult, OuterMatchingResult, CompositeProfileResult,
       PDEVerificationResult

abstract type AbstractRecoilProblem end

struct ProblemResult{P<:AbstractRecoilProblem,S}
    problem::P
    solution::S
    parameters::NamedTuple
    domain::NamedTuple
    mesh::NamedTuple
    diagnostics::NamedTuple
    provenance::NamedTuple
end

struct ConeSimilarityProblem <: AbstractRecoilProblem
    parameters::NamedTuple
    domain::NamedTuple
    solver::NamedTuple
    provenance::NamedTuple
end

struct OuterMatchingProblem <: AbstractRecoilProblem
    inner::Any
    parameters::NamedTuple
    domain::NamedTuple
    solver::NamedTuple
    provenance::NamedTuple
end

struct CompositeProfileProblem <: AbstractRecoilProblem
    inner::Any
    outer::Any
    parameters::NamedTuple
    domain::NamedTuple
    solver::NamedTuple
    provenance::NamedTuple
end

struct PDEVerificationProblem <: AbstractRecoilProblem
    parameters::NamedTuple
    domain::NamedTuple
    solver::NamedTuple
    provenance::NamedTuple
end

const ConeSimilarityResult = ProblemResult{ConeSimilarityProblem,InnerSolution}
const OuterMatchingResult = ProblemResult{OuterMatchingProblem,OuterSolution}
const CompositeProfileResult = ProblemResult{CompositeProfileProblem,CompositeSolution}
const PDEVerificationResult = ProblemResult{PDEVerificationProblem,PDESolution}

function _base_provenance(kind::Symbol; solver_settings::NamedTuple=(;),
                          benchmark_ids=(), artifacts=())
    metadata = default_recoil_provenance_metadata(
        kind; solver_settings=solver_settings, benchmark_ids=benchmark_ids,
        artifacts=artifacts)
    metadata_fields = _provenance_metadata_fields(metadata)
    (problem_kind=kind,
     source_status=_DEFAULT_SOURCE_STATUS,
     source_ledger=_SOURCE_LEDGER_PATH,
     canonical_source_doi=_CANONICAL_CONE_DOI,
     status_note="Current primitive-variable equations, boundary data, matching, and numerical constants are local reconstructed implementation data unless the source ledger marks a fact as C2001 or C2008-meta.",
     implementation_status="local reconstruction pending Decent-King 2008 article body",
     metadata=metadata,
     metadata_fields...)
end

_stored_user_provenance(provenance::NamedTuple) =
    haskey(provenance, :user_provenance) ? provenance.user_provenance : provenance

_provenance_metadata_fields(metadata::ProvenanceMetadata) =
    (source_citations=metadata.source_citations,
     source_ids=metadata.source_ids,
     assumptions=metadata.assumptions,
     benchmark_ids=metadata.benchmark_ids,
     solver_settings=metadata.solver_settings,
     artifacts=metadata.artifacts,
     package=metadata.package)

function _merge_provenance(kind::Symbol, provenance::NamedTuple;
                           solver_settings::NamedTuple=(;),
                           benchmark_ids=(), artifacts=())
    user_provenance = _stored_user_provenance(provenance)
    user_provenance isa NamedTuple ||
        throw(ArgumentError("provenance.user_provenance must be a NamedTuple; got $(typeof(user_provenance))"))
    merged = merge(_base_provenance(kind; solver_settings=solver_settings,
                                    benchmark_ids=benchmark_ids,
                                    artifacts=artifacts),
                   user_provenance,
                   (user_provenance=user_provenance,))
    if haskey(user_provenance, :metadata)
        metadata = user_provenance.metadata
        metadata isa ProvenanceMetadata ||
            throw(ArgumentError("provenance.metadata must be a ProvenanceMetadata; got $(typeof(metadata))"))
        merged = merge(merged, _provenance_metadata_fields(metadata))
    end
    merged
end

_merge_provenance(kind::Symbol, provenance; kwargs...) =
    throw(ArgumentError("provenance must be a NamedTuple; got $(typeof(provenance))"))

function _require_positive_finite(name::AbstractString, x)
    isfinite(x) && x > 0 ||
        throw(ArgumentError("$name must be positive and finite; got $x"))
    nothing
end

function _require_nonnegative_finite(name::AbstractString, x)
    isfinite(x) && x >= 0 ||
        throw(ArgumentError("$name must be nonnegative and finite; got $x"))
    nothing
end

function _require_finite(name::AbstractString, x)
    isfinite(x) ||
        throw(ArgumentError("$name must be finite; got $x"))
    nothing
end

function _require_positive_int(name::AbstractString, n::Int)
    n > 0 || throw(ArgumentError("$name must be positive; got $n"))
    nothing
end

_payload(x) = x isa ProblemResult ? x.solution : x

function _require_inner_solution_payload(inner; context::AbstractString)
    _payload(inner) isa InnerSolution ||
        throw(ArgumentError("$context requires an InnerSolution or ConeSimilarityResult"))
    nothing
end

function _require_outer_solution_payload(outer; context::AbstractString)
    _payload(outer) isa OuterSolution ||
        throw(ArgumentError("$context requires an OuterSolution or OuterMatchingResult"))
    nothing
end

function _query_grid_or_nothing(ξ_grid)
    ξ_grid === nothing && return nothing
    grid = Float64.(collect(ξ_grid))
    _validate_query_grid(grid; context="CompositeProfileProblem ξ_grid")
    grid
end

"""
    ConeSimilarityProblem(; ε=0.1, ξ₀=2.79, S₀=0.28, Sξξ₀=0.57,
                            ξ_max=60.0, newton_iters=30,
                            newton_tol=1e-4, ode_maxiters=2_000_000,
                            throw_on_failure=false, provenance=(;))

Public problem object for the current reconstructed cone-similarity inner BVP.
The legacy `solve_inner_bvp(; ...)` call still returns `InnerSolution`; calling
`solve_inner_bvp(problem)` returns `ConeSimilarityResult`.
"""
function ConeSimilarityProblem(; ε::Real=0.1, epsilon::Union{Nothing,Real}=nothing,
                                 ξ₀::Real=2.79, S₀::Real=0.28,
                                 Sξξ₀::Real=0.57, ξ_max::Real=60.0,
                                 newton_iters::Int=30,
                                 newton_tol::Real=1e-4,
                                 ode_maxiters::Int=2_000_000,
                                 throw_on_failure::Bool=false,
                                 provenance::NamedTuple=(;))
    εv = Float64(something(epsilon, ε))
    ξ₀v = Float64(ξ₀)
    S₀v = Float64(S₀)
    Sξξ₀v = Float64(Sξξ₀)
    ξ_maxv = Float64(ξ_max)
    tol = Float64(newton_tol)
    _require_positive_finite("ε", εv)
    _require_positive_finite("ξ₀", ξ₀v)
    _require_positive_finite("S₀", S₀v)
    _require_finite("Sξξ₀", Sξξ₀v)
    _require_finite("ξ_max", ξ_maxv)
    ξ_maxv > ξ₀v ||
        throw(ArgumentError("ξ_max must be greater than ξ₀; got ξ_max=$ξ_maxv and ξ₀=$ξ₀v"))
    _require_positive_int("newton_iters", newton_iters)
    _require_positive_int("ode_maxiters", ode_maxiters)
    _require_positive_finite("newton_tol", tol)

    parameters = (ε=εv, ξ₀=ξ₀v, S₀=S₀v, Sξξ₀=Sξξ₀v)
    domain = (ξ_max=ξ_maxv,)
    solver = (newton_iters=newton_iters, newton_tol=tol,
              ode_maxiters=ode_maxiters,
              throw_on_failure=throw_on_failure)
    ConeSimilarityProblem(parameters, domain, solver,
                          _merge_provenance(
                              :cone_similarity, provenance;
                              solver_settings=merge(parameters, domain, solver)))
end

"""
    OuterMatchingProblem(inner; ξ_match=15.0, ξ_max=100.0,
                         maxiters=500_000, provenance=(;))

Public problem object for the current local outer matching construction seeded
from an inner solution. `inner` may be an `InnerSolution` or `ConeSimilarityResult`.
"""
function OuterMatchingProblem(inner; ξ_match::Real=15.0, ξ_max::Real=100.0,
                              maxiters::Int=500_000,
                              provenance::NamedTuple=(;))
    ξ_matchv = Float64(ξ_match)
    ξ_maxv = Float64(ξ_max)
    _require_finite("ξ_match", ξ_matchv)
    _require_finite("ξ_max", ξ_maxv)
    ξ_maxv > ξ_matchv ||
        throw(ArgumentError("ξ_max must be greater than ξ_match; got ξ_max=$ξ_maxv and ξ_match=$ξ_matchv"))
    _require_positive_int("maxiters", maxiters)
    _require_inner_solution_payload(inner; context="OuterMatchingProblem")

    parameters = (ε_source=:inferred_from_inner,)
    domain = (ξ_match=ξ_matchv, ξ_max=ξ_maxv)
    solver = (maxiters=maxiters,)
    OuterMatchingProblem(inner, parameters, domain, solver,
                         _merge_provenance(
                             :outer_matching, provenance;
                             solver_settings=merge(parameters, domain, solver)))
end

"""
    CompositeProfileProblem(inner, outer; ξ_grid=nothing, n_points=500,
                            ξ_match=nothing, provenance=(;))

Public problem object for the current additive composite profile construction.
`inner`/`outer` may be legacy solution structs or corresponding ProblemResult
wrappers.
"""
function CompositeProfileProblem(inner, outer; ξ_grid=nothing, n_points::Int=500,
                                 ξ_match::Union{Nothing,Real}=nothing,
                                 provenance::NamedTuple=(;))
    n_points >= 2 ||
        throw(ArgumentError("n_points must be at least 2; got $n_points"))
    grid = _query_grid_or_nothing(ξ_grid)
    match = ξ_match === nothing ? nothing : Float64(ξ_match)
    match === nothing || _require_finite("ξ_match", match)
    _require_inner_solution_payload(inner; context="CompositeProfileProblem inner")
    _require_outer_solution_payload(outer; context="CompositeProfileProblem outer")
    parameters = (ε_source=:from_outer_solution,)
    domain = (ξ_match=match,)
    solver = (ξ_grid=grid, n_points=n_points)
    CompositeProfileProblem(inner, outer, parameters, domain, solver,
                            _merge_provenance(
                                :composite_profile, provenance;
                                solver_settings=merge(parameters, domain, solver)))
end

"""
    PDEVerificationProblem(; ε=0.1, N=200, z_min=0.01, z_max=10.0,
                           t_end=1.0, n_snapshots=10,
                           maxiters=1_000_000, provenance=(;))

Public problem object for the current method-of-lines PDE verification solve.
"""
function PDEVerificationProblem(; ε::Real=0.1, epsilon::Union{Nothing,Real}=nothing,
                                  N::Int=200, z_min::Real=0.01,
                                  z_max::Real=10.0, t_end::Real=1.0,
                                  n_snapshots::Int=10,
                                  maxiters::Int=1_000_000,
                                  provenance::NamedTuple=(;))
    εv = Float64(something(epsilon, ε))
    z_minv = Float64(z_min)
    z_maxv = Float64(z_max)
    t_endv = Float64(t_end)
    _require_positive_finite("ε", εv)
    N >= 3 || throw(ArgumentError("N must be at least 3; got $N"))
    _require_finite("z_min", z_minv)
    _require_finite("z_max", z_maxv)
    z_minv < z_maxv ||
        throw(ArgumentError("z_min must be less than z_max; got z_min=$z_minv and z_max=$z_maxv"))
    _require_nonnegative_finite("t_end", t_endv)
    _require_positive_int("n_snapshots", n_snapshots)
    _require_positive_int("maxiters", maxiters)

    parameters = (ε=εv,)
    domain = (z_min=z_minv, z_max=z_maxv, t_end=t_endv)
    solver = (N=N, n_snapshots=n_snapshots,
              maxiters=maxiters)
    PDEVerificationProblem(parameters, domain, solver,
                           _merge_provenance(
                               :pde_verification, provenance;
                               solver_settings=merge(parameters, domain, solver)))
end

function _inner_result_specific_diagnostics(sol::InnerSolution)
    (final_residual=sol.final_residual,
     final_residual_norm=sol.final_residual_norm,
     iterations=sol.iterations,
     converged=sol.converged,
     termination_reason=sol.termination_reason)
end

_result_diagnostics(sol; problem_kind::Symbol, existing::NamedTuple=(;)) =
    merge(existing, _normalized_diagnostic_fields(sol; problem_kind=problem_kind))

function _cone_result_parameters(problem::ConeSimilarityProblem,
                                 solution::InnerSolution)
    (ε=problem.parameters.ε,
     initial_ξ₀=problem.parameters.ξ₀,
     initial_S₀=problem.parameters.S₀,
     initial_Sξξ₀=problem.parameters.Sξξ₀,
     ξ₀=solution.ξ₀,
     S₀=solution.S₀,
     Sξξ₀=solution.Sξξ₀)
end

function ProblemResult(problem::ConeSimilarityProblem, solution::InnerSolution)
    domain = (ξ_min=first(solution.ξ), ξ_max=last(solution.ξ),
              ξ₀=solution.ξ₀)
    mesh = (ξ=solution.ξ,)
    ProblemResult{ConeSimilarityProblem,InnerSolution}(
        problem, solution, _cone_result_parameters(problem, solution), domain, mesh,
        _result_diagnostics(solution; problem_kind=:cone_similarity,
                            existing=_inner_result_specific_diagnostics(solution)),
        problem.provenance)
end

function ProblemResult(problem::OuterMatchingProblem, solution::OuterSolution)
    domain = (ξ_min=first(solution.ξ), ξ_max=last(solution.ξ),
              ξ_match=problem.domain.ξ_match)
    mesh = (ξ=solution.ξ,)
    diagnostics = _result_diagnostics(solution; problem_kind=:outer_matching,
                                      existing=solution.diagnostics)
    ProblemResult{OuterMatchingProblem,OuterSolution}(
        problem, solution, (ε=solution.ε,), domain, mesh, diagnostics,
        problem.provenance)
end

function ProblemResult(problem::CompositeProfileProblem, solution::CompositeSolution)
    domain = (ξ_min=first(solution.ξ), ξ_max=last(solution.ξ),
              ξ_match=solution.diagnostics.ξ_match)
    mesh = (ξ=solution.ξ,)
    diagnostics = _result_diagnostics(solution; problem_kind=:composite_profile,
                                      existing=solution.diagnostics)
    ProblemResult{CompositeProfileProblem,CompositeSolution}(
        problem, solution, (ε=solution.ε,), domain, mesh, diagnostics,
        problem.provenance)
end

function ProblemResult(problem::PDEVerificationProblem, solution::PDESolution)
    domain = (z_min=first(solution.z), z_max=last(solution.z),
              t_start=first(solution.t_snapshots),
              t_end=last(solution.t_snapshots))
    mesh = (z=solution.z, t=solution.t_snapshots)
    diagnostics = _result_diagnostics(solution; problem_kind=:pde_verification,
                                      existing=solution.diagnostics)
    ProblemResult{PDEVerificationProblem,PDESolution}(
        problem, solution, (ε=solution.ε,), domain, mesh, diagnostics,
        problem.provenance)
end

function _effective_problem(problem::ConeSimilarityProblem, overrides::NamedTuple)
    opts = merge(problem.parameters, problem.domain, problem.solver,
                 (provenance=problem.provenance,), overrides)
    ConeSimilarityProblem(; opts...)
end

function _effective_problem(problem::OuterMatchingProblem, overrides::NamedTuple)
    opts = merge(problem.domain, problem.solver, (provenance=problem.provenance,),
                 overrides)
    OuterMatchingProblem(problem.inner; opts...)
end

function _effective_problem(problem::CompositeProfileProblem,
                            overrides::NamedTuple)
    opts = merge(problem.domain, problem.solver, (provenance=problem.provenance,),
                 overrides)
    CompositeProfileProblem(problem.inner, problem.outer; opts...)
end

function _effective_problem(problem::PDEVerificationProblem,
                            overrides::NamedTuple)
    opts = merge(problem.parameters, problem.domain, problem.solver,
                 (provenance=problem.provenance,), overrides)
    PDEVerificationProblem(; opts...)
end

function solve_inner_bvp(problem::ConeSimilarityProblem; kwargs...)
    effective = _effective_problem(problem, (; kwargs...))
    call_kwargs = merge(effective.parameters, effective.domain, effective.solver)
    ProblemResult(effective, solve_inner_bvp(; call_kwargs...))
end

function solve_outer_matched(problem::OuterMatchingProblem; kwargs...)
    effective = _effective_problem(problem, (; kwargs...))
    inner = _payload(effective.inner)
    call_kwargs = merge(effective.domain, effective.solver)
    ProblemResult(effective, solve_outer_matched(inner; call_kwargs...))
end

function composite_solution(problem::CompositeProfileProblem; kwargs...)
    effective = _effective_problem(problem, (; kwargs...))
    inner = _payload(effective.inner)
    outer = _payload(effective.outer)
    call_kwargs = merge(effective.domain, effective.solver)
    ProblemResult(effective, composite_solution(inner, outer; call_kwargs...))
end

function solve_pde(problem::PDEVerificationProblem; kwargs...)
    effective = _effective_problem(problem, (; kwargs...))
    call_kwargs = merge(effective.parameters, effective.domain, effective.solver)
    ProblemResult(effective, solve_pde(; call_kwargs...))
end
