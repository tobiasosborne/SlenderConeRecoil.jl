# Lightweight diagnostic summaries shared by legacy solver results and the
# public ProblemResult API.

export DiagnosticSummary, residual_norm, mesh_summary, domain_summary,
       diagnostic_summary, diagnostics_succeeded

struct DiagnosticSummary
    problem_kind::Symbol
    successful::Bool
    residual_norm::Float64
    final_residual_norm::Float64
    retcode::Any
    mesh::NamedTuple
    domain::NamedTuple
    checks::NamedTuple
end

function DiagnosticSummary(problem_kind::Symbol; successful::Bool=false,
                           residual_norm=NaN,
                           final_residual_norm=residual_norm,
                           retcode=nothing, mesh::NamedTuple=(;),
                           domain::NamedTuple=(;), checks::NamedTuple=(;))
    DiagnosticSummary(problem_kind, successful, _float_or_nan(residual_norm),
                      _float_or_nan(final_residual_norm), retcode,
                      mesh, domain, checks)
end

_float_or_nan(x) = x === nothing ? NaN : Float64(x)

function _finite_vector(xs)
    all(x -> x isa Real && isfinite(x), xs)
end

function _max_abs_or_nan(xs)
    isempty(xs) && return NaN
    maximum(abs, xs)
end

function _min_or_nan(xs)
    isempty(xs) && return NaN
    minimum(xs)
end

function _max_or_nan(xs)
    isempty(xs) && return NaN
    maximum(xs)
end

function _vector_norm_or_nan(xs)
    isempty(xs) && return NaN
    values = Float64.(collect(xs))
    _finite_vector(values) || return NaN
    sqrt(sum(abs2, values))
end

function _retcode_successful(retcode)
    retcode === nothing && return false
    string(retcode) in ("Success", "Terminated")
end

function _endpoint_successful(diagnostics::NamedTuple)
    (haskey(diagnostics, :endpoint) &&
     haskey(diagnostics, :requested_endpoint)) || return true
    endpoint = Float64(diagnostics.endpoint)
    requested = Float64(diagnostics.requested_endpoint)
    isfinite(endpoint) || return false
    isfinite(requested) || return true
    abs(endpoint - requested) <= _endpoint_tolerance(requested)
end

function residual_norm(diagnostics::NamedTuple)
    haskey(diagnostics, :final_residual_norm) &&
        return _float_or_nan(diagnostics.final_residual_norm)
    haskey(diagnostics, :residual_norm) &&
        return _float_or_nan(diagnostics.residual_norm)
    haskey(diagnostics, :final_residual) &&
        return _vector_norm_or_nan(diagnostics.final_residual)
    NaN
end

residual_norm(summary::DiagnosticSummary) = summary.residual_norm
residual_norm(sol::InnerSolution) = _float_or_nan(sol.final_residual_norm)

function residual_norm(x)
    hasproperty(x, :diagnostics) && return residual_norm(getproperty(x, :diagnostics))
    NaN
end

function mesh_summary(xs; variable::Symbol=:x)
    n = length(xs)
    spacing_min = NaN
    spacing_max = NaN
    spacing_ratio = NaN
    strictly_increasing = n >= 2
    if n >= 2
        spacings = [Float64(xs[i] - xs[i-1]) for i in 2:n]
        finite_spacings = _finite_vector(spacings)
        strictly_increasing = finite_spacings && all(>(0), spacings)
        if finite_spacings
            spacing_min = minimum(spacings)
            spacing_max = maximum(spacings)
            spacing_ratio = spacing_min > 0 ? spacing_max / spacing_min : Inf
        end
    end
    (mesh_variable=variable,
     mesh_points=n,
     mesh_spacing_min=spacing_min,
     mesh_spacing_max=spacing_max,
     mesh_spacing_ratio=spacing_ratio,
     mesh_strictly_increasing=strictly_increasing)
end

function mesh_summary(mesh::NamedTuple)
    haskey(mesh, :ξ) && return mesh_summary(mesh.ξ; variable=:ξ)
    haskey(mesh, :z) && return mesh_summary(mesh.z; variable=:z)
    haskey(mesh, :x) && return mesh_summary(mesh.x; variable=:x)
    throw(ArgumentError("mesh_summary requires a NamedTuple with ξ, z, or x"))
end

mesh_summary(sol::InnerSolution) = mesh_summary(sol.ξ; variable=:ξ)
mesh_summary(sol::OuterSolution) = mesh_summary(sol.ξ; variable=:ξ)
mesh_summary(sol::CompositeSolution) = mesh_summary(sol.ξ; variable=:ξ)
mesh_summary(sol::HierarchySolution) = mesh_summary(sol.ξ; variable=:ξ)
mesh_summary(sol::PDESolution) = mesh_summary(sol.z; variable=:z)

function _time_domain_summary(times)
    n = length(times)
    start = n == 0 ? NaN : Float64(first(times))
    finish = n == 0 ? NaN : Float64(last(times))
    (time_points=n,
     time_start=start,
     time_end=finish,
     time_span=finish - start)
end

_empty_time_domain_summary() =
    (time_points=0, time_start=NaN, time_end=NaN, time_span=NaN)

function domain_summary(xs; variable::Symbol=:x, time=nothing)
    n = length(xs)
    start = n == 0 ? NaN : Float64(first(xs))
    finish = n == 0 ? NaN : Float64(last(xs))
    time_fields = time === nothing ? _empty_time_domain_summary() :
                  _time_domain_summary(time)
    merge((domain_variable=variable,
           domain_start=start,
           domain_end=finish,
           domain_span=finish - start),
          time_fields)
end

function domain_summary(domain::NamedTuple)
    if haskey(domain, :ξ_min) && haskey(domain, :ξ_max)
        return merge((domain_variable=:ξ,
                      domain_start=Float64(domain.ξ_min),
                      domain_end=Float64(domain.ξ_max),
                      domain_span=Float64(domain.ξ_max - domain.ξ_min)),
                     _empty_time_domain_summary())
    elseif haskey(domain, :z_min) && haskey(domain, :z_max)
        time_fields =
            haskey(domain, :t_start) && haskey(domain, :t_end) ?
            (time_points=0,
             time_start=Float64(domain.t_start),
             time_end=Float64(domain.t_end),
             time_span=Float64(domain.t_end - domain.t_start)) :
            _empty_time_domain_summary()
        return merge((domain_variable=:z,
                      domain_start=Float64(domain.z_min),
                      domain_end=Float64(domain.z_max),
                      domain_span=Float64(domain.z_max - domain.z_min)),
                     time_fields)
    end
    throw(ArgumentError("domain_summary requires ξ_min/ξ_max or z_min/z_max"))
end

domain_summary(sol::InnerSolution) = domain_summary(sol.ξ; variable=:ξ)
domain_summary(sol::OuterSolution) = domain_summary(sol.ξ; variable=:ξ)
domain_summary(sol::CompositeSolution) = domain_summary(sol.ξ; variable=:ξ)
domain_summary(sol::HierarchySolution) = domain_summary(sol.ξ; variable=:ξ)
domain_summary(sol::PDESolution) =
    domain_summary(sol.z; variable=:z, time=sol.t_snapshots)

function diagnostics_succeeded(diagnostics::NamedTuple)
    if haskey(diagnostics, :successful)
        Bool(diagnostics.successful) || return false
        if haskey(diagnostics, :retcode) && diagnostics.retcode !== nothing
            return _retcode_successful(diagnostics.retcode) &&
                   _endpoint_successful(diagnostics)
        end
        return true
    end
    if haskey(diagnostics, :retcode) && diagnostics.retcode !== nothing
        return _retcode_successful(diagnostics.retcode) &&
               _endpoint_successful(diagnostics)
    end
    haskey(diagnostics, :converged) && return Bool(diagnostics.converged)
    false
end

diagnostics_succeeded(summary::DiagnosticSummary) = summary.successful

function diagnostics_succeeded(sol::InnerSolution)
    !isempty(sol.ξ) && _finite_vector(sol.ξ) && _finite_vector(sol.S) &&
        _finite_vector(sol.Sξ) && _finite_vector(sol.Sξξ) &&
        _finite_vector(sol.U) && isfinite(sol.final_residual_norm)
end

function diagnostics_succeeded(sol::OuterSolution)
    diagnostics_succeeded(sol.diagnostics) && _finite_vector(sol.ξ) &&
        _finite_vector(sol.s₁) && _finite_vector(sol.s₁ξ) &&
        _finite_vector(sol.s₁ξξ) && _finite_vector(sol.u₁)
end

function diagnostics_succeeded(sol::CompositeSolution)
    haskey(sol.diagnostics, :fit_points) && sol.diagnostics.fit_points > 0 &&
        haskey(sol.diagnostics, :overlap_slope) &&
        isfinite(sol.diagnostics.overlap_slope) &&
        haskey(sol.diagnostics, :overlap_intercept) &&
        isfinite(sol.diagnostics.overlap_intercept) &&
        _finite_vector(sol.ξ) && _finite_vector(sol.S) && _finite_vector(sol.U)
end

function diagnostics_succeeded(sol::HierarchySolution)
    diagnostics_succeeded(sol.diagnostics) && _finite_vector(sol.ξ) &&
        _finite_vector(sol.S) && _finite_vector(sol.Sξ) &&
        _finite_vector(sol.Sξξ) && _finite_vector(sol.U)
end

function diagnostics_succeeded(sol::PDESolution)
    diagnostics_succeeded(sol.diagnostics) || return false
    _pde_solution_data_valid(sol) || return false

    diagnostics = sol.diagnostics
    for field in (:pde_data_valid, :finite_state, :positive_radius,
                  :grid_strictly_increasing, :grid_finite, :time_finite,
                  :time_strictly_increasing)
        if haskey(diagnostics, field)
            value = getfield(diagnostics, field)
            value isa Bool && value || return false
        end
    end
    if haskey(diagnostics, :retcode_successful)
        value = diagnostics.retcode_successful
        if value !== missing && !(value isa Bool && value)
            return false
        end
    end
    if haskey(diagnostics, :endpoint_reached)
        value = diagnostics.endpoint_reached
        if value !== missing && !(value isa Bool && value)
            return false
        end
    end
    if haskey(diagnostics, :minimum_radius) &&
       !(isfinite(diagnostics.minimum_radius) && diagnostics.minimum_radius > 0)
        return false
    end
    if haskey(diagnostics, :radius_positivity_margin) &&
       !(isfinite(diagnostics.radius_positivity_margin) &&
         diagnostics.radius_positivity_margin > 0)
        return false
    end
    if haskey(diagnostics, :saved_points_match_reported) &&
       diagnostics.saved_points_match_reported !== missing &&
       !(diagnostics.saved_points_match_reported isa Bool &&
         diagnostics.saved_points_match_reported)
        return false
    end
    true
end

function _inner_diagnostic_checks(sol::InnerSolution)
    (converged=sol.converged,
     iterations=sol.iterations,
     termination_reason=sol.termination_reason,
     residual_length=length(sol.final_residual))
end

function _outer_diagnostic_checks(sol::OuterSolution)
    (maximum_abs_perturbation=_max_abs_or_nan(sol.s₁),
     maximum_abs_velocity=_max_abs_or_nan(sol.u₁))
end

function _composite_diagnostic_checks(sol::CompositeSolution)
    (minimum_S=_min_or_nan(sol.S),
     maximum_S=_max_or_nan(sol.S),
     maximum_abs_U=_max_abs_or_nan(sol.U))
end

function _hierarchy_diagnostic_checks(sol::HierarchySolution)
    (minimum_S=_min_or_nan(sol.S),
     maximum_S=_max_or_nan(sol.S),
     maximum_abs_U=_max_abs_or_nan(sol.U))
end

function _pde_diagnostic_checks(sol::PDESolution)
    radii = isempty(sol.R) ? Float64[] : reduce(vcat, sol.R)
    velocities = isempty(sol.u) ? Float64[] : reduce(vcat, sol.u)
    diagnostics = sol.diagnostics
    (minimum_radius=get(diagnostics, :minimum_radius, _min_or_nan(radii)),
     radius_positivity_margin=
         get(diagnostics, :radius_positivity_margin, _min_or_nan(radii)),
     maximum_radius=get(diagnostics, :maximum_radius, _max_or_nan(radii)),
     maximum_abs_velocity=
         get(diagnostics, :maximum_abs_velocity, _max_abs_or_nan(velocities)),
     grid_spacing_min=get(diagnostics, :grid_spacing_min, NaN),
     grid_spacing_max=get(diagnostics, :grid_spacing_max, NaN),
     grid_spacing_ratio=get(diagnostics, :grid_spacing_ratio, NaN),
     grid_strictly_increasing=
         get(diagnostics, :grid_strictly_increasing, false),
     saved_time_points=get(diagnostics, :saved_time_points,
                           length(sol.t_snapshots)),
     state_snapshots=get(diagnostics, :state_snapshots, length(sol.R)),
     retcode_string=get(diagnostics, :retcode_string,
                        string(get(diagnostics, :retcode, nothing))),
     retcode_successful=get(diagnostics, :retcode_successful, false),
     endpoint_reached=get(diagnostics, :endpoint_reached, false),
     initial_area_mass=get(diagnostics, :initial_area_mass, NaN),
     final_area_mass=get(diagnostics, :final_area_mass, NaN),
     final_area_mass_drift=get(diagnostics, :final_area_mass_drift, NaN),
     final_relative_area_mass_drift=
         get(diagnostics, :final_relative_area_mass_drift, NaN),
     max_abs_area_mass_balance_residual=
         get(diagnostics, :max_abs_area_mass_balance_residual, NaN),
     max_abs_left_boundary_area_flux=
         get(diagnostics, :max_abs_left_boundary_area_flux, NaN),
     max_abs_right_boundary_area_flux=
         get(diagnostics, :max_abs_right_boundary_area_flux, NaN),
     similarity_collapse_status=
         get(diagnostics, :similarity_collapse_status, :not_computed),
     similarity_collapse_successful=
         get(diagnostics, :similarity_collapse_successful, false),
     similarity_collapse_score=
         get(diagnostics, :similarity_collapse_score, NaN),
     similarity_collapse_xi_min=
         get(diagnostics, :similarity_collapse_xi_min, NaN),
     similarity_collapse_xi_max=
         get(diagnostics, :similarity_collapse_xi_max, NaN),
     similarity_collapse_grid_points=
         get(diagnostics, :similarity_collapse_grid_points, 0),
     similarity_collapse_included_snapshots=
         get(diagnostics, :similarity_collapse_included_snapshots, 0),
     similarity_collapse_excluded_snapshots=
         get(diagnostics, :similarity_collapse_excluded_snapshots, 0),
     similarity_collapse_profile_relative_rms=
         get(diagnostics, :similarity_collapse_profile_relative_rms, NaN),
     similarity_collapse_slope_relative_rms=
         get(diagnostics, :similarity_collapse_slope_relative_rms, NaN),
     similarity_collapse_curvature_relative_rms=
         get(diagnostics, :similarity_collapse_curvature_relative_rms, NaN),
     similarity_collapse_velocity_relative_rms=
         get(diagnostics, :similarity_collapse_velocity_relative_rms, NaN),
     similarity_collapse_wave_phase_status=
         get(diagnostics, :similarity_collapse_wave_phase_status, :not_computed),
     similarity_collapse_wave_phase_relative_rms=
         get(diagnostics, :similarity_collapse_wave_phase_relative_rms, NaN))
end

function diagnostic_summary(sol::InnerSolution; problem_kind::Symbol=:cone_similarity)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_inner_diagnostic_checks(sol))
end

function diagnostic_summary(sol::OuterSolution; problem_kind::Symbol=:outer_matching)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      retcode=get(sol.diagnostics, :retcode, nothing),
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_outer_diagnostic_checks(sol))
end

function diagnostic_summary(sol::CompositeSolution; problem_kind::Symbol=:composite_profile)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_composite_diagnostic_checks(sol))
end

function diagnostic_summary(sol::HierarchySolution; problem_kind::Symbol=:outer_hierarchy)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      retcode=get(sol.diagnostics, :retcode, nothing),
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_hierarchy_diagnostic_checks(sol))
end

function diagnostic_summary(sol::PDESolution; problem_kind::Symbol=:pde_verification)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      retcode=get(sol.diagnostics, :retcode, nothing),
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_pde_diagnostic_checks(sol))
end

function diagnostic_summary(diagnostics::NamedTuple; problem_kind::Symbol=:unknown)
    rn = residual_norm(diagnostics)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(diagnostics),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      retcode=get(diagnostics, :retcode, nothing))
end

function _problem_kind_from_result(x, fallback::Symbol)
    if hasproperty(x, :problem)
        problem = getproperty(x, :problem)
        if hasproperty(problem, :provenance) &&
           haskey(getproperty(problem, :provenance), :problem_kind)
            return getproperty(problem, :provenance).problem_kind
        end
    end
    fallback
end

function diagnostic_summary(x; problem_kind::Symbol=:unknown)
    if hasproperty(x, :solution)
        kind = _problem_kind_from_result(x, problem_kind)
        return diagnostic_summary(getproperty(x, :solution); problem_kind=kind)
    elseif hasproperty(x, :diagnostics)
        return diagnostic_summary(getproperty(x, :diagnostics);
                                  problem_kind=problem_kind)
    end
    throw(ArgumentError("diagnostic_summary does not know how to summarize $(typeof(x))"))
end

function as_namedtuple(summary::DiagnosticSummary)
    merge((problem_kind=summary.problem_kind,
           successful=summary.successful,
           residual_norm=summary.residual_norm,
           final_residual_norm=summary.final_residual_norm,
           retcode=summary.retcode),
          summary.mesh,
          summary.domain,
          summary.checks)
end

_normalized_diagnostic_fields(x; problem_kind::Symbol) =
    as_namedtuple(diagnostic_summary(x; problem_kind=problem_kind))
