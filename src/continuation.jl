# Natural-parameter continuation helpers for the current reconstructed inner BVP.
#
# Ledger status: IMPL-inferred. These helpers orchestrate existing local
# shooting solves and do not upgrade the source fidelity of the underlying BVP.

export InnerContinuationResult, continue_inner_bvp_domain, continue_inner_bvp

struct InnerContinuationResult
    steps::Vector{NamedTuple}
    solutions::Vector{InnerSolution}
    parameters::NamedTuple
    diagnostics::NamedTuple
    provenance::NamedTuple
end

_unknowns_tuple(u::InnerBVPUnknowns) =
    (ξ₀=u.ξ₀, S₀=u.S₀, Sξξ₀=u.Sξξ₀)

_unknowns_from_solution(sol::InnerSolution) =
    InnerBVPUnknowns(sol.ξ₀, sol.S₀, sol.Sξξ₀)

function _continuation_values(primary, alias, primary_name::AbstractString,
                              alias_name::AbstractString)
    if primary !== nothing && alias !== nothing
        throw(ArgumentError("provide either $primary_name or $alias_name, not both"))
    end
    raw = primary === nothing ? alias : primary
    raw === nothing &&
        throw(ArgumentError("$primary_name is required"))
    values = raw isa Real ? [Float64(raw)] : Float64.(collect(raw))
    isempty(values) &&
        throw(ArgumentError("$primary_name must contain at least one value"))
    all(isfinite, values) ||
        throw(ArgumentError("$primary_name must contain only finite values"))
    values
end

function _normalise_continuation_axis(axis)
    symbol = axis isa Symbol ? axis : Symbol(axis)
    symbol in (:ξ_max, :xi_max) ||
        throw(ArgumentError("only axis=:ξ_max is currently supported; got $axis"))
    :ξ_max
end

function _continuation_step_diagnostics(index::Int, axis::Symbol,
                                        requested_value::Float64,
                                        sol::InnerSolution)
    finite_solution = diagnostics_succeeded(sol)
    successful = sol.converged && finite_solution
    merge(
        mesh_summary(sol),
        domain_summary(sol),
        (problem_kind=:inner_bvp_continuation_step,
         step=index,
         axis=axis,
         requested_value=requested_value,
         converged=sol.converged,
         diagnostic_successful=finite_solution,
         successful=successful,
         residual_norm=sol.final_residual_norm,
         final_residual_norm=sol.final_residual_norm,
         final_residual=sol.final_residual,
         iterations=sol.iterations,
         termination_reason=sol.termination_reason,
         source_status=_DEFAULT_SOURCE_STATUS))
end

function _continuation_step(index::Int, axis::Symbol,
                            requested_value::Float64,
                            initial::InnerBVPUnknowns,
                            sol::InnerSolution)
    diagnostics = _continuation_step_diagnostics(index, axis, requested_value, sol)
    final = _unknowns_from_solution(sol)
    (index=index,
     axis=axis,
     requested_value=requested_value,
     initial_guess=_unknowns_tuple(initial),
     final_unknowns=_unknowns_tuple(final),
     residual_norm=sol.final_residual_norm,
     converged=sol.converged,
     successful=diagnostics.successful,
     diagnostic_successful=diagnostics.diagnostic_successful,
     termination_reason=sol.termination_reason,
     diagnostics=diagnostics,
     solution=sol)
end

function _continuation_error_step(index::Int, axis::Symbol,
                                  requested_value::Float64,
                                  initial::InnerBVPUnknowns, err)
    message = sprint(showerror, err)
    diagnostics =
        (problem_kind=:inner_bvp_continuation_step,
         step=index,
         axis=axis,
         requested_value=requested_value,
         converged=false,
         diagnostic_successful=false,
         successful=false,
         residual_norm=NaN,
         final_residual_norm=NaN,
         final_residual=Float64[],
         iterations=0,
         termination_reason=message,
         error_type=string(typeof(err)),
         source_status=_DEFAULT_SOURCE_STATUS)
    (index=index,
     axis=axis,
     requested_value=requested_value,
     initial_guess=_unknowns_tuple(initial),
     final_unknowns=nothing,
     residual_norm=NaN,
     converged=false,
     successful=false,
     diagnostic_successful=false,
     termination_reason=message,
     diagnostics=diagnostics,
     solution=nothing)
end

function _continuation_diagnostics(axis::Symbol, values::Vector{Float64},
                                   steps::Vector{NamedTuple})
    failed_steps = [step.index for step in steps if !step.successful]
    completed = count(step -> step.successful, steps)
    final_residual = isempty(steps) ? NaN : steps[end].residual_norm
    successful = completed == length(values)
    reason = isempty(steps) ? "no continuation steps attempted" :
             successful ? "completed all continuation steps" :
             steps[end].termination_reason
    (problem_kind=:inner_bvp_continuation,
     axis=axis,
     parameter_sequence=copy(values),
     steps_requested=length(values),
     steps_attempted=length(steps),
     steps_completed=completed,
     failed_steps=failed_steps,
     successful=successful,
     residual_norm=final_residual,
     final_residual_norm=final_residual,
     termination_reason=reason,
     step_diagnostics=Tuple(step.diagnostics for step in steps),
     final_diagnostics=isempty(steps) ? (;) : steps[end].diagnostics,
     source_status=_DEFAULT_SOURCE_STATUS)
end

function _continuation_provenance(parameters::NamedTuple,
                                  diagnostics::NamedTuple,
                                  provenance::NamedTuple)
    solver_settings = merge(parameters,
                            (steps_attempted=diagnostics.steps_attempted,
                             steps_completed=diagnostics.steps_completed,
                             failed_steps=diagnostics.failed_steps))
    _merge_provenance(:inner_bvp_continuation, provenance;
                      solver_settings=solver_settings)
end

function _continuation_failure_message(step)
    "inner BVP continuation failed at step $(step.index) for " *
        "$(step.axis)=$(step.requested_value): $(step.termination_reason)"
end

"""
    continue_inner_bvp_domain(; ξ_max_values, ε=0.1, ξ₀=2.79,
                                S₀=0.28, Sξξ₀=0.57, ...)

Run a natural-parameter continuation ladder over `ξ_max`. Each successful step
reuses the final `(ξ₀, S₀, Sξξ₀)` from the previous inner shooting solve as the
initial guess for the next requested domain length.

If a step fails and `throw_on_failure=false`, the returned
`InnerContinuationResult` records the failed step and stops. With
`throw_on_failure=true`, the same failure is reported as an error.

Ledger status: IMPL-inferred local reconstruction pending Decent-King 2008
article-body verification.
"""
function continue_inner_bvp_domain(; ξ_max_values=nothing,
                                    xi_max_values=nothing,
                                    ε::Real=0.1,
                                    epsilon::Union{Nothing,Real}=nothing,
                                    ξ₀::Real=2.79,
                                    S₀::Real=0.28,
                                    Sξξ₀::Real=0.57,
                                    newton_iters::Int=30,
                                    newton_tol::Real=1e-4,
                                    ode_maxiters::Int=2_000_000,
                                    throw_on_failure::Bool=false,
                                    provenance::NamedTuple=(;))
    values = _continuation_values(ξ_max_values, xi_max_values,
                                  "ξ_max_values", "xi_max_values")
    εv = Float64(something(epsilon, ε))
    tol = Float64(newton_tol)
    _require_positive_finite("ε", εv)
    _require_positive_finite("newton_tol", tol)
    _require_positive_int("ode_maxiters", ode_maxiters)
    newton_iters >= 0 ||
        throw(ArgumentError("newton_iters must be nonnegative; got $newton_iters"))

    initial = InnerBVPUnknowns(ξ₀, S₀, Sξξ₀)
    current = initial
    steps = Vector{NamedTuple}()
    solutions = InnerSolution[]

    for (index, ξ_max) in pairs(values)
        step = try
            sol = solve_inner_bvp(; ε=εv, ξ₀=current.ξ₀, S₀=current.S₀,
                                  Sξξ₀=current.Sξξ₀, ξ_max=ξ_max,
                                  newton_iters=newton_iters,
                                  newton_tol=tol,
                                  ode_maxiters=ode_maxiters,
                                  throw_on_failure=false)
            push!(solutions, sol)
            _continuation_step(index, :ξ_max, ξ_max, current, sol)
        catch err
            _continuation_error_step(index, :ξ_max, ξ_max, current, err)
        end

        push!(steps, step)
        if !step.successful
            throw_on_failure && error(_continuation_failure_message(step))
            break
        end
        current = _unknowns_from_solution(step.solution)
    end

    parameters =
        (axis=:ξ_max,
         ε=εv,
         initial_unknowns=_unknowns_tuple(initial),
         ξ_max_values=copy(values),
         newton_iters=newton_iters,
         newton_tol=tol,
         ode_maxiters=ode_maxiters,
         throw_on_failure=throw_on_failure)
    diagnostics = _continuation_diagnostics(:ξ_max, values, steps)
    InnerContinuationResult(steps, solutions, parameters, diagnostics,
                            _continuation_provenance(parameters, diagnostics,
                                                     provenance))
end

function continue_inner_bvp_domain(problem::ConeSimilarityProblem; kwargs...)
    opts = merge(problem.parameters, problem.solver,
                 (provenance=problem.provenance,), (; kwargs...))
    continue_inner_bvp_domain(; opts...)
end

"""
    continue_inner_bvp(; axis=:ξ_max, values, kwargs...)

Generic entry point for inner-BVP continuation. Currently `axis=:ξ_max` is the
supported natural-parameter path and dispatches to `continue_inner_bvp_domain`.
"""
function continue_inner_bvp(; axis=:ξ_max, values=nothing, kwargs...)
    normalized_axis = _normalise_continuation_axis(axis)
    normalized_axis == :ξ_max || error("unreachable continuation axis")
    opts = (; kwargs...)
    if values !== nothing
        if haskey(opts, :ξ_max_values) || haskey(opts, :xi_max_values)
            throw(ArgumentError("provide either values or ξ_max_values/xi_max_values, not both"))
        end
        opts = merge((ξ_max_values=values,), opts)
    end
    continue_inner_bvp_domain(; opts...)
end

function continue_inner_bvp(problem::ConeSimilarityProblem; axis=:ξ_max,
                            values=nothing, kwargs...)
    normalized_axis = _normalise_continuation_axis(axis)
    normalized_axis == :ξ_max || error("unreachable continuation axis")
    opts = (; kwargs...)
    if values !== nothing
        if haskey(opts, :ξ_max_values) || haskey(opts, :xi_max_values)
            throw(ArgumentError("provide either values or ξ_max_values/xi_max_values, not both"))
        end
        opts = merge((ξ_max_values=values,), opts)
    end
    continue_inner_bvp_domain(problem; opts...)
end
