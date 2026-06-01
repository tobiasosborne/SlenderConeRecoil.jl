# Shared checks for SciML solve objects.

using DifferentialEquations

const _SOLVER_ENDPOINT_RTOL = 1e-7
const _SOLVER_ENDPOINT_ATOL = 1e-8

function _endpoint_tolerance(expected;
                             rtol::Float64=_SOLVER_ENDPOINT_RTOL,
                             atol::Float64=_SOLVER_ENDPOINT_ATOL)
    max(atol, rtol * max(1.0, abs(Float64(expected))))
end

function _solution_diagnostics(context::AbstractString, endpoint, requested_endpoint;
                               retcode=nothing, saved_points::Int=0)
    successful = retcode !== nothing &&
                 string(retcode) in ("Success", "Terminated") &&
                 isfinite(Float64(endpoint))
    if successful && isfinite(Float64(requested_endpoint))
        tol = _endpoint_tolerance(requested_endpoint)
        successful = abs(Float64(endpoint) - Float64(requested_endpoint)) <= tol
    end
    (context=String(context),
     retcode=retcode,
     successful=successful,
     endpoint=Float64(endpoint),
     requested_endpoint=Float64(requested_endpoint),
     saved_points=saved_points)
end

function _manual_solution_diagnostics(context::AbstractString, xs)
    endpoint = isempty(xs) ? NaN : xs[end]
    _solution_diagnostics(context, endpoint, NaN)
end

function _require_successful_solution(sol, expected_endpoint;
                                      context::AbstractString,
                                      expected_times=nothing,
                                      rtol::Float64=_SOLVER_ENDPOINT_RTOL,
                                      atol::Float64=_SOLVER_ENDPOINT_ATOL)
    if !SciMLBase.successful_retcode(sol)
        error("$context failed: SciML solver returned retcode=$(sol.retcode)")
    end

    if isempty(sol.t) || isempty(sol.u)
        error("$context failed: SciML solver returned no saved states")
    end

    endpoint = Float64(sol.t[end])
    endpoint_tol = _endpoint_tolerance(expected_endpoint; rtol=rtol, atol=atol)
    if !isfinite(endpoint) || abs(endpoint - Float64(expected_endpoint)) > endpoint_tol
        error("$context failed: reached endpoint $endpoint, expected $(Float64(expected_endpoint)) within $endpoint_tol")
    end

    if expected_times !== nothing
        if length(sol.t) != length(expected_times)
            error("$context failed: saved $(length(sol.t)) time point(s), expected $(length(expected_times))")
        end
        for i in eachindex(expected_times)
            expected_t = Float64(expected_times[i])
            actual_t = Float64(sol.t[i])
            time_tol = _endpoint_tolerance(expected_t; rtol=rtol, atol=atol)
            if !isfinite(actual_t) || abs(actual_t - expected_t) > time_tol
                error("$context failed: saved time $i was $actual_t, expected $expected_t within $time_tol")
            end
        end
    end

    for i in eachindex(sol.u)
        if !all(isfinite, sol.u[i])
            error("$context failed: non-finite state at saved index $i")
        end
    end

    _solution_diagnostics(context, endpoint, expected_endpoint;
                          retcode=sol.retcode, saved_points=length(sol.t))
end
