# Midpoint-collocation solver for the current reconstructed inner BVP.
#
# Ledger status: IMPL-inferred. This solves the same local primitive-variable
# BVP exposed by `inner_bvp_residual`; it is not a Decent-King 2008
# source-transcribed collocation formulation.

using LinearAlgebra: norm

export solve_inner_bvp_collocation

function _collocation_rhs(u, ξ)
    S, Sp, Spp, U = u
    if S <= 0
        return zeros(promote_type(eltype(u), typeof(ξ)), 4)
    end
    v = U - ξ
    Up = -2 - 2 * Sp * v / S
    Sppp = -(2/9) * U + (4/9) * v * Up - Sp / S^2
    [Sp, Spp, Sppp, Up]
end

function _validate_inner_collocation_args(ε::Float64, ξ₀::Float64, S₀::Float64,
                                          Sξξ₀::Float64, ξ_max::Float64,
                                          nodes::Int, max_nodes::Int,
                                          maxiters::Int,
                                          residual_tol::Float64,
                                          fd_step::Float64,
                                          damping_steps::Int)
    _require_positive_finite("ε", ε)
    _require_positive_finite("ξ₀", ξ₀)
    _require_positive_finite("S₀", S₀)
    _require_finite("Sξξ₀", Sξξ₀)
    _require_finite("ξ_max", ξ_max)
    ξ_max > ξ₀ ||
        throw(ArgumentError("ξ_max must be greater than ξ₀; got ξ_max=$ξ_max and ξ₀=$ξ₀"))
    nodes >= 3 || throw(ArgumentError("nodes must be at least 3; got $nodes"))
    max_nodes >= nodes ||
        throw(ArgumentError("max_nodes must be at least nodes; got max_nodes=$max_nodes and nodes=$nodes"))
    _require_positive_int("maxiters", maxiters)
    _require_positive_finite("residual_tol", residual_tol)
    _require_positive_finite("fd_step", fd_step)
    _require_positive_int("damping_steps", damping_steps)
    nothing
end

function _collocation_seed_from_vector(y::AbstractVector, old_nodes::Int,
                                       new_nodes::Int)
    old_states, unknowns = _collocation_unpack(y, old_nodes)
    old_η = collect(range(0.0, 1.0; length=old_nodes))
    new_η = collect(range(0.0, 1.0; length=new_nodes))
    new_states = zeros(Float64, 4, new_nodes)
    for row in 1:4, (j, ηj) in pairs(new_η)
        new_states[row, j] = _interp_linear(old_η, view(old_states, row, :), ηj)
    end
    _collocation_pack(new_states, unknowns)
end

function _interp_linear(xs, ys, x)
    i = searchsortedlast(xs, x)
    i = clamp(i, 1, length(xs) - 1)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    ys[i] + t * (ys[i+1] - ys[i])
end

function _collocation_pack(states::AbstractMatrix, unknowns::InnerBVPUnknowns)
    vcat(vec(states), [unknowns.ξ₀, unknowns.S₀, unknowns.Sξξ₀])
end

function _collocation_unpack(y::AbstractVector, nodes::Int)
    states = reshape(view(y, 1:4*nodes), 4, nodes)
    unknowns = InnerBVPUnknowns(y[4*nodes+1], y[4*nodes+2], y[4*nodes+3])
    states, unknowns
end

function _collocation_seed_from_solution(seed::InnerSolution, nodes::Int,
                                         ξ_max::Float64)
    η = collect(range(0.0, 1.0; length=nodes))
    states = zeros(Float64, 4, nodes)
    for (j, ηj) in pairs(η)
        ξ = seed.ξ₀ + ηj * (ξ_max - seed.ξ₀)
        states[:, j] .= (
            _interp_linear(seed.ξ, seed.S, ξ),
            _interp_linear(seed.ξ, seed.Sξ, ξ),
            _interp_linear(seed.ξ, seed.Sξξ, ξ),
            _interp_linear(seed.ξ, seed.U, ξ),
        )
    end
    _collocation_pack(states, InnerBVPUnknowns(seed.ξ₀, seed.S₀, seed.Sξξ₀))
end

function _collocation_residual_vector(y::AbstractVector, nodes::Int,
                                      ε::Float64, ξ_max::Float64)
    states, unknowns = _collocation_unpack(y, nodes)
    ξ₀ = unknowns.ξ₀
    scale = ξ_max - ξ₀
    η = range(0.0, 1.0; length=nodes)
    residual = Vector{Float64}(undef, 4 * (nodes - 1) + 7)
    k = 1
    for i in 1:(nodes - 1)
        h = Float64(η[i+1] - η[i])
        ηmid = Float64((η[i] + η[i+1]) / 2)
        ξmid = ξ₀ + ηmid * scale
        umid = (states[:, i] .+ states[:, i+1]) ./ 2
        defect = states[:, i+1] .- states[:, i] .-
                 h .* scale .* _collocation_rhs(umid, ξmid)
        residual[k:k+3] .= defect
        k += 4
    end

    left = states[:, 1]
    right = states[:, end]
    residual[k] = left[1] - unknowns.S₀
    residual[k+1] = left[2]
    residual[k+2] = left[3] - unknowns.Sξξ₀
    residual[k+3] = left[4] - 4/5 * unknowns.ξ₀
    residual[k+4] = right[1] / ξ_max - ε
    residual[k+5] = right[4]
    residual[k+6] = right[3]
    residual
end

function _collocation_newton(y0::Vector{Float64}, nodes::Int, ε::Float64,
                             ξ_max::Float64; maxiters::Int,
                             residual_tol::Float64, fd_step::Float64,
                             damping_steps::Int)
    y = copy(y0)
    F = _collocation_residual_vector(y, nodes, ε, ξ_max)
    F_norm = norm(F)
    converged = F_norm < residual_tol
    iterations = 0
    reason = converged ? "residual tolerance reached" : "iteration limit reached"

    for iter in 1:maxiters
        converged && break
        J = Matrix{Float64}(undef, length(F), length(y))
        for j in eachindex(y)
            h = fd_step * max(1.0, abs(y[j]))
            yp = copy(y)
            yp[j] += h
            J[:, j] .= (_collocation_residual_vector(yp, nodes, ε, ξ_max) .- F) ./ h
        end
        if !all(isfinite, J)
            reason = "non-finite finite-difference Jacobian"
            break
        end
        step = try
            J \ F
        catch err
            reason = "linearized collocation solve failed: $(sprint(showerror, err))"
            break
        end
        if !all(isfinite, step)
            reason = "non-finite Newton step"
            break
        end

        accepted = false
        for damp in 0:(damping_steps - 1)
            α = 0.5^damp
            trial = y .- α .* step
            trial_F = try
                _collocation_residual_vector(trial, nodes, ε, ξ_max)
            catch
                continue
            end
            trial_norm = norm(trial_F)
            if all(isfinite, trial_F) && trial_norm < F_norm
                y = trial
                F = trial_F
                F_norm = trial_norm
                accepted = true
                break
            end
        end
        iterations = iter
        if F_norm < residual_tol
            converged = true
            reason = "residual tolerance reached"
            break
        elseif !accepted
            reason = "line search failed to reduce collocation residual"
            break
        end
    end

    (y=y, residual=F, residual_norm=F_norm, iterations=iterations,
     converged=converged, termination_reason=reason)
end

function _collocation_solution_from_vector(y::AbstractVector, nodes::Int,
                                           ε::Float64, ξ_max::Float64,
                                           final_residual::Vector{Float64},
                                           final_residual_norm::Float64,
                                           iterations::Int,
                                           converged::Bool,
                                           termination_reason::AbstractString)
    states, unknowns = _collocation_unpack(y, nodes)
    ξ = collect(range(unknowns.ξ₀, ξ_max; length=nodes))
    InnerSolution(ξ, collect(states[1, :]), collect(states[2, :]),
                  collect(states[3, :]), collect(states[4, :]),
                  unknowns.ξ₀, unknowns.S₀, unknowns.Sξξ₀,
                  final_residual, final_residual_norm, iterations, converged,
                  String(termination_reason))
end

function _collocation_problem(ε::Float64, initial::InnerBVPUnknowns,
                              ξ_max::Float64, solver::NamedTuple,
                              provenance::NamedTuple)
    parameters = (ε=ε, ξ₀=initial.ξ₀, S₀=initial.S₀, Sξξ₀=initial.Sξξ₀)
    domain = (ξ_max=ξ_max,)
    ConeSimilarityProblem(parameters, domain, solver,
                          _merge_provenance(:inner_bvp_collocation, provenance;
                                            solver_settings=merge(parameters,
                                                                  domain,
                                                                  solver)))
end

function _collocation_result(problem::ConeSimilarityProblem,
                             solution::InnerSolution,
                             collocation_diagnostics::NamedTuple)
    result = ProblemResult(problem, solution)
    diagnostics = merge(result.diagnostics, collocation_diagnostics)
    ProblemResult{ConeSimilarityProblem,InnerSolution}(
        result.problem, result.solution, result.parameters, result.domain,
        result.mesh, diagnostics, result.provenance)
end

"""
    solve_inner_bvp_collocation(; kwargs...)

Solve the current reconstructed inner BVP with a finite-dimensional midpoint
collocation system over `nodes` points in the scaled coordinate `η ∈ [0, 1]`.
The unknown tip parameters `(ξ₀, S₀, Sξξ₀)` are part of the nonlinear system.

The default initial guess is an explicit shooting-seeded guess; this is recorded
in diagnostics and is not a silent fallback. If the collocation Newton solve
does not meet `residual_tol`, the function throws by default.

Ledger status: IMPL-inferred local reconstruction pending Decent-King 2008
article-body verification.
"""
function solve_inner_bvp_collocation(; ε::Real=0.1,
                                      epsilon::Union{Nothing,Real}=nothing,
                                      ξ₀::Real=2.79, S₀::Real=0.28,
                                      Sξξ₀::Real=0.57,
                                      ξ_max::Real=60.0,
                                      nodes::Int=16,
                                      max_nodes::Int=nodes,
                                      maxiters::Int=8,
                                      residual_tol::Real=1e-6,
                                      fd_step::Real=1e-6,
                                      damping_steps::Int=10,
                                      seed_newton_iters::Int=30,
                                      seed_newton_tol::Real=1e-4,
                                      seed_ode_maxiters::Int=2_000_000,
                                      throw_on_failure::Bool=true,
                                      provenance::NamedTuple=(;))
    εv = Float64(something(epsilon, ε))
    ξ₀v = Float64(ξ₀)
    S₀v = Float64(S₀)
    Sξξ₀v = Float64(Sξξ₀)
    ξ_maxv = Float64(ξ_max)
    tol = Float64(residual_tol)
    fd = Float64(fd_step)
    _validate_inner_collocation_args(εv, ξ₀v, S₀v, Sξξ₀v, ξ_maxv, nodes,
                                     max_nodes, maxiters, tol, fd,
                                     damping_steps)
    _require_nonnegative_finite("seed_newton_tol", Float64(seed_newton_tol))
    _require_positive_int("seed_ode_maxiters", seed_ode_maxiters)
    seed_newton_iters >= 0 ||
        throw(ArgumentError("seed_newton_iters must be nonnegative; got $seed_newton_iters"))

    seed = solve_inner_bvp(; ε=εv, ξ₀=ξ₀v, S₀=S₀v, Sξξ₀=Sξξ₀v,
                           ξ_max=ξ_maxv,
                           newton_iters=seed_newton_iters,
                           newton_tol=Float64(seed_newton_tol),
                           ode_maxiters=seed_ode_maxiters,
                           throw_on_failure=false)
    current_nodes = nodes
    refinement_steps = 0
    y0 = _collocation_seed_from_solution(seed, current_nodes, ξ_maxv)
    solve = _collocation_newton(y0, current_nodes, εv, ξ_maxv;
                                maxiters=maxiters,
                                residual_tol=tol,
                                fd_step=fd,
                                damping_steps=damping_steps)
    while !solve.converged && current_nodes < max_nodes
        next_nodes = min(max_nodes, 2 * current_nodes - 1)
        y0 = _collocation_seed_from_vector(solve.y, current_nodes, next_nodes)
        current_nodes = next_nodes
        refinement_steps += 1
        solve = _collocation_newton(y0, current_nodes, εv, ξ_maxv;
                                    maxiters=maxiters,
                                    residual_tol=tol,
                                    fd_step=fd,
                                    damping_steps=damping_steps)
    end

    states, unknowns = _collocation_unpack(solve.y, current_nodes)
    right = states[:, end]
    far_field_residual = Float64[right[1] / ξ_maxv - εv, right[4], right[3]]
    far_field_norm = norm(far_field_residual)
    solution = _collocation_solution_from_vector(
        solve.y, current_nodes, εv, ξ_maxv, far_field_residual, far_field_norm,
        solve.iterations, solve.converged, solve.termination_reason)

    if !solve.converged && throw_on_failure
        error("inner collocation solve failed: $(solve.termination_reason) after $(solve.iterations) iteration(s); collocation_residual_norm=$(solve.residual_norm)")
    end

    solver = (method=:midpoint_collocation,
              seed_source=:shooting,
              seed_newton_iters=seed_newton_iters,
              seed_newton_tol=Float64(seed_newton_tol),
              seed_ode_maxiters=seed_ode_maxiters,
              initial_nodes=nodes,
              nodes=current_nodes,
              max_nodes=max_nodes,
              maxiters=maxiters,
              residual_tol=tol,
              fd_step=fd,
              damping_steps=damping_steps,
              throw_on_failure=throw_on_failure)
    problem = _collocation_problem(εv, InnerBVPUnknowns(ξ₀v, S₀v, Sξξ₀v),
                                   ξ_maxv, solver, provenance)
    collocation_diagnostics =
        (problem_kind=:inner_bvp_collocation,
         solver_backend=:midpoint_collocation,
         seed_source=:shooting,
         seed_converged=seed.converged,
         seed_residual_norm=seed.final_residual_norm,
         collocation_converged=solve.converged,
         collocation_iterations=solve.iterations,
         collocation_refinement_steps=refinement_steps,
         collocation_residual_norm=solve.residual_norm,
         collocation_residual_length=length(solve.residual),
         collocation_unknowns=length(solve.y),
         collocation_initial_nodes=nodes,
         collocation_nodes=current_nodes,
         collocation_max_nodes=max_nodes,
         termination_reason=solve.termination_reason,
         successful=solve.converged && isfinite(solve.residual_norm),
         source_status=_DEFAULT_SOURCE_STATUS)
    _collocation_result(problem, solution, collocation_diagnostics)
end

function solve_inner_bvp_collocation(problem::ConeSimilarityProblem; kwargs...)
    seed_options = (seed_newton_iters=problem.solver.newton_iters,
                    seed_newton_tol=problem.solver.newton_tol,
                    seed_ode_maxiters=problem.solver.ode_maxiters)
    opts = merge(problem.parameters, problem.domain, seed_options,
                 (provenance=problem.provenance,), (; kwargs...))
    solve_inner_bvp_collocation(; opts...)
end
