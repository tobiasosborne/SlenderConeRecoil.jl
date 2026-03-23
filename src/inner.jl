# Inner BVP solver: solve the nonlinear similarity ODE near the tip.
#
# The similarity ODE system (correct signs from Bernoulli):
#   2S + 2S'(U - ξ) + S·U' = 0              [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0    [momentum]
#
# Solved as an IVP by shooting from the tip ξ₀:
#   BCs at tip: S'(ξ₀) = 0 (smooth rounded tip),
#               S(ξ₀) > 0 (finite radius)
#   Far-field: S(ξ) ~ εξ as ξ → ∞
#
# Shooting parameters: (ξ₀, S₀) where S₀ = S(ξ₀).
# From S'(ξ₀) = 0 and the ODEs, we can derive U(ξ₀) and U'(ξ₀).

using DifferentialEquations

export solve_inner_bvp, InnerSolution

struct InnerSolution
    ξ::Vector{Float64}
    S::Vector{Float64}
    U::Vector{Float64}
    ξ₀::Float64
    S₀::Float64
end

# ── ODE right-hand side ────────────────────────────────────────────────
"""
Resolve the similarity ODEs for (S', U') given (S, U) at ξ.

From mass:      U' = -2 - 2S'(U-ξ)/S
From momentum:  S'/S² = -(2/9)U + (4/9)(U-ξ)U'   [note: sign flip vs wrong version]
                S' = S²[-(2/9)U + (4/9)(U-ξ)U']

Substitute mass into momentum:
  S' = S²[-(2/9)U + (4/9)(U-ξ)(-2 - 2S'(U-ξ)/S)]
  S' = S²[-(2/9)U - (8/9)(U-ξ) - (8/9)(U-ξ)²S'/S]
  S' + (8/9)S(U-ξ)²S' = S²[-(2/9)U - (8/9)(U-ξ)]
  S'[1 + (8/9)S(U-ξ)²] = -S²[(2/9)U + (8/9)(U-ξ)]
  S' = -S²[(2/9)U + (8/9)(U-ξ)] / [1 + (8/9)S(U-ξ)²]

Note: denominator 1 + (8/9)S(U-ξ)² > 0 always (no finite-ξ singularity).
"""
function inner_rhs!(du, u, p, ξ)
    S, U = u
    if S ≤ 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    v = U - ξ  # velocity in similarity frame

    # S' from resolved system (correct sign)
    numer = -S^2 * (2/9 * U + 8/9 * v)
    denom = 1 + 8/9 * S * v^2    # always positive

    Sξ = numer / denom
    Uξ = -2 - 2 * Sξ * v / S

    du[1] = Sξ
    du[2] = Uξ
end

# ── Initial conditions at the tip ──────────────────────────────────────
"""
    tip_initial_conditions(ξ₀, S₀)

Compute (S, U) at the tip ξ = ξ₀.
At the tip: S'(ξ₀) = 0 (symmetry).

From mass with S' = 0:   2S + S·U' = 0  →  U' = -2
From momentum with S' = 0:
  -(2/9)U + (4/9)(U - ξ₀)(-2) = 0
  -(2/9)U - (8/9)(U - ξ₀) = 0
  -(10/9)U + (8/9)ξ₀ = 0
  U = (4/5)ξ₀
"""
function tip_initial_conditions(ξ₀::Float64, S₀::Float64)
    U₀ = 4/5 * ξ₀
    (S₀, U₀)
end

# ── Shooting solver with Newton iteration ──────────────────────────────
"""
    solve_inner_bvp(; ξ₀=0.0, S₀=1.0, ξ_max=50.0, ε=0.1,
                      newton_iters=20, newton_tol=1e-8)

Solve the inner problem by shooting from the tip with Newton iteration
on S₀ to match the far-field condition S(ξ)/ξ → ε.

Parameters:
- ξ₀: tip position
- S₀: initial guess for tip radius S(ξ₀)
- ξ_max: integration endpoint
- ε: cone half-angle (far-field target S ~ εξ)
- newton_iters: max Newton iterations for shooting
- newton_tol: convergence tolerance

Returns an InnerSolution.
"""
function solve_inner_bvp(; ξ₀::Float64=0.0, S₀::Float64=1.0,
                           ξ_max::Float64=50.0, ε::Float64=0.1,
                           newton_iters::Int=20, newton_tol::Float64=1e-8)
    function shoot(S₀_trial)
        S_val, U_val = tip_initial_conditions(ξ₀, S₀_trial)
        u0 = [S_val, U_val]
        tspan = (ξ₀ + 1e-6, ξ_max)
        prob = ODEProblem(inner_rhs!, u0, tspan)
        sol = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-12,
                    maxiters=1_000_000)
        ξ_vals = sol.t
        S_vals = [sol.u[i][1] for i in eachindex(sol.u)]
        U_vals = [sol.u[i][2] for i in eachindex(sol.u)]
        # Far-field slope: S(ξ_end)/ξ_end should approach ε
        slope = S_vals[end] / ξ_vals[end]
        (ξ_vals, S_vals, U_vals, slope)
    end

    # Newton iteration on S₀ to match far-field slope = ε
    S₀_curr = S₀
    δ = 1e-6
    for iter in 1:newton_iters
        ξv, Sv, Uv, slope = shoot(S₀_curr)
        residual = slope - ε
        if abs(residual) < newton_tol
            return InnerSolution(ξv, Sv, Uv, ξ₀, S₀_curr)
        end
        # Finite difference Jacobian
        _, _, _, slope_p = shoot(S₀_curr + δ)
        dres_dS₀ = (slope_p - slope) / δ
        if abs(dres_dS₀) < 1e-15
            break
        end
        S₀_curr -= residual / dres_dS₀
        if S₀_curr ≤ 0
            S₀_curr = δ  # keep positive
        end
    end

    # Return best solution even if Newton didn't converge
    ξv, Sv, Uv, _ = shoot(S₀_curr)
    InnerSolution(ξv, Sv, Uv, ξ₀, S₀_curr)
end

# ── Far-field diagnostics ──────────────────────────────────────────────
"""
    far_field_slope(sol::InnerSolution)

Compute S(ξ)/ξ at the last few points. For a correct solution,
this should approach ε.
"""
function far_field_slope(sol::InnerSolution)
    n = length(sol.ξ)
    idx = max(1, n-10):n
    [sol.S[i] / sol.ξ[i] for i in idx if sol.ξ[i] > 0]
end
