# Inner BVP solver: solve the nonlinear similarity ODE near the tip.
#
# The similarity ODE system:
#   2S + 2S'(U - ξ) + S·U' = 0        [mass]
#   -(2/9)U + (4/9)(U - ξ)U' + S'/S² = 0   [momentum]
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

From the mass equation:
  S·U' = -2S - 2S'(U - ξ)
  U' = -2 - 2S'(U - ξ)/S

From the momentum equation:
  S'/S² = (2/9)U - (4/9)(U - ξ)U'
  S' = S²·((2/9)U - (4/9)(U - ξ)U')

Strategy: first get U' from mass (requires S'), but S' depends on U' via
momentum. So we solve the coupled system.

From momentum:  S' = S²[(2/9)U - (4/9)(U-ξ)U']
From mass:      U' = -2 - 2S'(U-ξ)/S

Substitute mass into momentum:
  S' = S²[(2/9)U - (4/9)(U-ξ)(-2 - 2S'(U-ξ)/S)]
  S' = S²[(2/9)U + (8/9)(U-ξ) + (8/9)(U-ξ)²S'/S]
  S' - (8/9)S(U-ξ)²S' = S²[(2/9)U + (8/9)(U-ξ)]
  S'[1 - (8/9)S(U-ξ)²] = S²[(2/9)U + (8/9)(U-ξ)]
  S' = S²[(2/9)U + (8/9)(U-ξ)] / [1 - (8/9)S(U-ξ)²]

Then U' from mass.
"""
function inner_rhs!(du, u, p, ξ)
    S, U = u
    if S ≤ 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    v = U - ξ  # velocity in similarity frame

    # S' from resolved system
    numer = S^2 * (2/9 * U + 8/9 * v)
    denom = 1 - 8/9 * S * v^2

    if abs(denom) < 1e-12
        du[1] = 0.0
        du[2] = 0.0
        return
    end

    Sξ = numer / denom
    Uξ = -2 - 2 * Sξ * v / S

    du[1] = Sξ
    du[2] = Uξ
end

# ── Initial conditions at the tip ──────────────────────────────────────
"""
    tip_initial_conditions(ξ₀, S₀)

Compute (S, U, S', U') at the tip ξ = ξ₀.
At the tip: S'(ξ₀) = 0 (symmetry).
From the mass ODE with S' = 0:
  2S + S·U' = 0  →  U' = -2

From the momentum ODE with S' = 0:
  -(2/9)U + (4/9)(U - ξ₀)U' = 0
  -(2/9)U + (4/9)(U - ξ₀)(-2) = 0
  -(2/9)U - (8/9)(U - ξ₀) = 0
  -(2/9)U - (8/9)U + (8/9)ξ₀ = 0
  -(10/9)U + (8/9)ξ₀ = 0
  U = (4/5)ξ₀
"""
function tip_initial_conditions(ξ₀::Float64, S₀::Float64)
    U₀ = 4/5 * ξ₀
    (S₀, U₀)
end

# ── Shooting solver ────────────────────────────────────────────────────
"""
    solve_inner_bvp(; ξ₀=0.0, S₀=1.0, ξ_max=50.0, ε=0.1)

Solve the inner problem by shooting from the tip.

Parameters:
- ξ₀: tip position
- S₀: tip radius S(ξ₀)
- ξ_max: integration endpoint
- ε: cone half-angle (for far-field matching target S ~ εξ)

Returns an InnerSolution.
"""
function solve_inner_bvp(; ξ₀::Float64=0.0, S₀::Float64=1.0,
                           ξ_max::Float64=50.0, ε::Float64=0.1)
    S₀_val, U₀_val = tip_initial_conditions(ξ₀, S₀)
    u0 = [S₀_val, U₀_val]

    # Integrate from ξ₀ to ξ_max
    tspan = (ξ₀ + 1e-6, ξ_max)  # small offset to avoid singularity at tip
    prob = ODEProblem(inner_rhs!, u0, tspan)

    # Use a stiff-aware solver since the system can be stiff
    sol = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-12,
                maxiters=100_000, dtmin=1e-15)

    ξ_vals = sol.t
    S_vals = [sol.u[i][1] for i in eachindex(sol.u)]
    U_vals = [sol.u[i][2] for i in eachindex(sol.u)]

    InnerSolution(ξ_vals, S_vals, U_vals, ξ₀, S₀)
end

# ── Far-field diagnostics ──────────────────────────────────────────────
"""
    far_field_slope(sol::InnerSolution)

Compute S(ξ)/ξ at the last few points. For a correct solution matching
the far field, this should approach ε.
"""
function far_field_slope(sol::InnerSolution)
    n = length(sol.ξ)
    idx = max(1, n-10):n
    [sol.S[i] / sol.ξ[i] for i in idx if sol.ξ[i] > 0]
end
