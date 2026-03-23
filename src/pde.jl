# Time-dependent PDE solver: solve the full 1D slender model directly.
#
# 1D slender model (nondimensional, γ/ρ = 1):
#   ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   ∂u/∂t + u·∂u/∂z = -∂/∂z(1/R) = Rz/R²
#
# Spatial discretisation: 2nd-order finite differences on a stretched grid.
# Time integration: Tsit5 from DifferentialEquations.jl.
# IC: R(z,0) = εz (undisturbed cone), u(z,0) = 0.
# BC at z = L: outflow (zero-gradient).
# BC at z = z_min: u = 0 (no flow through truncated tip), dR/dt from interior.

using DifferentialEquations

export solve_pde, PDESolution, rescale_to_similarity

struct PDESolution
    z::Vector{Float64}
    t_snapshots::Vector{Float64}
    R::Vector{Vector{Float64}}   # R[i] = R(z, t_snapshots[i])
    u::Vector{Vector{Float64}}   # u[i] = u(z, t_snapshots[i])
    ε::Float64
end

# ── Grid construction ──────────────────────────────────────────────────
"""
    stretched_grid(N, z_min, z_max; β=2.0)

Create a stretched grid with N points, clustered near z_min (the tip).
Uses tanh-based stretching for finer resolution near the left boundary.
"""
function stretched_grid(N::Int, z_min::Float64, z_max::Float64; β::Float64=2.0)
    s = range(0, 1, length=N) |> collect
    # Cluster near s=0 (z_min) using tanh stretching
    z_min .+ (z_max - z_min) .* (1 .- tanh.(β .* (1 .- s)) ./ tanh(β))
end

# ── Finite difference on non-uniform grid (2nd-order) ──────────────────
"""
Compute df/dz on a non-uniform grid using 2nd-order finite differences.
Uses proper non-uniform stencils: central differences with variable spacing.
"""
function ddz!(df, f, z)
    N = length(f)
    # 2nd-order forward difference at left boundary (non-uniform grid)
    h1 = z[2] - z[1]
    h2 = z[3] - z[2]
    df[1] = (-h2*(2h1+h2)*f[1] + (h1+h2)^2*f[2] - h1^2*f[3]) / (h1*h2*(h1+h2))

    # Central differences for interior (non-uniform spacing)
    for i in 2:N-1
        hm = z[i] - z[i-1]
        hp = z[i+1] - z[i]
        df[i] = (f[i+1]*hm^2 + f[i]*(hp^2 - hm^2) - f[i-1]*hp^2) /
                (hp * hm * (hp + hm))
    end

    # 2nd-order backward difference at right boundary (non-uniform grid)
    hm1 = z[N] - z[N-1]
    hm2 = z[N-1] - z[N-2]
    df[N] = (hm2*(2hm1+hm2)*f[N] - (hm1+hm2)^2*f[N-1] + hm1^2*f[N-2]) /
            (hm1*hm2*(hm1+hm2))
end

# ── PDE right-hand side ───────────────────────────────────────────────
"""
The 1D slender model in primitive variables (R, u):
  Rt = -u·Rz - (R/2)·uz          [from mass, dividing by 2R]
  ut = -u·uz - ∂/∂z(1/R)         [momentum, correct sign]
     = -u·uz + Rz/R²
"""
function pde_rhs!(dw, w, p, t)
    N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf = p
    R_orig = @view w[1:N]
    u_orig = @view w[N+1:2N]
    dR = @view dw[1:N]
    du = @view dw[N+1:2N]

    # Use R directly (no clamping needed if IC is positive and solver is stable)

    # Compute spatial derivatives
    ddz!(Rz_buf, R_orig, z)
    ddz!(uz_buf, u_orig, z)

    # Compute ∂/∂z(1/R) — no allocation
    invR_buf .= 1.0 ./ R_orig
    ddz!(invR_z_buf, invR_buf, z)

    for i in 1:N
        dR[i] = -u_orig[i] * Rz_buf[i] - 0.5 * R_orig[i] * uz_buf[i]
        du[i] = -u_orig[i] * uz_buf[i] - invR_z_buf[i]  # correct sign: -∂/∂z(1/R) = Rz/R²
    end

    # Left BC at truncated tip: u = 0, R fixed (no unphysical shrinkage)
    du[1] = 0.0
    dR[1] = 0.0

    # Outflow BC at right boundary: zero gradient
    dR[N] = dR[N-1]
    du[N] = du[N-1]
end

# ── Solver ─────────────────────────────────────────────────────────────
"""
    solve_pde(; ε=0.1, N=200, z_min=0.01, z_max=10.0,
               t_end=1.0, n_snapshots=10)

Solve the 1D slender model PDE using method of lines.
"""
function solve_pde(; ε::Float64=0.1, N::Int=200, z_min::Float64=0.01,
                    z_max::Float64=10.0, t_end::Float64=1.0,
                    n_snapshots::Int=10)
    z = stretched_grid(N, z_min, z_max)

    # Initial conditions: R = εz, u = 0
    R0 = ε .* z
    u0 = zeros(N)
    w0 = vcat(R0, u0)

    # Buffers for spatial derivatives (no allocation in RHS)
    Rz_buf = zeros(N)
    uz_buf = zeros(N)
    invR_z_buf = zeros(N)
    invR_buf = zeros(N)

    p = (N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf)
    tspan = (0.0, t_end)

    # Save at specified times
    t_save = range(0.0, t_end, length=n_snapshots+1) |> collect

    prob = ODEProblem(pde_rhs!, w0, tspan, p)
    # Use implicit solver — capillary pressure makes the system stiff.
    # Use finite-diff Jacobian since RHS uses pre-allocated Float64 buffers.
    sol = solve(prob, FBDF(autodiff=false); reltol=1e-6, abstol=1e-8,
                saveat=t_save, maxiters=1_000_000)

    # Extract snapshots
    t_out = sol.t
    R_snapshots = [sol.u[i][1:N] for i in eachindex(sol.u)]
    u_snapshots = [sol.u[i][N+1:2N] for i in eachindex(sol.u)]

    PDESolution(z, t_out, R_snapshots, u_snapshots, ε)
end

# ── Rescale to similarity variables ────────────────────────────────────
"""
    rescale_to_similarity(pde::PDESolution, t_idx::Int)

Rescale the PDE snapshot at time index t_idx to similarity variables:
  ξ = z/t^{2/3}, S = R/t^{2/3}

Returns (ξ, S) vectors.
"""
function rescale_to_similarity(pde::PDESolution, t_idx::Int)
    t = pde.t_snapshots[t_idx]
    if t ≤ 0
        return (pde.z, pde.R[t_idx])  # can't rescale at t=0
    end
    ℓ = t^(2/3)
    ξ = pde.z ./ ℓ
    S = pde.R[t_idx] ./ ℓ
    (ξ, S)
end
