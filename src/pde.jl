# Time-dependent PDE solver: solve the full 1D slender model directly.
#
# 1D slender model (nondimensional, γ/ρ = 1):
#   ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   ∂u/∂t + u·∂u/∂z = ∂/∂z(1/R)
#
# Spatial discretisation: 4th-order finite differences on a stretched grid.
# Time integration: Tsit5 from DifferentialEquations.jl.
# IC: R(z,0) = εz (undisturbed cone), u(z,0) = 0.
# BC at z = L: outflow (zero-gradient).

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
Uses algebraic stretching: z = z_min + (z_max-z_min)·(η/β)·(1+β)/(1+η)
where η ∈ [0, β].
"""
function stretched_grid(N::Int, z_min::Float64, z_max::Float64; β::Float64=2.0)
    η = range(0, β, length=N) |> collect
    z_min .+ (z_max - z_min) .* (η ./ β) .* (1 + β) ./ (1 .+ η)
end

# ── Finite difference operators (4th-order interior, 2nd-order boundary)
function ddz!(df, f, z)
    N = length(f)
    for i in 3:N-2
        h = z[i+1] - z[i]  # assume locally uniform for FD stencil
        df[i] = (-f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2]) / (12h)
    end
    # 2nd-order at boundaries
    if N ≥ 3
        h = z[2] - z[1]
        df[1] = (-3f[1] + 4f[2] - f[3]) / (2h)
        df[2] = (f[3] - f[1]) / (z[3] - z[1])
        h = z[N] - z[N-1]
        df[N] = (3f[N] - 4f[N-1] + f[N-2]) / (2h)
        df[N-1] = (f[N] - f[N-2]) / (z[N] - z[N-2])
    end
end

# ── PDE right-hand side ───────────────────────────────────────────────
"""
The 1D slender model in primitive variables (R, u):
  Rt = -u·Rz - (R/2)·uz          [from mass, dividing by 2R]
  ut = -u·uz + ∂/∂z(1/R)         [momentum]
     = -u·uz - Rz/R²
"""
function pde_rhs!(dw, w, p, t)
    N, z, Rz_buf, uz_buf, invR_z_buf = p
    R = @view w[1:N]
    u = @view w[N+1:2N]
    dR = @view dw[1:N]
    du = @view dw[N+1:2N]

    # Enforce positivity
    for i in 1:N
        if R[i] < 1e-15
            R[i] = 1e-15
        end
    end

    # Compute spatial derivatives
    ddz!(Rz_buf, R, z)
    ddz!(uz_buf, u, z)

    # Compute ∂/∂z(1/R)
    invR = [1.0 / R[i] for i in 1:N]
    ddz!(invR_z_buf, invR, z)

    for i in 1:N
        dR[i] = -u[i] * Rz_buf[i] - 0.5 * R[i] * uz_buf[i]
        du[i] = -u[i] * uz_buf[i] + invR_z_buf[i]
    end

    # Outflow BC at right boundary: zero gradient
    dR[N] = dR[N-1]
    du[N] = du[N-1]

    # Symmetry-like BC at left boundary: u=0 maintained if tip is at z=0
    # Actually, for the cone tip we set dR[1] from the interior
    # and du[1] = 0 (no flow through the tip)
    # du[1] = 0.0  # optional: helps stability at the tip
end

# ── Solver ─────────────────────────────────────────────────────────────
"""
    solve_pde(; ε=0.1, N=200, z_min=0.01, z_max=10.0,
               t_end=1.0, n_snapshots=10)

Solve the 1D slender model PDE using method of lines.

Parameters:
- ε: cone half-angle
- N: number of grid points
- z_min: left boundary (near tip, > 0 to avoid singularity)
- z_max: right boundary
- t_end: final time
- n_snapshots: number of output snapshots

Returns a PDESolution.
"""
function solve_pde(; ε::Float64=0.1, N::Int=200, z_min::Float64=0.01,
                    z_max::Float64=10.0, t_end::Float64=1.0,
                    n_snapshots::Int=10)
    z = stretched_grid(N, z_min, z_max)

    # Initial conditions: R = εz, u = 0
    R0 = ε .* z
    u0 = zeros(N)
    w0 = vcat(R0, u0)

    # Buffers for spatial derivatives
    Rz_buf = zeros(N)
    uz_buf = zeros(N)
    invR_z_buf = zeros(N)

    p = (N, z, Rz_buf, uz_buf, invR_z_buf)
    tspan = (0.0, t_end)

    # Save at specified times
    t_save = range(0.0, t_end, length=n_snapshots+1) |> collect

    prob = ODEProblem(pde_rhs!, w0, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-8,
                saveat=t_save, maxiters=1_000_000, dtmin=1e-15)

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
