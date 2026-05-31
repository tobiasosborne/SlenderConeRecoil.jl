# Time-dependent PDE solver: solve the full 1D slender model directly.
#
# 1D slender model (nondimensional, γ/ρ = 1):
#   ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   ∂u/∂t + u·∂u/∂z = -∂/∂z(1/R) = Rz/R²
#
# Spatial discretisation: 2nd-order finite differences on a stretched grid.
# Time integration: FBDF (implicit BDF) from DifferentialEquations.jl.
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
    diagnostics::NamedTuple
end

PDESolution(z, t_snapshots, R, u, ε) =
    PDESolution(z, t_snapshots, R, u, ε,
                _manual_solution_diagnostics("PDESolution constructed directly", t_snapshots))

function _require_finite_scalar(name::AbstractString, x)
    if !isfinite(x)
        throw(ArgumentError("$name must be finite; got $x"))
    end
end

function _validate_stretched_grid_args(N::Int, z_min::Float64, z_max::Float64, β::Float64)
    N >= 3 || throw(ArgumentError("stretched_grid requires N >= 3; got N=$N"))
    _require_finite_scalar("z_min", z_min)
    _require_finite_scalar("z_max", z_max)
    z_min < z_max || throw(ArgumentError("stretched_grid requires z_min < z_max; got z_min=$z_min, z_max=$z_max"))
    _require_finite_scalar("β", β)
    β > 0 || throw(ArgumentError("stretched_grid requires β > 0; got β=$β"))
end

function _validate_grid_vector(z; context::AbstractString)
    N = length(z)
    N >= 3 || throw(ArgumentError("$context requires at least 3 grid points; got $N"))
    for i in eachindex(z)
        if !isfinite(z[i])
            throw(ArgumentError("$context requires finite grid values; got z[$i]=$(z[i])"))
        end
    end
    for i in 2:N
        if !(z[i] > z[i-1])
            throw(ArgumentError("$context requires a strictly increasing grid; got z[$(i-1)]=$(z[i-1]), z[$i]=$(z[i])"))
        end
    end
    N
end

function _validate_ddz_inputs(df, f, z)
    N = length(f)
    length(df) == N || throw(ArgumentError("ddz! requires df and f to have the same length; got length(df)=$(length(df)), length(f)=$N"))
    length(z) == N || throw(ArgumentError("ddz! requires f and z to have the same length; got length(f)=$N, length(z)=$(length(z))"))
    _validate_grid_vector(z; context="ddz!")
    N
end

function _validate_solve_pde_args(ε::Float64, N::Int, z_min::Float64, z_max::Float64,
                                  t_end::Float64, n_snapshots::Int)
    _require_finite_scalar("ε", ε)
    ε > 0 || throw(ArgumentError("solve_pde requires ε > 0; got ε=$ε"))
    _validate_stretched_grid_args(N, z_min, z_max, 2.0)
    _require_finite_scalar("t_end", t_end)
    t_end >= 0 || throw(ArgumentError("solve_pde requires t_end >= 0; got t_end=$t_end"))
    n_snapshots >= 1 || throw(ArgumentError("solve_pde requires n_snapshots >= 1; got n_snapshots=$n_snapshots"))
end

function _validate_positive_initial_radius(R0, z, ε)
    for i in eachindex(R0)
        if !(isfinite(R0[i]) && R0[i] > 0)
            throw(ArgumentError("solve_pde requires positive finite initial radius; got R0[$i]=$(R0[i]) from ε=$ε and z[$i]=$(z[i])"))
        end
    end
end

# ── Grid construction ──────────────────────────────────────────────────
"""
    stretched_grid(N, z_min, z_max; β=2.0)

Create a stretched grid with N points, clustered near z_min (the tip).
Uses tanh-based stretching for finer resolution near the left boundary.
"""
function stretched_grid(N::Int, z_min::Float64, z_max::Float64; β::Float64=2.0)
    _validate_stretched_grid_args(N, z_min, z_max, β)
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
    _validate_ddz_inputs(df, f, z)
    _ddz_unchecked!(df, f, z)
end

function _ddz_unchecked!(df, f, z)
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

function _make_pde_rhs_parameters(z)
    N = _validate_grid_vector(z; context="pde_rhs! grid")
    # Serial-only mutable scratch cache. Each ODEProblem created by solve_pde owns
    # a fresh tuple, and callers must not share this p object across concurrent
    # pde_rhs! evaluations.
    (N, z, zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))
end

function _unpack_pde_rhs_parameters(p)
    if !(p isa Tuple) || length(p) != 8
        throw(ArgumentError("pde_rhs! expects parameter tuple (N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf, Rzz_buf, Rzzz_buf)"))
    end
    N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf, Rzz_buf, Rzzz_buf = p
    N isa Integer || throw(ArgumentError("pde_rhs! parameter N must be an integer; got $(typeof(N))"))
    N = Int(N)
    _validate_grid_vector(z; context="pde_rhs! grid") == N ||
        throw(ArgumentError("pde_rhs! parameter N=$N does not match length(z)=$(length(z))"))
    for (name, buf) in ((:Rz_buf, Rz_buf), (:uz_buf, uz_buf), (:invR_z_buf, invR_z_buf),
                        (:invR_buf, invR_buf), (:Rzz_buf, Rzz_buf), (:Rzzz_buf, Rzzz_buf))
        length(buf) == N || throw(ArgumentError("pde_rhs! requires $name to have length $N; got $(length(buf))"))
    end
    (N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf, Rzz_buf, Rzzz_buf)
end

function _validate_pde_state(w, N::Int)
    length(w) == 2N || throw(ArgumentError("pde_rhs! requires state length 2N=$((2N)); got $(length(w))"))
    for i in eachindex(w)
        if !isfinite(w[i])
            throw(DomainError(w[i], "pde_rhs! requires a finite state; got w[$i]=$(w[i])"))
        end
    end
    for i in 1:N
        if !(w[i] > 0)
            throw(DomainError(w[i], "pde_rhs! requires positive radius values; got R[$i]=$(w[i])"))
        end
    end
end

function _validate_pde_derivative_buffer(dw, N::Int)
    length(dw) == 2N || throw(ArgumentError("pde_rhs! requires derivative buffer length 2N=$((2N)); got $(length(dw))"))
end

# ── PDE right-hand side ───────────────────────────────────────────────
"""
    pde_rhs!(dw, w, p, t)

The 1D slender model with axial curvature, in primitive variables (R, u):
  Rt = -u·Rz - (R/2)·uz                    [mass]
  ut = -u·uz + Rz/R² + Rzzz                [momentum with κ = 1/R - Rzz]

The Rzzz term is dispersive and produces capillary waves.

`p` is a serial-only mutable cache tuple owned by one ODEProblem. Do not share
the same `p` object across concurrent RHS evaluations. The radius state must
remain finite and strictly positive; otherwise this RHS throws `DomainError`
before forming `1/R`.
"""
function pde_rhs!(dw, w, p, t)
    N, z, Rz_buf, uz_buf, invR_z_buf, invR_buf, Rzz_buf, Rzzz_buf =
        _unpack_pde_rhs_parameters(p)
    _validate_pde_derivative_buffer(dw, N)
    _validate_pde_state(w, N)

    R_orig = @view w[1:N]
    u_orig = @view w[N+1:2N]
    dR = @view dw[1:N]
    du = @view dw[N+1:2N]

    # Compute spatial derivatives
    _ddz_unchecked!(Rz_buf, R_orig, z)
    _ddz_unchecked!(uz_buf, u_orig, z)

    # ∂/∂z(1/R)
    invR_buf .= 1.0 ./ R_orig
    _ddz_unchecked!(invR_z_buf, invR_buf, z)

    # Rzzz = d³R/dz³ (apply ddz! three times)
    _ddz_unchecked!(Rzz_buf, Rz_buf, z)
    _ddz_unchecked!(Rzzz_buf, Rzz_buf, z)

    for i in 1:N
        dR[i] = -u_orig[i] * Rz_buf[i] - 0.5 * R_orig[i] * uz_buf[i]
        du[i] = -u_orig[i] * uz_buf[i] - invR_z_buf[i] + Rzzz_buf[i]
    end

    # Left BC: u = 0, R fixed, Rzzz = 0 (symmetry)
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
                    n_snapshots::Int=10, maxiters::Int=1_000_000)
    _validate_solve_pde_args(ε, N, z_min, z_max, t_end, n_snapshots)
    z = stretched_grid(N, z_min, z_max)

    # Initial conditions: R = εz, u = 0
    R0 = ε .* z
    _validate_positive_initial_radius(R0, z, ε)
    u0 = zeros(N)
    w0 = vcat(R0, u0)

    if t_end == 0.0
        diagnostics = _solution_diagnostics("PDE method-of-lines solve", 0.0, 0.0;
                                            retcode=:Success, saved_points=1)
        return PDESolution(z, [0.0], [copy(R0)], [copy(u0)], ε, diagnostics)
    end

    p = _make_pde_rhs_parameters(z)
    tspan = (0.0, t_end)

    # Save at specified times
    t_save = range(0.0, t_end, length=n_snapshots+1) |> collect

    prob = ODEProblem(pde_rhs!, w0, tspan, p)
    # Use implicit solver — capillary pressure makes the system stiff.
    # Use finite-diff Jacobian since RHS uses pre-allocated Float64 buffers.
    sol = solve(prob, FBDF(autodiff=false); reltol=1e-6, abstol=1e-8,
                saveat=t_save, maxiters=maxiters)
    diagnostics = _require_successful_solution(sol, t_end;
                                               context="PDE method-of-lines solve",
                                               expected_times=t_save)

    # Extract snapshots
    t_out = sol.t
    R_snapshots = [sol.u[i][1:N] for i in eachindex(sol.u)]
    u_snapshots = [sol.u[i][N+1:2N] for i in eachindex(sol.u)]

    PDESolution(z, t_out, R_snapshots, u_snapshots, ε, diagnostics)
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
