# Local outer hierarchy experiment: CAS-derived algebra + nonlinear solver.
#
# Ledger status, per
# docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md:
# the 2001 precursor has an outer expansion with oscillatory harmonics and
# half-integer powers in different variables. The hierarchy below is a local
# reconstructed CAS experiment for the package's primitive S,U equations, not a
# source-confirmed Decent-King 2008 outer hierarchy.
#
# Supported CAS grammar is intentionally narrower than the source-specific
# matched-asymptotic structure: finite integer-power/Laurent ε series with
# ε-free symbolic coefficients. Half-integer and multiple-scale oscillatory
# source hierarchies are not implemented by derive_outer_equations.
#
# The CAS substitutes the ansatz S = εξ + ε³σ₁ + ε⁵σ₂, U = ε²ω₁ + ε⁴ω₂
# into the similarity ODEs, expands in ε, and collects at each order.
# This verifies the first-order equations and derives second-order forcing.
#
# The numerical solver goes beyond linearisation by integrating the FULL
# nonlinear similarity ODE outward from a matching point where the inner
# solution is well-resolved. This is a local diagnostic extension of the
# reconstructed S,U system, not a source-backed summation of the paper's
# asymptotic hierarchy.

using DifferentialEquations

export derive_outer_equations, eval_sexpr, solve_outer_full, solve_outer_linearised, HierarchySolution

# ── Numerical evaluation of symbolic expressions ─────────────────────
"""Evaluate a symbolic expression tree given numerical bindings."""
function eval_sexpr(e::SExpr, b::Dict{Symbol,Float64})::Float64
    if e isa Num; return Float64(e.val)
    elseif e isa Sym; return b[e.name]
    elseif e isa Add; return sum(eval_sexpr(t, b) for t in e.terms)
    elseif e isa Mul; return Float64(e.coeff) * prod(eval_sexpr(f, b) for f in e.factors)
    elseif e isa Pow; return eval_sexpr(e.base, b) ^ eval_sexpr(e.exp, b)
    elseif e isa Func
        a = eval_sexpr(e.args[1], b)
        e.name == :sin && return sin(a)
        e.name == :cos && return cos(a)
        e.name == :sqrt && return sqrt(a)
        error("unknown function $(e.name)")
    end
    error("eval_sexpr: unknown type")
end

# ── CAS derivation of the outer ODE hierarchy ────────────────────────
"""
    derive_outer_equations(; order=5)

Substitute the outer ansatz into the similarity ODEs, expand in ε,
and collect at each order. Returns Dict{Int, (mass, momentum)} of
symbolic equations for the local reconstructed S,U system.

The supported expansion family is
`S = εξ + ε³σ₁ + ε⁵σ₂`, `U = ε²ω₁ + ε⁴ω₂`, evaluated by the
integer-power/Laurent grammar in `expand_in`. This checks internal algebra
only; the ε-power structure and coefficient equations are local/internal and
are not source-confirmed Decent-King formulae. The half-integer and
multiple-scale oscillatory hierarchy described by the available 2001 precursor
is not implemented here.
"""
function derive_outer_equations(; order::Int=5)
    order >= 0 ||
        throw(ArgumentError("outer hierarchy derivation order must be nonnegative; got $order"))

    ε = Sym(:ε); ξ = Sym(:ξ)

    # Perturbation unknowns (scaled: s₁ ~ O(1), not O(ε³))
    σ₁  = Sym(:σ1);  σ₁p  = Sym(:σ1p);  σ₁ppp  = Sym(:σ1ppp)
    ω₁  = Sym(:ω1);  ω₁p  = Sym(:ω1p)
    σ₂  = Sym(:σ2);  σ₂p  = Sym(:σ2p);  σ₂ppp  = Sym(:σ2ppp)
    ω₂  = Sym(:ω2);  ω₂p  = Sym(:ω2p)

    # Outer ansatz: S = εξ + ε³σ₁ + ε⁵σ₂, U = ε²ω₁ + ε⁴ω₂
    ε3 = pow(ε, Num(3)); ε5 = pow(ε, Num(5))
    ε2 = pow(ε, Num(2)); ε4 = pow(ε, Num(4))

    S    = add(mul(ε, ξ), mul(ε3, σ₁), mul(ε5, σ₂))
    Sp   = add(ε, mul(ε3, σ₁p), mul(ε5, σ₂p))
    Sppp = add(mul(ε3, σ₁ppp), mul(ε5, σ₂ppp))
    U    = add(mul(ε2, ω₁), mul(ε4, ω₂))
    Up   = add(mul(ε2, ω₁p), mul(ε4, ω₂p))

    # Mass: 2S + 2S'(U - ξ) + SU'
    mass_raw = add(mul(Num(2), S),
                   mul(Num(2), Sp, add(U, neg(ξ))),
                   mul(S, Up))
    mass_exp = expand_in(mass_raw, ε, order)

    # Momentum: -(2/9)U + (4/9)(U-ξ)U' - S'/S² - S'''
    mom_raw = add(mul(Num(-2//9), U),
                  mul(Num(4//9), add(U, neg(ξ)), Up),
                  neg(mul(Sp, pow(S, Num(-2)))),
                  neg(Sppp))
    mom_exp = expand_in(mom_raw, ε, order)

    # Collect at each ε-order
    eqs = Dict{Int, @NamedTuple{mass::SExpr, momentum::SExpr}}()
    for k in -2:order
        mk = collect_order(mass_exp, ε, k)
        pk = collect_order(mom_exp, ε, k)
        z(e) = e isa Num && iszero(e.val)
        if !z(mk) || !z(pk)
            eqs[k] = (mass=mk, momentum=pk)
        end
    end
    eqs
end

# ── Full nonlinear outer solver ──────────────────────────────────────

struct HierarchySolution
    ξ::Vector{Float64}
    S::Vector{Float64}
    Sξ::Vector{Float64}
    Sξξ::Vector{Float64}
    U::Vector{Float64}
    ε::Float64
    diagnostics::NamedTuple
end

HierarchySolution(ξ, S, Sξ, Sξξ, U, ε) =
    HierarchySolution(ξ, S, Sξ, Sξξ, U, ε,
                      _manual_solution_diagnostics("HierarchySolution constructed directly", ξ))

"""
    solve_outer_full(inner; ξ_match=15.0, ξ_max=80.0)

Solve the FULL nonlinear similarity ODE outward from ξ_match, seeded
by the inner solution's state. This local diagnostic solve captures nonlinear
interactions in the reconstructed S,U system, but it is not a source-backed
higher-order outer hierarchy.

The result extends the inner solution into the far-field without the
drift that accumulates when shooting from the tip.
"""
function solve_outer_full(inner; ξ_match::Float64=15.0, ξ_max::Float64=80.0,
                          maxiters::Int=1_000_000)
    _require_matching_interval(inner, ξ_match, ξ_max; context="outer full solve")
    ε = inner.S[end] / inner.ξ[end]
    _require_positive_finite_epsilon(ε; context="outer full solve inferred ε")
    maxiters > 0 ||
        throw(ArgumentError("outer full solve requires maxiters > 0; got $maxiters"))

    # Extract inner state at matching point
    S_m  = _interp_val_strict(inner.ξ, inner.S, ξ_match; context="outer full solve S")
    Sp_m = _interp_val_strict(inner.ξ, inner.Sξ, ξ_match; context="outer full solve Sξ")
    Spp_m = _interp_val_strict(inner.ξ, inner.Sξξ, ξ_match; context="outer full solve Sξξ")
    U_m  = _interp_val_strict(inner.ξ, inner.U, ξ_match; context="outer full solve U")

    y0 = [S_m, Sp_m, Spp_m, U_m]
    prob = ODEProblem(inner_rhs!, y0, (ξ_match, ξ_max))
    sol = solve(prob, Rodas5P(); reltol=1e-10, abstol=1e-12, maxiters=maxiters)
    diagnostics = _require_successful_solution(sol, ξ_max; context="outer full solve")

    n = length(sol.t)
    HierarchySolution(
        sol.t,
        [sol.u[i][1] for i in 1:n],
        [sol.u[i][2] for i in 1:n],
        [sol.u[i][3] for i in 1:n],
        [sol.u[i][4] for i in 1:n],
        ε,
        diagnostics
    )
end

"""
    solve_outer_linearised(inner; ξ_match=15.0, ξ_max=80.0)

First-order linearised outer solver (for comparison with full nonlinear).
Solves the reconstructed linearised ODE from ξ_match outward. The comparison
is local/exploratory; it is not yet tied to source-confirmed matching constants
or Decent-King coefficient equations.
"""
function solve_outer_linearised(inner; ξ_match::Float64=15.0, ξ_max::Float64=80.0,
                                maxiters::Int=500_000)
    _require_matching_interval(inner, ξ_match, ξ_max; context="outer linearised hierarchy solve")
    ε = inner.S[end] / inner.ξ[end]
    _require_positive_finite_epsilon(ε; context="outer linearised hierarchy solve inferred ε")
    maxiters > 0 ||
        throw(ArgumentError("outer linearised hierarchy solve requires maxiters > 0; got $maxiters"))

    s_m  = _interp_val_strict(inner.ξ, inner.S, ξ_match; context="outer linearised hierarchy solve S") - ε * ξ_match
    sp_m = _interp_val_strict(inner.ξ, inner.Sξ, ξ_match; context="outer linearised hierarchy solve Sξ") - ε
    spp_m = _interp_val_strict(inner.ξ, inner.Sξξ, ξ_match; context="outer linearised hierarchy solve Sξξ")
    u_m  = _interp_val_strict(inner.ξ, inner.U, ξ_match; context="outer linearised hierarchy solve U")

    prob = ODEProblem(outer_rhs!, [s_m, sp_m, spp_m, u_m], (ξ_match, ξ_max), [ε])
    sol = solve(prob, Rodas5P(); reltol=1e-8, abstol=1e-10, maxiters=maxiters)
    diagnostics = _require_successful_solution(sol, ξ_max; context="outer linearised hierarchy solve")

    n = length(sol.t)
    S_lin = ε .* sol.t .+ [sol.u[i][1] for i in 1:n]
    Sp_lin = ε .+ [sol.u[i][2] for i in 1:n]
    Spp_lin = [sol.u[i][3] for i in 1:n]
    U_lin = [sol.u[i][4] for i in 1:n]

    HierarchySolution(sol.t, S_lin, Sp_lin, Spp_lin, U_lin, ε, diagnostics)
end
