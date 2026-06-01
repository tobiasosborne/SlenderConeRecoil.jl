# Matching and composite solution: combine inner and outer solutions.
#
# Ledger status, per
# docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md:
# - C2001: inner/outer regions and an intermediate-coordinate matching idea are
#   source-backed in the precursor, but in different A,R variables with
#   OCR-limited formulae and constants.
# - IMPL-inferred/BLOCKED-2008: the additive S_inner + S_outer - S_overlap
#   composite and the fitted linear common part below are local diagnostic
#   constructions, not source-confirmed Decent-King 2008 composite formulae.
#
# Inner solution: valid near the tip, S_inner(ξ), U_inner(ξ)
# Outer solution: S_outer = εξ + s₁(ξ), U_outer = u₁(ξ)
#
# Matching: in the overlap region, inner far-field ↔ outer near-field.
# Additive composite: S_comp = S_inner + S_outer - S_overlap
# where S_overlap is the common asymptotic form in the overlap region.

export composite_solution, CompositeSolution, overlap_residual

struct CompositeSolution
    ξ::Vector{Float64}
    S::Vector{Float64}
    U::Vector{Float64}
    ε::Float64
    diagnostics::NamedTuple
end

CompositeSolution(ξ, S, U, ε) =
    CompositeSolution(ξ, S, U, ε,
                      (overlap_slope=NaN, overlap_intercept=NaN,
                       ξ_match=NaN, ξ_min=NaN, ξ_max=NaN,
                       fit_points=0))

# ── Overlap extraction ─────────────────────────────────────────────────
"""
    inner_far_field(sol::InnerSolution, ξ_match)

Extract the asymptotic behavior of the inner solution for large ξ.
In the far field, S_inner ~ c₁·ξ (linear growth).
Returns (slope, intercept) from a linear fit to S(ξ) for ξ > ξ_match.
This is a local diagnostic common-part estimate, not a source-backed matching
constant.
"""
function inner_far_field(sol::InnerSolution, ξ_match::Float64)
    fit = _inner_far_field_fit(sol, ξ_match)
    (fit.slope, fit.intercept)
end

function _inner_far_field_fit(sol::InnerSolution, ξ_match::Float64)
    _validate_interpolation_table(sol.ξ, sol.S; context="inner far-field fit")
    _require_in_interpolation_domain(sol.ξ, ξ_match; context="inner far-field fit",
                                     x_name="ξ_match")
    idx = findall(ξ -> ξ > ξ_match, sol.ξ)
    if length(idx) < 3
        throw(ArgumentError("inner far-field fit requires at least three points with ξ > ξ_match; got $(length(idx)) for ξ_match=$ξ_match"))
    end
    ξs = sol.ξ[idx]
    Ss = sol.S[idx]
    # Linear fit: S ≈ a·ξ + b
    n = length(ξs)
    ξ_mean = sum(ξs) / n
    S_mean = sum(Ss) / n
    num = sum((ξs[i] - ξ_mean) * (Ss[i] - S_mean) for i in 1:n)
    den = sum((ξs[i] - ξ_mean)^2 for i in 1:n)
    slope = den > 0 ? num / den : 0.0
    intercept = S_mean - slope * ξ_mean
    (slope=slope, intercept=intercept, ξ_match=ξ_match, fit_points=n)
end

# ── Composite construction ─────────────────────────────────────────────
"""
    composite_solution(inner::InnerSolution, outer::OuterSolution;
                       ξ_grid=nothing, n_points=500, ξ_match=nothing)

Construct the additive composite solution:
  S_comp(ξ) = S_inner(ξ) + S_outer(ξ) - S_overlap(ξ)

where S_overlap = slope·ξ + intercept. This fitted common part is a local
diagnostic overlap estimate; the source-confirmed Decent-King composite and
matching constants remain blocked pending the 2008 article body.
When ξ_match is not supplied, the fit starts at the lower edge of the
common inner/outer overlap.

Both solutions are interpolated onto a common ξ grid.
"""
function composite_solution(inner::InnerSolution, outer::OuterSolution;
                            ξ_grid::Union{Nothing,Vector{Float64}}=nothing,
                            n_points::Int=500,
                            ξ_match::Union{Nothing,Float64}=nothing)
    ε = outer.ε
    _require_positive_finite_epsilon(ε; context="composite construction outer ε")

    # Determine overlap region: where both solutions are valid
    ξ_min, ξ_max = _common_overlap_interval(inner, outer; context="composite construction")

    if ξ_grid === nothing
        n_points >= 2 ||
            throw(ArgumentError("composite construction requires n_points >= 2; got $n_points"))
        ξ_grid = range(ξ_min, ξ_max, length=n_points) |> collect
    else
        _validate_query_grid(ξ_grid; context="composite construction ξ_grid")
        ξ_grid[1] >= ξ_min && ξ_grid[end] <= ξ_max ||
            throw(ArgumentError("composite construction ξ_grid must lie in common overlap [$ξ_min, $ξ_max]; got [$(ξ_grid[1]), $(ξ_grid[end])]"))
    end

    # Interpolate solutions onto common grid
    S_inner_interp = _interp_strict(inner.ξ, inner.S, ξ_grid; context="composite inner S")
    U_inner_interp = _interp_strict(inner.ξ, inner.U, ξ_grid; context="composite inner U")
    S_outer_interp = [ε * ξ + _interp_scalar_strict(outer.ξ, outer.s₁, ξ; context="composite outer s₁")
                      for ξ in ξ_grid]
    U_outer_interp = _interp_strict(outer.ξ, outer.u₁, ξ_grid; context="composite outer u₁")

    # Overlap: fitted common part from the inner far field.
    fit_ξ_match = something(ξ_match, ξ_min)
    _require_overlap_match_point(fit_ξ_match, ξ_min, ξ_max;
                                 context="composite construction")
    fit = _inner_far_field_fit(inner, fit_ξ_match)
    S_overlap = [fit.slope * ξ + fit.intercept for ξ in ξ_grid]
    U_overlap = zeros(length(ξ_grid))

    # Additive composite
    S_comp = S_inner_interp .+ S_outer_interp .- S_overlap
    U_comp = U_inner_interp .+ U_outer_interp .- U_overlap

    diagnostics = (overlap_slope=fit.slope,
                   overlap_intercept=fit.intercept,
                   ξ_match=fit.ξ_match,
                   ξ_min=ξ_min,
                   ξ_max=ξ_max,
                   fit_points=fit.fit_points)
    CompositeSolution(ξ_grid, S_comp, U_comp, ε, diagnostics)
end

# ── Simple linear interpolation ────────────────────────────────────────
function _interp(xs::Vector{Float64}, ys::Vector{Float64}, xq::Vector{Float64})
    [_interp_scalar(xs, ys, x) for x in xq]
end

function _interp_scalar(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    _interp_val(xs, ys, x)
end

function _interp_strict(xs::Vector{Float64}, ys::Vector{Float64}, xq::Vector{Float64};
                        context::AbstractString="strict interpolation")
    [_interp_scalar_strict(xs, ys, x; context=context) for x in xq]
end

function _interp_scalar_strict(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64;
                               context::AbstractString="strict interpolation")
    _interp_val_strict(xs, ys, x; context=context)
end

function _validate_query_grid(xs::Vector{Float64}; context::AbstractString)
    length(xs) >= 2 ||
        throw(ArgumentError("$context requires at least two points; got $(length(xs))"))
    all(isfinite, xs) ||
        throw(ArgumentError("$context requires finite coordinates"))
    for i in 2:length(xs)
        xs[i] > xs[i-1] ||
            throw(ArgumentError("$context requires a strictly increasing grid; found xs[$(i-1)]=$(xs[i-1]) and xs[$i]=$(xs[i])"))
    end
    nothing
end

function _common_overlap_interval(inner::InnerSolution, outer::OuterSolution;
                                  context::AbstractString)
    _validate_interpolation_table(inner.ξ, inner.S; context="$context inner S")
    _validate_interpolation_table(inner.ξ, inner.U; context="$context inner U")
    _validate_interpolation_table(outer.ξ, outer.s₁; context="$context outer s₁")
    _validate_interpolation_table(outer.ξ, outer.u₁; context="$context outer u₁")
    ξ_min = max(inner.ξ[1], outer.ξ[1])
    ξ_max = min(inner.ξ[end], outer.ξ[end])
    ξ_min < ξ_max ||
        throw(ArgumentError("$context requires a non-empty common ξ overlap; got [$ξ_min, $ξ_max]"))
    (ξ_min, ξ_max)
end

function _require_overlap_match_point(ξ_match::Float64, ξ_min::Float64, ξ_max::Float64;
                                      context::AbstractString)
    isfinite(ξ_match) ||
        throw(ArgumentError("$context requires finite ξ_match; got $ξ_match"))
    ξ_min <= ξ_match < ξ_max ||
        throw(ArgumentError("$context requires ξ_match in common overlap [$ξ_min, $ξ_max); got $ξ_match"))
    nothing
end

# ── Overlap residual diagnostic ────────────────────────────────────────
"""
    overlap_residual(inner::InnerSolution, outer::OuterSolution;
                     ξ_start=nothing, ξ_end=nothing)

Compute the mismatch between inner and outer solutions in the overlap region.
Returns the max |S_inner(ξ) - S_outer(ξ)| relative to S_outer. This is an
implementation diagnostic for the current reconstructed equations, not a
published Decent-King error norm.
"""
function overlap_residual(inner::InnerSolution, outer::OuterSolution;
                          ξ_start::Union{Nothing,Float64}=nothing,
                          ξ_end::Union{Nothing,Float64}=nothing)
    ε = outer.ε
    _require_positive_finite_epsilon(ε; context="overlap residual outer ε")
    ξ_min, ξ_max = _common_overlap_interval(inner, outer; context="overlap residual")
    ξ_lo = something(ξ_start, ξ_min)
    ξ_hi = something(ξ_end, ξ_max)
    isfinite(ξ_lo) ||
        throw(ArgumentError("overlap residual requires finite ξ_start; got $ξ_lo"))
    isfinite(ξ_hi) ||
        throw(ArgumentError("overlap residual requires finite ξ_end; got $ξ_hi"))
    ξ_lo >= ξ_min && ξ_hi <= ξ_max ||
        throw(ArgumentError("overlap residual interval [$ξ_lo, $ξ_hi] must lie in common overlap [$ξ_min, $ξ_max]"))
    ξ_lo < ξ_hi ||
        throw(ArgumentError("overlap residual requires ξ_start < ξ_end; got $ξ_lo and $ξ_hi"))

    ξ_test = range(ξ_lo, ξ_hi, length=50) |> collect

    max_res = 0.0
    for ξ in ξ_test
        Si = _interp_scalar_strict(inner.ξ, inner.S, ξ; context="overlap residual inner S")
        So = ε * ξ + _interp_scalar_strict(outer.ξ, outer.s₁, ξ; context="overlap residual outer s₁")
        ref = max(abs(So), 1e-10)
        res = abs(Si - So) / ref
        max_res = max(max_res, res)
    end
    max_res
end
