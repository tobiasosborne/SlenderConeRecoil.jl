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

export composite_solution, CompositeSolution, overlap_residual,
       overlap_window_diagnostics

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
                       fit_points=0, overlap_window=(;),
                       mismatch_norm=NaN, fit_coefficients=(;),
                       sensitivity=(;), truncation_sensitivity=(;),
                       successful=false))

function _composite_parts(inner::InnerSolution, outer::OuterSolution,
                          ξ_grid::Vector{Float64},
                          S_inner_interp::Vector{Float64},
                          U_inner_interp::Vector{Float64},
                          S_outer_interp::Vector{Float64},
                          U_outer_interp::Vector{Float64},
                          S_overlap::Vector{Float64},
                          fit, ξ_min::Float64, ξ_max::Float64)
    _composite_parts(; ξ=ξ_grid, S_inner=S_inner_interp,
                     U_inner=U_inner_interp, S_outer=S_outer_interp,
                     U_outer=U_outer_interp, S_common=S_overlap,
                     U_common=zeros(length(ξ_grid)), ε=outer.ε, fit=fit,
                     ξ_min=ξ_min, ξ_max=ξ_max)
end

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
    overlap_diagnostics = overlap_window_diagnostics(inner, outer;
                                                     ξ_start=fit_ξ_match,
                                                     ξ_end=ξ_max)
    fit = (slope=overlap_diagnostics.fit_coefficients.slope,
           intercept=overlap_diagnostics.fit_coefficients.intercept,
           ξ_match=fit_ξ_match,
           fit_ξ_min=overlap_diagnostics.window.ξ_min,
           fit_ξ_max=overlap_diagnostics.window.ξ_max,
           fit_points=overlap_diagnostics.window.fit_points)
    S_overlap = [fit.slope * ξ + fit.intercept for ξ in ξ_grid]
    U_overlap = zeros(length(ξ_grid))
    parts = _composite_parts(; ξ=ξ_grid, S_inner=S_inner_interp,
                             U_inner=U_inner_interp, S_outer=S_outer_interp,
                             U_outer=U_outer_interp, S_common=S_overlap,
                             U_common=U_overlap, ε=ε, fit=fit,
                             ξ_min=ξ_min, ξ_max=ξ_max)

    # Additive composite
    S_comp = S_inner_interp .+ S_outer_interp .- S_overlap
    U_comp = U_inner_interp .+ U_outer_interp .- U_overlap

    diagnostics = (overlap_slope=fit.slope,
                   overlap_intercept=fit.intercept,
                   ξ_match=fit.ξ_match,
                   ξ_min=ξ_min,
                   ξ_max=ξ_max,
                   fit_points=fit.fit_points,
                   overlap_window=overlap_diagnostics.window,
                   mismatch_norm=overlap_diagnostics.mismatch_norm,
                   mismatch_rms=overlap_diagnostics.mismatch_rms,
                   mismatch_max_relative=overlap_diagnostics.mismatch_max_relative,
                   fit_coefficients=overlap_diagnostics.fit_coefficients,
                   fit_residual_rms=overlap_diagnostics.fit_residual_rms,
                   sensitivity=overlap_diagnostics.sensitivity,
                   truncation_sensitivity=overlap_diagnostics.truncation_sensitivity,
                   overlap_stable=overlap_diagnostics.stable,
                   overlap_classification=overlap_diagnostics.classification,
                   overlap_diagnostics=overlap_diagnostics,
                   common_part=parts.common,
                   parts=parts,
                   source_status=_DEFAULT_SOURCE_STATUS,
                   successful=true)
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

function _require_overlap_window(ξ_lo::Float64, ξ_hi::Float64,
                                 ξ_min::Float64, ξ_max::Float64;
                                 context::AbstractString)
    isfinite(ξ_lo) ||
        throw(ArgumentError("$context requires finite ξ_start; got $ξ_lo"))
    isfinite(ξ_hi) ||
        throw(ArgumentError("$context requires finite ξ_end; got $ξ_hi"))
    ξ_lo >= ξ_min && ξ_hi <= ξ_max ||
        throw(ArgumentError("$context window [$ξ_lo, $ξ_hi] must lie in common overlap [$ξ_min, $ξ_max]"))
    ξ_lo < ξ_hi ||
        throw(ArgumentError("$context requires ξ_start < ξ_end; got $ξ_lo and $ξ_hi"))
    nothing
end

function _fit_inner_window(sol::InnerSolution, ξ_lo::Float64, ξ_hi::Float64;
                           min_points::Int, context::AbstractString)
    fit = _fit_inner_polynomial_window(sol, ξ_lo, ξ_hi; order=1,
                                       min_points=min_points,
                                       context=context)
    (slope=fit.slope, intercept=fit.intercept, ξ_match=ξ_lo,
     fit_ξ_min=ξ_lo, fit_ξ_max=ξ_hi, fit_points=fit.fit_points,
     residual_rms=fit.residual_rms, coefficients=fit.coefficients)
end

function _outer_radius_at(outer::OuterSolution, ξ::Float64; context::AbstractString)
    outer.ε * ξ + _interp_scalar_strict(outer.ξ, outer.s₁, ξ; context=context)
end

function _fit_inner_polynomial_window(sol::InnerSolution, ξ_lo::Float64,
                                      ξ_hi::Float64; order::Int,
                                      min_points::Int,
                                      context::AbstractString)
    _validate_interpolation_table(sol.ξ, sol.S; context="$context inner S")
    order >= 0 ||
        throw(ArgumentError("$context requires nonnegative polynomial order; got $order"))
    min_points >= 3 ||
        throw(ArgumentError("$context requires min_points >= 3; got $min_points"))
    required_points = max(min_points, order + 1)
    idx = findall(i -> sol.ξ[i] > ξ_lo && sol.ξ[i] <= ξ_hi, eachindex(sol.ξ))
    if length(idx) < required_points
        idx = findall(i -> sol.ξ[i] >= ξ_lo && sol.ξ[i] <= ξ_hi,
                      eachindex(sol.ξ))
    end
    length(idx) >= required_points ||
        throw(ArgumentError("$context requires at least $required_points inner points in [$ξ_lo, $ξ_hi] for order $order; got $(length(idx))"))

    ξs = sol.ξ[idx]
    Ss = sol.S[idx]
    n = length(ξs)
    A = [ξs[i]^(j - 1) for i in 1:n, j in 1:(order + 1)]
    coeffs = Float64.(A \ Ss)
    residuals = Ss .- A * coeffs
    residual_rms = sqrt(sum(abs2, residuals) / n)
    intercept = coeffs[1]
    slope = order >= 1 ? coeffs[2] : 0.0
    (order=order,
     coefficients=Tuple(coeffs),
     slope=slope,
     intercept=intercept,
     ξ_match=ξ_lo,
     fit_ξ_min=ξ_lo,
     fit_ξ_max=ξ_hi,
     fit_points=n,
     residual_rms=residual_rms)
end

function _fit_order_tuple(fit_orders, context::AbstractString)
    values = fit_orders isa Integer ? (fit_orders,) : Tuple(fit_orders)
    !isempty(values) ||
        throw(ArgumentError("$context requires at least one fit order"))
    orders = Int[]
    for order in values
        order isa Integer ||
            throw(ArgumentError("$context fit_orders must contain integers; got $(repr(order))"))
        order_int = Int(order)
        order_int >= 0 ||
            throw(ArgumentError("$context fit_orders must be nonnegative; got $order_int"))
        push!(orders, order_int)
    end
    length(unique(orders)) == length(orders) ||
        throw(ArgumentError("$context fit_orders must not contain duplicates; got $orders"))
    Tuple(orders)
end

function _polyval(coefficients, ξ::Float64)
    value = 0.0
    power = 1.0
    for coefficient in coefficients
        value += coefficient * power
        power *= ξ
    end
    value
end

function _common_outer_mismatch(coefficients, outer::OuterSolution,
                                ξ_lo::Float64, ξ_hi::Float64,
                                sample_points::Int;
                                context::AbstractString)
    sample_points >= 2 ||
        throw(ArgumentError("$context requires sample_points >= 2; got $sample_points"))
    ξs = range(ξ_lo, ξ_hi, length=sample_points)
    sq = 0.0
    refsq = 0.0
    max_relative = 0.0
    for ξ in ξs
        ξ_value = Float64(ξ)
        common = _polyval(coefficients, ξ_value)
        outer_radius = _outer_radius_at(outer, ξ_value;
                                        context="$context outer s₁")
        δ = common - outer_radius
        sq += δ^2
        refsq += outer_radius^2
        max_relative = max(max_relative,
                           abs(δ) / max(abs(outer_radius), 1e-10))
    end
    rms = sqrt(sq / sample_points)
    relative_rms = rms / max(sqrt(refsq / sample_points), 1e-10)
    (rms=rms, relative_rms=relative_rms, max_relative=max_relative,
     sample_points=sample_points)
end

function _window_mismatch(inner::InnerSolution, outer::OuterSolution,
                          ξ_lo::Float64, ξ_hi::Float64, sample_points::Int;
                          context::AbstractString)
    sample_points >= 2 ||
        throw(ArgumentError("$context requires sample_points >= 2; got $sample_points"))
    ξs = range(ξ_lo, ξ_hi, length=sample_points)
    sq = 0.0
    refsq = 0.0
    max_relative = 0.0
    for ξ in ξs
        Si = _interp_scalar_strict(inner.ξ, inner.S, Float64(ξ);
                                   context="$context inner S")
        So = _outer_radius_at(outer, Float64(ξ); context="$context outer s₁")
        δ = Si - So
        sq += δ^2
        refsq += So^2
        max_relative = max(max_relative, abs(δ) / max(abs(So), 1e-10))
    end
    rms = sqrt(sq / sample_points)
    relative_rms = rms / max(sqrt(refsq / sample_points), 1e-10)
    (rms=rms, relative_rms=relative_rms, max_relative=max_relative,
     sample_points=sample_points)
end

function _relative_span(values)
    isempty(values) && return NaN
    lo = minimum(values)
    hi = maximum(values)
    scale = max(maximum(abs, values), 1e-12)
    (hi - lo) / scale
end

function _truncation_order_sensitivity(inner::InnerSolution,
                                       outer::OuterSolution,
                                       ξ_lo::Float64, ξ_hi::Float64,
                                       fit_orders, min_points::Int,
                                       sample_points::Int,
                                       tolerance::Float64)
    order_results = map(fit_orders) do order
        fit = _fit_inner_polynomial_window(
            inner, ξ_lo, ξ_hi; order=order, min_points=min_points,
            context="overlap window diagnostics truncation-order fit")
        mismatch = _common_outer_mismatch(
            fit.coefficients, outer, ξ_lo, ξ_hi, sample_points;
            context="overlap window diagnostics truncation-order")
        (order=order,
         coefficients=fit.coefficients,
         fit_points=fit.fit_points,
         fit_residual_rms=fit.residual_rms,
         common_outer_mismatch_norm=mismatch.relative_rms,
         common_outer_mismatch_rms=mismatch.rms,
         common_outer_mismatch_max_relative=mismatch.max_relative)
    end
    fit_residual_relative_range =
        _relative_span([result.fit_residual_rms for result in order_results])
    mismatch_relative_range =
        _relative_span([result.common_outer_mismatch_norm for result in order_results])
    sensitivity_norm = max(fit_residual_relative_range,
                           mismatch_relative_range)
    stable = isfinite(sensitivity_norm) && sensitivity_norm <= tolerance
    (orders=Tuple(fit_orders),
     results=Tuple(order_results),
     fit_residual_relative_range=fit_residual_relative_range,
     common_outer_mismatch_relative_range=mismatch_relative_range,
     sensitivity_norm=sensitivity_norm,
     tolerance=tolerance,
     stable=stable)
end

function _overlap_sensitivity(windows, tolerance::Float64, truncation)
    slopes = [w.fit_coefficients.slope for w in windows]
    intercepts = [w.fit_coefficients.intercept for w in windows]
    mismatches = [w.mismatch_norm for w in windows]
    slope_relative_range = _relative_span(slopes)
    intercept_relative_range = _relative_span(intercepts)
    coefficient_relative_range = max(slope_relative_range,
                                     intercept_relative_range)
    mismatch_relative_range = _relative_span(mismatches)
    window_sensitivity_norm = max(coefficient_relative_range,
                                  mismatch_relative_range)
    truncation_norm = get(truncation, :sensitivity_norm, NaN)
    sensitivity_norm = max(window_sensitivity_norm, truncation_norm)
    stable = isfinite(sensitivity_norm) && sensitivity_norm <= tolerance
    (slope_relative_range=slope_relative_range,
     intercept_relative_range=intercept_relative_range,
     coefficient_relative_range=coefficient_relative_range,
     mismatch_relative_range=mismatch_relative_range,
     window_sensitivity_norm=window_sensitivity_norm,
     truncation=truncation,
     sensitivity_norm=sensitivity_norm,
     tolerance=tolerance,
     stable=stable)
end

"""
    overlap_window_diagnostics(inner::InnerSolution, outer::OuterSolution;
                               ξ_start=nothing, ξ_end=nothing,
                               window_fraction=0.5, window_count=5,
                               fit_orders=(1, 2),
                               min_points=3, sample_points=50,
                               stability_tolerance=0.1)

Report local overlap-window diagnostics for the current reconstructed matching
workflow. The selected window defaults to the common inner/outer domain; nested
trailing subwindows are used to estimate sensitivity of the fitted linear
common part and mismatch norm. Polynomial common-part fits over `fit_orders`
report sensitivity to truncation order. These diagnostics are IMPL-inferred/local
until the Decent-King 2008 article-body matching formulae are available.
"""
function overlap_window_diagnostics(inner::InnerSolution, outer::OuterSolution;
                                    ξ_start::Union{Nothing,Real}=nothing,
                                    ξ_end::Union{Nothing,Real}=nothing,
                                    window_fraction::Real=0.5,
                                    window_count::Int=5,
                                    fit_orders=(1, 2),
                                    min_points::Int=3,
                                    sample_points::Int=50,
                                    stability_tolerance::Real=0.1)
    ε = outer.ε
    _require_positive_finite_epsilon(ε; context="overlap window diagnostics outer ε")
    ξ_min, ξ_max = _common_overlap_interval(inner, outer;
                                            context="overlap window diagnostics")
    ξ_lo = ξ_start === nothing ? ξ_min : Float64(ξ_start)
    ξ_hi = ξ_end === nothing ? ξ_max : Float64(ξ_end)
    _require_overlap_window(ξ_lo, ξ_hi, ξ_min, ξ_max;
                            context="overlap window diagnostics")

    window_fraction_value = Float64(window_fraction)
    isfinite(window_fraction_value) && 0 < window_fraction_value < 1 ||
        throw(ArgumentError("overlap window diagnostics requires 0 < window_fraction < 1; got $window_fraction"))
    window_count >= 2 ||
        throw(ArgumentError("overlap window diagnostics requires window_count >= 2; got $window_count"))
    min_points >= 3 ||
        throw(ArgumentError("overlap window diagnostics requires min_points >= 3; got $min_points"))
    sample_points >= 2 ||
        throw(ArgumentError("overlap window diagnostics requires sample_points >= 2; got $sample_points"))
    tolerance = Float64(stability_tolerance)
    isfinite(tolerance) && tolerance >= 0 ||
        throw(ArgumentError("overlap window diagnostics requires nonnegative finite stability_tolerance; got $stability_tolerance"))
    fit_order_values = _fit_order_tuple(fit_orders,
                                        "overlap window diagnostics")

    fit = _fit_inner_window(inner, ξ_lo, ξ_hi; min_points=min_points,
                            context="overlap window diagnostics selected fit")
    mismatch = _window_mismatch(inner, outer, ξ_lo, ξ_hi, sample_points;
                                context="overlap window diagnostics selected")

    width = ξ_hi - ξ_lo
    narrowest_lower = ξ_hi - window_fraction_value * width
    candidate_lowers = collect(range(ξ_lo, narrowest_lower,
                                     length=window_count))
    windows = NamedTuple[]
    skipped_window_count = 0
    for lower in candidate_lowers
        lower64 = Float64(lower)
        try
            candidate_fit = _fit_inner_window(
                inner, lower64, ξ_hi; min_points=min_points,
                context="overlap window diagnostics sensitivity fit")
            candidate_mismatch = _window_mismatch(
                inner, outer, lower64, ξ_hi, sample_points;
                context="overlap window diagnostics sensitivity")
            push!(windows,
                  (ξ_min=lower64,
                   ξ_max=ξ_hi,
                   width=ξ_hi - lower64,
                   fit_points=candidate_fit.fit_points,
                   fit_residual_rms=candidate_fit.residual_rms,
                   fit_coefficients=(slope=candidate_fit.slope,
                                     intercept=candidate_fit.intercept),
                   mismatch_norm=candidate_mismatch.relative_rms,
                   mismatch_rms=candidate_mismatch.rms,
                   mismatch_max_relative=candidate_mismatch.max_relative))
        catch err
            err isa ArgumentError || rethrow()
            skipped_window_count += 1
        end
    end
    if isempty(windows)
        push!(windows,
              (ξ_min=ξ_lo,
               ξ_max=ξ_hi,
               width=width,
               fit_points=fit.fit_points,
               fit_residual_rms=fit.residual_rms,
               fit_coefficients=(slope=fit.slope,
                                 intercept=fit.intercept),
               mismatch_norm=mismatch.relative_rms,
               mismatch_rms=mismatch.rms,
               mismatch_max_relative=mismatch.max_relative))
    end
    truncation = _truncation_order_sensitivity(inner, outer, ξ_lo, ξ_hi,
                                               fit_order_values, min_points,
                                               sample_points, tolerance)
    sensitivity = _overlap_sensitivity(windows, tolerance, truncation)
    classification = sensitivity.stable ? :stable_overlap : :unstable_overlap

    (window=(ξ_min=ξ_lo,
             ξ_max=ξ_hi,
             width=width,
             overlap_ξ_min=ξ_min,
             overlap_ξ_max=ξ_max,
             fit_points=fit.fit_points,
             sample_points=sample_points,
             window_fraction=window_fraction_value,
             window_count=window_count,
             evaluated_window_count=length(windows),
             skipped_window_count=skipped_window_count,
             min_points=min_points,
             fit_orders=fit_order_values),
     mismatch_norm=mismatch.relative_rms,
     mismatch_rms=mismatch.rms,
     mismatch_max_relative=mismatch.max_relative,
     fit_coefficients=(slope=fit.slope, intercept=fit.intercept),
     fit_residual_rms=fit.residual_rms,
     sensitivity=sensitivity,
     truncation_sensitivity=truncation,
     fit_orders=fit_order_values,
     stable=sensitivity.stable,
     classification=classification,
     windows=Tuple(windows),
     source_status=_DEFAULT_SOURCE_STATUS,
     successful=true)
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
