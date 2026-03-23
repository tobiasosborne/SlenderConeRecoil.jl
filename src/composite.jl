# Matching and composite solution: combine inner and outer solutions.
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
end

# ── Overlap extraction ─────────────────────────────────────────────────
"""
    inner_far_field(sol::InnerSolution, ξ_match)

Extract the asymptotic behavior of the inner solution for large ξ.
In the far field, S_inner ~ c₁·ξ (linear growth).
Returns (slope, intercept) from a linear fit to S(ξ) for ξ > ξ_match.
"""
function inner_far_field(sol::InnerSolution, ξ_match::Float64)
    idx = findall(ξ -> ξ > ξ_match, sol.ξ)
    if length(idx) < 3
        # Not enough points; use last available
        idx = max(1, length(sol.ξ)-5):length(sol.ξ)
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
    (slope, intercept)
end

# ── Composite construction ─────────────────────────────────────────────
"""
    composite_solution(inner::InnerSolution, outer::OuterSolution;
                       ξ_grid=nothing, n_points=500)

Construct the additive composite solution:
  S_comp(ξ) = S_inner(ξ) + S_outer(ξ) - S_overlap(ξ)

where S_overlap = slope·ξ + intercept (the common asymptotic form).

Both solutions are interpolated onto a common ξ grid.
"""
function composite_solution(inner::InnerSolution, outer::OuterSolution;
                            ξ_grid::Union{Nothing,Vector{Float64}}=nothing,
                            n_points::Int=500)
    ε = outer.ε

    # Determine overlap region: where both solutions are valid
    ξ_min = max(inner.ξ[1], outer.ξ[1])
    ξ_max = min(inner.ξ[end], outer.ξ[end])

    if ξ_grid === nothing
        ξ_grid = range(ξ_min, ξ_max, length=n_points) |> collect
    end

    # Get inner far-field asymptotic form
    ξ_match = (ξ_min + ξ_max) / 2  # midpoint of overlap
    slope, intercept = inner_far_field(inner, ξ_match)

    # Interpolate solutions onto common grid
    S_inner_interp = _interp(inner.ξ, inner.S, ξ_grid)
    U_inner_interp = _interp(inner.ξ, inner.U, ξ_grid)
    S_outer_interp = [ε * ξ + _interp_scalar(outer.ξ, outer.s₁, ξ) for ξ in ξ_grid]
    U_outer_interp = _interp(outer.ξ, outer.u₁, ξ_grid)

    # Overlap: the common asymptotic form
    S_overlap = [slope * ξ + intercept for ξ in ξ_grid]
    U_overlap = zeros(length(ξ_grid))  # U → 0 in the far field of inner

    # Additive composite
    S_comp = S_inner_interp .+ S_outer_interp .- S_overlap
    U_comp = U_inner_interp .+ U_outer_interp .- U_overlap

    CompositeSolution(ξ_grid, S_comp, U_comp, ε)
end

# ── Simple linear interpolation ────────────────────────────────────────
function _interp(xs::Vector{Float64}, ys::Vector{Float64}, xq::Vector{Float64})
    [_interp_scalar(xs, ys, x) for x in xq]
end

function _interp_scalar(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    if x ≤ xs[1]
        return ys[1]
    elseif x ≥ xs[end]
        return ys[end]
    end
    # Binary search for interval
    lo, hi = 1, length(xs)
    while hi - lo > 1
        mid = div(lo + hi, 2)
        if xs[mid] ≤ x
            lo = mid
        else
            hi = mid
        end
    end
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    ys[lo] + t * (ys[hi] - ys[lo])
end

# ── Overlap residual diagnostic ────────────────────────────────────────
"""
    overlap_residual(inner::InnerSolution, outer::OuterSolution;
                     ξ_start=nothing, ξ_end=nothing)

Compute the mismatch between inner and outer solutions in the overlap region.
Returns the max |S_inner(ξ) - S_outer(ξ)| relative to S_outer.
"""
function overlap_residual(inner::InnerSolution, outer::OuterSolution;
                          ξ_start::Union{Nothing,Float64}=nothing,
                          ξ_end::Union{Nothing,Float64}=nothing)
    ε = outer.ε
    ξ_lo = something(ξ_start, max(inner.ξ[1], outer.ξ[1]))
    ξ_hi = something(ξ_end, min(inner.ξ[end], outer.ξ[end]))

    ξ_test = range(ξ_lo, ξ_hi, length=50) |> collect

    max_res = 0.0
    for ξ in ξ_test
        Si = _interp_scalar(inner.ξ, inner.S, ξ)
        So = ε * ξ + _interp_scalar(outer.ξ, outer.s₁, ξ)
        ref = max(abs(So), 1e-10)
        res = abs(Si - So) / ref
        max_res = max(max_res, res)
    end
    max_res
end
