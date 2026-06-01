# Independent finite-difference checks for the time-dependent PDE operator.
#
# This backend intentionally does not call ddz! or _ddz_unchecked!. It builds
# local polynomial finite-difference weights directly and compares the resulting
# primitive RHS against the package method-of-lines implementation.

export independent_pde_rhs, independent_pde_residual_diagnostics,
       pde_discretization_comparison

const _INDEPENDENT_PDE_BACKEND_NAME = "local-polynomial-finite-difference"
const _INDEPENDENT_PDE_BACKEND_DESCRIPTION =
    "direct local polynomial finite differences assembled from Vandermonde moment equations; does not call ddz! or _ddz_unchecked!"

function _independent_pde_stencil_width(width::Int)
    width >= 5 ||
        throw(ArgumentError("independent PDE verifier requires stencil_width >= 5 for third derivatives; got $width"))
    isodd(width) ||
        throw(ArgumentError("independent PDE verifier requires an odd stencil_width; got $width"))
    width
end

function _validate_independent_pde_state(z, R, u; context::AbstractString)
    z_values = _as_finite_float_vector(z; name="$context grid")
    N = _validate_grid_vector(z_values; context="$context grid")
    R_values = _as_finite_float_vector(R; name="$context radius")
    u_values = _as_finite_float_vector(u; name="$context velocity")
    length(R_values) == N ||
        throw(ArgumentError("$context requires radius length $N; got $(length(R_values))"))
    length(u_values) == N ||
        throw(ArgumentError("$context requires velocity length $N; got $(length(u_values))"))
    for i in 1:N
        R_values[i] > 0 ||
            throw(DomainError(R_values[i], "$context requires positive radius values; got R[$i]=$(R_values[i])"))
    end
    (z=z_values, R=R_values, u=u_values)
end

function _independent_float_time(t; context::AbstractString)
    value = Float64(t)
    isfinite(value) ||
        throw(ArgumentError("$context requires a finite time; got $t"))
    value
end

function _local_stencil_indices(i::Int, N::Int, width::Int)
    width <= N ||
        throw(ArgumentError("finite-difference stencil_width=$width exceeds grid points N=$N"))
    half_width = fld(width, 2)
    first_index = clamp(i - half_width, 1, N - width + 1)
    first_index:(first_index + width - 1)
end

function _solve_small_dense(A::Matrix{Float64}, b::Vector{Float64})
    n = length(b)
    size(A, 1) == n && size(A, 2) == n ||
        throw(ArgumentError("small dense solve requires a square matrix matching b"))
    M = copy(A)
    rhs = copy(b)

    for col in 1:n
        pivot = col
        pivot_abs = abs(M[col, col])
        for row in (col + 1):n
            candidate = abs(M[row, col])
            if candidate > pivot_abs
                pivot = row
                pivot_abs = candidate
            end
        end
        pivot_abs > eps(Float64) ||
            throw(ArgumentError("finite-difference stencil produced a singular moment system"))
        if pivot != col
            for j in col:n
                M[col, j], M[pivot, j] = M[pivot, j], M[col, j]
            end
            rhs[col], rhs[pivot] = rhs[pivot], rhs[col]
        end
        for row in (col + 1):n
            factor = M[row, col] / M[col, col]
            M[row, col] = 0.0
            for j in (col + 1):n
                M[row, j] -= factor * M[col, j]
            end
            rhs[row] -= factor * rhs[col]
        end
    end

    x = zeros(n)
    for row in n:-1:1
        total = rhs[row]
        for j in (row + 1):n
            total -= M[row, j] * x[j]
        end
        x[row] = total / M[row, row]
    end
    x
end

function _finite_difference_weights(x0::Float64, nodes::AbstractVector{Float64},
                                    derivative_order::Int)
    derivative_order >= 0 ||
        throw(ArgumentError("derivative_order must be nonnegative; got $derivative_order"))
    n = length(nodes)
    derivative_order < n ||
        throw(ArgumentError("derivative_order=$derivative_order requires more than $derivative_order stencil nodes; got $n"))
    A = zeros(n, n)
    b = zeros(n)
    offsets = Float64.(nodes .- x0)
    for power in 0:(n - 1)
        row = power + 1
        for j in 1:n
            A[row, j] = offsets[j]^power
        end
        b[row] = power == derivative_order ? factorial(derivative_order) : 0.0
    end
    _solve_small_dense(A, b)
end

function _independent_derivative(values::Vector{Float64},
                                 z::Vector{Float64},
                                 derivative_order::Int;
                                 stencil_width::Int)
    N = length(z)
    width = _independent_pde_stencil_width(stencil_width)
    width <= N ||
        throw(ArgumentError("independent derivative requires at least stencil_width=$width grid points; got N=$N"))
    derivative = zeros(N)
    for i in 1:N
        indices = _local_stencil_indices(i, N, width)
        weights = _finite_difference_weights(z[i], @view(z[indices]),
                                             derivative_order)
        value = 0.0
        k = 1
        for j in indices
            value += weights[k] * values[j]
            k += 1
        end
        derivative[i] = value
    end
    derivative
end

"""
    independent_pde_rhs(z, R, u; time=0, stencil_width=5)

Evaluate the primitive slender-cone PDE RHS with an independent finite-
difference backend. The backend assembles local polynomial weights directly
from moment equations on the supplied nonuniform grid; it does not call
`ddz!` or `_ddz_unchecked!`.

The returned vector is ordered like `pde_rhs!`: first `dR/dt`, then `du/dt`.
This is an implementation-level check of the package PDE operator, not a
Decent-King benchmark value.
"""
function independent_pde_rhs(z, R, u; time::Real=0,
                             stencil_width::Int=5)
    data = _validate_independent_pde_state(z, R, u;
                                           context="independent_pde_rhs")
    _independent_float_time(time; context="independent_pde_rhs")
    width = _independent_pde_stencil_width(stencil_width)
    N = length(data.z)
    width <= N ||
        throw(ArgumentError("independent_pde_rhs requires at least stencil_width=$width grid points; got N=$N"))

    Rz = _independent_derivative(data.R, data.z, 1;
                                 stencil_width=width)
    uz = _independent_derivative(data.u, data.z, 1;
                                 stencil_width=width)
    invR = 1.0 ./ data.R
    invR_z = _independent_derivative(invR, data.z, 1;
                                     stencil_width=width)
    Rzzz = _independent_derivative(data.R, data.z, 3;
                                   stencil_width=width)

    dR = zeros(N)
    du = zeros(N)
    for i in 1:N
        dR[i] = -data.u[i] * Rz[i] - 0.5 * data.R[i] * uz[i]
        du[i] = -data.u[i] * uz[i] - invR_z[i] + Rzzz[i]
    end

    dR[1] = 0.0
    du[1] = 0.0
    dR[N] = dR[N - 1]
    du[N] = du[N - 1]
    vcat(dR, du)
end

function _rms(values::AbstractVector{Float64})
    isempty(values) && return NaN
    sqrt(sum(abs2, values) / length(values))
end

function _max_abs(values::AbstractVector{Float64})
    isempty(values) && return NaN
    maximum(abs, values)
end

function _independent_margin(margin::Int, N::Int)
    margin >= 0 ||
        throw(ArgumentError("comparison margin must be nonnegative; got $margin"))
    2margin < N ||
        throw(ArgumentError("comparison margin=$margin leaves no interior points for N=$N"))
    margin
end

function _component_error_fields(package_values, independent_values, indices,
                                 atol::Float64, rtol::Float64)
    diff_values = independent_values[indices] .- package_values[indices]
    reference_max = _max_abs(package_values[indices])
    reference_rms = _rms(package_values[indices])
    max_abs_error = _max_abs(diff_values)
    rms_error = _rms(diff_values)
    denominator_max = max(reference_max, eps(Float64))
    denominator_rms = max(reference_rms, eps(Float64))
    tolerance = atol + rtol * reference_max
    (max_abs_error=max_abs_error,
     rms_error=rms_error,
     reference_max_abs=reference_max,
     reference_rms=reference_rms,
     relative_max_abs_error=max_abs_error / denominator_max,
     relative_rms_error=rms_error / denominator_rms,
     tolerance=tolerance,
     agrees=max_abs_error <= tolerance)
end

function _validate_comparison_tolerances(atol, rtol)
    atol_value = Float64(atol)
    rtol_value = Float64(rtol)
    isfinite(atol_value) && atol_value >= 0 ||
        throw(ArgumentError("absolute comparison tolerance must be finite and nonnegative; got $atol"))
    isfinite(rtol_value) && rtol_value >= 0 ||
        throw(ArgumentError("relative comparison tolerance must be finite and nonnegative; got $rtol"))
    (atol=atol_value, rtol=rtol_value)
end

function _expected_rhs_from_parts(expected_Rt, expected_ut, N::Int)
    if expected_Rt === nothing && expected_ut === nothing
        return nothing
    end
    expected_Rt !== nothing && expected_ut !== nothing ||
        throw(ArgumentError("expected_Rt and expected_ut must be supplied together"))
    Rt = _as_finite_float_vector(expected_Rt; name="expected_Rt")
    ut = _as_finite_float_vector(expected_ut; name="expected_ut")
    length(Rt) == N ||
        throw(ArgumentError("expected_Rt length must be $N; got $(length(Rt))"))
    length(ut) == N ||
        throw(ArgumentError("expected_ut length must be $N; got $(length(ut))"))
    vcat(Rt, ut)
end

function _expected_rhs_vector(expected_rhs, expected_Rt, expected_ut, N::Int)
    parts = _expected_rhs_from_parts(expected_Rt, expected_ut, N)
    if expected_rhs === nothing
        return parts
    end
    parts === nothing ||
        throw(ArgumentError("provide either expected_rhs or expected_Rt/expected_ut, not both"))
    if expected_rhs isa NamedTuple
        haskey(expected_rhs, :Rt) && haskey(expected_rhs, :ut) ||
            throw(ArgumentError("expected_rhs NamedTuple requires Rt and ut fields"))
        return _expected_rhs_from_parts(expected_rhs.Rt, expected_rhs.ut, N)
    end
    values = _as_finite_float_vector(expected_rhs; name="expected_rhs")
    length(values) == 2N ||
        throw(ArgumentError("expected_rhs length must be 2N=$((2N)); got $(length(values))"))
    values
end

function _mms_residual_fields(package_rhs::Vector{Float64},
                              independent_rhs::Vector{Float64},
                              expected_rhs,
                              radius_indices,
                              velocity_indices)
    expected_rhs === nothing &&
        return (mms_expected_provided=false,)
    package_radius = package_rhs[radius_indices] .- expected_rhs[radius_indices]
    package_velocity = package_rhs[velocity_indices] .- expected_rhs[velocity_indices]
    independent_radius =
        independent_rhs[radius_indices] .- expected_rhs[radius_indices]
    independent_velocity =
        independent_rhs[velocity_indices] .- expected_rhs[velocity_indices]
    (mms_expected_provided=true,
     package_mms_radius_max_abs_residual=_max_abs(package_radius),
     package_mms_radius_rms_residual=_rms(package_radius),
     package_mms_velocity_max_abs_residual=_max_abs(package_velocity),
     package_mms_velocity_rms_residual=_rms(package_velocity),
     independent_mms_radius_max_abs_residual=_max_abs(independent_radius),
     independent_mms_radius_rms_residual=_rms(independent_radius),
     independent_mms_velocity_max_abs_residual=_max_abs(independent_velocity),
     independent_mms_velocity_rms_residual=_rms(independent_velocity))
end

"""
    pde_discretization_comparison(z, R, u; kwargs...)

Compare `pde_rhs!` against `independent_pde_rhs` on one snapshot. The default
comparison excludes four points at each boundary so the diagnostic measures
interior operator agreement rather than the deliberately imposed boundary
closures. The default tolerances are implementation-level cheap-test
tolerances: `atol=5e-3`, `rtol=5e-2`.

Optional manufactured-solution residuals can be recorded with either
`expected_rhs=vcat(Rt, ut)` or `expected_Rt=Rt, expected_ut=ut`.
"""
function pde_discretization_comparison(z, R, u; time::Real=0,
                                       margin::Int=4,
                                       atol::Real=5e-3,
                                       rtol::Real=5e-2,
                                       stencil_width::Int=5,
                                       expected_rhs=nothing,
                                       expected_Rt=nothing,
                                       expected_ut=nothing)
    data = _validate_independent_pde_state(z, R, u;
                                           context="pde_discretization_comparison")
    time_value = _independent_float_time(time;
                                         context="pde_discretization_comparison")
    N = length(data.z)
    width = _independent_pde_stencil_width(stencil_width)
    width <= N ||
        throw(ArgumentError("pde_discretization_comparison requires at least stencil_width=$width grid points; got N=$N"))
    margin_value = _independent_margin(margin, N)
    tolerances = _validate_comparison_tolerances(atol, rtol)
    expected = _expected_rhs_vector(expected_rhs, expected_Rt, expected_ut, N)

    state = vcat(data.R, data.u)
    package_rhs = zeros(2N)
    pde_rhs!(package_rhs, state, _make_pde_rhs_parameters(data.z), time_value)
    independent_rhs = independent_pde_rhs(
        data.z, data.R, data.u; time=time_value, stencil_width=width)

    radius_indices = (margin_value + 1):(N - margin_value)
    velocity_indices = (N + margin_value + 1):(2N - margin_value)
    radius_errors = _component_error_fields(
        package_rhs, independent_rhs, radius_indices,
        tolerances.atol, tolerances.rtol)
    velocity_errors = _component_error_fields(
        package_rhs, independent_rhs, velocity_indices,
        tolerances.atol, tolerances.rtol)
    agrees = radius_errors.agrees && velocity_errors.agrees
    grid_stats = _spacing_stats(data.z)
    mms_fields = _mms_residual_fields(package_rhs, independent_rhs, expected,
                                      radius_indices, velocity_indices)
    residual_norm_value = max(radius_errors.rms_error,
                              velocity_errors.rms_error)

    merge(
        (status=agrees ? :ok : :disagreement,
         successful=agrees,
         source_status=_DEFAULT_SOURCE_STATUS,
         diagnostic_basis="implementation-level independent finite-difference comparison; not Decent-King benchmark data",
         backend_name=_INDEPENDENT_PDE_BACKEND_NAME,
         backend_description=_INDEPENDENT_PDE_BACKEND_DESCRIPTION,
         package_backend="pde_rhs! with package finite-difference cache",
         calls_package_derivative_operator=false,
         time=time_value,
         residual_norm=residual_norm_value,
         final_residual_norm=residual_norm_value,
         tolerances=(absolute=tolerances.atol,
                     relative=tolerances.rtol,
                     basis="cheap smooth-state agreement between independent local polynomial FD and package pde_rhs!"),
         grid=(variable=:z,
               points=N,
               spacing_min=grid_stats.min,
               spacing_max=grid_stats.max,
               spacing_ratio=grid_stats.ratio,
               strictly_increasing=grid_stats.strictly_increasing),
         comparison_window=(margin=margin_value,
                            first_index=first(radius_indices),
                            last_index=last(radius_indices),
                            points=length(radius_indices)),
         independent_stencil=(method=:vandermonde_moment_weights,
                              width=width,
                              derivative_orders=(1, 3)),
         radius_max_abs_error=radius_errors.max_abs_error,
         radius_rms_error=radius_errors.rms_error,
         radius_reference_max_abs=radius_errors.reference_max_abs,
         radius_reference_rms=radius_errors.reference_rms,
         radius_relative_max_abs_error=
             radius_errors.relative_max_abs_error,
         radius_relative_rms_error=radius_errors.relative_rms_error,
         radius_tolerance=radius_errors.tolerance,
         radius_agrees=radius_errors.agrees,
         velocity_max_abs_error=velocity_errors.max_abs_error,
         velocity_rms_error=velocity_errors.rms_error,
         velocity_reference_max_abs=velocity_errors.reference_max_abs,
         velocity_reference_rms=velocity_errors.reference_rms,
         velocity_relative_max_abs_error=
             velocity_errors.relative_max_abs_error,
         velocity_relative_rms_error=velocity_errors.relative_rms_error,
         velocity_tolerance=velocity_errors.tolerance,
         velocity_agrees=velocity_errors.agrees,
         max_abs_component_error=
             max(radius_errors.max_abs_error,
                 velocity_errors.max_abs_error),
         rms_component_error=residual_norm_value),
        mms_fields)
end

function pde_discretization_comparison(sol::PDESolution,
                                       snapshot::Integer=lastindex(sol.R);
                                       kwargs...)
    index = Int(snapshot)
    checkbounds(sol.R, index)
    checkbounds(sol.u, index)
    checkbounds(sol.t_snapshots, index)
    pde_discretization_comparison(sol.z, sol.R[index], sol.u[index];
                                  time=sol.t_snapshots[index], kwargs...)
end

"""
    independent_pde_residual_diagnostics(z, R, u; kwargs...)

Alias for `pde_discretization_comparison` that emphasizes manufactured-
solution or reference-RHS residual checks. Use `expected_rhs` or
`expected_Rt`/`expected_ut` to include MMS residual fields.
"""
independent_pde_residual_diagnostics(z, R, u; kwargs...) =
    pde_discretization_comparison(z, R, u; kwargs...)

independent_pde_residual_diagnostics(sol::PDESolution,
                                     snapshot::Integer=lastindex(sol.R);
                                     kwargs...) =
    pde_discretization_comparison(sol, snapshot; kwargs...)
