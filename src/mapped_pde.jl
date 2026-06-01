# Similarity-frame PDE verification utilities.
#
# This is a diagnostic verifier, not a full transformed-PDE solver. It uses a
# prescribed stationary similarity-frame profile on a fixed xi grid, maps it
# into physical coordinates, and reuses the existing PDE conservation and
# collapse diagnostics on a common physical reconstruction window.

export MappedPDESolution, mapped_coordinate_transform, solve_mapped_pde

struct MappedPDESolution
    xi::Vector{Float64}
    t_snapshots::Vector{Float64}
    length_scales::Vector{Float64}
    z::Vector{Vector{Float64}}
    S::Vector{Vector{Float64}}
    U::Vector{Vector{Float64}}
    R::Vector{Vector{Float64}}
    u::Vector{Vector{Float64}}
    epsilon::Float64
    exponent::Float64
    time_offset::Float64
    diagnostics::NamedTuple
end

function _mapped_finite_float(name::AbstractString, value)
    x = Float64(value)
    isfinite(x) || throw(ArgumentError("$name must be finite; got $value"))
    x
end

function _mapped_positive_float(name::AbstractString, value)
    x = _mapped_finite_float(name, value)
    x > 0 || throw(ArgumentError("$name must be positive; got $value"))
    x
end

function _mapped_length_scale(time::Float64; exponent::Float64,
                              time_offset::Float64)
    effective_time = time - time_offset
    effective_time > 0 ||
        throw(ArgumentError("mapped verifier requires positive effective time; got t=$time and time_offset=$time_offset"))
    effective_time^exponent
end

function _mapped_velocity_scale(time::Float64; exponent::Float64,
                                time_offset::Float64)
    effective_time = time - time_offset
    effective_time > 0 ||
        throw(ArgumentError("mapped verifier requires positive effective time; got t=$time and time_offset=$time_offset"))
    exponent * effective_time^(exponent - 1)
end

function _validate_mapped_xi(xi; context::AbstractString)
    xi_values = _as_finite_float_vector(xi; name="$context xi")
    _validate_grid_vector(xi_values; context="$context xi")
    first(xi_values) > 0 ||
        throw(ArgumentError("$context requires xi values to be positive; got xi[1]=$(first(xi_values))"))
    xi_values
end

function _validate_mapped_state(xi_values::Vector{Float64}, S, U;
                                context::AbstractString)
    S_values = _as_finite_float_vector(S; name="$context S")
    U_values = _as_finite_float_vector(U; name="$context U")
    length(S_values) == length(xi_values) ||
        throw(ArgumentError("$context requires S to have length $(length(xi_values)); got $(length(S_values))"))
    length(U_values) == length(xi_values) ||
        throw(ArgumentError("$context requires U to have length $(length(xi_values)); got $(length(U_values))"))
    for i in eachindex(S_values)
        S_values[i] > 0 ||
            throw(DomainError(S_values[i], "$context requires positive similarity profile values; got S[$i]=$(S_values[i])"))
    end
    (S=S_values, U=U_values)
end

"""
    mapped_coordinate_transform(xi, S, U, t; exponent=2/3, time_offset=0)

Map a similarity-frame snapshot to physical coordinates using
`z = ell(t) * xi`, `R = ell(t) * S`, and
`u = ell_dot(t) * U`, where
`ell(t) = (t - time_offset)^exponent` and
`ell_dot(t) = exponent * (t - time_offset)^(exponent - 1)`.
"""
function mapped_coordinate_transform(xi, S, U, t; exponent::Real=2 / 3,
                                     time_offset::Real=0.0)
    exponent_value = _mapped_positive_float("mapped exponent", exponent)
    time_offset_value = _mapped_finite_float("time_offset", time_offset)
    time_value = _mapped_finite_float("mapped time", t)
    xi_values = _validate_mapped_xi(xi; context="mapped coordinate transform")
    state = _validate_mapped_state(xi_values, S, U;
                                   context="mapped coordinate transform")
    ell = _mapped_length_scale(time_value; exponent=exponent_value,
                               time_offset=time_offset_value)
    velocity_scale = _mapped_velocity_scale(time_value;
                                            exponent=exponent_value,
                                            time_offset=time_offset_value)
    (xi=xi_values,
     z=ell .* xi_values,
     S=state.S,
     U=state.U,
     R=ell .* state.S,
     u=velocity_scale .* state.U,
     time=time_value,
     effective_time=time_value - time_offset_value,
     length_scale=ell,
     velocity_scale=velocity_scale,
     exponent=exponent_value,
     time_offset=time_offset_value)
end

function _validate_solve_mapped_pde_args(epsilon::Float64, N::Int,
                                         xi_min::Float64, xi_max::Float64,
                                         effective_time_start::Float64,
                                         effective_time_end::Float64,
                                         n_snapshots::Int,
                                         exponent::Float64,
                                         time_offset::Float64)
    epsilon > 0 && isfinite(epsilon) ||
        throw(ArgumentError("solve_mapped_pde requires positive finite epsilon; got $epsilon"))
    N >= 5 ||
        throw(ArgumentError("solve_mapped_pde requires N >= 5; got N=$N"))
    isfinite(xi_min) && isfinite(xi_max) && 0 < xi_min < xi_max ||
        throw(ArgumentError("solve_mapped_pde requires 0 < xi_min < xi_max; got xi_min=$xi_min, xi_max=$xi_max"))
    effective_time_start > 0 && isfinite(effective_time_start) ||
        throw(ArgumentError("solve_mapped_pde requires positive finite effective_time_start; got $effective_time_start"))
    effective_time_end > effective_time_start && isfinite(effective_time_end) ||
        throw(ArgumentError("solve_mapped_pde requires effective_time_end > effective_time_start; got $effective_time_end <= $effective_time_start"))
    n_snapshots >= 2 ||
        throw(ArgumentError("solve_mapped_pde requires n_snapshots >= 2; got $n_snapshots"))
    exponent > 0 && isfinite(exponent) ||
        throw(ArgumentError("solve_mapped_pde requires positive finite exponent; got $exponent"))
    isfinite(time_offset) ||
        throw(ArgumentError("solve_mapped_pde requires finite time_offset; got $time_offset"))
    nothing
end

function _mapped_reference_profile(xi::Vector{Float64}, epsilon::Float64)
    span = last(xi) - first(xi)
    center = 0.5 * (first(xi) + last(xi))
    width = max(0.25 * span, eps(Float64))
    bump = @. exp(-((xi - center) / width)^2)
    S = @. epsilon * xi * (1 + 0.05 * bump)
    U = @. xi / last(xi) + 0.02 * sin(2pi * (xi - first(xi)) / span)
    (S=S, U=U)
end

function _mapped_stationarity_diagnostics(t_snapshots::Vector{Float64},
                                          S, U)
    reference_S = S[end]
    reference_U = U[end]
    profile_scale = max(sqrt(sum(abs2, reference_S) / length(reference_S)),
                        eps(Float64))
    velocity_scale = max(sqrt(sum(abs2, reference_U) / length(reference_U)),
                         eps(Float64))
    profile_mse = 0.0
    velocity_mse = 0.0
    max_abs_profile = 0.0
    max_abs_velocity = 0.0
    temporal_profile = 0.0
    temporal_velocity = 0.0
    temporal_count = 0

    per_snapshot = NamedTuple[]
    for k in eachindex(S)
        dS = S[k] .- reference_S
        dU = U[k] .- reference_U
        profile_rms = sqrt(sum(abs2, dS) / length(dS))
        velocity_rms = sqrt(sum(abs2, dU) / length(dU))
        profile_mse += sum(abs2, dS) / length(dS)
        velocity_mse += sum(abs2, dU) / length(dU)
        max_abs_profile = max(max_abs_profile, maximum(abs, dS))
        max_abs_velocity = max(max_abs_velocity, maximum(abs, dU))
        push!(per_snapshot,
              (index=k,
               time=t_snapshots[k],
               profile_rms=profile_rms,
               profile_relative_rms=profile_rms / profile_scale,
               velocity_rms=velocity_rms,
               velocity_relative_rms=velocity_rms / velocity_scale))
        if k > 1
            dt = t_snapshots[k] - t_snapshots[k - 1]
            dSdt = (S[k] .- S[k - 1]) ./ dt
            dUdt = (U[k] .- U[k - 1]) ./ dt
            temporal_profile += sum(abs2, dSdt) / length(dSdt)
            temporal_velocity += sum(abs2, dUdt) / length(dUdt)
            temporal_count += 1
        end
    end

    n = length(S)
    profile_rms = sqrt(profile_mse / n)
    velocity_rms = sqrt(velocity_mse / n)
    profile_relative = profile_rms / profile_scale
    velocity_relative = velocity_rms / velocity_scale
    aggregate = sqrt(0.5 * (profile_relative^2 + velocity_relative^2))
    temporal_profile_rms = temporal_count == 0 ? NaN :
                           sqrt(temporal_profile / temporal_count)
    temporal_velocity_rms = temporal_count == 0 ? NaN :
                            sqrt(temporal_velocity / temporal_count)
    (status=:ok,
     successful=isfinite(aggregate),
     basis="temporal stationarity of prescribed similarity-frame S and U",
     profile_rms=profile_rms,
     profile_relative_rms=profile_relative,
     velocity_rms=velocity_rms,
     velocity_relative_rms=velocity_relative,
     aggregate_relative_rms=aggregate,
     max_abs_profile_residual=max_abs_profile,
     max_abs_velocity_residual=max_abs_velocity,
     temporal_profile_rms=temporal_profile_rms,
     temporal_velocity_rms=temporal_velocity_rms,
     per_snapshot=per_snapshot)
end

function _mapped_physical_reconstruction(transforms, n_grid::Int)
    z_min = maximum(first(transform.z) for transform in transforms)
    z_max = minimum(last(transform.z) for transform in transforms)
    if !(isfinite(z_min) && isfinite(z_max) && z_min < z_max)
        throw(ArgumentError("mapped verifier requires an overlapping physical reconstruction window; got z_min=$z_min and z_max=$z_max"))
    end
    z_left = nextfloat(z_min)
    z_right = prevfloat(z_max)
    z_left < z_right ||
        throw(ArgumentError("mapped verifier common physical window is too small after endpoint guard"))
    z_common = collect(range(z_left, z_right, length=n_grid))

    R_common = Vector{Vector{Float64}}()
    u_common = Vector{Vector{Float64}}()
    for transform in transforms
        xi_query = z_common ./ transform.length_scale
        S_interp = _interp_strict(transform.xi, transform.S, xi_query;
                                  context="mapped verifier S reconstruction")
        U_interp = _interp_strict(transform.xi, transform.U, xi_query;
                                  context="mapped verifier U reconstruction")
        push!(R_common, transform.length_scale .* S_interp)
        push!(u_common, transform.velocity_scale .* U_interp)
    end
    (z=z_common, R=R_common, u=u_common)
end

function _mapped_longer_time_reference(transforms, xi_max::Float64,
                                       exponent::Float64,
                                       fixed_reference_z_max)
    reference_z_max = fixed_reference_z_max === nothing ?
                      last(first(transforms).z) :
                      _mapped_positive_float("fixed_reference_z_max",
                                             fixed_reference_z_max)
    limit = (reference_z_max / xi_max)^(1 / exponent)
    final_effective_time = last(transforms).effective_time
    (fixed_reference_z_max=reference_z_max,
     fixed_grid_effective_time_limit=limit,
     final_effective_time=final_effective_time,
     effective_time_gain=final_effective_time - limit,
     effective_time_ratio=final_effective_time / limit,
     reaches_longer_effective_similarity_time=
         final_effective_time > limit * (1 + sqrt(eps(Float64))))
end

function _mapped_collapse_fields(z, t_snapshots, R, u; epsilon::Float64,
                                 exponent::Float64,
                                 time_offset::Float64)
    diagnostic = try
        similarity_collapse_diagnostics(
            z, t_snapshots, R, u; epsilon=epsilon, exponent=exponent,
            time_offset=time_offset, n_grid=max(16, min(128, length(z))))
    catch err
        (status=:failed,
         successful=false,
         message=sprint(showerror, err),
         source_status=_DEFAULT_SOURCE_STATUS,
         aggregate_score=NaN,
         xi_window=(min=NaN, max=NaN, source=:not_available),
         interpolation_grid=(points=0,),
         included_snapshots=NamedTuple[],
         excluded_snapshots=NamedTuple[],
         norms=(profile=(relative_rms=NaN,),
                slope=(relative_rms=NaN,),
                curvature=(relative_rms=NaN,),
                velocity=(relative_rms=NaN,),
                wave_phase=(status=:not_computed,
                            aggregate_relative_rms=NaN)))
    end
    (similarity_collapse=diagnostic,
     similarity_collapse_status=diagnostic.status,
     similarity_collapse_successful=diagnostic.successful,
     similarity_collapse_score=diagnostic.aggregate_score,
     similarity_collapse_xi_min=diagnostic.xi_window.min,
     similarity_collapse_xi_max=diagnostic.xi_window.max,
     similarity_collapse_grid_points=diagnostic.interpolation_grid.points,
     similarity_collapse_included_snapshots=length(diagnostic.included_snapshots),
     similarity_collapse_excluded_snapshots=length(diagnostic.excluded_snapshots),
     similarity_collapse_profile_relative_rms=
         diagnostic.norms.profile.relative_rms,
     similarity_collapse_slope_relative_rms=
         diagnostic.norms.slope.relative_rms,
     similarity_collapse_curvature_relative_rms=
         diagnostic.norms.curvature.relative_rms,
     similarity_collapse_velocity_relative_rms=
         diagnostic.norms.velocity.relative_rms,
     similarity_collapse_wave_phase_status=
         diagnostic.norms.wave_phase.status,
     similarity_collapse_wave_phase_relative_rms=
         diagnostic.norms.wave_phase.aggregate_relative_rms)
end

function _mapped_solution_diagnostics(xi::Vector{Float64},
                                      t_snapshots::Vector{Float64},
                                      transforms, S, U,
                                      reconstruction; epsilon::Float64,
                                      exponent::Float64,
                                      time_offset::Float64,
                                      fixed_reference_z_max)
    conservation = pde_conservation_diagnostics(
        reconstruction.z, t_snapshots, reconstruction.R, reconstruction.u;
        retcode=:Success, endpoint=last(t_snapshots),
        requested_endpoint=last(t_snapshots), saved_points=length(t_snapshots))
    collapse = _mapped_collapse_fields(
        reconstruction.z, t_snapshots, reconstruction.R, reconstruction.u;
        epsilon=epsilon, exponent=exponent, time_offset=time_offset)
    stationarity = _mapped_stationarity_diagnostics(t_snapshots, S, U)
    reference = _mapped_longer_time_reference(
        transforms, last(xi), exponent, fixed_reference_z_max)
    physical_ranges = [(time=transform.time,
                        effective_time=transform.effective_time,
                        length_scale=transform.length_scale,
                        velocity_scale=transform.velocity_scale,
                        z_min=first(transform.z),
                        z_max=last(transform.z))
                       for transform in transforms]
    successful = conservation.pde_data_valid &&
                 conservation.finite_state &&
                 conservation.positive_radius &&
                 conservation.grid_strictly_increasing &&
                 conservation.time_strictly_increasing &&
                 conservation.retcode_successful &&
                 conservation.endpoint_reached &&
                 collapse.similarity_collapse_successful &&
                 stationarity.successful &&
                 reference.reaches_longer_effective_similarity_time
    base = (context="mapped-coordinate PDE verifier",
            retcode=:Success,
            successful=successful,
            endpoint=last(t_snapshots),
            requested_endpoint=last(t_snapshots),
            saved_points=length(t_snapshots),
            source_status=_DEFAULT_SOURCE_STATUS,
            diagnostic_basis="prescribed similarity-frame verifier; not a full transformed PDE solve and not Decent-King benchmark data",
            coordinate_system=:similarity_mapped,
            coordinate_transform="z = ell(t) * xi; R = ell(t) * S; u = ell_dot(t) * U",
            length_scale_convention="ell(t) = (t - time_offset)^exponent",
            velocity_scaling_convention="ell_dot(t) = exponent * (t - time_offset)^(exponent - 1)",
            exponent=exponent,
            time_offset=time_offset,
            xi_domain=(min=first(xi), max=last(xi), points=length(xi)),
            physical_z_ranges=physical_ranges,
            length_scales=[transform.length_scale for transform in transforms],
            velocity_scales=[transform.velocity_scale for transform in transforms],
            common_physical_grid=(z_min=first(reconstruction.z),
                                  z_max=last(reconstruction.z),
                                  points=length(reconstruction.z)),
            fixed_grid_reference=reference,
            reaches_longer_effective_similarity_time=
                reference.reaches_longer_effective_similarity_time,
            mapped_stationarity=stationarity,
            mapped_stationarity_residual_norm=
                stationarity.aggregate_relative_rms,
            residual_norm=stationarity.aggregate_relative_rms,
            final_residual_norm=stationarity.aggregate_relative_rms,
            conservation_diagnostics=conservation)
    merge(base, conservation, collapse, (successful=successful,))
end

"""
    solve_mapped_pde(; kwargs...)

Construct a cheap similarity-frame PDE verifier on a fixed `xi` grid.

The current implementation prescribes a stationary similarity-frame profile,
maps snapshots to physical coordinates, and evaluates conservation,
positivity, and collapse diagnostics on their common physical reconstruction
window. It is intended as an implementation-level mapped-coordinate verifier,
not as a source-backed transformed-PDE benchmark.
"""
function solve_mapped_pde(; ε::Real=0.1,
                          epsilon::Union{Nothing,Real}=nothing,
                          N::Int=64,
                          xi_min::Real=1.0, xi_max::Real=8.0,
                          effective_time_start::Real=0.25,
                          effective_time_end::Real=1.0,
                          n_snapshots::Int=5,
                          exponent::Real=2 / 3,
                          time_offset::Real=0.0,
                          fixed_reference_z_max=nothing)
    epsilon_value = Float64(something(epsilon, ε))
    xi_min_value = Float64(xi_min)
    xi_max_value = Float64(xi_max)
    tau_start = Float64(effective_time_start)
    tau_end = Float64(effective_time_end)
    exponent_value = Float64(exponent)
    time_offset_value = Float64(time_offset)
    _validate_solve_mapped_pde_args(
        epsilon_value, N, xi_min_value, xi_max_value, tau_start, tau_end,
        n_snapshots, exponent_value, time_offset_value)

    xi = collect(range(xi_min_value, xi_max_value, length=N))
    profile = _mapped_reference_profile(xi, epsilon_value)
    t_snapshots = time_offset_value .+
                  collect(range(tau_start, tau_end, length=n_snapshots))
    S = [copy(profile.S) for _ in t_snapshots]
    U = [copy(profile.U) for _ in t_snapshots]
    transforms = [mapped_coordinate_transform(
                      xi, S[k], U[k], t_snapshots[k];
                      exponent=exponent_value,
                      time_offset=time_offset_value)
                  for k in eachindex(t_snapshots)]
    reconstruction = _mapped_physical_reconstruction(transforms, N)
    diagnostics = _mapped_solution_diagnostics(
        xi, t_snapshots, transforms, S, U, reconstruction;
        epsilon=epsilon_value, exponent=exponent_value,
        time_offset=time_offset_value,
        fixed_reference_z_max=fixed_reference_z_max)
    MappedPDESolution(
        xi,
        t_snapshots,
        [transform.length_scale for transform in transforms],
        [transform.z for transform in transforms],
        S,
        U,
        [transform.R for transform in transforms],
        [transform.u for transform in transforms],
        epsilon_value,
        exponent_value,
        time_offset_value,
        diagnostics)
end

function _mapped_nested_finite(vectors)
    all(vector -> _finite_vector(vector), vectors)
end

function _mapped_nested_positive(vectors)
    all(vector -> _finite_vector(vector) && all(>(0), vector), vectors)
end

function diagnostics_succeeded(sol::MappedPDESolution)
    diagnostics_succeeded(sol.diagnostics) || return false
    _finite_vector(sol.xi) && all(>(0), diff(sol.xi)) || return false
    _finite_vector(sol.t_snapshots) && all(>(0), diff(sol.t_snapshots)) ||
        return false
    _finite_vector(sol.length_scales) && all(>(0), sol.length_scales) ||
        return false
    _mapped_nested_finite(sol.z) || return false
    _mapped_nested_positive(sol.S) || return false
    _mapped_nested_finite(sol.U) || return false
    _mapped_nested_positive(sol.R) || return false
    _mapped_nested_finite(sol.u) || return false
    diagnostics = sol.diagnostics
    get(diagnostics, :reaches_longer_effective_similarity_time, false) ||
        return false
    isfinite(get(diagnostics, :mapped_stationarity_residual_norm, NaN)) ||
        return false
    get(diagnostics, :similarity_collapse_successful, false) || return false
    true
end

mesh_summary(sol::MappedPDESolution) = mesh_summary(sol.xi; variable=:xi)

function domain_summary(sol::MappedPDESolution)
    time_fields = _time_domain_summary(sol.t_snapshots)
    z_min = minimum(first(snapshot) for snapshot in sol.z)
    z_max = maximum(last(snapshot) for snapshot in sol.z)
    merge((domain_variable=:xi,
           domain_start=first(sol.xi),
           domain_end=last(sol.xi),
           domain_span=last(sol.xi) - first(sol.xi),
           physical_z_min=z_min,
           physical_z_max=z_max,
           physical_z_span=z_max - z_min),
          time_fields)
end

function _mapped_diagnostic_checks(sol::MappedPDESolution)
    diagnostics = sol.diagnostics
    (coordinate_system=get(diagnostics, :coordinate_system,
                           :similarity_mapped),
     exponent=get(diagnostics, :exponent, sol.exponent),
     time_offset=get(diagnostics, :time_offset, sol.time_offset),
     reaches_longer_effective_similarity_time=
         get(diagnostics, :reaches_longer_effective_similarity_time, false),
     fixed_grid_effective_time_limit=
         get(get(diagnostics, :fixed_grid_reference, (;)),
             :fixed_grid_effective_time_limit, NaN),
     final_effective_time=
         get(get(diagnostics, :fixed_grid_reference, (;)),
             :final_effective_time, NaN),
     mapped_stationarity_residual_norm=
         get(diagnostics, :mapped_stationarity_residual_norm, NaN),
     minimum_radius=get(diagnostics, :minimum_radius, NaN),
     radius_positivity_margin=
         get(diagnostics, :radius_positivity_margin, NaN),
     max_abs_area_mass_balance_residual=
         get(diagnostics, :max_abs_area_mass_balance_residual, NaN),
     similarity_collapse_status=
         get(diagnostics, :similarity_collapse_status, :not_computed),
     similarity_collapse_successful=
         get(diagnostics, :similarity_collapse_successful, false),
     similarity_collapse_score=
         get(diagnostics, :similarity_collapse_score, NaN))
end

function diagnostic_summary(sol::MappedPDESolution;
                            problem_kind::Symbol=:mapped_pde_verification)
    rn = residual_norm(sol)
    DiagnosticSummary(problem_kind;
                      successful=diagnostics_succeeded(sol),
                      residual_norm=rn,
                      final_residual_norm=rn,
                      retcode=get(sol.diagnostics, :retcode, nothing),
                      mesh=mesh_summary(sol),
                      domain=domain_summary(sol),
                      checks=_mapped_diagnostic_checks(sol))
end
