# Quantitative similarity-collapse diagnostics for PDE snapshots.
#
# These metrics are implementation-level regression diagnostics. They check
# whether rescaled PDE snapshots agree on a common similarity window; they are
# not Decent-King benchmark values unless a caller supplies source-backed data.

export similarity_collapse_diagnostics

const _DEFAULT_COLLAPSE_COMPONENT_WEIGHTS =
    (profile=1.0, slope=0.5, curvature=0.25, velocity=0.5, wave_phase=0.0)

function _collapse_finite_float(name::AbstractString, value)
    x = Float64(value)
    isfinite(x) || throw(ArgumentError("$name must be finite; got $value"))
    x
end

function _collapse_positive_float(name::AbstractString, value)
    x = _collapse_finite_float(name, value)
    x > 0 || throw(ArgumentError("$name must be positive; got $value"))
    x
end

function _collapse_nonnegative_float(name::AbstractString, value)
    x = _collapse_finite_float(name, value)
    x >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    x
end

function _collapse_component_weights(weights)
    merged = merge(_DEFAULT_COLLAPSE_COMPONENT_WEIGHTS, weights)
    for name in keys(_DEFAULT_COLLAPSE_COMPONENT_WEIGHTS)
        value = _collapse_nonnegative_float("collapse component weight $name",
                                            getfield(merged, name))
        merged = merge(merged, NamedTuple{(name,)}((value,)))
    end
    sum(getfield(merged, name) for name in keys(_DEFAULT_COLLAPSE_COMPONENT_WEIGHTS)) > 0 ||
        throw(ArgumentError("collapse component weights must contain at least one positive weight"))
    merged
end

function _collapse_window_values(xi_window)
    xi_window === nothing && return nothing
    if xi_window isa NamedTuple
        lo = if haskey(xi_window, :min)
            xi_window.min
        elseif haskey(xi_window, :xi_min)
            xi_window.xi_min
        else
            throw(ArgumentError("xi_window NamedTuple requires min/xi_min"))
        end
        hi = if haskey(xi_window, :max)
            xi_window.max
        elseif haskey(xi_window, :xi_max)
            xi_window.xi_max
        else
            throw(ArgumentError("xi_window NamedTuple requires max/xi_max"))
        end
        return (_collapse_finite_float("xi_window minimum", lo),
                _collapse_finite_float("xi_window maximum", hi))
    end
    values = Float64.(collect(xi_window))
    length(values) == 2 ||
        throw(ArgumentError("xi_window requires exactly two values; got $(length(values))"))
    all(isfinite, values) ||
        throw(ArgumentError("xi_window requires finite values"))
    (values[1], values[2])
end

function _collapse_failure(status::Symbol, message::AbstractString; exponent,
                           time_offset, n_grid, component_weights,
                           xi_window=nothing,
                           xi_window_source=:not_available,
                           xi_grid=Float64[],
                           included_snapshots=NamedTuple[],
                           excluded_snapshots=NamedTuple[],
                           excluded_regions=NamedTuple[],
                           snapshot_weights_raw=Float64[],
                           snapshot_weights=Float64[],
                           wave_phase=(status=:not_computed,
                                       successful=false,
                                       aggregate_rms=NaN,
                                       aggregate_relative_rms=NaN,
                                       max_abs=NaN,
                                       per_snapshot=NamedTuple[]))
    window = xi_window === nothing ?
             (min=NaN, max=NaN, source=:not_available) :
             (min=xi_window[1], max=xi_window[2],
              source=xi_window_source)
    (status=status,
     successful=false,
     message=String(message),
     source_status=_DEFAULT_SOURCE_STATUS,
     diagnostic_basis="implementation similarity-collapse metric; not Decent-King benchmark data",
     exponent=Float64(exponent),
     time_offset=Float64(time_offset),
     length_scale_convention="ell(t) = (t - time_offset)^exponent",
     velocity_scaling_convention="U = u / (exponent * (t - time_offset)^(exponent - 1))",
     xi_window=window,
     xi_grid=xi_grid,
     interpolation_grid=(variable=:xi, points=length(xi_grid),
                         min=isempty(xi_grid) ? NaN : first(xi_grid),
                         max=isempty(xi_grid) ? NaN : last(xi_grid)),
     component_weights=component_weights,
     snapshot_weights=snapshot_weights,
     snapshot_weights_raw=snapshot_weights_raw,
     included_snapshots=included_snapshots,
     excluded_snapshots=excluded_snapshots,
     excluded_regions=excluded_regions,
     reference_snapshot_index=missing,
     reference_time=NaN,
     reference_included_position=missing,
     norms=(profile=(status=:not_computed, rms=NaN,
                     relative_rms=NaN, max_abs=NaN,
                     reference_rms=NaN, per_snapshot=NamedTuple[]),
            slope=(status=:not_computed, rms=NaN,
                   relative_rms=NaN, max_abs=NaN,
                   reference_rms=NaN, per_snapshot=NamedTuple[]),
            curvature=(status=:not_computed, rms=NaN,
                       relative_rms=NaN, max_abs=NaN,
                       reference_rms=NaN, per_snapshot=NamedTuple[]),
            velocity=(status=:not_computed, rms=NaN,
                      relative_rms=NaN, max_abs=NaN,
                      reference_rms=NaN, per_snapshot=NamedTuple[]),
            wave_phase=wave_phase),
     per_snapshot_scores=NamedTuple[],
     aggregate_score=NaN,
     aggregate_components=NamedTuple[])
end

function _collapse_derivatives(xi::Vector{Float64}, S::Vector{Float64})
    slope = zeros(length(S))
    curvature = zeros(length(S))
    ddz!(slope, S, xi)
    ddz!(curvature, slope, xi)
    (slope=slope, curvature=curvature)
end

function _rescaled_collapse_snapshots(data; exponent::Float64,
                                      time_offset::Float64,
                                      edge_fraction::Float64)
    included = NamedTuple[]
    excluded = NamedTuple[]
    for k in eachindex(data.t)
        t = data.t[k]
        effective_time = t - time_offset
        if !(isfinite(effective_time) && effective_time > 0)
            push!(excluded, (index=k, time=t, effective_time=effective_time,
                             reason=:nonpositive_effective_time))
            continue
        end

        ell = effective_time^exponent
        velocity_scale = exponent * effective_time^(exponent - 1)
        if !(isfinite(ell) && ell > 0)
            push!(excluded, (index=k, time=t, effective_time=effective_time,
                             reason=:invalid_length_scale))
            continue
        elseif !(isfinite(velocity_scale) && velocity_scale != 0)
            push!(excluded, (index=k, time=t, effective_time=effective_time,
                             reason=:invalid_velocity_scale))
            continue
        end

        xi = data.z ./ ell
        S = data.R[k] ./ ell
        U = data.u[k] ./ velocity_scale
        _validate_interpolation_table(xi, S; context="similarity collapse profile snapshot $k")
        _validate_interpolation_table(xi, U; context="similarity collapse velocity snapshot $k")
        derivatives = _collapse_derivatives(xi, S)
        span = last(xi) - first(xi)
        trim = edge_fraction * span
        trusted_min = first(xi) + trim
        trusted_max = last(xi) - trim
        if !(trusted_min < trusted_max)
            push!(excluded, (index=k, time=t, effective_time=effective_time,
                             reason=:empty_trusted_window_after_edge_trim))
            continue
        end
        push!(included, (index=k, time=t, effective_time=effective_time,
                         ell=ell, velocity_scale=velocity_scale,
                         xi=xi, S=S, U=U,
                         slope=derivatives.slope,
                         curvature=derivatives.curvature,
                         trusted_min=trusted_min,
                         trusted_max=trusted_max))
    end
    (included=included, excluded=excluded)
end

function _collapse_reference_position(reference, included)
    reference === :latest &&
        return argmax([snapshot.time for snapshot in included])
    reference === :first && return 1
    reference isa Integer || throw(ArgumentError("reference must be :latest, :first, or an included snapshot index"))
    for i in eachindex(included)
        included[i].index == Int(reference) && return i
    end
    throw(ArgumentError("reference snapshot index $reference was not included in the collapse window"))
end

function _collapse_snapshot_weights(weights, included_count::Int)
    if weights === nothing
        raw = fill(1.0, included_count)
    else
        raw = Float64.(collect(weights))
        length(raw) == included_count ||
            throw(ArgumentError("snapshot_weights length must match included snapshots; got $(length(raw)) and $included_count"))
        all(isfinite, raw) ||
            throw(ArgumentError("snapshot_weights requires finite values"))
        all(>=(0), raw) ||
            throw(ArgumentError("snapshot_weights requires nonnegative values"))
        sum(raw) > 0 ||
            throw(ArgumentError("snapshot_weights requires at least one positive value"))
    end
    total = sum(raw)
    (raw=raw, normalized=raw ./ total)
end

function _collapse_component_norm(name::Symbol, values, included,
                                  reference_position::Int,
                                  snapshot_weights::Vector{Float64})
    ref = values[reference_position]
    n = length(ref)
    reference_rms = sqrt(sum(abs2, ref) / n)
    denom = max(reference_rms, 1.0)
    per = NamedTuple[]
    weighted_mse = 0.0
    max_abs_error = 0.0
    for i in eachindex(values)
        diff = values[i] .- ref
        mse = sum(abs2, diff) / n
        rms = sqrt(mse)
        maxerr = maximum(abs, diff)
        weighted_mse += snapshot_weights[i] * mse
        max_abs_error = max(max_abs_error, maxerr)
        push!(per, (index=included[i].index,
                    time=included[i].time,
                    rms=rms,
                    relative_rms=rms / denom,
                    max_abs=maxerr))
    end
    rms = sqrt(weighted_mse)
    (status=:ok,
     rms=rms,
     relative_rms=rms / denom,
     max_abs=max_abs_error,
     reference_rms=reference_rms,
     relative_denominator=denom,
     per_snapshot=per)
end

function _collapse_wave_phase_norm(xi_grid::Vector{Float64}, profiles,
                                   included, reference_position::Int,
                                   snapshot_weights::Vector{Float64},
                                   epsilon)
    epsilon === nothing &&
        return (status=:epsilon_not_available,
                successful=false,
                aggregate_rms=NaN,
                aggregate_relative_rms=NaN,
                max_abs=NaN,
                per_snapshot=NamedTuple[],
                diagnostics=NamedTuple[])

    per_diagnostics = NamedTuple[]
    for i in eachindex(profiles)
        diagnostic = try
            wave_diagnostics(xi_grid, profiles[i]; epsilon=epsilon)
        catch err
            push!(per_diagnostics,
                  (index=included[i].index,
                   time=included[i].time,
                   status=:failed,
                   successful=false,
                   message=sprint(showerror, err),
                   zero_crossing_count=0,
                   crest_count=0,
                   trough_count=0,
                   zero_crossing_xi=Float64[]))
            continue
        end
        push!(per_diagnostics,
              (index=included[i].index,
               time=included[i].time,
               status=diagnostic.status,
               successful=diagnostic.successful,
               message="",
               zero_crossing_count=diagnostic.zero_crossings.count,
               crest_count=diagnostic.extrema.crest_count,
               trough_count=diagnostic.extrema.trough_count,
               zero_crossing_xi=diagnostic.zero_crossings.xi))
    end

    ref = per_diagnostics[reference_position]
    if length(ref.zero_crossing_xi) < 2
        return (status=:insufficient_wave_train,
                successful=false,
                aggregate_rms=NaN,
                aggregate_relative_rms=NaN,
                max_abs=NaN,
                per_snapshot=per_diagnostics,
                diagnostics=per_diagnostics)
    end

    per = NamedTuple[]
    weighted_mse = 0.0
    max_abs_error = 0.0
    used = 0
    ref_scale = max(sqrt(sum(abs2, ref.zero_crossing_xi) /
                         length(ref.zero_crossing_xi)), eps(Float64))
    for i in eachindex(per_diagnostics)
        current = per_diagnostics[i]
        n = min(length(ref.zero_crossing_xi), length(current.zero_crossing_xi))
        if n < 2
            push!(per, merge(current, (phase_rms=NaN,
                                       phase_relative_rms=NaN,
                                       phase_max_abs=NaN,
                                       compared_zero_crossings=n)))
            continue
        end
        diff = current.zero_crossing_xi[1:n] .- ref.zero_crossing_xi[1:n]
        mse = sum(abs2, diff) / n
        rms = sqrt(mse)
        maxerr = maximum(abs, diff)
        weighted_mse += snapshot_weights[i] * mse
        max_abs_error = max(max_abs_error, maxerr)
        used += 1
        push!(per, merge(current, (phase_rms=rms,
                                   phase_relative_rms=rms / ref_scale,
                                   phase_max_abs=maxerr,
                                   compared_zero_crossings=n)))
    end

    used < 2 &&
        return (status=:insufficient_wave_train,
                successful=false,
                aggregate_rms=NaN,
                aggregate_relative_rms=NaN,
                max_abs=NaN,
                per_snapshot=per,
                diagnostics=per_diagnostics)
    rms = sqrt(weighted_mse)
    (status=:ok,
     successful=true,
     aggregate_rms=rms,
     aggregate_relative_rms=rms / ref_scale,
     max_abs=max_abs_error,
     per_snapshot=per,
     diagnostics=per_diagnostics)
end

function _collapse_score_components(norms, component_weights)
    components = NamedTuple[]
    for name in (:profile, :slope, :curvature, :velocity)
        weight = getfield(component_weights, name)
        norm = getfield(norms, name)
        if weight > 0 && norm.status == :ok && isfinite(norm.relative_rms)
            push!(components, (component=name, weight=weight,
                               relative_rms=norm.relative_rms))
        end
    end
    wave_weight = component_weights.wave_phase
    wave = norms.wave_phase
    if wave_weight > 0 && wave.status == :ok &&
       isfinite(wave.aggregate_relative_rms)
        push!(components, (component=:wave_phase, weight=wave_weight,
                           relative_rms=wave.aggregate_relative_rms))
    end
    isempty(components) && return (score=NaN, components=components)
    total = sum(component.weight for component in components)
    score = sqrt(sum(component.weight * component.relative_rms^2
                     for component in components) / total)
    (score=score, components=components)
end

function _collapse_per_snapshot_scores(norms, included, component_weights)
    scores = NamedTuple[]
    component_names = (:profile, :slope, :curvature, :velocity)
    for i in eachindex(included)
        terms = NamedTuple[]
        for name in component_names
            weight = getfield(component_weights, name)
            norm = getfield(norms, name)
            if weight > 0 && norm.status == :ok
                rel = norm.per_snapshot[i].relative_rms
                isfinite(rel) &&
                    push!(terms, (component=name, weight=weight,
                                  relative_rms=rel))
            end
        end
        total = sum(term.weight for term in terms; init=0.0)
        score = total > 0 ?
                sqrt(sum(term.weight * term.relative_rms^2
                         for term in terms) / total) :
                NaN
        push!(scores, (index=included[i].index,
                       time=included[i].time,
                       effective_time=included[i].effective_time,
                       score=score,
                       components=terms))
    end
    scores
end

"""
    similarity_collapse_diagnostics(z, t_snapshots, R, u; kwargs...)
    similarity_collapse_diagnostics(pde::PDESolution; kwargs...)

Quantify collapse of PDE snapshots in the similarity variables
`xi = z / ell(t)`, `S = R / ell(t)`, and
`U = u / (exponent * (t - time_offset)^(exponent - 1))`, where
`ell(t) = (t - time_offset)^exponent`.

The returned `NamedTuple` records the trusted `xi` window, interpolation grid,
snapshot/component weights, profile/slope/curvature/velocity norms, wave-phase
status when `epsilon` is available, excluded snapshots and boundary regions,
and an aggregate weighted score. Nonpositive effective-time snapshots are
excluded rather than treated as hard failures.
"""
function similarity_collapse_diagnostics(z, t_snapshots, R, u;
                                         epsilon=nothing,
                                         exponent::Real=2 / 3,
                                         time_offset::Real=0.0,
                                         xi_window=nothing,
                                         n_grid::Int=128,
                                         reference=:latest,
                                         component_weights::NamedTuple=(;),
                                         snapshot_weights=nothing,
                                         edge_fraction::Real=0.05,
                                         min_snapshots::Int=2)
    exponent_value = _collapse_positive_float("similarity exponent", exponent)
    time_offset_value = _collapse_finite_float("time_offset", time_offset)
    edge_fraction_value = _collapse_nonnegative_float("edge_fraction", edge_fraction)
    edge_fraction_value < 0.5 ||
        throw(ArgumentError("edge_fraction must be less than 0.5; got $edge_fraction"))
    n_grid >= 5 ||
        throw(ArgumentError("similarity collapse requires n_grid >= 5; got $n_grid"))
    min_snapshots >= 2 ||
        throw(ArgumentError("similarity collapse requires min_snapshots >= 2; got $min_snapshots"))
    weights = _collapse_component_weights(component_weights)
    supplied_window = _collapse_window_values(xi_window)
    if supplied_window !== nothing && !(supplied_window[1] < supplied_window[2])
        throw(ArgumentError("xi_window requires min < max; got $(supplied_window[1]) and $(supplied_window[2])"))
    end

    data = _validate_pde_solution_data(z, t_snapshots, R, u;
                                       context="similarity collapse diagnostics")
    rescaled = _rescaled_collapse_snapshots(
        data; exponent=exponent_value,
        time_offset=time_offset_value,
        edge_fraction=edge_fraction_value)
    included = rescaled.included
    excluded = rescaled.excluded

    if length(included) < min_snapshots
        return _collapse_failure(
            :insufficient_snapshots,
            "similarity collapse requires at least $min_snapshots positive-effective-time snapshots";
            exponent=exponent_value, time_offset=time_offset_value,
            n_grid=n_grid, component_weights=weights,
            included_snapshots=[(index=s.index, time=s.time,
                                 effective_time=s.effective_time)
                                for s in included],
            excluded_snapshots=excluded)
    end

    window_source = supplied_window === nothing ? :intersection_of_trusted_windows :
                    :provided
    xi_min, xi_max = if supplied_window === nothing
        (maximum(snapshot.trusted_min for snapshot in included),
         minimum(snapshot.trusted_max for snapshot in included))
    else
        supplied_window
    end

    if supplied_window !== nothing
        retained = NamedTuple[]
        for snapshot in included
            if snapshot.trusted_min <= xi_min && xi_max <= snapshot.trusted_max
                push!(retained, snapshot)
            else
                push!(excluded,
                      (index=snapshot.index,
                       time=snapshot.time,
                       effective_time=snapshot.effective_time,
                       reason=:outside_supplied_xi_window))
            end
        end
        included = retained
    end

    window = (min=xi_min, max=xi_max, source=window_source)
    if length(included) < min_snapshots
        return _collapse_failure(
            :insufficient_snapshots,
            "fewer than $min_snapshots snapshots cover the requested similarity window";
            exponent=exponent_value, time_offset=time_offset_value,
            n_grid=n_grid, component_weights=weights,
            xi_window_source=window_source,
            xi_window=(xi_min, xi_max),
            included_snapshots=[(index=s.index, time=s.time,
                                 effective_time=s.effective_time)
                                for s in included],
            excluded_snapshots=excluded)
    end

    if !(isfinite(xi_min) && isfinite(xi_max) && xi_min < xi_max)
        return _collapse_failure(
            :empty_xi_overlap,
            "included snapshots do not share a non-empty trusted similarity window";
            exponent=exponent_value, time_offset=time_offset_value,
            n_grid=n_grid, component_weights=weights,
            xi_window_source=window_source,
            xi_window=(xi_min, xi_max),
            included_snapshots=[(index=s.index, time=s.time,
                                 effective_time=s.effective_time)
                                for s in included],
            excluded_snapshots=excluded)
    end

    xi_grid = collect(range(xi_min, xi_max, length=n_grid))
    _validate_query_grid(xi_grid; context="similarity collapse interpolation grid")
    excluded_regions =
        [(index=snapshot.index,
          time=snapshot.time,
          effective_time=snapshot.effective_time,
          left_boundary=(min=first(snapshot.xi), max=xi_min,
                         reason=:truncated_tip_or_pretrusted_region),
          right_boundary=(min=xi_max, max=last(snapshot.xi),
                          reason=:far_outflow_or_posttrusted_region))
         for snapshot in included]
    snapshot_weight_data = _collapse_snapshot_weights(snapshot_weights,
                                                      length(included))
    normalized_weights = snapshot_weight_data.normalized

    profiles = [_interp_strict(snapshot.xi, snapshot.S, xi_grid;
                               context="similarity collapse profile")
                for snapshot in included]
    slopes = [_interp_strict(snapshot.xi, snapshot.slope, xi_grid;
                             context="similarity collapse slope")
              for snapshot in included]
    curvatures = [_interp_strict(snapshot.xi, snapshot.curvature, xi_grid;
                                 context="similarity collapse curvature")
                  for snapshot in included]
    velocities = [_interp_strict(snapshot.xi, snapshot.U, xi_grid;
                                 context="similarity collapse velocity")
                  for snapshot in included]

    reference_position = _collapse_reference_position(reference, included)
    norms = (profile=_collapse_component_norm(:profile, profiles, included,
                                              reference_position,
                                              normalized_weights),
             slope=_collapse_component_norm(:slope, slopes, included,
                                            reference_position,
                                            normalized_weights),
             curvature=_collapse_component_norm(:curvature, curvatures,
                                                included, reference_position,
                                                normalized_weights),
             velocity=_collapse_component_norm(:velocity, velocities,
                                               included, reference_position,
                                               normalized_weights),
             wave_phase=_collapse_wave_phase_norm(xi_grid, profiles, included,
                                                  reference_position,
                                                  normalized_weights,
                                                  epsilon))
    score = _collapse_score_components(norms, weights)
    per_snapshot_scores = _collapse_per_snapshot_scores(norms, included,
                                                        weights)
    score_successful = isfinite(score.score)

    (status=score_successful ? :ok : :no_weighted_components,
     successful=score_successful,
     message=score_successful ? "similarity collapse diagnostics computed" :
             "no finite weighted collapse components were available",
     source_status=_DEFAULT_SOURCE_STATUS,
     diagnostic_basis="implementation similarity-collapse metric; not Decent-King benchmark data",
     exponent=exponent_value,
     time_offset=time_offset_value,
     length_scale_convention="ell(t) = (t - time_offset)^exponent",
     velocity_scaling_convention="U = u / (exponent * (t - time_offset)^(exponent - 1))",
     xi_window=window,
     xi_grid=xi_grid,
     interpolation_grid=(variable=:xi, points=length(xi_grid),
                         min=first(xi_grid), max=last(xi_grid),
                         spacing=length(xi_grid) > 1 ?
                                 xi_grid[2] - xi_grid[1] : NaN),
     component_weights=weights,
     snapshot_weights=normalized_weights,
     snapshot_weights_raw=snapshot_weight_data.raw,
     included_snapshots=[(index=s.index, time=s.time,
                          effective_time=s.effective_time,
                          ell=s.ell,
                          velocity_scale=s.velocity_scale,
                          xi_min=first(s.xi),
                          xi_max=last(s.xi),
                          trusted_min=s.trusted_min,
                          trusted_max=s.trusted_max)
                         for s in included],
     excluded_snapshots=excluded,
     excluded_regions=excluded_regions,
     reference_snapshot_index=included[reference_position].index,
     reference_time=included[reference_position].time,
     reference_included_position=reference_position,
     norms=norms,
     per_snapshot_scores=per_snapshot_scores,
     aggregate_score=score.score,
     aggregate_components=score.components)
end

function similarity_collapse_diagnostics(pde::PDESolution; kwargs...)
    similarity_collapse_diagnostics(pde.z, pde.t_snapshots, pde.R, pde.u;
                                    epsilon=pde.ε, kwargs...)
end
