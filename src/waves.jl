# Capillary-wave diagnostics for excess profiles S(xi) - epsilon*xi.
#
# Ledger status:
# - C2001 backs only the qualitative observation that Figure 1 has
#   high-frequency oscillations modulated on a longer scale.
# - All numerical phase, wavelength, and envelope quantities reported here are
#   IMPL-inferred local diagnostics for the current reconstructed profiles.

export wave_diagnostics

const _WAVE_C2001_QUALITATIVE_ID = "DK2001-figure1-oscillation-qualitative"

function _wave_epsilon(epsilon; context::AbstractString)
    epsilon === nothing &&
        throw(ArgumentError("$context requires positive finite epsilon"))
    eps_value = Float64(epsilon)
    isfinite(eps_value) && eps_value > 0 ||
        throw(ArgumentError("$context requires positive finite epsilon; got $epsilon"))
    eps_value
end

function _wave_input_vectors(xi, S, epsilon; context::AbstractString,
                             min_samples::Int)
    min_samples >= 5 ||
        throw(ArgumentError("$context requires min_samples >= 5; got $min_samples"))
    xi_values = Float64.(collect(xi))
    S_values = Float64.(collect(S))
    length(xi_values) == length(S_values) ||
        throw(ArgumentError("$context requires xi and S to have equal lengths; got $(length(xi_values)) and $(length(S_values))"))
    length(xi_values) >= min_samples ||
        throw(ArgumentError("$context requires at least $min_samples samples; got $(length(xi_values))"))
    _validate_interpolation_table(xi_values, S_values; context=context)
    eps_value = _wave_epsilon(epsilon; context=context)
    excess = S_values .- eps_value .* xi_values
    all(isfinite, excess) ||
        throw(ArgumentError("$context produced non-finite excess values"))
    (xi=xi_values, S=S_values, epsilon=eps_value, excess=excess)
end

function _wave_sign(value::Float64, deadband::Float64)
    value > deadband && return 1
    value < -deadband && return -1
    0
end

function _zero_location(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    if y1 == y2
        return (x1 + x2) / 2
    end
    x = x1 - y1 * (x2 - x1) / (y2 - y1)
    clamp(x, min(x1, x2), max(x1, x2))
end

function _zero_crossings(xi::Vector{Float64}, values::Vector{Float64},
                         deadband::Float64)
    locations = Float64[]
    brackets = Tuple{Int,Int}[]
    last_index = 0
    last_sign = 0
    for i in eachindex(values)
        current_sign = _wave_sign(values[i], deadband)
        current_sign == 0 && continue
        if last_sign != 0 && current_sign != last_sign
            push!(locations,
                  _zero_location(xi[last_index], values[last_index],
                                 xi[i], values[i]))
            push!(brackets, (last_index, i))
        end
        last_index = i
        last_sign = current_sign
    end
    (count=length(locations), xi=locations, brackets=brackets)
end

function _quadratic_turning_point(x0::Float64, y0::Float64,
                                  x1::Float64, y1::Float64,
                                  x2::Float64, y2::Float64)
    h0 = x0 - x1
    h2 = x2 - x1
    if h0 == 0 || h2 == 0 || h0 == h2
        return (xi=x1, value=y1)
    end
    A = [h0^2 h0 1.0;
         0.0  0.0 1.0;
         h2^2 h2 1.0]
    a, b, c = A \ [y0, y1, y2]
    if !isfinite(a) || !isfinite(b) || abs(a) <= eps(Float64)
        return (xi=x1, value=y1)
    end
    h_vertex = -b / (2a)
    if min(h0, h2) <= h_vertex <= max(h0, h2)
        value = a * h_vertex^2 + b * h_vertex + c
        return (xi=x1 + h_vertex, value=value)
    end
    (xi=x1, value=y1)
end

function _turning_points(xi::Vector{Float64}, values::Vector{Float64},
                         deadband::Float64)
    crest_xi = Float64[]
    crest_values = Float64[]
    crest_indices = Int[]
    trough_xi = Float64[]
    trough_values = Float64[]
    trough_indices = Int[]
    for i in 2:(length(values) - 1)
        abs(values[i]) > deadband || continue
        left = values[i] - values[i - 1]
        right = values[i + 1] - values[i]
        if left > 0 && right < 0
            refined = _quadratic_turning_point(xi[i - 1], values[i - 1],
                                               xi[i], values[i],
                                               xi[i + 1], values[i + 1])
            push!(crest_xi, refined.xi)
            push!(crest_values, refined.value)
            push!(crest_indices, i)
        elseif left < 0 && right > 0
            refined = _quadratic_turning_point(xi[i - 1], values[i - 1],
                                               xi[i], values[i],
                                               xi[i + 1], values[i + 1])
            push!(trough_xi, refined.xi)
            push!(trough_values, refined.value)
            push!(trough_indices, i)
        end
    end
    (crest_count=length(crest_xi),
     trough_count=length(trough_xi),
     count=length(crest_xi) + length(trough_xi),
     crests=(xi=crest_xi, value=crest_values, indices=crest_indices),
     troughs=(xi=trough_xi, value=trough_values, indices=trough_indices))
end

function _median_or_nan(values)
    isempty(values) && return NaN
    sorted = sort(Float64.(collect(values)))
    n = length(sorted)
    isodd(n) ? sorted[(n + 1) ÷ 2] : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2
end

function _mean_or_nan(values)
    isempty(values) && return NaN
    sum(values) / length(values)
end

function _wavelengths_from_locations(locations::Vector{Float64})
    if length(locations) < 2
        return (xi=Float64[], value=Float64[], count=0)
    end
    xi_mid = [(locations[i] + locations[i + 1]) / 2
              for i in 1:(length(locations) - 1)]
    values = [2 * (locations[i + 1] - locations[i])
              for i in 1:(length(locations) - 1)]
    (xi=xi_mid, value=values, count=length(values))
end

function _same_kind_wavelengths(locations::Vector{Float64})
    if length(locations) < 2
        return (xi=Float64[], value=Float64[], count=0)
    end
    xi_mid = [(locations[i] + locations[i + 1]) / 2
              for i in 1:(length(locations) - 1)]
    values = [locations[i + 1] - locations[i]
              for i in 1:(length(locations) - 1)]
    (xi=xi_mid, value=values, count=length(values))
end

function _local_max_spacing(xi::Vector{Float64}, center::Float64,
                            width::Float64)
    spacings = Float64[]
    lo = center - width / 2
    hi = center + width / 2
    for i in 1:(length(xi) - 1)
        midpoint = (xi[i] + xi[i + 1]) / 2
        if lo <= midpoint <= hi
            push!(spacings, xi[i + 1] - xi[i])
        end
    end
    isempty(spacings) ? maximum(diff(xi)) : maximum(spacings)
end

function _resolution_summary(xi::Vector{Float64}, wavelengths;
                             min_samples_per_wavelength::Float64)
    if isempty(wavelengths.value)
        return (status=:insufficient_wave_train,
                threshold=min_samples_per_wavelength,
                count=0,
                samples_per_wavelength=Float64[],
                min_samples_per_wavelength=NaN,
                median_samples_per_wavelength=NaN)
    end
    samples = Float64[]
    for (center, wavelength) in zip(wavelengths.xi, wavelengths.value)
        spacing = _local_max_spacing(xi, center, wavelength)
        push!(samples, wavelength / spacing)
    end
    min_samples = minimum(samples)
    median_samples = _median_or_nan(samples)
    status = min_samples >= min_samples_per_wavelength ?
             :sufficient_resolution : :insufficient_resolution
    (status=status,
     threshold=min_samples_per_wavelength,
     count=length(samples),
     samples_per_wavelength=samples,
     min_samples_per_wavelength=min_samples,
     median_samples_per_wavelength=median_samples)
end

function _linear_fit(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    xmean = sum(xs) / n
    ymean = sum(ys) / n
    den = sum((x - xmean)^2 for x in xs)
    if den <= 0
        return (slope=NaN, intercept=NaN, residual_rms=NaN)
    end
    slope = sum((xs[i] - xmean) * (ys[i] - ymean) for i in 1:n) / den
    intercept = ymean - slope * xmean
    residuals = [ys[i] - (intercept + slope * xs[i]) for i in 1:n]
    (slope=slope,
     intercept=intercept,
     residual_rms=sqrt(sum(abs2, residuals) / n))
end

function _monotone_decay_fraction(amplitudes::Vector{Float64})
    length(amplitudes) < 2 && return NaN
    decreases = count(i -> amplitudes[i + 1] <= amplitudes[i],
                      1:(length(amplitudes) - 1))
    decreases / (length(amplitudes) - 1)
end

function _envelope_summary(extrema, deadband::Float64, min_points::Int)
    points = [(xi=x, amplitude=abs(v), value=v, kind=:crest)
              for (x, v) in zip(extrema.crests.xi, extrema.crests.value)]
    append!(points, [(xi=x, amplitude=abs(v), value=v, kind=:trough)
                     for (x, v) in zip(extrema.troughs.xi,
                                       extrema.troughs.value)])
    sort!(points, by=p -> p.xi)
    filtered = [p for p in points if p.amplitude > deadband]
    if length(filtered) < min_points
        return (status=:insufficient_extrema,
                points=length(filtered),
                xi=[p.xi for p in filtered],
                amplitude=[p.amplitude for p in filtered],
                value=[p.value for p in filtered],
                kind=[p.kind for p in filtered],
                log_fit_slope=NaN,
                log_fit_intercept=NaN,
                fit_rms=NaN,
                decay_rate=NaN,
                decay_length=NaN,
                amplitude_ratio=NaN,
                monotone_decay_fraction=NaN)
    end
    envelope_xi = [p.xi for p in filtered]
    amplitudes = [p.amplitude for p in filtered]
    logs = log.(amplitudes)
    fit = _linear_fit(envelope_xi, logs)
    decay_rate = -fit.slope
    decay_length = decay_rate > 0 ? 1 / decay_rate : Inf
    (status=:ok,
     points=length(filtered),
     xi=envelope_xi,
     amplitude=amplitudes,
     value=[p.value for p in filtered],
     kind=[p.kind for p in filtered],
     log_fit_slope=fit.slope,
     log_fit_intercept=fit.intercept,
     fit_rms=fit.residual_rms,
     decay_rate=decay_rate,
     decay_length=decay_length,
     amplitude_ratio=last(amplitudes) / first(amplitudes),
     monotone_decay_fraction=_monotone_decay_fraction(amplitudes))
end

function _phase_summary(zero_crossings, wavelengths)
    phases = [Float64(i - 1) * pi for i in 1:zero_crossings.count]
    local_phase = if isempty(wavelengths.xi)
        Float64[]
    else
        [(phases[i] + phases[i + 1]) / 2 for i in 1:(length(phases) - 1)]
    end
    (status=zero_crossings.count >= 2 ? :ok : :insufficient_zero_crossings,
     origin=zero_crossings.count == 0 ? NaN : first(zero_crossings.xi),
     convention="successive zero crossings advance phase by pi; absolute phase is local",
     zero_crossing_xi=zero_crossings.xi,
     zero_crossing_phase=phases,
     local_phase=(xi=wavelengths.xi,
                  phase=local_phase,
                  wavelength=wavelengths.value,
                  wavenumber=[2pi / w for w in wavelengths.value]))
end

function _wave_source_basis()
    (qualitative=(source_id=_WAVE_C2001_QUALITATIVE_ID,
                  source_status="C2001",
                  fact_type=:qualitative,
                  claim="high-frequency oscillations modulated on a longer scale"),
     quantitative=(source_status=_DEFAULT_SOURCE_STATUS,
                   note="phase, wavelength, and envelope numbers are local implementation diagnostics, not Decent-King 2008 laws"))
end

function _wave_status(zero_crossings, extrema, envelope, resolution;
                      min_zero_crossings::Int, min_extrema::Int)
    if zero_crossings.count < min_zero_crossings || extrema.count < min_extrema
        return :insufficient_wave_train
    elseif resolution.status != :sufficient_resolution
        return :insufficient_resolution
    elseif envelope.status != :ok
        return :insufficient_envelope
    end
    :ok
end

"""
    wave_diagnostics(xi, S; epsilon, ...)

Measure local capillary-wave diagnostics for the excess profile
`S(xi) - epsilon*xi`. The returned `NamedTuple` reports zero crossings,
crests/troughs, local wavelength and phase estimates, envelope decay, and
resolution status. Quantitative values are local implementation diagnostics;
only the qualitative existence of high-frequency oscillations modulated on a
longer scale is C2001-backed in the current source ledger.
"""
function wave_diagnostics(xi, S; epsilon=nothing,
                          epsilon_source::Symbol=:provided,
                          profile_kind::Symbol=:sampled_radius,
                          min_samples::Int=5,
                          min_zero_crossings::Int=4,
                          min_extrema::Int=3,
                          min_envelope_points::Int=3,
                          min_samples_per_wavelength::Real=8.0,
                          absolute_deadband::Real=0.0,
                          relative_deadband::Real=1e-6,
                          min_modulation_scale_ratio::Real=2.0)
    context = "wave diagnostics"
    min_zero_crossings >= 2 ||
        throw(ArgumentError("$context requires min_zero_crossings >= 2; got $min_zero_crossings"))
    min_extrema >= 1 ||
        throw(ArgumentError("$context requires min_extrema >= 1; got $min_extrema"))
    min_envelope_points >= 2 ||
        throw(ArgumentError("$context requires min_envelope_points >= 2; got $min_envelope_points"))
    samples_per_wavelength_threshold = Float64(min_samples_per_wavelength)
    isfinite(samples_per_wavelength_threshold) &&
        samples_per_wavelength_threshold > 0 ||
        throw(ArgumentError("$context requires positive finite min_samples_per_wavelength; got $min_samples_per_wavelength"))
    abs_deadband = Float64(absolute_deadband)
    rel_deadband = Float64(relative_deadband)
    isfinite(abs_deadband) && abs_deadband >= 0 ||
        throw(ArgumentError("$context requires nonnegative finite absolute_deadband; got $absolute_deadband"))
    isfinite(rel_deadband) && rel_deadband >= 0 ||
        throw(ArgumentError("$context requires nonnegative finite relative_deadband; got $relative_deadband"))
    scale_ratio_threshold = Float64(min_modulation_scale_ratio)
    isfinite(scale_ratio_threshold) && scale_ratio_threshold >= 0 ||
        throw(ArgumentError("$context requires nonnegative finite min_modulation_scale_ratio; got $min_modulation_scale_ratio"))

    data = _wave_input_vectors(xi, S, epsilon; context=context,
                               min_samples=min_samples)
    max_abs_excess = maximum(abs, data.excess)
    deadband = max(abs_deadband, rel_deadband * max_abs_excess)

    zero_crossings = _zero_crossings(data.xi, data.excess, deadband)
    extrema = _turning_points(data.xi, data.excess, deadband)
    zero_wavelengths = _wavelengths_from_locations(zero_crossings.xi)
    crest_wavelengths = _same_kind_wavelengths(extrema.crests.xi)
    trough_wavelengths = _same_kind_wavelengths(extrema.troughs.xi)
    wavelength_values = zero_wavelengths.value
    resolution = _resolution_summary(
        data.xi, zero_wavelengths;
        min_samples_per_wavelength=samples_per_wavelength_threshold)
    envelope = _envelope_summary(extrema, deadband, min_envelope_points)
    phase = _phase_summary(zero_crossings, zero_wavelengths)

    median_wavelength = _median_or_nan(wavelength_values)
    mean_wavelength = _mean_or_nan(wavelength_values)
    modulation_scale_ratio =
        isfinite(envelope.decay_length) && isfinite(median_wavelength) &&
        median_wavelength > 0 ? envelope.decay_length / median_wavelength : Inf
    status = _wave_status(zero_crossings, extrema, envelope, resolution;
                          min_zero_crossings=min_zero_crossings,
                          min_extrema=min_extrema)
    oscillatory = zero_crossings.count >= min_zero_crossings &&
                  extrema.count >= min_extrema
    longer_scale = oscillatory && envelope.status == :ok &&
                   modulation_scale_ratio >= scale_ratio_threshold
    qualitative_consistency = oscillatory && longer_scale ?
                              :consistent_with_C2001_qualitative :
                              :not_established_by_current_profile

    (status=status,
     successful=status == :ok,
     profile_kind=profile_kind,
     epsilon=data.epsilon,
     epsilon_source=epsilon_source,
     source_status=_DEFAULT_SOURCE_STATUS,
     quantitative_source_status=_DEFAULT_SOURCE_STATUS,
     source_basis=_wave_source_basis(),
     samples=(count=length(data.xi),
              xi_min=first(data.xi),
              xi_max=last(data.xi),
              xi_span=last(data.xi) - first(data.xi),
              max_abs_excess=max_abs_excess,
              deadband=deadband),
     zero_crossings=merge(zero_crossings,
                          (status=zero_crossings.count >= min_zero_crossings ?
                                  :ok : :insufficient_zero_crossings,)),
     extrema=merge(extrema,
                   (status=extrema.count >= min_extrema ?
                           :ok : :insufficient_extrema,)),
     wavelength=(status=zero_wavelengths.count > 0 ? :ok :
                        :insufficient_zero_crossings,
                 zero_crossing=merge(zero_wavelengths,
                                     (mean=mean_wavelength,
                                      median=median_wavelength)),
                 crest_to_crest=crest_wavelengths,
                 trough_to_trough=trough_wavelengths),
     phase=phase,
     envelope=envelope,
     resolution=resolution,
     qualitative=(oscillatory=oscillatory,
                  modulated_on_longer_scale=longer_scale,
                  modulation_scale_ratio=modulation_scale_ratio,
                  min_modulation_scale_ratio=scale_ratio_threshold,
                  source_consistency=qualitative_consistency))
end

function wave_diagnostics(sol::InnerSolution; epsilon=nothing, kwargs...)
    eps_value = if epsilon === nothing
        _validate_interpolation_table(sol.ξ, sol.S;
                                      context="inner wave diagnostics")
        sol.S[end] / sol.ξ[end]
    else
        epsilon
    end
    eps_source = epsilon === nothing ? :inferred_endpoint_slope : :provided
    wave_diagnostics(sol.ξ, sol.S; epsilon=eps_value,
                     epsilon_source=eps_source,
                     profile_kind=:inner_solution, kwargs...)
end

function wave_diagnostics(sol::OuterSolution; epsilon=nothing, kwargs...)
    eps_value = epsilon === nothing ? sol.ε : epsilon
    radius = sol.ε .* sol.ξ .+ sol.s₁
    eps_source = epsilon === nothing ? :outer_solution : :provided
    wave_diagnostics(sol.ξ, radius; epsilon=eps_value,
                     epsilon_source=eps_source,
                     profile_kind=:outer_solution, kwargs...)
end

function wave_diagnostics(sol::CompositeSolution; epsilon=nothing, kwargs...)
    eps_value = epsilon === nothing ? sol.ε : epsilon
    eps_source = epsilon === nothing ? :composite_solution : :provided
    wave_diagnostics(sol.ξ, sol.S; epsilon=eps_value,
                     epsilon_source=eps_source,
                     profile_kind=:composite_solution, kwargs...)
end
