# Solver-independent residuals for the current reconstructed inner BVP.
#
# Source status: IMPL-inferred local reconstruction. These residual components
# expose the existing shooting boundary data without upgrading any claim to
# Decent-King 2008 article-body support.

using LinearAlgebra: norm

export InnerBVPUnknowns, InnerBVPResidual,
       inner_bvp_residual, inner_bvp_residual_components,
       inner_bvp_residual_vector

const _INNER_BVP_RESIDUAL_COMPONENT_NAMES =
    (:far_field_slope, :far_field_velocity, :far_field_curvature)
const _INNER_BVP_RESIDUAL_COMPONENT_STATUS =
    (far_field_slope=_DEFAULT_SOURCE_STATUS,
     far_field_velocity=_DEFAULT_SOURCE_STATUS,
     far_field_curvature=_DEFAULT_SOURCE_STATUS)
const _INNER_BVP_BOUNDARY_STATUS =
    (tip_slope=_DEFAULT_SOURCE_STATUS,
     tip_velocity=_DEFAULT_SOURCE_STATUS,
     tip_third_derivative=_DEFAULT_SOURCE_STATUS,
     far_field_slope=_DEFAULT_SOURCE_STATUS,
     far_field_velocity=_DEFAULT_SOURCE_STATUS,
     far_field_curvature=_DEFAULT_SOURCE_STATUS)

struct InnerBVPUnknowns
    ξ₀::Float64
    S₀::Float64
    Sξξ₀::Float64

    function InnerBVPUnknowns(ξ₀::Real, S₀::Real, Sξξ₀::Real)
        ξ₀v = Float64(ξ₀)
        S₀v = Float64(S₀)
        Sξξ₀v = Float64(Sξξ₀)
        _require_inner_bvp_nonnegative_finite("ξ₀", ξ₀v)
        _require_inner_bvp_positive_finite("S₀", S₀v)
        _require_inner_bvp_finite("Sξξ₀", Sξξ₀v)
        new(ξ₀v, S₀v, Sξξ₀v)
    end
end

function InnerBVPUnknowns(xs)
    (xs isa Tuple || xs isa AbstractVector) ||
        throw(ArgumentError("InnerBVPUnknowns requires a 3-tuple or vector; got $(typeof(xs))"))
    length(xs) == 3 ||
        throw(ArgumentError("InnerBVPUnknowns requires exactly 3 values: (ξ₀, S₀, Sξξ₀); got $(length(xs))"))
    InnerBVPUnknowns(xs[1], xs[2], xs[3])
end

InnerBVPUnknowns(u::InnerBVPUnknowns) = u
InnerBVPUnknowns(; ξ₀, S₀, Sξξ₀) = InnerBVPUnknowns(ξ₀, S₀, Sξξ₀)

struct InnerBVPResidual
    unknowns::InnerBVPUnknowns
    components::NamedTuple
    parameters::NamedTuple
    domain::NamedTuple
    mesh::NamedTuple
    diagnostics::NamedTuple
    provenance::NamedTuple
end

function _require_inner_bvp_positive_finite(name::AbstractString, x)
    isfinite(x) && x > 0 ||
        throw(ArgumentError("$name must be positive and finite; got $x"))
    nothing
end

function _require_inner_bvp_nonnegative_finite(name::AbstractString, x)
    isfinite(x) && x >= 0 ||
        throw(ArgumentError("$name must be nonnegative and finite; got $x"))
    nothing
end

function _require_inner_bvp_finite(name::AbstractString, x)
    isfinite(x) ||
        throw(ArgumentError("$name must be finite; got $x"))
    nothing
end

function _validate_inner_bvp_residual_inputs(unknowns::InnerBVPUnknowns,
                                             ε::Float64,
                                             ξ_max::Float64,
                                             maxiters::Int)
    _require_inner_bvp_positive_finite("ε", ε)
    _require_inner_bvp_finite("ξ_max", ξ_max)
    ξ_max > unknowns.ξ₀ ||
        throw(ArgumentError("ξ_max must be greater than ξ₀; got ξ_max=$ξ_max and ξ₀=$(unknowns.ξ₀)"))
    maxiters > 0 ||
        throw(ArgumentError("maxiters must be positive; got $maxiters"))
    nothing
end

function _inner_bvp_components(ξ_vals, S_vals, Spp_vals, U_vals, ε::Float64)
    (far_field_slope=Float64(S_vals[end] / ξ_vals[end] - ε),
     far_field_velocity=Float64(U_vals[end]),
     far_field_curvature=Float64(Spp_vals[end]))
end

_inner_bvp_residual_vector(components::NamedTuple) =
    Float64[getfield(components, name) for name in _INNER_BVP_RESIDUAL_COMPONENT_NAMES]

function _inner_bvp_evaluate(unknowns::InnerBVPUnknowns;
                             ε::Real=0.1,
                             epsilon::Union{Nothing,Real}=nothing,
                             ξ_max::Real=60.0,
                             maxiters::Int=2_000_000,
                             context::AbstractString="inner BVP residual")
    εv = Float64(something(epsilon, ε))
    ξ_maxv = Float64(ξ_max)
    _validate_inner_bvp_residual_inputs(unknowns, εv, ξ_maxv, maxiters)

    ξ_vals, S_vals, Sp_vals, Spp_vals, U_vals, shoot_diagnostics =
        _shoot_with_diagnostics(unknowns.ξ₀, unknowns.S₀, unknowns.Sξξ₀,
                                ξ_maxv; maxiters=maxiters, context=context)
    components = _inner_bvp_components(ξ_vals, S_vals, Spp_vals, U_vals, εv)
    vector = _inner_bvp_residual_vector(components)
    residual_norm_value = norm(vector)
    diagnostics = merge(
        shoot_diagnostics,
        (successful=Bool(get(shoot_diagnostics, :successful, false)) &&
                    all(isfinite, vector),
         residual_norm=residual_norm_value,
         final_residual_norm=residual_norm_value,
         component_names=_INNER_BVP_RESIDUAL_COMPONENT_NAMES,
         component_status=_INNER_BVP_RESIDUAL_COMPONENT_STATUS,
         endpoint_state=(S=Float64(S_vals[end]),
                         Sξ=Float64(Sp_vals[end]),
                         Sξξ=Float64(Spp_vals[end]),
                         U=Float64(U_vals[end])),
         source_status=_DEFAULT_SOURCE_STATUS))
    (components=components,
     vector=vector,
     ξ=ξ_vals,
     diagnostics=diagnostics,
     ε=εv,
     ξ_max=ξ_maxv,
     maxiters=maxiters)
end

function _inner_bvp_parameters(ε::Float64)
    (ε=ε,
     far_field_targets=(slope=ε, velocity=0.0, curvature=0.0),
     tip_conditions=(slope=0.0, velocity_factor=4/5, third_derivative=0.0),
     component_names=_INNER_BVP_RESIDUAL_COMPONENT_NAMES,
     boundary_status=_INNER_BVP_BOUNDARY_STATUS)
end

function _inner_bvp_domain(unknowns::InnerBVPUnknowns, ξ_vals, ξ_max::Float64)
    (independent_variable=:ξ,
     ξ₀=unknowns.ξ₀,
     ξ_start=Float64(first(ξ_vals)),
     ξ_max=ξ_max,
     integration_offset=Float64(first(ξ_vals) - unknowns.ξ₀))
end

function _inner_bvp_mesh(ξ_vals)
    merge((ξ=ξ_vals,), mesh_summary(ξ_vals; variable=:ξ))
end

function _inner_bvp_provenance(settings::NamedTuple)
    metadata = default_recoil_provenance_metadata(:inner_bvp_residual;
                                                  solver_settings=settings)
    (problem_kind=:inner_bvp_residual,
     source_status=_DEFAULT_SOURCE_STATUS,
     source_ledger=_SOURCE_LEDGER_PATH,
     canonical_source_doi=_CANONICAL_CONE_DOI,
     status_note="Current inner BVP residual components, tip regularity, and far-field boundary data are local reconstructed implementation data pending Decent-King 2008 article-body verification.",
     component_status=_INNER_BVP_RESIDUAL_COMPONENT_STATUS,
     boundary_status=_INNER_BVP_BOUNDARY_STATUS,
     metadata=metadata,
     source_citations=metadata.source_citations,
     source_ids=metadata.source_ids,
     assumptions=metadata.assumptions,
     benchmark_ids=metadata.benchmark_ids,
     solver_settings=metadata.solver_settings,
     artifacts=metadata.artifacts,
     package=metadata.package)
end

"""
    inner_bvp_residual(unknowns; ε=0.1, ξ_max=60.0,
                       maxiters=2_000_000)

Evaluate the current reconstructed inner BVP endpoint residual independently of
the Newton shooting solver. `unknowns` are explicit solve variables
`(ξ₀, S₀, Sξξ₀)`. Components are labeled as `:far_field_slope`,
`:far_field_velocity`, and `:far_field_curvature`.

Ledger status: IMPL-inferred local reconstruction pending Decent-King 2008
article-body verification.
"""
function inner_bvp_residual(unknowns::InnerBVPUnknowns; kwargs...)
    eval = _inner_bvp_evaluate(unknowns; kwargs...)
    parameters = _inner_bvp_parameters(eval.ε)
    domain = _inner_bvp_domain(unknowns, eval.ξ, eval.ξ_max)
    solver = (maxiters=eval.maxiters,
              algorithm="Rodas5P endpoint shooting residual")
    settings = merge(parameters, domain, solver)
    InnerBVPResidual(unknowns, eval.components, parameters, domain,
                     _inner_bvp_mesh(eval.ξ), eval.diagnostics,
                     _inner_bvp_provenance(settings))
end

inner_bvp_residual(xs; kwargs...) =
    inner_bvp_residual(InnerBVPUnknowns(xs); kwargs...)

inner_bvp_residual_components(residual::InnerBVPResidual) = residual.components
inner_bvp_residual_components(unknowns::InnerBVPUnknowns; kwargs...) =
    inner_bvp_residual_components(inner_bvp_residual(unknowns; kwargs...))
inner_bvp_residual_components(xs; kwargs...) =
    inner_bvp_residual_components(InnerBVPUnknowns(xs); kwargs...)

inner_bvp_residual_vector(residual::InnerBVPResidual) =
    _inner_bvp_residual_vector(residual.components)

function inner_bvp_residual_vector(unknowns::InnerBVPUnknowns; kwargs...)
    eval = _inner_bvp_evaluate(unknowns; kwargs...)
    eval.vector
end

inner_bvp_residual_vector(xs; kwargs...) =
    inner_bvp_residual_vector(InnerBVPUnknowns(xs); kwargs...)
