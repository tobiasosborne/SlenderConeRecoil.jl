# Explicit region and common-part metadata for the current reconstructed
# matched-asymptotic workflow.
#
# Ledger status: IMPL-inferred unless a caller provides stronger source
# metadata. These objects describe the implementation's diagnostic
# inner/outer/common-part data; they do not claim Decent-King 2008 article-body
# fidelity.

export AsymptoticRegion, CommonPart, CompositeParts,
       evaluate_common_part, common_part_values

struct AsymptoticRegion
    name::Symbol
    variables::NamedTuple
    transform::NamedTuple
    assumptions::Tuple{Vararg{Assumption}}
    truncation::NamedTuple
    source_status::String
    source_ids::Tuple{Vararg{SourceID}}
    references::Tuple{Vararg{SourceCitation}}
    provenance::NamedTuple
    mesh::NamedTuple
    data::NamedTuple
    diagnostics::NamedTuple
end

struct CommonPart
    name::Symbol
    expression::Symbol
    variables::NamedTuple
    transform::NamedTuple
    assumptions::Tuple{Vararg{Assumption}}
    truncation::NamedTuple
    source_status::String
    source_ids::Tuple{Vararg{SourceID}}
    references::Tuple{Vararg{SourceCitation}}
    provenance::NamedTuple
    slope::Float64
    intercept::Float64
    ξ_match::Float64
    fit_window::NamedTuple
    mesh::NamedTuple
    data::NamedTuple
    diagnostics::NamedTuple
end

struct CompositeParts
    inner::AsymptoticRegion
    outer::AsymptoticRegion
    common::CommonPart
    assumptions::Tuple{Vararg{String}}
    diagnostics::NamedTuple
end

function _typed_tuple(::Type{T}, xs, field::AbstractString) where {T}
    xs === nothing && return ()
    xs isa T && return (xs,)
    xs isa AbstractString &&
        throw(ArgumentError("$field must contain $(T) objects, not a string"))
    values = Tuple(xs)
    all(x -> x isa T, values) ||
        throw(ArgumentError("$field must contain only $(T) objects"))
    values
end

function _source_id_tuple(ids)
    ids === nothing && return ()
    ids isa SourceID && return (ids,)
    if ids isa AbstractString || ids isa Symbol
        return (SourceID(string(ids);
                         kind=:implementation,
                         status=_DEFAULT_SOURCE_STATUS,
                         ledger_path=_SOURCE_LEDGER_PATH),)
    end
    values = Tuple(ids)
    map(values) do id
        id isa SourceID && return id
        (id isa AbstractString || id isa Symbol) &&
            return SourceID(string(id);
                            kind=:implementation,
                            status=_DEFAULT_SOURCE_STATUS,
                            ledger_path=_SOURCE_LEDGER_PATH)
        throw(ArgumentError("source_ids must contain SourceID, string, or symbol values"))
    end
end

function _citation_tuple(refs)
    refs === nothing && return ()
    refs isa SourceCitation && return (refs,)
    refs isa AbstractString &&
        return (SourceCitation(refs; status="unverified",
                               ledger_path=_SOURCE_LEDGER_PATH),)
    values = Tuple(refs)
    map(values) do ref
        ref isa SourceCitation && return ref
        ref isa AbstractString &&
            return SourceCitation(ref; status="unverified",
                                  ledger_path=_SOURCE_LEDGER_PATH)
        throw(ArgumentError("references must contain SourceCitation or string values"))
    end
end

_canonical_reference_tuple() = (_canonical_cone_citation(),)

_region_provenance(kind::Symbol; source_status::AbstractString=_DEFAULT_SOURCE_STATUS) =
    (source_status=String(source_status),
     source_ledger=_SOURCE_LEDGER_PATH,
     canonical_source_doi=_CANONICAL_CONE_DOI,
     implementation_status="local reconstructed asymptotic-region metadata pending Decent-King 2008 article body",
     problem_kind=kind)

function _region_source_ids(label::AbstractString)
    (SourceID("DK2008-canonical-metadata";
              kind=:source,
              description="Canonical cone-recoil paper metadata; article-body matching and composite formula extraction remains blocked.",
              citation_id="DK2008-cone",
              status="C2008-meta",
              ledger_path=_SOURCE_LEDGER_PATH),
     SourceID("IMPL-inferred-$(label)";
              kind=:implementation,
              description="Current $(label) asymptotic-region/common-part object is reconstructed from implementation data.",
              citation_id="DK2008-cone",
              status=_DEFAULT_SOURCE_STATUS,
              ledger_path=_SOURCE_LEDGER_PATH))
end

function _region_assumption(label::AbstractString, description::AbstractString)
    (Assumption("current-$(label)-local-reconstruction";
                description=description,
                status=_DEFAULT_SOURCE_STATUS,
                source_ids=("IMPL-inferred-$(label)",)),)
end

function AsymptoticRegion(name::Symbol; variables::NamedTuple,
                          transform::NamedTuple=(;),
                          assumptions=(), truncation::NamedTuple=(;),
                          source_status::AbstractString=_DEFAULT_SOURCE_STATUS,
                          source_ids=(), references=(),
                          provenance::NamedTuple=(;),
                          mesh::NamedTuple=(;), data::NamedTuple=(;),
                          diagnostics::NamedTuple=(;))
    assumption_tuple = _typed_tuple(Assumption, assumptions, "assumptions")
    AsymptoticRegion(name, variables, transform,
                     assumption_tuple, truncation, String(source_status),
                     _source_id_tuple(source_ids), _citation_tuple(references),
                     provenance, mesh, data, diagnostics)
end

function CommonPart(name::Symbol=:linear_inner_far_field;
                    expression::Symbol=:slope_times_xi_plus_intercept,
                    variables::NamedTuple=(independent=:ξ, dependent=:S),
                    transform::NamedTuple=(coordinate=:overlap_ξ,
                                           representation=:linear_fit),
                    assumptions=(), truncation::NamedTuple=(order=:linear,),
                    source_status::AbstractString=_DEFAULT_SOURCE_STATUS,
                    source_ids=(), references=(),
                    provenance::NamedTuple=(;),
                    slope::Real=NaN, intercept::Real=NaN,
                    ξ_match::Real=NaN, fit_window::NamedTuple=(;),
                    mesh::NamedTuple=(;), data::NamedTuple=(;),
                    diagnostics::NamedTuple=(;))
    assumption_tuple = _typed_tuple(Assumption, assumptions, "assumptions")
    CommonPart(name, expression, variables, transform, assumption_tuple,
               truncation, String(source_status), _source_id_tuple(source_ids),
               _citation_tuple(references), provenance, Float64(slope),
               Float64(intercept), Float64(ξ_match), fit_window, mesh, data,
               diagnostics)
end

function evaluate_common_part(common::CommonPart, ξ)
    common.slope .* ξ .+ common.intercept
end

common_part_values(common::CommonPart, ξs) = evaluate_common_part(common, ξs)

function Base.getproperty(common::CommonPart, name::Symbol)
    name === :coefficients &&
        return (slope=getfield(common, :slope),
                intercept=getfield(common, :intercept))
    getfield(common, name)
end

function Base.propertynames(common::CommonPart, private::Bool=false)
    names = fieldnames(CommonPart)
    private ? (:coefficients, names...) : (:coefficients, names...)
end

function _region_summary(region::AsymptoticRegion)
    (name=region.name,
     source_status=region.source_status,
     variables=region.variables,
     transform=region.transform,
     truncation=region.truncation,
     mesh_points=haskey(region.mesh, :ξ) ? length(region.mesh.ξ) : 0,
     data_keys=Tuple(keys(region.data)))
end

function _common_part_summary(common::CommonPart)
    (name=common.name,
     expression=common.expression,
     source_status=common.source_status,
     slope=common.slope,
     intercept=common.intercept,
     ξ_match=common.ξ_match,
     fit_window=common.fit_window,
     mesh_points=haskey(common.mesh, :ξ) ? length(common.mesh.ξ) : 0,
     data_keys=Tuple(keys(common.data)))
end

function _composite_parts_summary(parts::CompositeParts)
    (inner=_region_summary(parts.inner),
     outer=_region_summary(parts.outer),
     common=_common_part_summary(parts.common),
     assumptions=parts.assumptions,
     diagnostics=parts.diagnostics)
end

function _inner_region(ξ, S, U)
    label = "composite-inner-region"
    AsymptoticRegion(
        :inner;
        variables=(coordinate=:ξ, radius=:S_inner, velocity=:U_inner),
        transform=(from=:inner_similarity_solution,
                   to=:composite_overlap_grid,
                   interpolation=:linear_strict),
        assumptions=_region_assumption(
            label,
            "Inner-region data are the current reconstructed inner solution interpolated onto the composite overlap grid."),
        truncation=(order=:leading_inner_profile,
                    common_coordinate=:ξ,
                    source_status=_DEFAULT_SOURCE_STATUS),
        source_ids=_region_source_ids(label),
        references=_canonical_reference_tuple(),
        provenance=_region_provenance(:composite_profile),
        mesh=(ξ=ξ,),
        data=(S=S, U=U),
        diagnostics=(mesh_points=length(ξ),
                     source_status=_DEFAULT_SOURCE_STATUS))
end

function _outer_region(ξ, S, U, ε)
    label = "composite-outer-region"
    AsymptoticRegion(
        :outer;
        variables=(coordinate=:ξ, radius=:S_outer,
                   perturbation=:s₁, velocity=:U_outer),
        transform=(from=:outer_solution,
                   to=:composite_overlap_grid,
                   reconstruction=:εξ_plus_s₁,
                   interpolation=:linear_strict),
        assumptions=_region_assumption(
            label,
            "Outer-region radius data use S_outer = εξ + s₁ on the current reconstructed outer solution."),
        truncation=(order=:leading_outer_profile,
                    common_coordinate=:ξ,
                    source_status=_DEFAULT_SOURCE_STATUS),
        source_ids=_region_source_ids(label),
        references=_canonical_reference_tuple(),
        provenance=merge(_region_provenance(:composite_profile), (ε=Float64(ε),)),
        mesh=(ξ=ξ,),
        data=(S=S, U=U),
        diagnostics=(mesh_points=length(ξ),
                     source_status=_DEFAULT_SOURCE_STATUS,
                     ε=Float64(ε)))
end

function _linear_common_part(ξ, S_common, U_common, fit;
                             ξ_min::Real, ξ_max::Real)
    label = "composite-common-part"
    ξ_match = Float64(fit.ξ_match)
    fit_window = (ξ_min=Float64(get(fit, :fit_ξ_min, ξ_match)),
                  ξ_max=Float64(get(fit, :fit_ξ_max, ξ_max)),
                  overlap_ξ_min=Float64(ξ_min),
                  overlap_ξ_max=Float64(ξ_max),
                  ξ_match=ξ_match,
                  fit_points=Int(fit.fit_points),
                  selection=:inner_ξ_greater_than_match)
    CommonPart(
        :linear_inner_far_field;
        expression=:slope_times_xi_plus_intercept,
        variables=(independent=:ξ, dependent=:S_common,
                   velocity=:U_common),
        transform=(from=:inner_far_field_fit,
                   to=:composite_overlap_grid,
                   representation=:linear_common_part),
        assumptions=_region_assumption(
            label,
            "The common part is a local linear fit to the current inner far field and is not source-confirmed against Decent-King 2008."),
        truncation=(order=:linear,
                    omitted_terms=:not_estimated,
                    source_status=_DEFAULT_SOURCE_STATUS),
        source_ids=_region_source_ids(label),
        references=_canonical_reference_tuple(),
        provenance=_region_provenance(:composite_profile),
        slope=fit.slope,
        intercept=fit.intercept,
        ξ_match=ξ_match,
        fit_window=fit_window,
        mesh=(ξ=ξ,),
        data=(S=S_common, U=U_common),
        diagnostics=(source_status=_DEFAULT_SOURCE_STATUS,
                     fit_points=Int(fit.fit_points),
                     sampled_values=length(ξ),
                     slope=Float64(fit.slope),
                     intercept=Float64(fit.intercept),
                     ξ_match=ξ_match))
end

function _composite_parts(; ξ, S_inner, U_inner, S_outer, U_outer,
                          S_common, U_common, ε, fit, ξ_min, ξ_max)
    inner = _inner_region(ξ, S_inner, U_inner)
    outer = _outer_region(ξ, S_outer, U_outer, ε)
    common = _linear_common_part(ξ, S_common, U_common, fit;
                                 ξ_min=ξ_min, ξ_max=ξ_max)
    assumptions = ("additive composite S_inner + S_outer - S_common",
                   "linear fitted common part is IMPL-inferred")
    diagnostics = (source_status=_DEFAULT_SOURCE_STATUS,
                   region_names=(inner.name, outer.name),
                   common_part=common.name,
                   mesh_points=length(ξ))
    CompositeParts(inner, outer, common, assumptions, diagnostics)
end

function _legacy_composite_parts(ξ, S, U, diagnostics::NamedTuple)
    ξ_min = Float64(get(diagnostics, :ξ_min, isempty(ξ) ? NaN : first(ξ)))
    ξ_max = Float64(get(diagnostics, :ξ_max, isempty(ξ) ? NaN : last(ξ)))
    ξ_match = Float64(get(diagnostics, :ξ_match, NaN))
    fit_points = Int(get(diagnostics, :fit_points, 0))
    fit = (slope=Float64(get(diagnostics, :overlap_slope, NaN)),
           intercept=Float64(get(diagnostics, :overlap_intercept, NaN)),
           ξ_match=ξ_match,
           fit_ξ_min=ξ_match,
           fit_ξ_max=ξ_max,
           fit_points=fit_points)
    empty_values = Float64[]
    _composite_parts(; ξ=ξ, S_inner=Float64.(S), U_inner=Float64.(U),
                     S_outer=empty_values, U_outer=empty_values,
                     S_common=empty_values, U_common=empty_values,
                     ε=NaN, fit=fit, ξ_min=ξ_min, ξ_max=ξ_max)
end

as_namedtuple(region::AsymptoticRegion) =
    (name=region.name,
     variables=region.variables,
     transform=region.transform,
     assumptions=as_namedtuple(region.assumptions),
     truncation=region.truncation,
     source_status=region.source_status,
     source_ids=as_namedtuple(region.source_ids),
     references=as_namedtuple(region.references),
     provenance=region.provenance,
     mesh=region.mesh,
     data=region.data,
     diagnostics=region.diagnostics)

as_namedtuple(common::CommonPart) =
    (name=common.name,
     expression=common.expression,
     variables=common.variables,
     transform=common.transform,
     assumptions=as_namedtuple(common.assumptions),
     truncation=common.truncation,
     source_status=common.source_status,
     source_ids=as_namedtuple(common.source_ids),
     references=as_namedtuple(common.references),
     provenance=common.provenance,
     slope=common.slope,
     intercept=common.intercept,
     ξ_match=common.ξ_match,
     fit_window=common.fit_window,
     mesh=common.mesh,
     data=common.data,
     diagnostics=common.diagnostics)

as_namedtuple(parts::CompositeParts) =
    (inner=as_namedtuple(parts.inner),
     outer=as_namedtuple(parts.outer),
     common=as_namedtuple(parts.common),
     assumptions=parts.assumptions,
     diagnostics=parts.diagnostics)
