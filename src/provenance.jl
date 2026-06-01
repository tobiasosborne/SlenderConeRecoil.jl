# Source-fidelity and reproducibility metadata.
#
# These types intentionally use only Base types so provenance can be attached to
# solver results, benchmark records, and generated artifacts without pulling
# plotting, artifact, or paper-fetch tooling into the core package.

export SourceCitation, SourceID, Assumption, BenchmarkID,
       SolverSettings, ArtifactMetadata, PackageMetadata, ProvenanceMetadata,
       as_namedtuple, package_metadata, default_recoil_provenance_metadata

const _SOURCE_LEDGER_PATH =
    "docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md"
const _CANONICAL_CONE_DOI = "10.1093/imamat/hxm043"
const _DEFAULT_SOURCE_STATUS = "IMPL-inferred"

_optional_string(x) = x === nothing ? nothing : string(x)

function _string_tuple(xs)
    xs === nothing && return ()
    (xs isa AbstractString || xs isa Symbol) && return (string(xs),)
    Tuple(string(x) for x in xs)
end

function _version_or_nothing(v)
    v === nothing && return nothing
    v isa VersionNumber && return v
    try
        VersionNumber(string(v))
    catch err
        throw(ArgumentError("version must be parseable as a VersionNumber; got $(repr(v)): $err"))
    end
end

"""
    SourceCitation(id; doi=nothing, title="", authors=(), year=nothing,
                   status="unverified", ledger_path=nothing)

Bibliographic citation metadata for a source used by a result, benchmark, or
artifact. `status` should use the source-fidelity labels from the source
ledger, for example `C2008-meta`, `C2001`, `IMPL-inferred`, or `BLOCKED-2008`.
"""
struct SourceCitation
    id::String
    doi::Union{Nothing,String}
    title::String
    authors::Tuple{Vararg{String}}
    year::Union{Nothing,Int}
    status::String
    ledger_path::Union{Nothing,String}
end

function SourceCitation(id; doi=nothing, title::AbstractString="",
                        authors=(), year=nothing,
                        status::AbstractString="unverified",
                        ledger_path=nothing)
    year_value = year === nothing ? nothing : Int(year)
    SourceCitation(string(id), _optional_string(doi), string(title),
                   _string_tuple(authors), year_value, string(status),
                   _optional_string(ledger_path))
end

"""
    SourceID(id; kind=:source, description="", citation_id=nothing,
             status="IMPL-inferred", ledger_path=nothing)

Stable identifier for a source-backed or implementation-inferred object such as
an equation, figure, table, benchmark record, or implementation convention.
"""
struct SourceID
    id::String
    kind::Symbol
    description::String
    citation_id::Union{Nothing,String}
    status::String
    ledger_path::Union{Nothing,String}
end

function SourceID(id; kind::Symbol=:source, description::AbstractString="",
                  citation_id=nothing,
                  status::AbstractString=_DEFAULT_SOURCE_STATUS,
                  ledger_path=nothing)
    SourceID(string(id), kind, string(description), _optional_string(citation_id),
             string(status), _optional_string(ledger_path))
end

"""
    Assumption(id; description, status="IMPL-inferred", source_ids=())

Assumption metadata attached to a model, solve, benchmark, or generated
artifact.
"""
struct Assumption
    id::String
    description::String
    status::String
    source_ids::Tuple{Vararg{String}}
end

function Assumption(id; description::AbstractString,
                    status::AbstractString=_DEFAULT_SOURCE_STATUS,
                    source_ids=())
    Assumption(string(id), string(description), string(status),
               _string_tuple(source_ids))
end

"""
    BenchmarkID(id; description="", source_ids=(), status="unverified",
                quantity=nothing)

Identifier for benchmark data or reference quantities. This records where a
benchmark came from without storing the benchmark value itself.
"""
struct BenchmarkID
    id::String
    description::String
    source_ids::Tuple{Vararg{String}}
    status::String
    quantity::Union{Nothing,String}
end

function BenchmarkID(id; description::AbstractString="", source_ids=(),
                     status::AbstractString="unverified", quantity=nothing)
    BenchmarkID(string(id), string(description), _string_tuple(source_ids),
                string(status), _optional_string(quantity))
end

"""
    SolverSettings(problem_kind; algorithm="", settings=(;))

Solver configuration metadata. `settings` is a `NamedTuple` so callers can emit
or inspect settings without depending on package-internal solver structs.
"""
struct SolverSettings
    problem_kind::Symbol
    algorithm::String
    settings::NamedTuple
end

function SolverSettings(problem_kind::Symbol; algorithm::AbstractString="",
                        settings::NamedTuple=(;))
    SolverSettings(problem_kind, string(algorithm), settings)
end

"""
    ArtifactMetadata(path; checksum=nothing, checksum_algorithm="sha256",
                     checksum_status=nothing, status="not_generated",
                     source_ids=())

Metadata for a generated or local artifact. The checksum is recorded when a
caller already has it; this core type does not compute hashes itself.
"""
struct ArtifactMetadata
    path::String
    checksum::Union{Nothing,String}
    checksum_algorithm::String
    checksum_status::String
    status::String
    source_ids::Tuple{Vararg{String}}
end

function ArtifactMetadata(path; checksum=nothing,
                          checksum_algorithm::AbstractString="sha256",
                          checksum_status=nothing,
                          status::AbstractString="not_generated",
                          source_ids=())
    checksum_value = _optional_string(checksum)
    checksum_status_value =
        checksum_status === nothing ? (checksum_value === nothing ? "unchecked" : "recorded") :
        string(checksum_status)
    ArtifactMetadata(string(path), checksum_value, string(checksum_algorithm),
                     checksum_status_value, string(status),
                     _string_tuple(source_ids))
end

"""
    PackageMetadata(; name="SlenderConeRecoil", version=nothing, commit=nothing,
                    commit_status="unavailable", dirty=nothing, root=nothing)

Package version and Git commit metadata. `commit_status` is `"available"` when
a commit hash was detected and `"unavailable"` otherwise.
"""
struct PackageMetadata
    name::String
    version::Union{Nothing,VersionNumber}
    commit::Union{Nothing,String}
    commit_status::String
    dirty::Union{Nothing,Bool}
    root::Union{Nothing,String}
end

function PackageMetadata(; name::AbstractString="SlenderConeRecoil",
                         version=nothing, commit=nothing,
                         commit_status::AbstractString="unavailable",
                         dirty::Union{Nothing,Bool}=nothing,
                         root=nothing)
    PackageMetadata(string(name), _version_or_nothing(version),
                    _optional_string(commit), string(commit_status), dirty,
                    _optional_string(root))
end

"""
    ProvenanceMetadata(; source_citations=(), source_ids=(), assumptions=(),
                       benchmark_ids=(), solver_settings=SolverSettings(:unknown),
                       artifacts=(), package=package_metadata())

Structured provenance payload attached under `result.provenance.metadata`.
"""
struct ProvenanceMetadata
    source_citations::Tuple{Vararg{SourceCitation}}
    source_ids::Tuple{Vararg{SourceID}}
    assumptions::Tuple{Vararg{Assumption}}
    benchmark_ids::Tuple{Vararg{BenchmarkID}}
    solver_settings::SolverSettings
    artifacts::Tuple{Vararg{ArtifactMetadata}}
    package::PackageMetadata
end

function _metadata_tuple(::Type{T}, xs, field::AbstractString) where {T}
    xs === nothing && return ()
    xs isa T && return (xs,)
    xs isa AbstractString &&
        throw(ArgumentError("$field must contain $(T) objects, not a string"))
    values = Tuple(xs)
    all(x -> x isa T, values) ||
        throw(ArgumentError("$field must contain only $(T) objects"))
    values
end

function ProvenanceMetadata(; source_citations=(), source_ids=(),
                            assumptions=(), benchmark_ids=(),
                            solver_settings::SolverSettings=SolverSettings(:unknown),
                            artifacts=(), package::PackageMetadata=package_metadata())
    ProvenanceMetadata(
        _metadata_tuple(SourceCitation, source_citations, "source_citations"),
        _metadata_tuple(SourceID, source_ids, "source_ids"),
        _metadata_tuple(Assumption, assumptions, "assumptions"),
        _metadata_tuple(BenchmarkID, benchmark_ids, "benchmark_ids"),
        solver_settings,
        _metadata_tuple(ArtifactMetadata, artifacts, "artifacts"),
        package)
end

as_namedtuple(x) = x
as_namedtuple(x::NamedTuple) = x
as_namedtuple(xs::Tuple) = map(as_namedtuple, xs)
as_namedtuple(c::SourceCitation) =
    (id=c.id, doi=c.doi, title=c.title, authors=c.authors, year=c.year,
     status=c.status, ledger_path=c.ledger_path)
as_namedtuple(id::SourceID) =
    (id=id.id, kind=id.kind, description=id.description,
     citation_id=id.citation_id, status=id.status, ledger_path=id.ledger_path)
as_namedtuple(a::Assumption) =
    (id=a.id, description=a.description, status=a.status,
     source_ids=a.source_ids)
as_namedtuple(b::BenchmarkID) =
    (id=b.id, description=b.description, source_ids=b.source_ids,
     status=b.status, quantity=b.quantity)
as_namedtuple(s::SolverSettings) =
    (problem_kind=s.problem_kind, algorithm=s.algorithm, settings=s.settings)
as_namedtuple(a::ArtifactMetadata) =
    (path=a.path, checksum=a.checksum,
     checksum_algorithm=a.checksum_algorithm,
     checksum_status=a.checksum_status, status=a.status,
     source_ids=a.source_ids)
as_namedtuple(p::PackageMetadata) =
    (name=p.name, version=p.version, commit=p.commit,
     commit_status=p.commit_status, dirty=p.dirty, root=p.root)
as_namedtuple(m::ProvenanceMetadata) =
    (source_citations=as_namedtuple(m.source_citations),
     source_ids=as_namedtuple(m.source_ids),
     assumptions=as_namedtuple(m.assumptions),
     benchmark_ids=as_namedtuple(m.benchmark_ids),
     solver_settings=as_namedtuple(m.solver_settings),
     artifacts=as_namedtuple(m.artifacts),
     package=as_namedtuple(m.package))

function _project_root()
    normpath(joinpath(@__DIR__, ".."))
end

function _project_version_from_file(project_file::AbstractString)
    isfile(project_file) || return nothing
    for line in eachline(project_file)
        m = match(Regex("^\\s*version\\s*=\\s*\"([^\"]+)\""), line)
        if m !== nothing
            return _version_or_nothing(m.captures[1])
        end
    end
    nothing
end

function _git_metadata(root::AbstractString)
    marker = joinpath(root, ".git")
    (isdir(marker) || isfile(marker)) ||
        return (commit=nothing, commit_status="unavailable", dirty=nothing)

    commit = try
        readchomp(`git -C $root rev-parse --verify HEAD`)
    catch
        nothing
    end
    status_text = try
        readchomp(`git -C $root status --porcelain`)
    catch
        nothing
    end
    (commit=commit,
     commit_status=commit === nothing ? "unavailable" : "available",
     dirty=status_text === nothing ? nothing : !isempty(status_text))
end

"""
    package_metadata(; root=<package root>, name="SlenderConeRecoil")

Detect package version and Git commit metadata when available. If `.git` is not
present or Git is unavailable, commit fields are left empty with
`commit_status == "unavailable"`.
"""
function package_metadata(; root::AbstractString=_project_root(),
                          name::AbstractString="SlenderConeRecoil")
    project_file = joinpath(root, "Project.toml")
    git = _git_metadata(root)
    PackageMetadata(; name=name,
                    version=_project_version_from_file(project_file),
                    commit=git.commit,
                    commit_status=git.commit_status,
                    dirty=git.dirty,
                    root=root)
end

function _canonical_cone_citation()
    SourceCitation("DK2008-cone";
                   doi=_CANONICAL_CONE_DOI,
                   title="Surface-tension-driven flow in a slender cone",
                   authors=("S. P. Decent", "A. C. King"),
                   year=2008,
                   status="C2008-meta",
                   ledger_path=_SOURCE_LEDGER_PATH)
end

function _default_source_ids(kind::Symbol)
    implementation_id = "IMPL-inferred-$(kind)"
    (SourceID("DK2008-canonical-metadata";
              kind=:source,
              description="Canonical cone-recoil paper metadata; article-body equation and numerical extraction remains blocked.",
              citation_id="DK2008-cone",
              status="C2008-meta",
              ledger_path=_SOURCE_LEDGER_PATH),
     SourceID(implementation_id;
              kind=:implementation,
              description="Current $(kind) equations, boundary data, matching data, and numerical constants are local reconstructed implementation data unless the source ledger states otherwise.",
              citation_id="DK2008-cone",
              status=_DEFAULT_SOURCE_STATUS,
              ledger_path=_SOURCE_LEDGER_PATH))
end

function _default_assumptions(kind::Symbol)
    (Assumption("current-$(kind)-local-reconstruction";
                description="Use the current primitive-variable implementation as local reconstructed model data pending Decent-King 2008 article-body verification.",
                status=_DEFAULT_SOURCE_STATUS,
                source_ids=("IMPL-inferred-$(kind)",)),)
end

function _default_solver_algorithm(kind::Symbol)
    kind === :cone_similarity &&
        return "Newton shooting with Rodas5P IVP integrations"
    kind === :inner_bvp_residual &&
        return "Rodas5P endpoint shooting residual"
    kind === :inner_bvp_collocation &&
        return "shooting-seeded midpoint collocation Newton solve"
    kind === :outer_matching &&
        return "matched outward Rodas5P integration of reconstructed linearised outer ODE"
    kind === :composite_profile &&
        return "additive composite profile with fitted linear common part"
    kind === :pde_verification &&
        return "method-of-lines finite-difference verification"
    "unspecified"
end

"""
    default_recoil_provenance_metadata(kind; solver_settings=(;),
                                       benchmark_ids=(), artifacts=())

Default source-fidelity metadata for the current cone-recoil implementation.
The default status is intentionally `IMPL-inferred` while the canonical 2008
article body remains unavailable.
"""
function default_recoil_provenance_metadata(kind::Symbol;
                                            solver_settings::NamedTuple=(;),
                                            benchmark_ids=(),
                                            artifacts=(),
                                            package::PackageMetadata=package_metadata())
    ProvenanceMetadata(
        source_citations=(_canonical_cone_citation(),),
        source_ids=_default_source_ids(kind),
        assumptions=_default_assumptions(kind),
        benchmark_ids=benchmark_ids,
        solver_settings=SolverSettings(kind;
                                       algorithm=_default_solver_algorithm(kind),
                                       settings=solver_settings),
        artifacts=artifacts,
        package=package)
end
