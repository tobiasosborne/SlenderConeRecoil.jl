using Test
using TOML
using SHA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE_DATA_PATH = joinpath(@__DIR__, "reference", "decent_king_cone_reference.toml")
const REFERENCE_SCHEMA_PATH = joinpath(PROJECT_ROOT, "docs", "reference_data_schema.md")
const REFERENCE_DATA_ROOT = joinpath(PROJECT_ROOT, "data", "reference")
const GITIGNORE_PATH = joinpath(PROJECT_ROOT, ".gitignore")
const DK2001_DOI = "10.1007/978-94-010-0796-2_10"
const DK2008_DOI = "10.1093/imamat/hxm043"
const SHA256_TAG_RE = r"^sha256:[0-9a-f]{64}$"

function _record_by_id(records)
    Dict(record["id"] => record for record in records)
end

function _has_nonempty_string(record, key)
    haskey(record, key) && record[key] isa AbstractString && !isempty(strip(record[key]))
end

function _canonical_hash_value(value)
    if value isa AbstractDict
        keys_sorted = sort(collect(keys(value)))
        return "{" * join(("$(key)=" * _canonical_hash_value(value[key])
                            for key in keys_sorted), ";") * "}"
    elseif value isa AbstractVector
        return "[" * join((_canonical_hash_value(item) for item in value), ";") * "]"
    else
        return repr(value)
    end
end

function _reference_record_payload(record)
    keys_sorted = filter(!=("data_hash"), sort(collect(keys(record))))
    join(("$(key)=" * _canonical_hash_value(record[key])
          for key in keys_sorted), "\n")
end

function _reference_record_hash(record)
    payload = _reference_record_payload(record)
    "sha256:" * bytes2hex(sha256(Vector{UInt8}(codeunits(payload))))
end

_is_sha256_tag(value) = value isa AbstractString && occursin(SHA256_TAG_RE, value)

@testset "Decent-King cone reference data" begin
    data = TOML.parsefile(REFERENCE_DATA_PATH)
    dataset_metadata = data["metadata"]
    records = data["records"]
    by_id = _record_by_id(records)

    required_string_fields = (
        "id",
        "source",
        "status",
        "fact_type",
        "quantity",
        "doi",
        "page_equation_figure_reference",
        "extraction_method",
        "units_scaling",
        "tolerance_rationale",
        "source_status_note",
        "local_source_artifact_path",
        "local_source_checksum_algorithm",
        "local_source_checksum",
        "local_source_checksum_status",
        "source_license_access_status",
        "recorded_date",
        "data_hash",
    )
    required_table_fields = ("parameters", "solver_settings", "tolerances")
    valid_statuses = Set(["C2001", "C2008-meta", "BLOCKED-2008"])
    valid_checksum_statuses = Set(["recorded", "missing", "blocked", "not_applicable"])

    @test dataset_metadata["ledger_path"] ==
          "docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md"
    @test dataset_metadata["schema_path"] == "docs/reference_data_schema.md"
    @test dataset_metadata["canonical_reusable_data_location"] == "data/reference/"
    @test dataset_metadata["package_name"] == "SlenderConeRecoil"
    @test dataset_metadata["package_version"] == "0.1.0"
    @test length(records) >= 8
    @test length(keys(by_id)) == length(records)

    @testset "Repository schema and artifact policy" begin
        @test isfile(REFERENCE_SCHEMA_PATH)
        @test isdir(REFERENCE_DATA_ROOT)
        @test isfile(joinpath(REFERENCE_DATA_ROOT, "README.md"))

        schema_text = read(REFERENCE_SCHEMA_PATH, String)
        schema_terms = (
            "doi",
            "local_source_checksum",
            "extraction_method",
            "parameters",
            "solver_settings",
            "tolerances",
            "data_hash",
            "Artifacts.toml",
        )
        for term in schema_terms
            @test occursin(term, schema_text)
        end
        @test occursin("licensed pdfs", lowercase(schema_text))

        ignore_text = read(GITIGNORE_PATH, String)
        @test occursin("docs/papers/*.pdf", ignore_text)
        @test occursin("docs/papers/*.png", ignore_text)
    end

    @testset "Schema and source-status separation" begin
        for record in records
            for key in required_string_fields
                @test _has_nonempty_string(record, key)
            end
            for key in required_table_fields
                @test haskey(record, key)
                @test record[key] isa AbstractDict
                @test !isempty(record[key])
            end
            @test record["status"] in valid_statuses
            @test record["local_source_checksum_algorithm"] == "sha256"
            @test record["local_source_checksum_status"] in valid_checksum_statuses
            if record["local_source_checksum_status"] == "recorded"
                @test _is_sha256_tag(record["local_source_checksum"])
            else
                @test record["local_source_checksum"] in
                      ("pending-acquisition", "not-applicable")
            end
            @test _is_sha256_tag(record["data_hash"])
            @test record["data_hash"] == _reference_record_hash(record)

            if record["fact_type"] == "numeric"
                @test haskey(record, "value")
                @test haskey(record, "tolerance")
                @test record["value"] isa Number
                @test record["tolerance"] isa Number
                @test record["tolerance"] >= 0
                @test record["status"] != "BLOCKED-2008"
                @test record["tolerances"]["absolute"] == record["tolerance"]
                @test haskey(record["tolerances"], "rationale")
            end

            if record["status"] == "BLOCKED-2008"
                @test record["doi"] == DK2008_DOI
                @test record["fact_type"] == "blocked_quantitative"
                @test get(record, "no_numeric_benchmark", false) == true
                @test !haskey(record, "value")
                @test !haskey(record, "tolerance")
                @test record["parameters"]["status"] == "blocked"
                @test record["solver_settings"]["status"] == "not_applicable"
                @test record["tolerances"]["status"] == "unavailable"
            end
        end

        source_backed_values = filter(record -> haskey(record, "value"), records)
        @test all(record["doi"] == DK2001_DOI for record in source_backed_values)
        @test all(record["status"] == "C2001" for record in source_backed_values)
        @test !any(record -> occursin("xi0", record["quantity"]) ||
                            occursin("S0", record["quantity"]) ||
                            occursin("S''0", record["source_status_note"]),
                  source_backed_values)
    end

    @testset "2001 precursor numeric facts" begin
        figure_epsilon = by_id["DK2001-figure1-epsilon"]
        @test figure_epsilon["doi"] == DK2001_DOI
        @test figure_epsilon["value"] == 0.01
        @test figure_epsilon["tolerance"] == 0.0
        @test occursin("Figure 1", figure_epsilon["page_equation_figure_reference"])
        @test occursin("Dimensionless", figure_epsilon["units_scaling"])

        y0_limit = by_id["DK2001-direct-shooting-y0-limit"]
        @test y0_limit["value"] ≈ 1.50 atol=y0_limit["tolerance"]
        @test y0_limit["tolerance"] == 0.005
        @test "OCR-needs-check" in y0_limit["caveats"]
        @test occursin("Do not compare directly", y0_limit["source_status_note"])

        p0_limit = by_id["DK2001-direct-shooting-p0-limit"]
        @test p0_limit["value"] ≈ 2.65 atol=p0_limit["tolerance"]
        @test p0_limit["tolerance"] == 0.005
        @test "OCR-needs-check" in p0_limit["caveats"]
        @test occursin("Do not compare directly", p0_limit["source_status_note"])
    end

    @testset "2001 qualitative wave fact" begin
        oscillations = by_id["DK2001-figure1-oscillation-qualitative"]
        @test oscillations["doi"] == DK2001_DOI
        @test !haskey(oscillations, "value")
        @test occursin("high-frequency oscillations", oscillations["text_value"])
        @test occursin("longer scale", oscillations["text_value"])
        @test occursin("Figure 1", oscillations["page_equation_figure_reference"])
    end

    @testset "2008 article-body quantities remain blocked" begin
        metadata = by_id["DK2008-canonical-metadata"]
        @test metadata["doi"] == DK2008_DOI
        @test metadata["status"] == "C2008-meta"
        @test !haskey(metadata, "value")
        @test occursin("article-body numerical values remain blocked",
                       lowercase(metadata["source_status_note"]))

        blocked_records = filter(record -> record["status"] == "BLOCKED-2008", records)
        @test length(blocked_records) >= 3
        @test all(record -> record["doi"] == DK2008_DOI, blocked_records)
        @test all(record -> !haskey(record, "value"), blocked_records)
        @test all(record -> !haskey(record, "tolerance"), blocked_records)
        @test any(record -> record["quantity"] == "similarity_tip_constants", blocked_records)
        @test any(record -> record["quantity"] == "inner_outer_matching_constants", blocked_records)
        @test any(record -> record["quantity"] == "figure_derived_quantities", blocked_records)
    end
end
