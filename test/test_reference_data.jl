using Test
using TOML

const REFERENCE_DATA_PATH = joinpath(@__DIR__, "reference", "decent_king_cone_reference.toml")
const DK2001_DOI = "10.1007/978-94-010-0796-2_10"
const DK2008_DOI = "10.1093/imamat/hxm043"

function _record_by_id(records)
    Dict(record["id"] => record for record in records)
end

function _has_nonempty_string(record, key)
    haskey(record, key) && record[key] isa AbstractString && !isempty(strip(record[key]))
end

@testset "Decent-King cone reference data" begin
    data = TOML.parsefile(REFERENCE_DATA_PATH)
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
    )
    valid_statuses = Set(["C2001", "C2008-meta", "BLOCKED-2008"])

    @test data["metadata"]["ledger_path"] ==
          "docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md"
    @test length(records) >= 8
    @test length(keys(by_id)) == length(records)

    @testset "Schema and source-status separation" begin
        for record in records
            for key in required_string_fields
                @test _has_nonempty_string(record, key)
            end
            @test record["status"] in valid_statuses

            if record["fact_type"] == "numeric"
                @test haskey(record, "value")
                @test haskey(record, "tolerance")
                @test record["value"] isa Number
                @test record["tolerance"] isa Number
                @test record["tolerance"] >= 0
                @test record["status"] != "BLOCKED-2008"
            end

            if record["status"] == "BLOCKED-2008"
                @test record["doi"] == DK2008_DOI
                @test record["fact_type"] == "blocked_quantitative"
                @test get(record, "no_numeric_benchmark", false) == true
                @test !haskey(record, "value")
                @test !haskey(record, "tolerance")
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
