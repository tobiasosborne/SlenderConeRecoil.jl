using Test
using SlenderConeRecoil

@testset "Core provenance metadata structs are inspectable" begin
    citation = SourceCitation("DK2008";
                              doi="10.1093/imamat/hxm043",
                              title="Surface-tension-driven flow in a slender cone",
                              authors=("S. P. Decent", "A. C. King"),
                              year=2008,
                              status="C2008-meta",
                              ledger_path="ledger.md")
    source_id = SourceID("DK2008-eq-blocked";
                         kind=:equation,
                         description="Equation extraction blocked pending article body",
                         citation_id=citation.id,
                         status="BLOCKED-2008",
                         ledger_path="ledger.md")
    assumption = Assumption("local-reconstruction";
                            description="Use current primitive-variable solver as local reconstruction",
                            source_ids=(source_id.id,))
    benchmark = BenchmarkID("DK2008-tip-constant";
                            description="Placeholder for future source-backed tip constant",
                            source_ids=(source_id.id,),
                            status="BLOCKED-2008",
                            quantity="tip constant")
    solver = SolverSettings(:benchmark;
                            algorithm="fixture",
                            settings=(ε=0.1, maxiters=42))
    artifact = ArtifactMetadata("docs/papers/DecentKing2008.pdf";
                                checksum="abc123",
                                checksum_status="recorded",
                                status="missing",
                                source_ids=(source_id.id,))
    package = PackageMetadata(; version=v"0.1.0", commit="abc123",
                              commit_status="available", dirty=false)

    metadata = ProvenanceMetadata(; source_citations=(citation,),
                                  source_ids=(source_id,),
                                  assumptions=(assumption,),
                                  benchmark_ids=(benchmark,),
                                  solver_settings=solver,
                                  artifacts=(artifact,),
                                  package=package)

    nt = as_namedtuple(metadata)
    @test nt.source_citations[1].doi == "10.1093/imamat/hxm043"
    @test nt.source_ids[1].kind == :equation
    @test nt.assumptions[1].source_ids == ("DK2008-eq-blocked",)
    @test nt.benchmark_ids[1].status == "BLOCKED-2008"
    @test nt.solver_settings.settings.maxiters == 42
    @test nt.artifacts[1].checksum_status == "recorded"
    @test nt.package.commit == "abc123"
end

@testset "Default problem provenance stays honest and includes solver settings" begin
    problem = PDEVerificationProblem(; ε=0.2, N=7, z_min=0.01, z_max=0.5,
                                     t_end=0.0, n_snapshots=1,
                                     maxiters=123)

    @test problem.provenance.source_status == "IMPL-inferred"
    @test problem.provenance.canonical_source_doi == "10.1093/imamat/hxm043"
    @test problem.provenance.source_ledger ==
          "docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md"

    metadata = problem.provenance.metadata
    @test metadata isa ProvenanceMetadata
    @test metadata.source_citations[1].doi == "10.1093/imamat/hxm043"
    @test metadata.source_citations[1].status == "C2008-meta"
    @test any(id -> id.status == "IMPL-inferred", metadata.source_ids)
    @test metadata.assumptions[1].status == "IMPL-inferred"
    @test metadata.solver_settings.problem_kind == :pde_verification
    @test occursin("method-of-lines", metadata.solver_settings.algorithm)
    @test metadata.solver_settings.settings.ε == 0.2
    @test metadata.solver_settings.settings.N == 7
    @test metadata.solver_settings.settings.maxiters == 123
    @test metadata.package.name == "SlenderConeRecoil"
    @test metadata.package.version == v"0.1.0"
    @test metadata.package.commit_status in ("available", "unavailable")

    @test problem.provenance.solver_settings === metadata.solver_settings
    @test problem.provenance.package === metadata.package
end

@testset "Problem override provenance reports effective solver settings" begin
    problem = PDEVerificationProblem(; ε=0.1, N=3, z_min=0.01, z_max=0.1,
                                     t_end=0.0, n_snapshots=1,
                                     maxiters=10,
                                     provenance=(custom_tag=:kept,))
    result = solve_pde(problem; ε=0.3, N=4, z_max=0.2, maxiters=11)

    @test result.provenance.custom_tag == :kept
    @test result.provenance.user_provenance.custom_tag == :kept
    @test result.problem.parameters.ε == 0.3
    @test result.provenance.metadata.solver_settings.settings.ε == 0.3
    @test result.provenance.metadata.solver_settings.settings.N == 4
    @test result.provenance.metadata.solver_settings.settings.z_max == 0.2
    @test result.provenance.metadata.solver_settings.settings.maxiters == 11
end

@testset "Package metadata is robust without git metadata" begin
    mktempdir() do dir
        meta = package_metadata(; root=dir)
        @test meta.name == "SlenderConeRecoil"
        @test meta.version === nothing
        @test meta.commit === nothing
        @test meta.commit_status == "unavailable"
        @test meta.dirty === nothing
    end
end
