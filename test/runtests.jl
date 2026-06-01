using Test
using SlenderConeRecoil

const TEST_GROUP = lowercase(get(ENV, "SLENDER_RECOIL_TEST_GROUP", "fast"))
const FAST_TEST_FILES = (
    "test_public_api_wrappers.jl",
    "test_provenance_metadata.jl",
    "test_bead1.jl",
    "test_bead2.jl",
    "test_bead3.jl",
    "test_bead4.jl",
    "test_outer_hierarchy.jl",
    "test_reference_data.jl",
)
const SLOW_TEST_FILES = (
    "test_bead5.jl",
    "test_bead6.jl",
    "test_bead7.jl",
    "test_bead8.jl",
    "test_numerical_regressions.jl",
)
const VALID_TEST_GROUPS = ("fast", "slow", "all")

TEST_GROUP in VALID_TEST_GROUPS || error(
    "Invalid SLENDER_RECOIL_TEST_GROUP=$(repr(TEST_GROUP)); expected one of " *
    join(VALID_TEST_GROUPS, ", ")
)

function include_tests(files)
    for file in files
        include(joinpath(@__DIR__, file))
    end
end

@testset "Package load and public API" begin
    @test Sym(:x) isa SExpr
    @test solve_inner_bvp === SlenderConeRecoil.solve_inner_bvp
    @test solve_outer_matched === SlenderConeRecoil.solve_outer_matched
    @test solve_outer_linearised === SlenderConeRecoil.solve_outer_linearised

    public_api = (
        :SExpr, :Sym, :Num, :Add, :Mul, :Pow, :Func,
        :substitute, :differentiate, :walk, :add, :mul, :pow, :neg,
        :expand_in, :collect_order,
        :slender_mass_eq, :slender_momentum_eq, :slender_system,
        :similarity_ode_mass, :similarity_ode_momentum, :similarity_system,
        :solve_inner_bvp, :InnerSolution,
        :solve_outer, :solve_outer_driven, :solve_outer_matched, :OuterSolution,
        :composite_solution, :CompositeSolution, :overlap_residual,
        :derive_outer_equations, :eval_sexpr, :solve_outer_full,
        :solve_outer_linearised, :HierarchySolution,
        :solve_pde, :PDESolution, :rescale_to_similarity,
        :AbstractRecoilProblem,
        :ConeSimilarityProblem, :OuterMatchingProblem,
        :CompositeProfileProblem, :PDEVerificationProblem,
        :ProblemResult,
        :ConeSimilarityResult, :OuterMatchingResult,
        :CompositeProfileResult, :PDEVerificationResult,
        :SourceCitation, :SourceID, :Assumption, :BenchmarkID,
        :SolverSettings, :ArtifactMetadata, :PackageMetadata,
        :ProvenanceMetadata, :as_namedtuple, :package_metadata,
        :default_recoil_provenance_metadata,
    )

    exported = names(SlenderConeRecoil)
    for name in public_api
        @test isdefined(SlenderConeRecoil, name)
        @test name in exported
    end
end

@testset "SlenderConeRecoil $(TEST_GROUP) gate" begin
    if TEST_GROUP == "fast"
        include_tests(FAST_TEST_FILES)
    elseif TEST_GROUP == "slow"
        include_tests(SLOW_TEST_FILES)
    else
        include_tests(FAST_TEST_FILES)
        include_tests(SLOW_TEST_FILES)
    end
end
