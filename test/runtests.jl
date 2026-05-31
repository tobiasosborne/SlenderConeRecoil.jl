using Test
using SlenderConeRecoil

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
    )

    exported = names(SlenderConeRecoil)
    for name in public_api
        @test isdefined(SlenderConeRecoil, name)
        @test name in exported
    end
end

@testset "SlenderConeRecoil" begin
    include("test_bead1.jl")
    include("test_bead2.jl")
    include("test_bead3.jl")
    include("test_bead4.jl")
    include("test_bead5.jl")
    include("test_bead6.jl")
    include("test_outer_hierarchy.jl")
    include("test_bead7.jl")
    include("test_bead8.jl")
end
