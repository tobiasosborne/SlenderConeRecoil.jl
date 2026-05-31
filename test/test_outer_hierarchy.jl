using Test
using SlenderConeRecoil

@testset "Outer hierarchy CAS" begin
    @testset "derive_outer_equations expands Laurent coefficients cleanly" begin
        eqs = derive_outer_equations(order=5)
        rendered = sprint(show, eqs)

        @test Set(keys(eqs)) == Set([-1, 1, 2, 3, 4, 5])
        @test !occursin("1//0", rendered)
        @test !occursin("Inf", rendered)
        @test !occursin("NaN", rendered)
        @test !occursin("ε", rendered)
        @test eqs[-1].mass == Num(0)
        @test eval_sexpr(eqs[-1].momentum, Dict(:ξ => 2.0)) ≈ -0.25
    end
end
