using Test
using SlenderConeRecoil

@testset "Outer hierarchy CAS (local reconstructed S,U algebra)" begin
    @testset "derive_outer_equations expands Laurent coefficients cleanly" begin
        # Internal consistency only. These coefficients are not source-backed
        # Decent-King outer hierarchy equations.
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

    @testset "hierarchy grammar is explicit integer/Laurent local algebra" begin
        ε = Sym(:ε)
        eqs = derive_outer_equations(order=5)
        rendered = sprint(show, eqs)

        @test all(k -> k isa Int, keys(eqs))
        @test !occursin("1//2", rendered)
        @test !occursin("sqrt", rendered)
        @test_throws ArgumentError derive_outer_equations(order=-1)

        err = try
            expand_in(pow(ε, Num(1//2)), ε, 3)
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin("integer", sprint(showerror, err))
        @test occursin("Laurent", sprint(showerror, err))
    end
end
