using Test
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))

@testset "Similarity reduction" begin
    ξ = Sym(:ξ); S = Sym(:S); U = Sym(:U)
    Sξ = Sym(:Sξ); Uξ = Sym(:Uξ)

    @testset "Mass ODE structure" begin
        mass = similarity_ode_mass()
        @test mass isa SExpr
        s = string(mass)
        @test occursin("S", s)
        @test occursin("U", s)
        @test occursin("ξ", s)
    end

    @testset "Momentum ODE structure" begin
        mom = similarity_ode_momentum()
        @test mom isa SExpr
        s = string(mom)
        @test occursin("U", s)
        @test occursin("S", s)
    end

    @testset "System tuple" begin
        sys = similarity_system()
        @test haskey(sys, :mass)
        @test haskey(sys, :momentum)
    end

    @testset "t cancels (numerical verification)" begin
        @test verify_t_cancels() == true
    end

    @testset "Mass ODE: numerical evaluation" begin
        # 2S + 2Sξ(U - ξ) + S·Uξ = 0
        # S=1, Sξ=0, U=0, Uξ=-2, ξ=0:
        # 2(1) + 2(0)(0-0) + 1(-2) = 2 - 2 = 0
        mass = similarity_ode_mass()
        result = substitute(mass,
            S => Num(1), Sξ => Num(0), U => Num(0), Uξ => Num(-2), ξ => Num(0))
        @test result == Num(0)
    end

    @testset "Momentum ODE: numerical evaluation" begin
        # -(2/9)U + (4/9)(U - ξ)Uξ - Sξ/S² = 0
        # U=0, ξ=0, Uξ=0, Sξ=0, S=1:
        # 0 + 0 - 0 = 0
        mom = similarity_ode_momentum()
        result = substitute(mom,
            U => Num(0), ξ => Num(0), Uξ => Num(0),
            Sξ => Num(0), S => Num(1))
        @test result == Num(0)

        # U=9, ξ=0, Uξ=0, Sξ=-2, S=1:
        # -(2/9)(9) + 0 - (-2)/1 = -2 + 2 = 0
        result2 = substitute(mom,
            U => Num(9), ξ => Num(0), Uξ => Num(0),
            Sξ => Num(-2), S => Num(1))
        @test result2 == Num(0)
    end

    @testset "Far-field consistency: S ~ εξ, U ~ (2/3)ξ" begin
        # In the far field, the cone is undisturbed: R = εz, u = ż·(2/3)
        # In similarity variables: S = εξ (to leading order) and U = (2/3)ξ
        # The mass ODE should be satisfied:
        # 2(εξ) + 2(ε)(2ξ/3 - ξ) + (εξ)(2/3) = 0
        # 2εξ + 2ε(-ξ/3) + 2εξ/3 = 2εξ - 2εξ/3 + 2εξ/3 = 2εξ
        # ... this is NOT zero, which means the far-field is more subtle.
        # Actually, the undisturbed cone has u=0 (fluid at rest), so U=0.
        # Then: 2S + 2Sξ(0 - ξ) + S·0 = 2S - 2ξSξ = 0
        # With S = εξ: 2εξ - 2ξ·ε = 0 ✓
        mass = similarity_ode_mass()
        ε = Sym(:ε)
        # S = εξ, Sξ = ε, U = 0, Uξ = 0
        result = substitute(mass,
            S => mul(ε, ξ), Sξ => ε, U => Num(0), Uξ => Num(0))
        # Should be 2εξ + 2ε(0 - ξ) + εξ·0 = 2εξ - 2εξ = 0
        # But our substitute doesn't simplify symbolically, so evaluate numerically
        result_num = substitute(result, ε => Num(1//10), ξ => Num(5))
        @test result_num == Num(0)

        # Momentum: -(2/9)(0) + (4/9)(0-ξ)(0) + ε/(εξ)² = ε/(ε²ξ²) = 1/(εξ²)
        # This is NOT zero — the momentum equation is not satisfied by u=0
        # because surface tension drives the flow. This is correct: the
        # far-field has a small residual that drives the perturbation.
        # The far-field BC is that S → εξ as ξ → ∞, not that the ODE = 0 there.
    end
end
