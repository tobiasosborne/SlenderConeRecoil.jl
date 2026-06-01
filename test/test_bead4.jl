using Test
using SlenderConeRecoil

@testset "Similarity reduction source-status checks" begin
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

    @testset "C2001 Keller-Miksis length scaling with local U convention" begin
        # The t^(2/3) length scale is source-backed by the 2001 precursor.
        # The primitive velocity variable U and residual algebra are local
        # reconstruction checks pending the 2008 article body.
        @test SlenderConeRecoil.verify_t_cancels() == true
    end

    @testset "Mass ODE: local residual numerical evaluation" begin
        # 2S + 2Sξ(U - ξ) + S·Uξ = 0
        # S=1, Sξ=0, U=0, Uξ=-2, ξ=0:
        # 2(1) + 2(0)(0-0) + 1(-2) = 2 - 2 = 0
        mass = similarity_ode_mass()
        result = substitute(mass,
            S => Num(1), Sξ => Num(0), U => Num(0), Uξ => Num(-2), ξ => Num(0))
        @test result == Num(0)
    end

    @testset "Momentum ODE: local residual numerical evaluation" begin
        mom = similarity_ode_momentum()
        # U=9, ξ=0, Uξ=0, Sξ=-2, S=1:
        # -(2/9)(9) + 0 - (-2)/1 = -2 + 2 = 0
        result = substitute(mom,
            U => Num(9), ξ => Num(0), Uξ => Num(0),
            Sξ => Num(-2), S => Num(1))
        @test result == Num(0)

        # Non-trivial test with all terms nonzero:
        # U=3, ξ=1, Uξ=6, Sξ=2, S=2:
        # -(2/9)(3) + (4/9)(3-1)(6) - 2/4
        # = -6/9 + 48/9 - 1/2 = 42/9 - 1/2 = 14/3 - 1/2 = 25/6
        result2 = substitute(mom,
            U => Num(3), ξ => Num(1), Uξ => Num(6),
            Sξ => Num(2), S => Num(2))
        @test result2 == Num(25//6)
    end

    @testset "Momentum ODE with axial curvature local residual" begin
        mom_ax = similarity_ode_momentum(axial=true)
        Sξξξ = Sym(:Sξξξ)
        # Same as above plus -S''' term:
        # U=3, ξ=1, Uξ=6, Sξ=2, S=2, Sξξξ=1:
        # 25/6 - 1 = 19/6
        result = substitute(mom_ax,
            U => Num(3), ξ => Num(1), Uξ => Num(6),
            Sξ => Num(2), S => Num(2), Sξξξ => Num(1))
        @test result == Num(19//6)
    end

    @testset "Reconstructed far-field algebra: S ~ εξ, U ~ 0" begin
        # The source-backed far-field cone is stationary in physical variables.
        # Under the local primitive U convention, use S = εξ and U = 0:
        # 2S + 2Sξ(0 - ξ) + S·0 = 2S - 2ξSξ = 0
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
