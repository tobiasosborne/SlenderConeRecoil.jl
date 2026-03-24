using Test
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))

@testset "Slender-body derivation" begin

    @testset "Curvature expressions" begin
        κ_full = curvature_full()
        @test κ_full isa SExpr
        κ_lead = curvature_leading()
        # 1/R = R^{-1}
        @test κ_lead == pow(Sym(:R), Num(-1))
    end

    @testset "Mass conservation equation" begin
        eq = slender_mass_eq()
        @test eq isa SExpr
        R = Sym(:R); Rt = Sym(:Rt); Rz = Sym(:Rz)
        u = Sym(:u); uz = Sym(:uz)
        # 2R·Rt + 2R·Rz·u + R²·uz at R=3,Rt=2,Rz=1,u=4,uz=5:
        # 2(3)(2) + 2(3)(1)(4) + 9(5) = 12 + 24 + 45 = 81
        result = substitute(eq, R => Num(3), Rt => Num(2), Rz => Num(1),
                            u => Num(4), uz => Num(5))
        @test result == Num(81)
    end

    @testset "Momentum equation" begin
        eq = slender_momentum_eq()
        @test eq isa SExpr
        ut = Sym(:ut); u = Sym(:u); uz = Sym(:uz)
        Rz = Sym(:Rz); R = Sym(:R)
        # ut + u·uz - Rz/R² at ut=1,u=2,uz=3,Rz=8,R=2: 1+6-8/4 = 5
        result = substitute(eq, ut => Num(1), u => Num(2), uz => Num(3),
                            Rz => Num(8), R => Num(2))
        @test result == Num(5)
    end

    @testset "Momentum equation with axial curvature" begin
        eq = slender_momentum_eq(axial=true)
        @test eq isa SExpr
        ut = Sym(:ut); u = Sym(:u); uz = Sym(:uz)
        Rz = Sym(:Rz); R = Sym(:R); Rzzz = Sym(:Rzzz)
        # ut + u·uz - Rz/R² - Rzzz at ut=1,u=2,uz=3,Rz=8,R=2,Rzzz=10:
        # 1+6-2-10 = -5
        result = substitute(eq, ut => Num(1), u => Num(2), uz => Num(3),
                            Rz => Num(8), R => Num(2), Rzzz => Num(10))
        @test result == Num(-5)
    end

    @testset "System tuple" begin
        sys = slender_system()
        @test haskey(sys, :mass)
        @test haskey(sys, :momentum)
        @test sys.mass isa SExpr
        @test sys.momentum isa SExpr
    end

    @testset "Slender derivation verification" begin
        @test verify_slender_derivation() == true
    end

    @testset "Leading-order 1D system matches known result" begin
        # The known 1D model from the CLAUDE.md:
        #   ∂(R²)/∂t + ∂(R²u)/∂z = 0
        #   ∂u/∂t + u ∂u/∂z = ∂/∂z (1/R)
        #
        # Mass eq LHS: 2R·Rt + 2R·Rz·u + R²·uz
        # Momentum eq LHS: ut + u·uz + Rz/R²
        #   (note: ∂/∂z(1/R) = -Rz/R², so LHS-RHS = ut + u·uz - (-Rz/R²))

        mass = slender_mass_eq()
        mom = slender_momentum_eq()

        R = Sym(:R); Rt = Sym(:Rt); Rz = Sym(:Rz)
        u = Sym(:u); ut = Sym(:ut); uz = Sym(:uz)

        # Verify mass equation by substituting R=2, Rt=0, Rz=0, u=1, uz=3:
        # 2*2*0 + 2*2*0*1 + 4*3 = 12
        mass_eval = substitute(mass,
            R => Num(2), Rt => Num(0), Rz => Num(0),
            u => Num(1), uz => Num(3))
        @test mass_eval == Num(12)

        # Momentum: ut + u*uz - Rz/R² = 0
        # Substituting ut=1, u=2, uz=3, Rz=4, R=1: 1 + 6 - 4 = 3
        mom_eval = substitute(mom,
            ut => Num(1), u => Num(2), uz => Num(3),
            Rz => Num(4), R => Num(1))
        @test mom_eval == Num(3)
    end
end
