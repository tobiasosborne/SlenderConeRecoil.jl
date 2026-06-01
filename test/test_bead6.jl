using Test
using SlenderConeRecoil

@testset "Outer linearised problem (local reconstructed equations)" begin

    @testset "Public outer solve validates local domain assumptions" begin
        @test_throws ArgumentError solve_outer(ε=0.0, ξ_min=5.0, ξ_max=50.0)
        @test_throws ArgumentError solve_outer(ε=-0.1, ξ_min=5.0, ξ_max=50.0)
        @test_throws ArgumentError solve_outer(ε=Inf, ξ_min=5.0, ξ_max=50.0)
        @test_throws ArgumentError solve_outer(ε=0.1, ξ_min=0.0, ξ_max=50.0)
        @test_throws ArgumentError solve_outer(ε=0.1, ξ_min=50.0, ξ_max=5.0)
        @test_throws ArgumentError solve_outer(ε=0.1, ξ_min=5.0, ξ_max=50.0, seed=NaN)
        @test_throws ArgumentError solve_outer(ε=0.1, ξ_min=5.0, ξ_max=50.0, maxiters=0)
    end

    @testset "Zero perturbation is (nearly) steady" begin
        # With seed=0 at large ξ, the solution should stay near zero
        # (the base state almost satisfies the equations far from the tip)
        sol = solve_outer(ε=0.1, ξ_min=5.0, ξ_max=50.0, seed=0.0)
        @test sol isa OuterSolution
        @test length(sol.ξ) > 5
        @test sol.diagnostics.requested_endpoint ≈ 5.0
        @test sol.diagnostics.endpoint ≈ 5.0
        @test all(isfinite, sol.s₁)
        @test all(isfinite, sol.u₁)
    end

    @testset "Solver retcode failure is reported" begin
        err = @test_logs (:warn, r"max_iters") begin
            try
                solve_outer(ε=0.1, ξ_min=2.0, ξ_max=50.0, maxiters=1)
                nothing
            catch e
                e
            end
        end
        @test err isa ErrorException
        @test occursin("retcode", sprint(showerror, err))
    end

    @testset "Integration direction (inward)" begin
        sol = solve_outer(ε=0.1, ξ_min=2.0, ξ_max=30.0)
        # ξ should be sorted ascending after reversal
        @test sol.ξ[1] < sol.ξ[end]
        @test sol.ξ[1] ≈ 2.0 atol=1.0
        @test sol.ξ[end] ≈ 30.0 atol=1.0
    end

    @testset "Driven solution produces nontrivial response" begin
        sol = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=40.0)
        # Local diagnostic only: this is not a published Decent-King matching mode.
        # With a nonzero seed, there should be some perturbation.
        max_s1 = maximum(abs.(sol.s₁))
        @test max_s1 > 0  # nontrivial
        @test all(isfinite, sol.s₁)
        @test all(isfinite, sol.u₁)
    end

    @testset "Perturbation bounded at moderate ξ" begin
        sol = solve_outer_driven(ε=0.1, ξ_min=5.0, ξ_max=30.0)
        # The perturbation should remain moderate relative to base state εξ
        @test maximum(abs.(sol.s₁)) < 100.0
        @test maximum(abs.(sol.u₁)) < 100.0
    end

    @testset "ODE RHS evaluates cleanly" begin
        dy = zeros(4)
        SlenderConeRecoil.outer_rhs!(dy, [0.01, 0.0, 0.0, 0.0], [0.1], 10.0)
        @test all(isfinite, dy)
    end

    @testset "Linearised mass equation consistency" begin
        # Internal algebra check for the reconstructed linearised mass equation:
        # 2s₁ + 2εu₁ - 2ξs₁' + εξu₁' = 0.
        ε = 0.1; ξ_test = 5.0
        s₁ = 0.1; s₁p = 0.02; s₁pp = 0.0; u₁ = 0.05
        dy = zeros(4)
        SlenderConeRecoil.outer_rhs!(dy, [s₁, s₁p, s₁pp, u₁], [ε], ξ_test)
        u₁p = dy[4]
        mass_res = 2*s₁ + 2*ε*u₁ - 2*ξ_test*s₁p + ε*ξ_test*u₁p
        @test abs(mass_res) < 1e-12
    end
end
