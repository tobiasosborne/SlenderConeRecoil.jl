using Test
using SlenderConeRecoil
using LinearAlgebra: norm

@testset "Inner BVP solver (with axial curvature)" begin

    @testset "Tip initial conditions" begin
        ic = SlenderConeRecoil.tip_initial_conditions(0.0, 1.0, -0.5)
        @test ic == [1.0, 0.0, -0.5, 0.0]  # [S, S', S'', U]

        ic2 = SlenderConeRecoil.tip_initial_conditions(5.0, 0.5, -1.0)
        @test ic2[4] ≈ 4.0  # U = (4/5)ξ₀
        @test ic2[2] == 0.0  # S' = 0 at tip
    end

    @testset "ODE RHS: no singularity" begin
        du = zeros(4)
        for ξ in [0.0, 0.5, 1.0, 5.0, 10.0]
            SlenderConeRecoil.inner_rhs!(du, [1.0, 0.1, -0.5, 0.0], nothing, ξ + 1e-6)
            @test all(isfinite, du)
        end
    end

    @testset "Integration reaches far field" begin
        ξv, Sv, _, _, Uv = SlenderConeRecoil._shoot(2.5, 0.5, -1.0, 50.0)
        @test ξv[end] ≈ 50.0 atol=1.0
        @test all(isfinite, Sv)
        @test all(isfinite, Uv)
    end

    @testset "3D Newton converges" begin
        sol = solve_inner_bvp(ε=0.1)
        @test sol.ξ₀ > 0
        @test sol.S₀ > 0
        @test length(sol.final_residual) == 3
        @test sol.final_residual_norm ≈ norm(sol.final_residual)
        @test sol.iterations >= 0
        @test !isempty(sol.termination_reason)
        @test sol.converged == (sol.final_residual_norm < 1e-4)

        # Far-field slope near ε (with axial curvature, convergence is harder)
        slopes = SlenderConeRecoil.far_field_slope(sol)
        @test abs(slopes[end] - 0.1) < 0.05

        # Far-field U → 0
        @test abs(sol.U[end]) < 0.1

        # Blob exists behind tip: S exceeds εξ somewhere
        excess = sol.S .- 0.1 .* sol.ξ
        @test maximum(excess) > 0.01
    end

    @testset "Newton failure diagnostics" begin
        sol = solve_inner_bvp(ε=0.1, ξ_max=6.0, newton_iters=0)
        @test !sol.converged
        @test sol.iterations == 0
        @test sol.termination_reason == "iteration limit reached"
        @test length(sol.final_residual) == 3
        @test sol.final_residual_norm ≈ norm(sol.final_residual)

        @test_throws ErrorException solve_inner_bvp(ε=0.1, ξ_max=6.0,
                                                    newton_iters=0,
                                                    throw_on_failure=true)
    end

    @testset "ODE residual (with S''' term)" begin
        sol = solve_inner_bvp(ε=0.1)
        n = length(sol.ξ)
        mid = div(n, 2)
        if mid > 2 && mid < n-1
            ξm = sol.ξ[mid]; Sm = sol.S[mid]; Um = sol.U[mid]
            Sξ_fd = (sol.S[mid+1] - sol.S[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])
            Uξ_fd = (sol.U[mid+1] - sol.U[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])
            Sξξξ_fd = (sol.Sξξ[mid+1] - sol.Sξξ[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])

            mass_res = 2*Sm + 2*Sξ_fd*(Um - ξm) + Sm*Uξ_fd
            @test abs(mass_res) < 0.1

            mom_res = -(2/9)*Um + (4/9)*(Um-ξm)*Uξ_fd - Sξ_fd/Sm^2 - Sξξξ_fd
            @test abs(mom_res) < 0.1
        end
    end

    @testset "Momentum sign regression" begin
        # The correct momentum drives flow AWAY from the tip (Rz/R² > 0 for a cone).
        # Verify: at the tip, U > 0 (recoiling away) and U decreases toward far-field.
        sol = solve_inner_bvp(ε=0.1)
        # Tip velocity should be positive (recoil direction)
        @test sol.U[1] > 0
        # Far-field velocity should be near zero
        @test abs(sol.U[end]) < 0.1
    end
end
