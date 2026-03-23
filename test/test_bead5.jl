using Test
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))
include(joinpath(@__DIR__, "..", "src", "inner.jl"))

@testset "Inner BVP solver" begin

    @testset "Tip initial conditions" begin
        S₀, U₀ = tip_initial_conditions(0.0, 1.0)
        @test S₀ == 1.0
        @test U₀ ≈ 0.0

        S₀, U₀ = tip_initial_conditions(5.0, 0.5)
        @test S₀ == 0.5
        @test U₀ ≈ 4.0
    end

    @testset "ODE RHS: no singularity (correct denominator)" begin
        # With correct sign, denom = 1 + (8/9)S·v² > 0 always
        du = zeros(2)
        for ξ in [0.0, 0.5, 1.0, 5.0, 10.0]
            inner_rhs!(du, [1.0, 0.0], nothing, ξ + 1e-6)
            @test all(isfinite, du)
        end
    end

    @testset "Integration reaches far field" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=50.0, ε=0.1,
                              newton_iters=0)  # no Newton, just raw integration
        @test sol isa InnerSolution
        # With correct signs, should reach ξ_max without blowup
        @test sol.ξ[end] ≈ 50.0 atol=1.0
        @test all(isfinite, sol.S)
        @test all(isfinite, sol.U)
    end

    @testset "Solution positivity" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1,
                              newton_iters=0)
        @test all(s -> s > 0, sol.S)
    end

    @testset "Newton shooting converges" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=50.0, ε=0.1)
        slopes = far_field_slope(sol)
        @test length(slopes) > 0
        # After Newton iteration, far-field slope should match ε closely
        @test abs(slopes[end] - 0.1) < 0.001
    end

    @testset "ODE residual along solution" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1,
                              newton_iters=0)
        n = length(sol.ξ)
        mid = div(n, 2)
        if mid > 1 && mid < n
            Sm = sol.S[mid]; Um = sol.U[mid]; ξm = sol.ξ[mid]
            Sξ_fd = (sol.S[mid+1] - sol.S[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])
            Uξ_fd = (sol.U[mid+1] - sol.U[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])

            # Mass: 2S + 2Sξ(U-ξ) + S·Uξ
            mass_res = 2*Sm + 2*Sξ_fd*(Um - ξm) + Sm*Uξ_fd
            @test abs(mass_res) < 0.1

            # Momentum (correct sign): -(2/9)U + (4/9)(U-ξ)Uξ - Sξ/S²
            mom_res = -(2/9)*Um + (4/9)*(Um-ξm)*Uξ_fd - Sξ_fd/Sm^2
            @test abs(mom_res) < 0.1
        end
    end
end
