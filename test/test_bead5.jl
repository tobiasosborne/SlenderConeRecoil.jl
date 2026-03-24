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

    @testset "ODE RHS: no singularity" begin
        du = zeros(2)
        for ξ in [0.0, 0.5, 1.0, 5.0, 10.0]
            inner_rhs!(du, [1.0, 0.0], nothing, ξ + 1e-6)
            @test all(isfinite, du)
        end
    end

    @testset "Integration reaches far field" begin
        ξv, Sv, Uv = _shoot(2.5, 0.5, 50.0)
        @test ξv[end] ≈ 50.0 atol=1.0
        @test all(isfinite, Sv)
        @test all(isfinite, Uv)
    end

    @testset "2D Newton converges: both slope and U match" begin
        sol = solve_inner_bvp(ξ₀=2.5, S₀=0.5, ξ_max=100.0, ε=0.1)
        @test sol.ξ₀ > 0  # tip has recoiled
        @test sol.S₀ > 0  # positive tip radius

        # Far-field slope matches ε
        slopes = far_field_slope(sol)
        @test abs(slopes[end] - 0.1) < 0.001

        # Far-field U → 0
        @test abs(sol.U[end]) < 0.01

        # Blob: S₀ > ε·ξ₀ (tip is wider than undisturbed cone)
        @test sol.S₀ > 0.1 * sol.ξ₀
    end

    @testset "ODE residual along solution" begin
        sol = solve_inner_bvp(ξ₀=2.5, S₀=0.5, ξ_max=50.0, ε=0.1)
        n = length(sol.ξ)
        mid = div(n, 2)
        if mid > 1 && mid < n
            Sm = sol.S[mid]; Um = sol.U[mid]; ξm = sol.ξ[mid]
            Sξ_fd = (sol.S[mid+1] - sol.S[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])
            Uξ_fd = (sol.U[mid+1] - sol.U[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])

            mass_res = 2*Sm + 2*Sξ_fd*(Um - ξm) + Sm*Uξ_fd
            @test abs(mass_res) < 0.1

            mom_res = -(2/9)*Um + (4/9)*(Um-ξm)*Uξ_fd - Sξ_fd/Sm^2
            @test abs(mom_res) < 0.1
        end
    end
end
