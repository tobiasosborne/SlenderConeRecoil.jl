using Test
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))
include(joinpath(@__DIR__, "..", "src", "inner.jl"))

@testset "Inner BVP solver" begin

    @testset "Tip initial conditions" begin
        # At ξ₀=0: U₀ = (4/5)*0 = 0
        S₀, U₀ = tip_initial_conditions(0.0, 1.0)
        @test S₀ == 1.0
        @test U₀ ≈ 0.0

        # At ξ₀=5: U₀ = 4
        S₀, U₀ = tip_initial_conditions(5.0, 0.5)
        @test S₀ == 0.5
        @test U₀ ≈ 4.0
    end

    @testset "ODE RHS at tip (S'=0 check)" begin
        # At tip with S'=0, U'=-2: the RHS should give S'≈0
        ξ₀ = 0.0; S₀ = 1.0
        _, U₀ = tip_initial_conditions(ξ₀, S₀)
        du = zeros(2)
        inner_rhs!(du, [S₀, U₀], nothing, ξ₀ + 1e-6)
        # S' should be near 0 at the tip
        @test abs(du[1]) < 0.1
    end

    @testset "Integration runs without error" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        @test sol isa InnerSolution
        @test length(sol.ξ) > 10
        @test all(isfinite, sol.S)
        @test all(isfinite, sol.U)
    end

    @testset "Solution positivity" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        # S should remain positive (physical: radius > 0)
        @test all(s -> s > 0, sol.S)
    end

    @testset "Far-field slope diagnostic" begin
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=50.0, ε=0.1)
        slopes = far_field_slope(sol)
        @test length(slopes) > 0
        @test all(isfinite, slopes)
        # The slope S/ξ should be roughly constant in the far field
        # (approaching some value related to ε, though the exact
        # value depends on the shooting parameters ξ₀, S₀)
        if length(slopes) ≥ 3
            spread = maximum(slopes) - minimum(slopes)
            @test spread < 1.0  # not wildly varying
        end
    end

    @testset "ODE residual check" begin
        # Verify the ODE is satisfied along the solution
        sol = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        n = length(sol.ξ)
        # Check at a midpoint using finite differences for S', U'
        mid = div(n, 2)
        if mid > 1 && mid < n
            ξm = sol.ξ[mid]
            Sm = sol.S[mid]
            Um = sol.U[mid]
            h = sol.ξ[mid+1] - sol.ξ[mid]
            Sξ_fd = (sol.S[mid+1] - sol.S[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])
            Uξ_fd = (sol.U[mid+1] - sol.U[mid-1]) / (sol.ξ[mid+1] - sol.ξ[mid-1])

            # Mass residual: 2S + 2Sξ(U-ξ) + S·Uξ
            mass_res = 2*Sm + 2*Sξ_fd*(Um - ξm) + Sm*Uξ_fd
            @test abs(mass_res) < 0.1  # should be close to 0

            # Momentum residual: -(2/9)U + (4/9)(U-ξ)Uξ + Sξ/S²
            mom_res = -(2/9)*Um + (4/9)*(Um-ξm)*Uξ_fd + Sξ_fd/Sm^2
            @test abs(mom_res) < 0.1
        end
    end
end
