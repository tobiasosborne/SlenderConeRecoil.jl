using Test
using SlenderConeRecoil

@testset "Matching and composite" begin

    @testset "Linear interpolation" begin
        xs = [1.0, 2.0, 3.0, 4.0]
        ys = [10.0, 20.0, 30.0, 40.0]
        @test SlenderConeRecoil._interp_scalar(xs, ys, 1.5) ≈ 15.0
        @test SlenderConeRecoil._interp_scalar(xs, ys, 3.0) ≈ 30.0
        @test SlenderConeRecoil._interp_scalar(xs, ys, 0.5) ≈ 10.0  # clamp
        @test SlenderConeRecoil._interp_scalar(xs, ys, 5.0) ≈ 40.0  # clamp

        xq = [1.0, 2.5, 4.0]
        yq = SlenderConeRecoil._interp(xs, ys, xq)
        @test yq ≈ [10.0, 25.0, 40.0]
    end

    @testset "Inner far-field extraction" begin
        # Mock an inner solution with linear far-field S ≈ 0.1ξ
        ξ = collect(0.1:0.5:25.0)
        S = 0.1 .* ξ .+ 0.05  # slope=0.1, intercept=0.05
        Sξ = 0.1 .* ones(length(ξ))
        Sξξ = zeros(length(ξ))
        U = zeros(length(ξ))
        sol = InnerSolution(ξ, S, Sξ, Sξξ, U, 0.1, S[1], 0.0)

        slope, intercept = SlenderConeRecoil.inner_far_field(sol, 10.0)
        @test slope ≈ 0.1 atol=0.01
        @test intercept ≈ 0.05 atol=0.1
    end

    @testset "Composite construction runs" begin
        inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        outer = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=20.0)

        comp = composite_solution(inner, outer; n_points=100)
        @test comp isa CompositeSolution
        @test length(comp.ξ) == 100
        @test all(isfinite, comp.S)
        @test all(isfinite, comp.U)
    end

    @testset "Composite smoothness" begin
        inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        outer = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=20.0)
        comp = composite_solution(inner, outer; n_points=200)

        dS = diff(comp.S)
        dξ = diff(comp.ξ)
        slopes = dS ./ dξ
        @test all(isfinite, slopes)
        # Slopes should be bounded (O(ε) in far field, O(1) near tip)
        @test maximum(abs.(slopes)) < 100.0
    end

    @testset "Overlap residual is bounded" begin
        inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=15.0, ε=0.1)
        outer = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=15.0)

        res = overlap_residual(inner, outer)
        @test isfinite(res)
        @test res ≥ 0
        # Overlap mismatch should be moderate (not a perfect match due to Newton tolerance)
        @test res < 10.0
    end
end
