using Test
using SlenderConeRecoil

@testset "Matching and composite (local diagnostic construction)" begin

    @testset "Clamping interpolation helper" begin
        xs = [1.0, 2.0, 3.0, 4.0]
        ys = [10.0, 20.0, 30.0, 40.0]
        @test SlenderConeRecoil._interp_scalar(xs, ys, 1.5) ≈ 15.0
        @test SlenderConeRecoil._interp_scalar(xs, ys, 3.0) ≈ 30.0
        @test SlenderConeRecoil._interp_scalar(xs, ys, 0.5) ≈ 10.0
        @test SlenderConeRecoil._interp_scalar(xs, ys, 5.0) ≈ 40.0

        xq = [1.0, 2.5, 4.0]
        yq = SlenderConeRecoil._interp(xs, ys, xq)
        @test yq ≈ [10.0, 25.0, 40.0]

        @test_throws ArgumentError SlenderConeRecoil._interp_scalar([1.0, 1.0, 2.0], ys[1:3], 1.5)
        @test_throws ArgumentError SlenderConeRecoil._interp_scalar([1.0, 3.0, 2.0], ys[1:3], 1.5)
        @test_throws ArgumentError SlenderConeRecoil._interp_scalar(xs, ys[1:3], 1.5)
    end

    @testset "Inner far-field extraction" begin
        # Mock a local/exploratory common part with linear far-field S ≈ 0.1ξ.
        ξ = collect(0.1:0.5:25.0)
        S = 0.1 .* ξ .+ 0.05  # slope=0.1, intercept=0.05
        Sξ = 0.1 .* ones(length(ξ))
        Sξξ = zeros(length(ξ))
        U = zeros(length(ξ))
        sol = InnerSolution(ξ, S, Sξ, Sξξ, U, 0.1, S[1], 0.0)

        slope, intercept = SlenderConeRecoil.inner_far_field(sol, 10.0)
        @test slope ≈ 0.1 atol=0.01
        @test intercept ≈ 0.05 atol=0.1
        @test_throws ArgumentError SlenderConeRecoil.inner_far_field(sol, 24.5)
        @test_throws ArgumentError SlenderConeRecoil.inner_far_field(sol, 30.0)
    end

    @testset "Validated matching and overlap domains for local diagnostics" begin
        ξ = collect(0.0:1.0:10.0)
        S = 0.1 .* ξ .+ 0.05
        Sξ = 0.1 .* ones(length(ξ))
        Sξξ = zeros(length(ξ))
        U = zeros(length(ξ))
        inner = InnerSolution(ξ, S, Sξ, Sξξ, U, 0.0, S[1], 0.0)

        outer_ξ = collect(2.0:2.0:10.0)
        outer = OuterSolution(outer_ξ, zeros(length(outer_ξ)), zeros(length(outer_ξ)),
                              zeros(length(outer_ξ)), zeros(length(outer_ξ)), 0.1)

        @test_throws ArgumentError solve_outer_matched(inner; ξ_match=-1.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_matched(inner; ξ_match=11.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_matched(inner; ξ_match=5.0, ξ_max=5.0)
        @test_throws ArgumentError solve_outer_matched(inner; ξ_match=5.0, ξ_max=12.0, maxiters=0)
        @test_throws ArgumentError solve_outer_full(inner; ξ_match=-1.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_full(inner; ξ_match=5.0, ξ_max=12.0, maxiters=0)
        @test_throws ArgumentError solve_outer_linearised(inner; ξ_match=11.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_linearised(inner; ξ_match=5.0, ξ_max=12.0, maxiters=0)

        disjoint_outer = OuterSolution([20.0, 21.0], [0.0, 0.0], [0.0, 0.0],
                                       [0.0, 0.0], [0.0, 0.0], 0.1)
        touching_outer = OuterSolution([10.0, 11.0], [0.0, 0.0], [0.0, 0.0],
                                       [0.0, 0.0], [0.0, 0.0], 0.1)
        @test_throws ArgumentError composite_solution(inner, disjoint_outer; n_points=10)
        @test_throws ArgumentError composite_solution(inner, touching_outer; n_points=10)
        @test_throws ArgumentError overlap_residual(inner, disjoint_outer)
        @test_throws ArgumentError composite_solution(inner, outer; ξ_grid=[2.0, 1.0])
        @test_throws ArgumentError composite_solution(inner, outer; ξ_grid=[1.0, 3.0])
        @test_throws ArgumentError composite_solution(inner, outer; ξ_match=1.0)
        @test_throws ArgumentError composite_solution(inner, outer; ξ_match=10.0)
        @test_throws ArgumentError overlap_residual(inner, outer; ξ_start=1.0, ξ_end=3.0)

        invalid_epsilon_outer = OuterSolution(outer_ξ, zeros(length(outer_ξ)), zeros(length(outer_ξ)),
                                              zeros(length(outer_ξ)), zeros(length(outer_ξ)), 0.0)
        @test_throws ArgumentError composite_solution(inner, invalid_epsilon_outer; n_points=10)
        @test_throws ArgumentError overlap_residual(inner, invalid_epsilon_outer)

        zero_far_field = InnerSolution(ξ, zeros(length(ξ)), zeros(length(ξ)),
                                       zeros(length(ξ)), zeros(length(ξ)), 0.0, 0.0, 0.0)
        @test_throws ArgumentError solve_outer_matched(zero_far_field; ξ_match=5.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_full(zero_far_field; ξ_match=5.0, ξ_max=12.0)
        @test_throws ArgumentError solve_outer_linearised(zero_far_field; ξ_match=5.0, ξ_max=12.0)
    end

    @testset "Composite subtracts fitted overlap as a local common-part diagnostic" begin
        ξ = collect(0.0:1.0:10.0)
        S = 0.1 .* ξ .+ 0.05
        Sξ = 0.1 .* ones(length(ξ))
        Sξξ = zeros(length(ξ))
        U = zeros(length(ξ))
        inner = InnerSolution(ξ, S, Sξ, Sξξ, U, 0.0, S[1], 0.0)

        outer_ξ = collect(2.0:2.0:10.0)
        outer = OuterSolution(outer_ξ, zeros(length(outer_ξ)), zeros(length(outer_ξ)),
                              zeros(length(outer_ξ)), zeros(length(outer_ξ)), 0.1)

        comp = composite_solution(inner, outer; n_points=25)
        @test comp.diagnostics.overlap_slope ≈ 0.1 atol=1e-12
        @test comp.diagnostics.overlap_intercept ≈ 0.05 atol=1e-12
        @test comp.diagnostics.ξ_match ≈ 2.0
        @test comp.diagnostics.fit_points == count(x -> x > 2.0, ξ)
        @test comp.S ≈ 0.1 .* comp.ξ atol=1e-12
    end

    @testset "Composite construction runs for reconstructed inner/outer solves" begin
        inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=20.0, ε=0.1)
        outer = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=20.0)

        comp = composite_solution(inner, outer; n_points=100)
        @test comp isa CompositeSolution
        @test length(comp.ξ) == 100
        @test all(isfinite, comp.S)
        @test all(isfinite, comp.U)
    end

    @testset "Composite smoothness for reconstructed inner/outer solves" begin
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

    @testset "Overlap residual is bounded for reconstructed inner/outer solves" begin
        inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=15.0, ε=0.1)
        outer = solve_outer_driven(ε=0.1, ξ_min=2.0, ξ_max=15.0)

        res = overlap_residual(inner, outer)
        @test isfinite(res)
        @test res ≥ 0
        # Overlap mismatch should be moderate (not a perfect match due to Newton tolerance)
        @test res < 10.0
    end
end
