using Test
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))
include(joinpath(@__DIR__, "..", "src", "outer.jl"))

@testset "Outer linearised problem" begin

    @testset "Zero perturbation is (nearly) steady" begin
        # With s₁=0, u₁=0 at large ξ, the solution should stay near zero
        # (the base state almost satisfies the equations far from the tip)
        sol = solve_outer(ε=0.1, ξ_min=5.0, ξ_max=50.0, s₁_init=0.0, u₁_init=0.0)
        @test sol isa OuterSolution
        @test length(sol.ξ) > 5
        @test all(isfinite, sol.s₁)
        @test all(isfinite, sol.u₁)
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
        # With a nonzero seed, there should be some perturbation
        max_s1 = maximum(abs.(sol.s₁))
        @test max_s1 > 0  # nontrivial
        @test all(isfinite, sol.s₁)
        @test all(isfinite, sol.u₁)
    end

    @testset "Perturbation bounded at moderate ξ" begin
        sol = solve_outer_driven(ε=0.1, ξ_min=5.0, ξ_max=30.0)
        # The perturbation should not blow up for moderate ξ range
        @test maximum(abs.(sol.s₁)) < 1e6
        @test maximum(abs.(sol.u₁)) < 1e6
    end

    @testset "ODE RHS evaluates cleanly" begin
        dy = zeros(2)
        outer_rhs!(dy, [0.01, 0.0], [0.1], 10.0)
        @test all(isfinite, dy)
    end
end
