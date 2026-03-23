using Test
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))
include(joinpath(@__DIR__, "..", "src", "inner.jl"))
include(joinpath(@__DIR__, "..", "src", "outer.jl"))
include(joinpath(@__DIR__, "..", "src", "composite.jl"))
include(joinpath(@__DIR__, "..", "src", "pde.jl"))

@testset "Time-dependent PDE" begin

    @testset "Stretched grid" begin
        z = stretched_grid(100, 0.01, 10.0)
        @test length(z) == 100
        @test z[1] ≈ 0.01
        @test z[end] ≈ 10.0
        # Grid should be monotonically increasing
        @test all(diff(z) .> 0)
        # Grid is stretched (non-uniform spacing)
        @test z[2] - z[1] != z[end] - z[end-1]
    end

    @testset "Finite differences" begin
        # Test ddz! on f(z) = z² → f'(z) = 2z
        z = collect(range(0.0, 5.0, length=50))
        f = z .^ 2
        df = zeros(50)
        ddz!(df, f, z)
        # Check interior points
        for i in 5:45
            @test abs(df[i] - 2*z[i]) < 0.01
        end
    end

    @testset "PDE solver runs" begin
        pde = solve_pde(ε=0.1, N=100, z_min=0.05, z_max=5.0,
                        t_end=0.1, n_snapshots=3)
        @test pde isa PDESolution
        @test length(pde.z) == 100
        @test length(pde.t_snapshots) ≥ 2
        @test length(pde.R) == length(pde.t_snapshots)
    end

    @testset "Initial condition preserved at t≈0" begin
        pde = solve_pde(ε=0.1, N=100, z_min=0.05, z_max=5.0,
                        t_end=0.01, n_snapshots=2)
        # At t≈0, R should still be close to εz
        R_init = pde.R[1]
        R_expected = 0.1 .* pde.z
        @test maximum(abs.(R_init .- R_expected)) < 1e-10
    end

    @testset "R remains positive" begin
        pde = solve_pde(ε=0.1, N=100, z_min=0.05, z_max=5.0,
                        t_end=0.1, n_snapshots=3)
        for R_snap in pde.R
            @test all(R_snap .> 0)
        end
    end

    @testset "Rescale to similarity variables" begin
        pde = solve_pde(ε=0.1, N=100, z_min=0.05, z_max=5.0,
                        t_end=0.5, n_snapshots=3)
        # Rescale last snapshot
        t_idx = length(pde.t_snapshots)
        ξ, S = rescale_to_similarity(pde, t_idx)
        @test length(ξ) == length(pde.z)
        @test all(isfinite, ξ)
        @test all(isfinite, S)
        # For the initial condition (εz), rescaled S = εξ
        ξ0, S0 = rescale_to_similarity(pde, 1)
        # S0 = R0/t^{2/3} = εz/t^{2/3} = ε(z/t^{2/3}) = εξ
        # But at t=0 we can't rescale, so this returns raw values
        @test length(ξ0) == length(pde.z)
    end
end
