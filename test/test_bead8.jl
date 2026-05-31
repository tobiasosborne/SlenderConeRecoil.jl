using Test
using SlenderConeRecoil

@testset "Time-dependent PDE" begin

    @testset "Stretched grid" begin
        z = SlenderConeRecoil.stretched_grid(100, 0.01, 10.0)
        @test length(z) == 100
        @test z[1] ≈ 0.01
        @test z[end] ≈ 10.0
        # Grid should be monotonically increasing
        @test all(diff(z) .> 0)
        # Grid should be finer near the left (tip)
        @test z[2] - z[1] < z[end] - z[end-1]
    end

    @testset "Finite differences (uniform grid)" begin
        # Test ddz! on f(z) = z² → f'(z) = 2z
        z = collect(range(0.0, 5.0, length=50))
        f = z .^ 2
        df = zeros(50)
        SlenderConeRecoil.ddz!(df, f, z)
        for i in 5:45
            @test abs(df[i] - 2*z[i]) < 0.01
        end
    end

    @testset "Finite differences (stretched grid, sin)" begin
        # Test ddz! on f(z) = sin(z) on a non-uniform stretched grid
        z = SlenderConeRecoil.stretched_grid(100, 0.1, 5.0)
        f = sin.(z)
        df = zeros(100)
        SlenderConeRecoil.ddz!(df, f, z)
        # Interior points: should match cos(z) to 2nd order
        for i in 5:95
            @test abs(df[i] - cos(z[i])) < 0.01
        end
    end

    @testset "PDE solver runs" begin
        # Use z_min=1.0, coarser grid to reduce stiffness
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.01, n_snapshots=2)
        @test pde isa PDESolution
        @test length(pde.z) == 50
        @test length(pde.t_snapshots) ≥ 2
        @test length(pde.R) == length(pde.t_snapshots)
    end

    @testset "Initial condition preserved at t≈0" begin
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.001, n_snapshots=2)
        R_init = pde.R[1]
        R_expected = 0.1 .* pde.z
        @test maximum(abs.(R_init .- R_expected)) < 1e-10
    end

    @testset "R remains positive" begin
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.01, n_snapshots=2)
        for R_snap in pde.R
            @test all(R_snap .> 0)
        end
    end

    @testset "Rescale to similarity variables" begin
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.01, n_snapshots=2)
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
