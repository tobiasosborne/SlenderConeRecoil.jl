using Test
using SlenderConeRecoil

_pde_rhs_params(z) = (length(z), z,
                      zeros(length(z)), zeros(length(z)), zeros(length(z)),
                      zeros(length(z)), zeros(length(z)), zeros(length(z)))

# Manufactured profile for operator verification only. This is not cone
# benchmark data and should not be used as a physical reference solution.
function _mms_pde_profile(z)
    R = @. 1.4 + 0.08*sin(1.7*z) + 0.05*cos(2.3*z)
    u = @. 0.2*cos(1.1*z) + 0.07*sin(2.9*z)
    Rz = @. 0.08*1.7*cos(1.7*z) - 0.05*2.3*sin(2.3*z)
    uz = @. -0.2*1.1*sin(1.1*z) + 0.07*2.9*cos(2.9*z)
    invR_z = @. -Rz / R^2
    Rzzz = @. -0.08*1.7^3*cos(1.7*z) + 0.05*2.3^3*sin(2.3*z)
    Rt = @. -u*Rz - 0.5*R*uz
    ut = @. -u*uz - invR_z + Rzzz
    (R=R, u=u, Rt=Rt, ut=ut)
end

function _mms_pde_rhs_error(z; margin=8)
    N = length(z)
    N > 2margin + 1 || throw(ArgumentError("MMS grid too small for margin=$margin"))

    mms = _mms_pde_profile(z)
    w = vcat(mms.R, mms.u)
    dw = fill(NaN, 2N)
    SlenderConeRecoil.pde_rhs!(dw, w, _pde_rhs_params(z), 0.0)

    interior = (margin+1):(N-margin)
    velocity_interior = (N + margin + 1):(2N - margin)
    radius_error = maximum(abs.(dw[interior] .- mms.Rt[interior]))
    velocity_error = maximum(abs.(dw[velocity_interior] .- mms.ut[interior]))
    (radius_error=radius_error,
     velocity_error=velocity_error,
     max_error=max(radius_error, velocity_error),
     h=maximum(diff(z)),
     left_radius_bc=dw[1],
     left_velocity_bc=dw[N+1],
     right_radius_bc_delta=dw[N] - dw[N-1],
     right_velocity_bc_delta=dw[2N] - dw[2N-1])
end

_observed_order(coarse, fine, component::Symbol) =
    log(getproperty(coarse, component) / getproperty(fine, component)) /
    log(coarse.h / fine.h)

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

    @testset "Stretched grid validation" begin
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(2, 0.01, 10.0)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, NaN, 10.0)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, 0.01, Inf)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, 1.0, 1.0)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, 2.0, 1.0)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, 0.01, 10.0; β=0.0)
        @test_throws ArgumentError SlenderConeRecoil.stretched_grid(3, 0.01, 10.0; β=Inf)
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

    @testset "Finite difference validation" begin
        z = [0.0, 1.0, 2.0]
        f = [0.0, 1.0, 4.0]
        df = zeros(3)

        @test_throws ArgumentError SlenderConeRecoil.ddz!(zeros(2), f, z)
        @test_throws ArgumentError SlenderConeRecoil.ddz!(df, zeros(2), z)
        @test_throws ArgumentError SlenderConeRecoil.ddz!(zeros(2), zeros(2), [0.0, 1.0])
        @test_throws ArgumentError SlenderConeRecoil.ddz!(df, f, [0.0, NaN, 2.0])
        @test_throws ArgumentError SlenderConeRecoil.ddz!(df, f, [0.0, 1.0, 1.0])
        @test_throws ArgumentError SlenderConeRecoil.ddz!(df, f, [0.0, 2.0, 1.0])
    end

    @testset "PDE RHS boundary conditions" begin
        z = collect(range(1.0, 2.0, length=6))
        N = length(z)
        R = 0.2 .+ 0.1 .* z
        u = sin.(z)
        w = vcat(R, u)
        dw = fill(NaN, 2N)

        SlenderConeRecoil.pde_rhs!(dw, w, _pde_rhs_params(z), 0.0)

        @test all(isfinite, dw)
        @test dw[1] == 0.0
        @test dw[N+1] == 0.0
        @test dw[N] == dw[N-1]
        @test dw[2N] == dw[2N-1]
    end

    @testset "PDE RHS manufactured-solution convergence" begin
        cases = (
            (:uniform, N -> collect(range(0.2, 2.6, length=N))),
            (:stretched, N -> SlenderConeRecoil.stretched_grid(N, 0.2, 2.6; β=1.7)),
        )

        for (name, make_grid) in cases
            @testset "$name grid" begin
                errors = [_mms_pde_rhs_error(make_grid(N)) for N in (65, 129, 257)]

                for err in errors
                    @test isfinite(err.radius_error)
                    @test isfinite(err.velocity_error)
                    @test err.left_radius_bc == 0.0
                    @test err.left_velocity_bc == 0.0
                    @test err.right_radius_bc_delta == 0.0
                    @test err.right_velocity_bc_delta == 0.0
                end

                radius_orders = [_observed_order(errors[i], errors[i+1], :radius_error)
                                 for i in 1:(length(errors)-1)]
                velocity_orders = [_observed_order(errors[i], errors[i+1], :velocity_error)
                                   for i in 1:(length(errors)-1)]

                @test minimum(radius_orders) > 1.75
                @test minimum(velocity_orders) > 1.85
                @test radius_orders[end] > 1.9
                @test velocity_orders[end] > 1.9
                @test errors[end].radius_error < 8e-5
                @test errors[end].velocity_error < 5e-4
            end
        end
    end

    @testset "PDE RHS validation" begin
        z = collect(range(1.0, 2.0, length=5))
        N = length(z)
        R = 0.2 .+ 0.1 .* z
        u = zeros(N)
        w = vcat(R, u)
        dw = zeros(2N)
        p = _pde_rhs_params(z)

        @test_throws ArgumentError SlenderConeRecoil.pde_rhs!(zeros(2N - 1), w, p, 0.0)
        @test_throws ArgumentError SlenderConeRecoil.pde_rhs!(dw, w[1:end-1], p, 0.0)

        w_bad_radius = copy(w)
        w_bad_radius[3] = 0.0
        @test_throws DomainError SlenderConeRecoil.pde_rhs!(dw, w_bad_radius, p, 0.0)

        w_bad_velocity = copy(w)
        w_bad_velocity[N+2] = NaN
        @test_throws DomainError SlenderConeRecoil.pde_rhs!(dw, w_bad_velocity, p, 0.0)

        p_bad_grid = _pde_rhs_params([1.0, 1.5, 1.5, 1.75, 2.0])
        @test_throws ArgumentError SlenderConeRecoil.pde_rhs!(dw, w, p_bad_grid, 0.0)
    end

    @testset "PDE solver input validation" begin
        @test_throws ArgumentError solve_pde(ε=0.0, N=10, z_min=1.0, z_max=2.0,
                                             t_end=0.001, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=NaN, N=10, z_min=1.0, z_max=2.0,
                                             t_end=0.001, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=0.1, N=2, z_min=1.0, z_max=2.0,
                                             t_end=0.001, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=0.1, N=10, z_min=2.0, z_max=1.0,
                                             t_end=0.001, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=0.1, N=10, z_min=1.0, z_max=2.0,
                                             t_end=-0.001, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=0.1, N=10, z_min=1.0, z_max=2.0,
                                             t_end=Inf, n_snapshots=1)
        @test_throws ArgumentError solve_pde(ε=0.1, N=10, z_min=1.0, z_max=2.0,
                                             t_end=0.001, n_snapshots=0)
        @test_throws ArgumentError solve_pde(ε=0.1, N=10, z_min=0.0, z_max=1.0,
                                             t_end=0.001, n_snapshots=1)
    end

    @testset "PDE solver runs" begin
        # Use z_min=1.0, coarser grid to reduce stiffness
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.01, n_snapshots=2)
        @test pde isa PDESolution
        @test length(pde.z) == 50
        @test length(pde.t_snapshots) ≥ 2
        @test length(pde.R) == length(pde.t_snapshots)
        @test pde.diagnostics.requested_endpoint ≈ 0.01
        @test pde.diagnostics.endpoint ≈ 0.01
        @test pde.diagnostics.saved_points == length(pde.t_snapshots)
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
