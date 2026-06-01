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

function _manufactured_collapse_data(; perturb::Bool=false,
                                     times=[0.0, 0.25, 0.6, 1.0],
                                     exponent=2 / 3,
                                     time_offset=0.0,
                                     z=collect(range(0.5, 6.0, length=121)))
    R = Vector{Vector{Float64}}()
    u = Vector{Vector{Float64}}()
    for t in times
        effective_time = t - time_offset
        if effective_time <= 0
            push!(R, 0.2 .+ 0.1 .* z)
            push!(u, zeros(length(z)))
            continue
        end
        ell = effective_time^exponent
        velocity_scale = exponent * effective_time^(exponent - 1)
        xi = z ./ ell
        S = @. 1.0 + 0.2 * xi
        U = @. 0.35 - 0.04 * xi
        if perturb && t == last(times)
            S = @. S + 0.025 * sin(1.7 * xi)
            U = @. U + 0.015 * cos(1.3 * xi)
        end
        push!(R, ell .* S)
        push!(u, velocity_scale .* U)
    end
    (z=z, t=Float64.(times), R=R, u=u)
end

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

    @testset "PDE conservation and domain diagnostics" begin
        pde = solve_pde(ε=0.1, N=50, z_min=1.0, z_max=10.0,
                        t_end=0.01, n_snapshots=2)
        diagnostics = pde.diagnostics

        @test diagnostics_succeeded(pde)
        @test diagnostics.pde_data_valid
        @test diagnostics.finite_state
        @test diagnostics.positive_radius
        @test diagnostics.minimum_radius > 0
        @test diagnostics.radius_positivity_margin == diagnostics.minimum_radius
        @test diagnostics.grid_finite
        @test diagnostics.grid_strictly_increasing
        @test diagnostics.grid_spacing_min > 0
        @test diagnostics.grid_spacing_max ≥ diagnostics.grid_spacing_min
        @test diagnostics.grid_spacing_ratio ≥ 1
        @test diagnostics.saved_time_points == length(pde.t_snapshots)
        @test diagnostics.state_snapshots == length(pde.R)
        @test diagnostics.saved_points_match_reported
        @test diagnostics.time_finite
        @test diagnostics.time_strictly_increasing
        @test diagnostics.time_step_min > 0
        @test diagnostics.retcode_string == "Success"
        @test diagnostics.retcode_successful
        @test diagnostics.endpoint_reached
        @test diagnostics.solver_steps === missing ||
              diagnostics.solver_steps ≥ 0

        @test length(diagnostics.area_mass) == length(pde.t_snapshots)
        @test diagnostics.initial_area_mass == first(diagnostics.area_mass)
        @test diagnostics.final_area_mass == last(diagnostics.area_mass)
        @test diagnostics.area_mass_drift[1] == 0.0
        @test diagnostics.relative_area_mass_drift[1] == 0.0
        @test isfinite(diagnostics.final_area_mass_drift)
        @test isfinite(diagnostics.final_relative_area_mass_drift)
        @test length(diagnostics.left_boundary_area_flux) == length(pde.t_snapshots)
        @test length(diagnostics.right_boundary_area_flux) == length(pde.t_snapshots)
        @test diagnostics.boundary_integrated_area_mass_change[1] == 0.0
        @test diagnostics.area_mass_balance_residual[1] == 0.0
        @test isfinite(diagnostics.final_area_mass_balance_residual)
        @test isfinite(diagnostics.max_abs_area_mass_balance_residual)
        @test diagnostics.pde_diagnostic_source_status == "IMPL-inferred"
        @test occursin("not Decent-King benchmark",
                       diagnostics.pde_diagnostic_basis)
        @test haskey(diagnostics, :similarity_collapse)
        @test diagnostics.similarity_collapse_status == :ok
        @test diagnostics.similarity_collapse_successful
        @test isfinite(diagnostics.similarity_collapse_score)
        @test diagnostics.similarity_collapse_grid_points == 128
        @test diagnostics.similarity_collapse_included_snapshots >= 2
    end

    @testset "PDE conservation diagnostics validation" begin
        z = [1.0, 1.5, 2.0]
        times = [0.0, 0.1]
        R = [[1.0, 1.0, 1.0], [1.1, 1.1, 1.1]]
        u = [[0.0, 0.1, 0.2], [0.0, 0.1, 0.2]]

        diagnostics = pde_conservation_diagnostics(
            z, times, R, u; retcode=:Success, endpoint=0.1,
            requested_endpoint=0.1, saved_points=2)
        @test diagnostics.initial_area_mass ≈ 1.0
        @test diagnostics.final_area_mass ≈ 1.21
        @test diagnostics.retcode_successful
        @test diagnostics.endpoint_reached
        @test diagnostics.saved_points_match_reported

        @test_throws ArgumentError pde_conservation_diagnostics(
            [1.0, 1.0, 2.0], times, R, u)
        @test_throws ArgumentError pde_conservation_diagnostics(
            z, [0.0, 0.0], R, u)

        R_bad = deepcopy(R)
        R_bad[2][2] = 0.0
        @test_throws DomainError pde_conservation_diagnostics(
            z, times, R_bad, u)

        R_nan = deepcopy(R)
        R_nan[1][2] = NaN
        @test_throws DomainError pde_conservation_diagnostics(
            z, times, R_nan, u)

        u_bad = deepcopy(u)
        u_bad[1][2] = Inf
        @test_throws DomainError pde_conservation_diagnostics(
            z, times, R, u_bad)
    end

    @testset "PDE diagnostic success hard checks" begin
        z = [1.0, 1.5, 2.0]
        times = [0.0, 0.1]
        R = [[1.0, 1.0, 1.0], [1.1, 1.1, 1.1]]
        u = [[0.0, 0.1, 0.2], [0.0, 0.1, 0.2]]
        ok = (context="manual PDE",
              retcode=:Success,
              successful=true,
              endpoint=0.1,
              requested_endpoint=0.1,
              saved_points=2)

        @test diagnostics_succeeded(PDESolution(z, times, R, u, 0.1, ok))

        failed_retcode = merge(ok, (retcode=:MaxIters, successful=false))
        @test !diagnostics_succeeded(
            PDESolution(z, times, R, u, 0.1, failed_retcode))

        missed_endpoint = merge(ok, (endpoint=0.05,))
        @test !diagnostics_succeeded(
            PDESolution(z, times, R, u, 0.1, missed_endpoint))

        R_bad = deepcopy(R)
        R_bad[1][2] = -0.1
        bad_radius = PDESolution(z, times, R_bad, u, 0.1, ok)
        @test !diagnostics_succeeded(bad_radius)
        @test !diagnostic_summary(bad_radius).successful

        u_bad = deepcopy(u)
        u_bad[2][3] = NaN
        @test !diagnostics_succeeded(PDESolution(z, times, R, u_bad, 0.1, ok))

        @test !diagnostics_succeeded(
            PDESolution([1.0, 1.0, 2.0], times, R, u, 0.1, ok))
    end

    @testset "Similarity collapse diagnostics" begin
        data = _manufactured_collapse_data()
        diagnostics = similarity_collapse_diagnostics(
            data.z, data.t, data.R, data.u; epsilon=0.1, n_grid=64)

        @test diagnostics.successful
        @test diagnostics.status == :ok
        @test diagnostics.source_status == "IMPL-inferred"
        @test diagnostics.exponent ≈ 2 / 3
        @test diagnostics.time_offset == 0.0
        @test diagnostics.velocity_scaling_convention ==
              "U = u / (exponent * (t - time_offset)^(exponent - 1))"
        @test diagnostics.xi_window.min < diagnostics.xi_window.max
        @test diagnostics.xi_window.source == :intersection_of_trusted_windows
        @test length(diagnostics.xi_grid) == 64
        @test diagnostics.interpolation_grid.points == 64
        @test diagnostics.component_weights.profile == 1.0
        @test sum(diagnostics.snapshot_weights) ≈ 1.0
        @test length(diagnostics.snapshot_weights) == 3
        @test length(diagnostics.included_snapshots) == 3
        @test length(diagnostics.excluded_snapshots) == 1
        @test diagnostics.excluded_snapshots[1].reason ==
              :nonpositive_effective_time
        @test length(diagnostics.excluded_regions) == 3
        @test diagnostics.reference_snapshot_index == 4
        @test diagnostics.reference_time == 1.0
        @test diagnostics.norms.profile.relative_rms < 1e-12
        @test diagnostics.norms.slope.relative_rms < 1e-12
        @test diagnostics.norms.curvature.rms < 1e-12
        @test diagnostics.norms.velocity.relative_rms < 1e-12
        @test diagnostics.aggregate_score < 1e-12
        @test diagnostics.norms.wave_phase.status ==
              :insufficient_wave_train

        pde = PDESolution(data.z, data.t, data.R, data.u, 0.1,
                          (context="manufactured",
                           retcode=:Success,
                           successful=true,
                           endpoint=1.0,
                           requested_endpoint=1.0,
                           saved_points=length(data.t)))
        pde_diagnostics = similarity_collapse_diagnostics(pde; n_grid=64)
        @test pde_diagnostics.aggregate_score < 1e-12
        @test pde_diagnostics.norms.wave_phase.status ==
              :insufficient_wave_train

        offset_data = _manufactured_collapse_data(
            times=[0.05, 0.35, 0.7, 1.1], time_offset=0.1)
        offset_diagnostics = similarity_collapse_diagnostics(
            offset_data.z, offset_data.t, offset_data.R, offset_data.u;
            epsilon=0.1, time_offset=0.1, n_grid=32)
        @test offset_diagnostics.successful
        @test offset_diagnostics.time_offset == 0.1
        @test offset_diagnostics.aggregate_score < 1e-12
        @test length(offset_diagnostics.excluded_snapshots) == 1
        @test offset_diagnostics.excluded_snapshots[1].reason ==
              :nonpositive_effective_time

        perturbed_data = _manufactured_collapse_data(perturb=true)
        perturbed = similarity_collapse_diagnostics(
            perturbed_data.z, perturbed_data.t, perturbed_data.R,
            perturbed_data.u; epsilon=0.1, n_grid=64)
        @test perturbed.successful
        @test perturbed.aggregate_score > diagnostics.aggregate_score
        @test perturbed.norms.profile.relative_rms > 1e-4
        @test perturbed.norms.velocity.relative_rms > 1e-4
        @test any(score -> score.score > 1e-4,
                  perturbed.per_snapshot_scores)
    end

    @testset "Similarity collapse invalid snapshot sets" begin
        data = _manufactured_collapse_data()
        insufficient = _manufactured_collapse_data(times=[-1.0, 0.0])
        insufficient_diagnostics = similarity_collapse_diagnostics(
            insufficient.z, insufficient.t, insufficient.R, insufficient.u;
            epsilon=0.1, n_grid=32)
        @test !insufficient_diagnostics.successful
        @test insufficient_diagnostics.status == :insufficient_snapshots
        @test length(insufficient_diagnostics.excluded_snapshots) == 2

        nonoverlap = _manufactured_collapse_data(times=[0.001, 1e6],
                                                 z=[1.0, 2.0, 3.0, 4.0, 5.0])
        nonoverlap_diagnostics = similarity_collapse_diagnostics(
            nonoverlap.z, nonoverlap.t, nonoverlap.R, nonoverlap.u;
            epsilon=0.1, n_grid=32)
        @test !nonoverlap_diagnostics.successful
        @test nonoverlap_diagnostics.status == :empty_xi_overlap
        @test nonoverlap_diagnostics.xi_window.min >=
              nonoverlap_diagnostics.xi_window.max

        @test_throws ArgumentError similarity_collapse_diagnostics(
            [1.0, 1.0, 2.0], data.t, data.R, data.u)
        @test_throws ArgumentError similarity_collapse_diagnostics(
            data.z, data.t, data.R, data.u; n_grid=4)
        @test_throws ArgumentError similarity_collapse_diagnostics(
            data.z, data.t, data.R, data.u; exponent=0.0)
        @test_throws ArgumentError similarity_collapse_diagnostics(
            data.z, data.t, data.R, data.u; xi_window=(2.0, 1.0))
    end

    @testset "Mapped coordinate transform consistency" begin
        xi = [1.0, 2.0, 4.0]
        S = [0.2, 0.3, 0.5]
        U = [1.0, 0.5, -0.25]
        transform = mapped_coordinate_transform(xi, S, U, 8.0;
                                                exponent=2 / 3,
                                                time_offset=0.0)

        @test transform.effective_time == 8.0
        @test transform.length_scale ≈ 4.0
        @test transform.velocity_scale ≈ 1 / 3
        @test transform.z ≈ 4.0 .* xi
        @test transform.R ≈ 4.0 .* S
        @test transform.u ≈ (1 / 3) .* U
        @test !haskey(transform, :coordinate_transform)

        @test_throws ArgumentError mapped_coordinate_transform(
            [1.0, 1.0, 2.0], S, U, 8.0)
        @test_throws DomainError mapped_coordinate_transform(
            xi, [0.2, -0.1, 0.5], U, 8.0)
        @test_throws ArgumentError mapped_coordinate_transform(
            xi, S, U, 0.0; time_offset=0.0)
        @test_throws ArgumentError mapped_coordinate_transform(
            xi, S, U, 8.0; exponent=0.0)
    end

    @testset "Mapped PDE verifier diagnostics" begin
        mapped = solve_mapped_pde(epsilon=0.1, N=24, xi_min=1.0,
                                  xi_max=8.0,
                                  effective_time_start=0.25,
                                  effective_time_end=1.0,
                                  n_snapshots=4)
        alias_mapped = solve_mapped_pde(ε=0.1, N=8, xi_min=1.0,
                                        xi_max=3.0,
                                        effective_time_start=0.25,
                                        effective_time_end=0.5,
                                        n_snapshots=2)
        diagnostics = mapped.diagnostics

        @test mapped isa MappedPDESolution
        @test diagnostics_succeeded(alias_mapped)
        @test diagnostics_succeeded(mapped)
        @test length(mapped.xi) == 24
        @test length(mapped.t_snapshots) == 4
        @test length(mapped.z) == 4
        @test length(mapped.R) == 4
        @test length(mapped.u) == 4
        @test diagnostics.source_status == "IMPL-inferred"
        @test occursin("not a full transformed PDE solve",
                       diagnostics.diagnostic_basis)
        @test diagnostics.coordinate_system == :similarity_mapped
        @test diagnostics.coordinate_transform ==
              "z = ell(t) * xi; R = ell(t) * S; u = ell_dot(t) * U"
        @test diagnostics.xi_domain.min ≈ 1.0
        @test diagnostics.xi_domain.max ≈ 8.0
        @test diagnostics.xi_domain.points == 24
        @test length(diagnostics.length_scales) == 4
        @test all(diff(diagnostics.length_scales) .> 0)
        @test length(diagnostics.physical_z_ranges) == 4
        @test diagnostics.common_physical_grid.points == 24
        @test diagnostics.common_physical_grid.z_min > 0

        reference = diagnostics.fixed_grid_reference
        @test diagnostics.reaches_longer_effective_similarity_time
        @test reference.reaches_longer_effective_similarity_time
        @test reference.fixed_grid_effective_time_limit ≈ 0.25
        @test reference.final_effective_time ≈ 1.0
        @test reference.effective_time_ratio ≈ 4.0

        @test diagnostics.mapped_stationarity.status == :ok
        @test diagnostics.mapped_stationarity.successful
        @test diagnostics.mapped_stationarity_residual_norm < 1e-12
        @test diagnostics.residual_norm < 1e-12
        @test diagnostics.final_residual_norm < 1e-12

        @test haskey(diagnostics, :conservation_diagnostics)
        @test diagnostics.conservation_diagnostics.pde_data_valid
        @test diagnostics.pde_data_valid
        @test diagnostics.positive_radius
        @test diagnostics.grid_strictly_increasing
        @test diagnostics.time_strictly_increasing
        @test diagnostics.retcode_successful
        @test diagnostics.endpoint_reached
        @test diagnostics.minimum_radius > 0
        @test diagnostics.radius_positivity_margin > 0
        @test isfinite(diagnostics.max_abs_area_mass_balance_residual)
        @test diagnostics.similarity_collapse_successful
        @test diagnostics.similarity_collapse_status == :ok
        @test isfinite(diagnostics.similarity_collapse_score)
        @test diagnostics.similarity_collapse_score < 1e-2

        summary = diagnostic_summary(mapped)
        @test summary.successful
        @test summary.problem_kind == :mapped_pde_verification
        summary_fields = as_namedtuple(summary)
        @test summary_fields.mesh_variable == :xi
        @test summary_fields.reaches_longer_effective_similarity_time
        @test summary_fields.minimum_radius > 0

        @test_throws ArgumentError solve_mapped_pde(N=4)
        @test_throws ArgumentError solve_mapped_pde(xi_min=0.0)
        @test_throws ArgumentError solve_mapped_pde(
            effective_time_start=1.0, effective_time_end=0.5)
        @test_throws ArgumentError solve_mapped_pde(
            fixed_reference_z_max=0.0)
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
