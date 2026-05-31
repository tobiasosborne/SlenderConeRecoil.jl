using Test
using SlenderConeRecoil
using LinearAlgebra: norm

const LOCAL_REGRESSION_EPSILON = 0.1

function _central_derivative(xs, ys, i)
    hm = xs[i] - xs[i - 1]
    hp = xs[i + 1] - xs[i]
    (ys[i + 1] * hm^2 + ys[i] * (hp^2 - hm^2) - ys[i - 1] * hp^2) /
        (hp * hm * (hp + hm))
end

function _interior_sample_indices(n; count=7)
    lo = ceil(Int, 0.15 * n)
    hi = floor(Int, 0.85 * n)
    idxs = unique(round.(Int, range(lo, hi, length=count)))
    [i for i in idxs if 2 <= i <= n - 1]
end

function _inner_residuals_at(sol, i)
    xi = sol.ξ[i]
    S = sol.S[i]
    U = sol.U[i]
    Sp = _central_derivative(sol.ξ, sol.S, i)
    Up = _central_derivative(sol.ξ, sol.U, i)
    Sppp = _central_derivative(sol.ξ, sol.Sξξ, i)

    mass = 2S + 2Sp * (U - xi) + S * Up
    momentum = -(2 / 9) * U + (4 / 9) * (U - xi) * Up - Sp / S^2 - Sppp
    (mass=mass, momentum=momentum)
end

function _deadband_sign_changes(values; deadband)
    changes = 0
    last_sign = 0
    for value in values
        current_sign = value > deadband ? 1 : value < -deadband ? -1 : 0
        if current_sign != 0
            if last_sign != 0 && current_sign != last_sign
                changes += 1
            end
            last_sign = current_sign
        end
    end
    changes
end

function _large_turning_points(xs, values; deadband, xmax)
    idxs = Int[]
    for i in 2:length(values)-1
        if xs[i] < xmax &&
           abs(values[i]) > deadband &&
           (values[i] - values[i - 1]) * (values[i + 1] - values[i]) < 0
            push!(idxs, i)
        end
    end
    idxs
end

@testset "Source-status-aware numerical regressions" begin
    # These references are locally blessed against the current implementation
    # and Manifest. They are not published Decent-King benchmark values.
    inner = solve_inner_bvp(ε=LOCAL_REGRESSION_EPSILON)

    @testset "Inner locally blessed reference values" begin
        @test inner.ξ₀ ≈ 2.764907207109112 rtol=2e-3
        @test inner.S₀ ≈ 0.2383219144025037 rtol=2e-3
        @test inner.Sξξ₀ ≈ 1.158105692501843 rtol=2e-3

        @test !inner.converged
        @test inner.iterations == 21
        @test inner.termination_reason == "line search failed to reduce residual"
        @test inner.final_residual ≈ [-0.0010014845492979374,
                                      -0.006346091067721614,
                                       0.0005583784462499321] rtol=5e-2 atol=2e-5
        @test inner.final_residual_norm ≈ norm(inner.final_residual) rtol=1e-12
        @test inner.final_residual_norm ≈ 0.006448847155231315 rtol=5e-2

        far_slopes = SlenderConeRecoil.far_field_slope(inner)
        @test far_slopes[end] ≈ 0.09899851545070207 rtol=5e-4
        @test inner.S[end] / inner.ξ[end] ≈ 0.09899851545070207 rtol=5e-4
        @test inner.U[1] ≈ 2.21192576568729 rtol=2e-3
        @test inner.U[end] ≈ -0.006346091067721614 rtol=5e-2 atol=2e-5
    end

    @testset "Inner ODE residuals at multiple interior points" begin
        idxs = _interior_sample_indices(length(inner.ξ))
        @test length(idxs) == 7

        residuals = [_inner_residuals_at(inner, i) for i in idxs]
        mass_residuals = [r.mass for r in residuals]
        momentum_residuals = [r.momentum for r in residuals]

        @test maximum(abs.(mass_residuals)) < 1e-3
        @test maximum(abs.(momentum_residuals)) < 1e-3
        @test any(<(0), mass_residuals)
        @test any(>(0), mass_residuals)
        @test any(<(0), momentum_residuals)
        @test any(>(0), momentum_residuals)
    end

    @testset "Axial-curvature capillary-wave shape" begin
        excess = inner.S .- LOCAL_REGRESSION_EPSILON .* inner.ξ

        @test maximum(excess) ≈ 0.047525952370282065 rtol=1e-1
        @test minimum(excess) ≈ -0.0600890729578758 rtol=1e-1
        @test maximum(excess) > 0.04
        @test minimum(excess) < -0.04
        @test _deadband_sign_changes(excess; deadband=0.005) >= 6

        early_turns = _large_turning_points(inner.ξ, excess; deadband=0.002, xmax=8.0)
        @test length(early_turns) >= 8
        @test excess[early_turns[1]] < 0
        @test excess[early_turns[2]] > 0
        @test inner.ξ[early_turns[2]] ≈ 3.5783760968520104 rtol=5e-2
        @test excess[early_turns[2]] ≈ maximum(excess) rtol=1e-2
    end

    @testset "Outer/composite locally blessed overlap references" begin
        outer = solve_outer_matched(inner; ξ_match=15.0, ξ_max=30.0)

        @test string(outer.diagnostics.retcode) == "Success"
        @test outer.diagnostics.endpoint ≈ 30.0
        @test outer.ε ≈ 0.09899851545070207 rtol=5e-4
        @test outer.ξ[1] ≈ 15.0
        @test outer.s₁[1] ≈ 0.0030878268255205565 rtol=5e-2
        @test outer.u₁[1] ≈ 0.09019954334382589 rtol=5e-2
        @test outer.s₁[end] ≈ 0.00053422156606065 rtol=1e-1 atol=2e-5
        @test outer.u₁[end] ≈ 0.0020764169710089093 rtol=1e-1 atol=2e-5

        overlap = overlap_residual(inner, outer; ξ_start=15.0, ξ_end=30.0)
        @test overlap ≈ 9.727719665192918e-5 rtol=2.5e-1 atol=2e-5
        @test overlap < 2e-4

        comp = composite_solution(inner, outer; n_points=64, ξ_match=15.0)
        @test comp.diagnostics.overlap_slope ≈ 0.09898068182848431 rtol=5e-4
        @test comp.diagnostics.overlap_intercept ≈ 0.0010305094504263757 rtol=2e-1 atol=2e-4
        @test comp.diagnostics.fit_points > 1_000
        @test comp.S[1] ≈ 1.4903903802944123 rtol=5e-4
        @test comp.S[end] ≈ 2.9703536737919958 rtol=5e-4
        @test comp.U[1] ≈ 0.18039908668765178 rtol=5e-2
        @test comp.U[end] ≈ 0.0017048922522462091 rtol=1e-1 atol=2e-5
    end

    @testset "Short-time PDE internal consistency locally blessed" begin
        pde = solve_pde(ε=LOCAL_REGRESSION_EPSILON, N=20, z_min=1.0, z_max=3.0,
                        t_end=1e-4, n_snapshots=2)

        @test string(pde.diagnostics.retcode) == "Success"
        @test pde.diagnostics.endpoint ≈ 1e-4
        @test length(pde.t_snapshots) == 3

        initial_R = pde.R[1]
        final_R = pde.R[end]
        final_u = pde.u[end]

        @test maximum(abs.(initial_R .- LOCAL_REGRESSION_EPSILON .* pde.z)) < 1e-12
        @test final_R[1] ≈ initial_R[1] atol=1e-8
        @test final_u[1] == 0.0
        @test all(final_R .> 0)

        max_radius_drift = maximum(abs.(final_R .- initial_R))
        max_velocity = maximum(abs.(final_u))
        @test max_radius_drift ≈ 8.213770129328335e-8 rtol=3.5e-1 atol=5e-8
        @test max_velocity ≈ 0.0009664176229987088 rtol=1e-1

        xi, S = rescale_to_similarity(pde, length(pde.t_snapshots))
        similarity_drift = maximum(abs.(S .- LOCAL_REGRESSION_EPSILON .* xi))
        @test similarity_drift ≈ 3.8124943714024084e-5 rtol=3.5e-1 atol=2e-5
        @test similarity_drift < 8e-5
    end
end
