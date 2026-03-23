#!/usr/bin/env julia
# Generate all figures for SlenderConeRecoil.
# Usage: julia --project scripts/figures.jl

using Plots
default(dpi=150, size=(600, 400), linewidth=1.5, legend=:topright)

# Load modules
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))
include(joinpath(@__DIR__, "..", "src", "slender.jl"))
include(joinpath(@__DIR__, "..", "src", "similarity.jl"))
include(joinpath(@__DIR__, "..", "src", "inner.jl"))
include(joinpath(@__DIR__, "..", "src", "outer.jl"))
include(joinpath(@__DIR__, "..", "src", "composite.jl"))
include(joinpath(@__DIR__, "..", "src", "pde.jl"))

FIGDIR = joinpath(@__DIR__, "..", "figures")
mkpath(FIGDIR)

# ── Figure 1: Similarity profiles S(ξ) for several ε ──────────────────
function figure1()
    println("Figure 1: Similarity profiles...")
    p = plot(xlabel="ξ", ylabel="S(ξ)", title="Inner similarity profiles")

    for (ε, S₀, color) in [(0.05, 1.0, :blue), (0.1, 1.0, :red), (0.2, 1.0, :green)]
        sol = solve_inner_bvp(ξ₀=0.0, S₀=S₀, ξ_max=15.0, ε=ε)
        plot!(p, sol.ξ, sol.S, label="ε = $ε", color=color)
        # Also plot the far-field asymptote S = εξ
        ξ_ff = range(0, 15, length=50)
        plot!(p, ξ_ff, ε .* ξ_ff, label="", color=color, linestyle=:dash, alpha=0.4)
    end

    savefig(p, joinpath(FIGDIR, "fig1_similarity_profiles.pdf"))
    println("  Saved fig1_similarity_profiles.pdf")
end

# ── Figure 2: Inner, outer, and composite overlaid ─────────────────────
function figure2()
    println("Figure 2: Inner/outer/composite...")
    ε = 0.1

    inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=15.0, ε=ε)
    outer = solve_outer_driven(ε=ε, ξ_min=2.0, ξ_max=15.0)
    comp = composite_solution(inner, outer; n_points=300)

    p = plot(xlabel="ξ", ylabel="S", title="Inner, outer, composite (ε=$ε)")
    plot!(p, inner.ξ, inner.S, label="Inner", color=:blue, linewidth=2)
    plot!(p, outer.ξ, ε .* outer.ξ .+ outer.s₁, label="Outer (εξ + s₁)",
          color=:red, linestyle=:dash)
    plot!(p, comp.ξ, comp.S, label="Composite", color=:black, linestyle=:dot, linewidth=2)

    savefig(p, joinpath(FIGDIR, "fig2_inner_outer_composite.pdf"))
    println("  Saved fig2_inner_outer_composite.pdf")
end

# ── Figure 3: PDE snapshots rescaled to similarity variables ───────────
function figure3()
    println("Figure 3: PDE convergence to similarity...")
    ε = 0.1
    pde = solve_pde(ε=ε, N=200, z_min=0.05, z_max=10.0,
                    t_end=1.0, n_snapshots=5)

    p = plot(xlabel="ξ = z/t^{2/3}", ylabel="S = R/t^{2/3}",
             title="PDE rescaled to similarity (ε=$ε)")

    # Plot rescaled snapshots (skip t=0)
    colors = [:lightblue, :blue, :darkblue, :purple, :red]
    for (i, t_idx) in enumerate(2:length(pde.t_snapshots))
        t = pde.t_snapshots[t_idx]
        ξ, S = rescale_to_similarity(pde, t_idx)
        ci = min(i, length(colors))
        plot!(p, ξ, S, label="t = $(round(t, digits=2))", color=colors[ci])
    end

    # Overlay the undisturbed cone S = εξ
    ξ_ref = range(0, maximum(pde.z), length=100)
    plot!(p, ξ_ref, ε .* ξ_ref, label="S = εξ (cone)", color=:gray,
          linestyle=:dash, alpha=0.5)

    savefig(p, joinpath(FIGDIR, "fig3_pde_similarity_convergence.pdf"))
    println("  Saved fig3_pde_similarity_convergence.pdf")
end

# ── Figure 4: Velocity profiles ────────────────────────────────────────
function figure4()
    println("Figure 4: Velocity profiles...")
    ε = 0.1

    inner = solve_inner_bvp(ξ₀=0.0, S₀=1.0, ξ_max=15.0, ε=ε)

    p = plot(xlabel="ξ", ylabel="U(ξ)", title="Similarity velocity profile (ε=$ε)")
    plot!(p, inner.ξ, inner.U, label="U(ξ) inner", color=:blue, linewidth=2)
    hline!(p, [0], label="", color=:gray, linestyle=:dash, alpha=0.3)

    savefig(p, joinpath(FIGDIR, "fig4_velocity_profile.pdf"))
    println("  Saved fig4_velocity_profile.pdf")
end

# ── Run all ────────────────────────────────────────────────────────────
println("Generating figures for SlenderConeRecoil...")
figure1()
figure2()
figure3()
figure4()
println("\nAll figures saved to figures/")
