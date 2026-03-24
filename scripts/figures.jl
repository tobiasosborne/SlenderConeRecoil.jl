#!/usr/bin/env julia
# Generate all figures for SlenderConeRecoil.
# Usage: julia --project scripts/figures.jl

using Plots
default(dpi=200, size=(700, 450), linewidth=2, legend=:topright,
        guidefontsize=12, tickfontsize=10, legendfontsize=9)

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

# ── Figure 1: S(ξ) near tip showing the blob ──────────────────────────
function figure1(sol)
    println("Figure 1: Similarity profile S(ξ) near tip...")
    ε = 0.1

    # Zoomed view near tip where the blob is visible
    p = plot(xlabel="ξ", ylabel="S(ξ)",
             title="Similarity profile near tip (ε=$ε)")

    idx = sol.ξ .< 15
    plot!(p, sol.ξ[idx], sol.S[idx], label="S(ξ)", color=:blue, linewidth=2.5)

    # Undisturbed cone (only for ξ > 0)
    ξ_ref = range(0, 15, length=200)
    plot!(p, ξ_ref, ε .* ξ_ref, label="εξ (undisturbed)", color=:red,
          linestyle=:dash, linewidth=1.5)

    scatter!(p, [sol.ξ₀], [sol.S₀], label="tip: ξ₀=$(round(sol.ξ₀,digits=2)), S₀=$(round(sol.S₀,digits=3))",
             color=:red, markersize=6)

    # Annotate the blob excess
    annotate!(p, sol.ξ₀ + 1.5, sol.S₀ + 0.08,
              text("blob excess\n≈ $(round(sol.S₀ - ε*sol.ξ₀, digits=3))", 8, :left))

    xlims!(p, 0, 15)
    ylims!(p, 0, 1.8)

    savefig(p, joinpath(FIGDIR, "fig1_similarity_profile.pdf"))
    println("  Saved.")
end

# ── Figure 2: Excess profile S(ξ) - εξ showing blob structure ─────────
function figure2(sol)
    println("Figure 2: Excess profile S - εξ (blob structure)...")
    ε = 0.1

    excess = sol.S .- ε .* sol.ξ

    p = plot(xlabel="ξ", ylabel="S(ξ) - εξ",
             title="Blob excess above undisturbed cone (ε=$ε)")

    plot!(p, sol.ξ, excess, label="S(ξ) - εξ", color=:blue, linewidth=2.5)
    hline!(p, [0], label="", color=:gray, linestyle=:dash, alpha=0.5)
    vline!(p, [sol.ξ₀], label="tip ξ₀", color=:red, linestyle=:dot)

    # The excess starts at S₀ - εξ₀ and decays to 0
    scatter!(p, [sol.ξ₀], [sol.S₀ - ε*sol.ξ₀], label="",
             color=:red, markersize=5)

    xlims!(p, sol.ξ₀ - 0.5, 40)
    ylims!(p, -0.005, maximum(excess)*1.2)

    savefig(p, joinpath(FIGDIR, "fig2_blob_excess.pdf"))
    println("  Saved.")
end

# ── Figure 3: Physical recoil — cone profile at multiple times ─────────
function figure3(sol)
    println("Figure 3: Physical cone recoil at multiple times...")
    ε = 0.1

    p = plot(xlabel="z", ylabel="R(z,t)",
             title="Cone tip recoil: profile at different times (ε=$ε)")

    times = [0.01, 0.05, 0.2, 0.5, 1.0, 2.0]
    colors = cgrad(:viridis, length(times), categorical=true)

    for (k, t) in enumerate(times)
        ℓ = t^(2/3)
        z_vals = sol.ξ .* ℓ
        R_vals = sol.S .* ℓ
        z_tip = sol.ξ₀ * ℓ
        R_tip = sol.S₀ * ℓ

        # Add hemispherical cap at tip (upper half: axis → profile)
        z_cap, r_cap = tip_cap_upper(z_tip, R_tip)
        z_full = vcat(z_cap, z_vals)
        R_full = vcat(r_cap, R_vals)

        idx = z_full .< 10
        plot!(p, z_full[idx], R_full[idx], label="t=$t", color=colors[k])
    end

    # Undisturbed cone for reference
    z_ref = range(0, 10, length=100)
    plot!(p, z_ref, ε .* z_ref, label="R=εz", color=:black, linestyle=:dash,
          linewidth=1, alpha=0.5)

    # Mark tip positions at each time
    for (k, t) in enumerate(times)
        ℓ = t^(2/3)
        z_tip = sol.ξ₀ * ℓ
        R_tip = sol.S₀ * ℓ
        if z_tip < 10
            scatter!(p, [z_tip], [R_tip], label="", color=colors[k],
                     markersize=4, markerstrokewidth=0)
        end
    end

    xlims!(p, 0, 10)
    ylims!(p, 0, 1.2)

    savefig(p, joinpath(FIGDIR, "fig3_recoil_profiles.pdf"))
    println("  Saved.")
end

# ── Helper: hemispherical cap at the tip ────────────────────────────────
"""
Left-facing semicircular cap connecting 1D profile (R=R_tip at z=z_tip)
to R=0 on the axis. The cap curves LEFT from z_tip to z_tip - R_tip.

Full cap (upper+lower): θ from π/2 → 0 → -π/2
  z = z_tip - R_tip·cos(θ),  r = R_tip·sin(θ)

Upper half only: θ from 0 → π/2 (r goes 0 → R_tip, z goes z_tip-R_tip → z_tip)
"""
function tip_cap(z_tip, R_tip; n=50)
    # Full cap: from top (r=+R_tip) around left to bottom (r=-R_tip)
    θ = range(π/2, -π/2, length=n)
    z_cap = z_tip .- R_tip .* cos.(θ)
    r_cap = R_tip .* sin.(θ)
    (z_cap, r_cap)
end

function tip_cap_upper(z_tip, R_tip; n=30)
    # Upper half only: from axis (r=0) to profile start (r=R_tip)
    θ = range(0, π/2, length=n)
    z_cap = z_tip .- R_tip .* cos.(θ)
    r_cap = R_tip .* sin.(θ)
    (z_cap, r_cap)
end

# ── Figure 4: 2D axisymmetric shape — zoomed on tip ───────────────────
function figure4(sol)
    println("Figure 4: Axisymmetric tip shape (zoomed)...")
    ε = 0.1
    t = 1.0
    ℓ = t^(2/3)
    z_vals = sol.ξ .* ℓ
    R_vals = sol.S .* ℓ

    z_tip = sol.ξ₀ * ℓ
    R_tip = sol.S₀ * ℓ

    # Zoom tightly on tip: from just before cap to a few R_tip beyond
    z_lo = z_tip - 2*R_tip
    z_hi = z_tip + 12*R_tip
    idx = (z_vals .>= z_tip) .& (z_vals .< z_hi)

    # Build closed profile: cap + upper + lower (reversed)
    z_cap, r_cap = tip_cap(z_tip, R_tip, n=80)
    z_upper = z_vals[idx]
    r_upper = R_vals[idx]
    z_lower = reverse(z_upper)
    r_lower = -reverse(r_upper)

    z_outline = vcat(z_cap, z_upper, z_lower)
    r_outline = vcat(r_cap, r_upper, r_lower)

    p = plot(xlabel="z", ylabel="r",
             title="Recoiled tip shape (t=$t, ε=$ε)",
             aspect_ratio=:equal, size=(700, 500))

    # Filled fluid body
    plot!(p, z_outline, r_outline, seriestype=:shape,
          label="fluid body", fillcolor=:blue, alpha=0.3,
          linecolor=:blue, linewidth=2)

    # Undisturbed cone
    z_cone = range(z_lo, z_hi, length=100)
    plot!(p, z_cone, ε .* z_cone, label="undisturbed cone",
          color=:red, linestyle=:dash, linewidth=1.5)
    plot!(p, z_cone, -ε .* z_cone, label="", color=:red, linestyle=:dash,
          linewidth=1.5)

    xlims!(p, z_lo, z_hi)
    R_max = maximum(r_upper) * 1.5
    ylims!(p, -R_max, R_max)

    savefig(p, joinpath(FIGDIR, "fig4_tip_shape.pdf"))
    println("  Saved.")
end

# ── Figure 5: Velocity profile ─────────────────────────────────────────
function figure5(sol)
    println("Figure 5: Velocity profile...")
    ε = 0.1

    p = plot(xlabel="ξ", ylabel="U(ξ)",
             title="Similarity velocity (ε=$ε): fluid pushed away from tip")
    plot!(p, sol.ξ, sol.U, label="U(ξ)", color=:blue, linewidth=2.5)
    hline!(p, [0], label="", color=:gray, linestyle=:dash, alpha=0.3)
    vline!(p, [sol.ξ₀], label="tip ξ₀", color=:red, linestyle=:dot)

    # U > 0 near tip (fluid moves in +ξ direction = recoil)
    # U → 0 far from tip
    xlims!(p, sol.ξ₀ - 0.5, 30)

    savefig(p, joinpath(FIGDIR, "fig5_velocity.pdf"))
    println("  Saved.")
end

# ── Run all ────────────────────────────────────────────────────────────
println("Generating figures for SlenderConeRecoil...")
println("Solving inner BVP (3D Newton with axial curvature)...")
sol = solve_inner_bvp(ε=0.1)
println("  Converged: ξ₀=$(round(sol.ξ₀, digits=3)), S₀=$(round(sol.S₀, digits=4))")
println()

figure1(sol)
figure2(sol)
figure3(sol)
figure4(sol)
figure5(sol)
println("\nAll figures saved to figures/")
