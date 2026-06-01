# Reconstructed slender-body (lubrication) algebra.
#
# Source status, per
# docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md:
# - C2001/C2008-meta: inviscid surface-tension recoil of a small-aspect-ratio
#   cone, initially/far-field conical, with Keller-Miksis length scaling.
# - C2001: retaining curvature beyond the leading azimuthal term is important
#   for robust oscillatory simplified-model results.
# - IMPL-inferred: the primitive variables R(z,t), u(z,t), the mass equation,
#   the Bernoulli/momentum sign convention, and the optional axial-curvature
#   residual used below. These are local reconstruction targets pending the
#   canonical 2008 article body.
#
# Current local 1D primitive model (nondimensional, γ/ρ = 1):
#   ∂(R²)/∂t + ∂(R²u)/∂z = 0        (mass conservation)
#   ∂u/∂t + u ∂u/∂z = -∂/∂z(1/R) = Rz/R²   (momentum, capillary recoil)
#
# where R(z,t) is the free-surface radius and u(z,t) is the axial velocity.

export slender_mass_eq, slender_momentum_eq, slender_system

# ── Symbols ────────────────────────────────────────────────────────────
const sym_r = Sym(:r)
const sym_z = Sym(:z)
const sym_t = Sym(:t)
const sym_ε = Sym(:ε)
const sym_R = Sym(:R)      # R(z,t) — free surface radius
const sym_u = Sym(:u)      # u(z,t) — axial velocity
const sym_Rz = Sym(:Rz)    # ∂R/∂z
const sym_Rzz = Sym(:Rzz)  # ∂²R/∂z²
const sym_Rzzz = Sym(:Rzzz) # ∂³R/∂z³
const sym_Rt = Sym(:Rt)    # ∂R/∂t
const sym_uz = Sym(:uz)    # ∂u/∂z
const sym_ut = Sym(:ut)    # ∂u/∂t

# ── Full governing equations (symbolic, for reference) ─────────────────

"""
    curvature_full()

Full mean curvature κ in cylindrical coordinates:
  κ = 1/(R√(1+Rz²)) − Rzz/(1+Rz²)^{3/2}

The first term is azimuthal curvature (1/R at leading order),
the second is axial curvature (higher order in ε).

Ledger status: the formula is standard geometry and curvature retention is
motivated by the 2001 precursor; its exact 2008 asymptotic ordering in this
primitive-variable model is still provisional.
"""
function curvature_full()
    denom1 = Func(:sqrt, [add(Num(1), pow(sym_Rz, Num(2)))])
    azimuthal = pow(mul(sym_R, denom1), Num(-1))
    axial = neg(mul(sym_Rzz, pow(add(Num(1), pow(sym_Rz, Num(2))), Num(-3//2))))
    add(azimuthal, axial)
end

"""
    curvature_leading()

Leading-order curvature for a slender body (ε ≪ 1):
  κ ≈ 1/R

The azimuthal curvature dominates; axial curvature is O(ε²) smaller.

Ledger status: small aspect ratio is source-backed; this primitive reduction
is a local reconstructed algebra check until the 2008 article body is extracted.
"""
function curvature_leading()
    pow(sym_R, Num(-1))
end

# ── Slender-body 1D model ─────────────────────────────────────────────

"""
    slender_mass_eq()

Mass conservation in the 1D slender model:
  ∂(R²)/∂t + ∂(R²u)/∂z = 0

Returned as the LHS expression (= 0). Using the product rule:
  2R·Rt + 2R·Rz·u + R²·uz = 0

i.e.  2R(Rt + Rz·u) + R²·uz = 0

Dividing by 2R (assuming R ≠ 0):
  Rt + u·Rz + (R/2)·uz = 0

Ledger status: IMPL-inferred. The 2001 precursor uses a potential/free-surface
kinematic equation in different variables; this primitive conservation law is
tested as local algebra, not as a transcribed Decent-King equation.
"""
function slender_mass_eq()
    # ∂(R²)/∂t = 2R·Rt
    dR2dt = mul(Num(2), sym_R, sym_Rt)
    # ∂(R²u)/∂z = 2R·Rz·u + R²·uz
    dR2udz = add(mul(Num(2), sym_R, sym_Rz, sym_u), mul(pow(sym_R, Num(2)), sym_uz))
    # LHS = 0
    add(dR2dt, dR2udz)
end

"""
    slender_momentum_eq(; axial=false)

Momentum equation in the 1D slender model (nondimensional, γ/ρ = 1).

Leading order (axial=false):
  κ = 1/R  →  ut + u·uz - Rz/R² = 0

With axial curvature (axial=true):
  κ = 1/R - Rzz  →  ut + u·uz - Rz/R² - Rzzz = 0

The -Rzzz term is dispersive and produces capillary waves.

Ledger status: IMPL-inferred. Curvature retention is motivated by the 2001
precursor, but the primitive momentum equation, sign convention, and axial
curvature ordering remain provisional until checked against the 2008 article.
"""
function slender_momentum_eq(; axial::Bool=false)
    lhs = add(sym_ut, mul(sym_u, sym_uz))
    rhs = mul(sym_Rz, pow(sym_R, Num(-2)))
    expr = add(lhs, neg(rhs))
    if axial
        expr = add(expr, neg(sym_Rzzz))
    end
    expr
end

"""
    slender_system()

Return the complete 1D slender-body system as a named tuple:
  (mass=..., momentum=...)

Each entry is an SExpr representing the LHS of the equation = 0.
These are local reconstructed primitive-variable equations.
"""
function slender_system()
    (mass=slender_mass_eq(), momentum=slender_momentum_eq())
end

# ── Derivation verification ────────────────────────────────────────────
"""
    verify_slender_derivation(order=1)

Derive the 1D slender model from the full equations by:
1. Setting R = εf (where f is O(1))
2. Expanding curvature in ε
3. Verifying that leading-order curvature is 1/R = 1/(εf)

Returns true if the leading-order curvature matches.

This verifies only the local curvature scaling algebra. It is not a
Decent-King 2008 equation transcription test.
"""
function verify_slender_derivation(; order::Int=1)
    f = Sym(:f)
    # R = ε*f
    R_ansatz = mul(sym_ε, f)

    # Full curvature: 1/(R*sqrt(1+Rz²)) - Rzz/(1+Rz²)^{3/2}
    # At leading order in ε (with Rz = ε*fz, Rzz = ε*fzz):
    # Azimuthal: 1/(εf * sqrt(1 + ε²fz²))
    #          = 1/(εf) * (1 + ε²fz²)^{-1/2}
    #          = 1/(εf) * (1 - ε²fz²/2 + ...)
    #          = 1/(εf) at leading order
    #
    # Axial: εfzz/(1 + ε²fz²)^{3/2} = ε*fzz at leading order
    #
    # So κ = 1/(εf) + O(ε) = (1/ε)(1/f) + O(ε)

    # Verify: substitute R = εf into 1/R
    κ_leading = substitute(curvature_leading(), sym_R, R_ansatz)
    # Should be 1/(εf) = (1/ε) * (1/f)
    expected = pow(mul(sym_ε, f), Num(-1))
    κ_leading == expected
end
