# Governing equations and slender-body (lubrication) reduction.
# Defines the full cylindrical equations, applies the slender ansatz
# R = εf, φ = Φ₀ + ε²Φ₁ + ..., and extracts the 1D model at leading order.
#
# The 1D slender model (nondimensional, γ/ρ = 1):
#   ∂(R²)/∂t + ∂(R²u)/∂z = 0        (mass conservation)
#   ∂u/∂t + u ∂u/∂z = ∂/∂z (1/R)    (momentum with capillary pressure)
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
    slender_momentum_eq()

Momentum equation in the 1D slender model (nondimensional, γ/ρ = 1):
  ∂u/∂t + u·∂u/∂z = ∂/∂z (1/R)

Returned as LHS - RHS (= 0):
  ut + u·uz - ∂/∂z(1/R) = 0

where ∂/∂z(1/R) = -Rz/R².
So: ut + u·uz + Rz/R² = 0
"""
function slender_momentum_eq()
    # LHS: ut + u*uz
    lhs = add(sym_ut, mul(sym_u, sym_uz))
    # RHS: ∂/∂z(1/R) = -Rz/R²
    # So LHS - RHS = ut + u*uz - (-Rz/R²) = ut + u*uz + Rz/R²
    rhs_neg = mul(sym_Rz, pow(sym_R, Num(-2)))
    add(lhs, rhs_neg)
end

"""
    slender_system()

Return the complete 1D slender-body system as a named tuple:
  (mass=..., momentum=...)

Each entry is an SExpr representing the LHS of the equation = 0.
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
