# Similarity reduction: apply Keller-Miksis t^{2/3} scaling to the 1D model.
#
# Keller-Miksis scaling (nondimensional, γ/ρ = 1):
#   ℓ(t) = t^{2/3}
#   z = ℓ(t)·ξ ,  R = ℓ(t)·S(ξ) ,  u = ℓ̇(t)·U(ξ)
#
# where ℓ̇ = (2/3)t^{-1/3}.
#
# The PDE system becomes an ODE system in (S(ξ), U(ξ)):
#   2S(U - (2/3)ξ)S' + S²(U' - 2/3) = 0                [mass]
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0                 [momentum]
#   (replacing the notation: S' = dS/dξ, U' = dU/dξ)

export similarity_ode_mass, similarity_ode_momentum, similarity_system

# ── Similarity variables ───────────────────────────────────────────────
const sym_ξ = Sym(:ξ)
const sym_S = Sym(:S)       # S(ξ) — similarity profile for R
const sym_U = Sym(:U)       # U(ξ) — similarity profile for u
const sym_Sξ = Sym(:Sξ)     # dS/dξ
const sym_Sξξ = Sym(:Sξξ)   # d²S/dξ²
const sym_Sξξξ = Sym(:Sξξξ) # d³S/dξ³
const sym_Uξ = Sym(:Uξ)     # dU/dξ

# ── Derivation of the similarity ODEs ──────────────────────────────────
#
# Starting from the 1D slender model:
#   (1) ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   (2) ∂u/∂t + u·∂u/∂z = -∂/∂z(1/R)   [from Bernoulli: capillary recoil]
#
# Substituting z = t^{2/3}·ξ, R = t^{2/3}·S(ξ), u = (2/3)t^{-1/3}·U(ξ):
#
# Time derivatives (at fixed z, so ξ = z·t^{-2/3} varies):
#   ∂ξ/∂t|_z = -(2/3)·ξ·t^{-1}
#   ∂R/∂t|_z = t^{2/3}·Ṡ = (2/3)t^{-1/3}·S + t^{2/3}·S'·(-(2/3)ξ·t^{-1})
#            = (2/3)t^{-1/3}·(S - ξ·S')
#   ∂u/∂t|_z = (2/3)·(-(1/3))·t^{-4/3}·U + (2/3)t^{-1/3}·U'·(-(2/3)ξ·t^{-1})
#            = -(2/9)t^{-4/3}·(U + 2ξ·U')
#
# Spatial derivatives (∂/∂z = t^{-2/3}·∂/∂ξ):
#   ∂R/∂z = S'
#   ∂u/∂z = (2/3)t^{-1}·U'
#   ∂/∂z(1/R) = -S'/(t^{2/3}·S²) · t^{-2/3} = wait, let me redo:
#   1/R = 1/(t^{2/3}·S) = t^{-2/3}/S
#   ∂/∂z(1/R) = t^{-2/3}·∂/∂ξ(t^{-2/3}/S) = t^{-2/3}·t^{-2/3}·(-S'/S²)
#             = -t^{-4/3}·S'/S²
#
# Mass equation: ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   R² = t^{4/3}·S²
#   ∂(R²)/∂t = (4/3)t^{1/3}·S² + t^{4/3}·2S·S'·(-(2/3)ξ·t^{-1})
#            = (4/3)t^{1/3}·S² - (4/3)t^{1/3}·ξ·S·S'
#            = (4/3)t^{1/3}·S(S - ξS')
#   R²u = t^{4/3}·S²·(2/3)t^{-1/3}·U = (2/3)t·S²U
#   ∂(R²u)/∂z = t^{-2/3}·∂/∂ξ((2/3)t·S²U)
#             = (2/3)t^{1/3}·(2S·S'·U + S²·U')
#
#   Sum = 0: (4/3)t^{1/3}·S(S - ξS') + (2/3)t^{1/3}·(2SS'U + S²U') = 0
#   Divide by (2/3)t^{1/3}·S (nonzero):
#     2(S - ξS')/S + (2S'U + SU')/1 = wait, let me factor more carefully.
#
#   (4/3)S(S - ξS') + (2/3)(2SS'U + S²U') = 0
#   Divide by (2/3)S:
#     2(S - ξS') + 2S'U + SU' = 0
#     2S - 2ξS' + 2S'U + SU' = 0
#     2S + 2S'(U - ξ) + SU' = 0
#
# That's equivalent to: 2S(U - ξ)S'/S + ... hmm. Let me just write the
# canonical form as an ODE for S' and U'.
#
# Momentum equation: ∂u/∂t + u·∂u/∂z = -∂/∂z(1/R)
#   LHS: -(2/9)t^{-4/3}(U + 2ξU') + (4/9)t^{-4/3}UU'
#   RHS: -∂/∂z(1/R) = t^{-4/3}·S'/S²
#
#   Dividing by t^{-4/3}:
#   -(2/9)(U + 2ξU') + (4/9)UU' = S'/S²
#   -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0

# ── The similarity ODE system ──────────────────────────────────────────
#
# After the similarity substitution and simplification, the ODE system is:
#
#   Mass:     2S'(U - ξ) + S(U' - 2) + 2S = 0
#             i.e. 2S + 2S'(U - ξ) + SU' - 2S = 0 ... wait
#
# Let me be more careful. From the derivation above:
#   2S + 2S'(U - ξ) + SU' = 0
#
# So: SU' = -2S - 2S'(U - ξ) = -2(S + S'(U - ξ))
#     U' = -2(1 + S'(U - ξ)/S)   ... if S ≠ 0
#     U' = -2 - 2S'(U - ξ)/S
#
# Momentum: -(2/9)U - (4/9)ξU' + (4/9)UU' + S'/S² = 0
#   = -(2/9)U + (4/9)(U - ξ)U' + S'/S² = 0
#
# Multiply everything by 9/4:
#   -(1/2)U + (U - ξ)U' + (9/4)S'/S² = 0
#
# Hmm, that doesn't simplify as nicely. Let me re-derive more carefully
# using the substitution u = (2/3)ℓ̇·U where ℓ̇ = dℓ/dt.
#
# Actually, the standard convention in Decent & King uses:
#   u = ż = (2/3)t^{-1/3} times a velocity variable.
#
# For the code, let me just define the system in the form that will be
# solved numerically. The key expressions:

"""
    similarity_ode_mass()

Mass conservation ODE in similarity variables:
  2S + 2S'(U - ξ) + SU' = 0

where S' = dS/dξ, U' = dU/dξ.
Returns the LHS expression (= 0).
"""
function similarity_ode_mass()
    # 2S + 2Sξ(U - ξ) + S·Uξ = 0
    add(
        mul(Num(2), sym_S),
        mul(Num(2), sym_Sξ, add(sym_U, neg(sym_ξ))),
        mul(sym_S, sym_Uξ)
    )
end

"""
    similarity_ode_momentum(; axial=false)

Momentum ODE in similarity variables.

Leading order:  -(2/9)U + (4/9)(U - ξ)U' - S'/S² = 0
With axial:     -(2/9)U + (4/9)(U - ξ)U' - S'/S² - S''' = 0

The -S''' term produces capillary wave dispersion.
"""
function similarity_ode_momentum(; axial::Bool=false)
    expr = add(
        mul(Num(-2//9), sym_U),
        mul(Num(4//9), add(sym_U, neg(sym_ξ)), sym_Uξ),
        neg(mul(sym_Sξ, pow(sym_S, Num(-2))))
    )
    if axial
        expr = add(expr, neg(sym_Sξξξ))
    end
    expr
end

"""
    similarity_system()

Return the similarity ODE system as a named tuple:
  (mass=..., momentum=...)
"""
function similarity_system()
    (mass=similarity_ode_mass(), momentum=similarity_ode_momentum())
end

# ── Verify: t cancels in the similarity substitution ───────────────────
"""
    verify_t_cancels()

Substitute the similarity ansatz into the 1D PDE and verify that all
time dependence cancels, leaving only ODEs in ξ.

This is done numerically: evaluate the PDE residual at several values
of t and check that the ratio is independent of t.
"""
function verify_t_cancels()
    # Numerically verify that the similarity substitution into the PDE
    # gives a residual proportional to a power of t (i.e., t cancels
    # from the ODE after division).
    #
    # Mass PDE residual = t^{1/3} · (2/3) · S · [mass ODE]
    # Momentum PDE residual = t^{-4/3} · [momentum ODE expression]
    #
    # So: (PDE residual at t₁) / (PDE residual at t₂) should equal
    #     (t₁/t₂)^{power} for any (S, U, S', U', ξ).

    S, Sξ, U, Uξ, ξ = 2.0, 0.5, 1.0, -0.3, 1.5

    function mass_pde_residual(t)
        # R = t^{2/3}S, Rt = (2/3)t^{-1/3}(S - ξSξ)
        # Rz = Sξ, u = (2/3)t^{-1/3}U, uz = (2/3)t^{-1}Uξ
        R = t^(2/3) * S
        Rt = (2/3)*t^(-1/3) * (S - ξ*Sξ)
        Rz = Sξ
        u = (2/3)*t^(-1/3) * U
        uz = (2/3)*t^(-1) * Uξ
        2R*Rt + 2R*Rz*u + R^2*uz
    end

    function mom_pde_residual(t)
        u_val = (2/3)*t^(-1/3) * U
        ut = -(2/9)*t^(-4/3) * (U + 2ξ*Uξ)
        uz = (2/3)*t^(-1) * Uξ
        # ∂/∂z(1/R) = -t^{-4/3} Sξ/S²
        dinvR_dz = -t^(-4/3) * Sξ / S^2
        # PDE: ut + u·uz = -∂/∂z(1/R), residual = ut + u·uz + ∂/∂z(1/R)
        ut + u_val*uz + dinvR_dz
    end

    t1, t2 = 1.0, 2.0
    # Mass: residual ∝ t^{1/3}
    r1 = mass_pde_residual(t1)
    r2 = mass_pde_residual(t2)
    mass_ok = abs(r1/r2 - (t1/t2)^(1/3)) < 1e-10

    # Momentum: residual ∝ t^{-4/3}
    r1m = mom_pde_residual(t1)
    r2m = mom_pde_residual(t2)
    mom_ok = abs(r1m/r2m - (t1/t2)^(-4/3)) < 1e-10

    mass_ok && mom_ok
end
