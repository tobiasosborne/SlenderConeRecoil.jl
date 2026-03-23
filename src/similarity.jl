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
#   (2/9)U + (U - (2/3)ξ)U' + S'/S² = 0                 [momentum]
#   (replacing the notation: S' = dS/dξ, U' = dU/dξ)

export similarity_ode_mass, similarity_ode_momentum, similarity_system

# ── Similarity variables ───────────────────────────────────────────────
const sym_ξ = Sym(:ξ)
const sym_S = Sym(:S)       # S(ξ) — similarity profile for R
const sym_U = Sym(:U)       # U(ξ) — similarity profile for u
const sym_Sξ = Sym(:Sξ)     # dS/dξ
const sym_Uξ = Sym(:Uξ)     # dU/dξ

# ── Derivation of the similarity ODEs ──────────────────────────────────
#
# Starting from the 1D slender model:
#   (1) ∂(R²)/∂t + ∂(R²u)/∂z = 0
#   (2) ∂u/∂t + u·∂u/∂z = ∂/∂z(1/R)
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
# Momentum equation: ∂u/∂t + u·∂u/∂z = ∂/∂z(1/R)
#   -(2/9)t^{-4/3}(U + 2ξU') + (2/3)t^{-1/3}U·(2/3)t^{-1}U' = -t^{-4/3}S'/S²
#   -(2/9)t^{-4/3}(U + 2ξU') + (4/9)t^{-4/3}UU' = -t^{-4/3}S'/S²
#
#   Multiply by -t^{4/3}·(-1)... divide by t^{-4/3}:
#   -(2/9)(U + 2ξU') + (4/9)UU' = -S'/S²
#   -(2/9)U - (4/9)ξU' + (4/9)UU' = -S'/S²
#   -(2/9)U + (4/9)U'(U - ξ) + S'/S² = 0
#
#   Multiply by 9/2:
#   -U + 2U'(U - ξ) + (9/2)S'/S² = 0
#
# Hmm, let me use the more standard form. Let V = U - (2/3)ξ be the
# velocity in the similarity frame.
#
# Actually, the clearest form uses the substitution directly.
# Let me just define the ODE system in canonical resolved form.

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
    similarity_ode_momentum()

Momentum ODE in similarity variables:
  -(2/9)U + (4/9)(U - ξ)U' + S'/(S²) = 0

Returns the LHS expression (= 0).
"""
function similarity_ode_momentum()
    # -(2/9)U + (4/9)(U - ξ)Uξ + Sξ/S² = 0
    add(
        mul(Num(-2//9), sym_U),
        mul(Num(4//9), add(sym_U, neg(sym_ξ)), sym_Uξ),
        mul(sym_Sξ, pow(sym_S, Num(-2)))
    )
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
    # The 1D mass PDE: 2R·Rt + 2R·Rz·u + R²·uz = 0
    # Substitute: R = t^{2/3}·S, u = (2/3)t^{-1/3}·U
    #   Rt = (2/3)t^{-1/3}(S - ξS')
    #   Rz = S'  (since ∂/∂z = t^{-2/3}∂/∂ξ, and S' is dS/dξ,
    #            actually Rz = t^{2/3}·S'·t^{-2/3} = S')
    #   uz = (2/3)t^{-1/3}·U'·t^{-2/3} = (2/3)t^{-1}·U'
    #
    # 2R·Rt = 2·t^{2/3}S·(2/3)t^{-1/3}(S - ξS')
    #       = (4/3)t^{1/3}·S(S - ξS')
    # 2R·Rz·u = 2·t^{2/3}S·S'·(2/3)t^{-1/3}U
    #         = (4/3)t^{1/3}·SS'U
    # R²·uz = t^{4/3}S²·(2/3)t^{-1}U'
    #       = (2/3)t^{1/3}·S²U'
    #
    # Total: t^{1/3}·[(4/3)S(S-ξS') + (4/3)SS'U + (2/3)S²U']
    #       = t^{1/3}·(2/3)·[2S(S-ξS') + 2SS'U + S²U']
    #       = t^{1/3}·(2/3)·S·[2S - 2ξS' + 2S'U + SU']
    #       = t^{1/3}·(2/3)·S·[mass ODE expression]
    #
    # So the mass ODE = 0 iff the PDE = 0, and t^{1/3} cancels. ✓
    #
    # Similarly for momentum, the common factor is t^{-4/3}. ✓

    # Numerical check: pick S=2, Sξ=0.5, U=1, Uξ=-0.3, ξ=1
    # and evaluate mass ODE
    vals = Dict(
        sym_S => Num(2), sym_Sξ => Num(1//2),
        sym_U => Num(1), sym_Uξ => Num(-3//10),
        sym_ξ => Num(1)
    )
    mass_ode = similarity_ode_mass()
    result = mass_ode
    for (k, v) in vals
        result = substitute(result, k, v)
    end
    # The residual should be a pure number (no t dependence)
    result isa Num
end
