# Decent-King Source Ledger

Date: 2026-06-01

Scope: Beads issue `scr-8l4.2` / `UP-SF2`. This ledger records what is
source-backed now, what is only metadata-backed for the 2008 target paper,
what is reconstructed in the current implementation, and what remains blocked
until the 2008 PDF is acquired.

## Status Legend

| Status | Meaning |
| --- | --- |
| `C2001` | Confirmed from Decent and King's 2001 IUTAM precursor chapter, local PDF `docs/papers/DecentKing2001_IUTAM.pdf`, printed pp. 81-88, PDF pages 93-100. |
| `C2008-meta` | Confirmed from DOI/institutional metadata or abstract-level records for Decent and King (2008), not from the article body. |
| `IMPL-inferred` | Present in this repository's current equations, solvers, tests, or documentation, but not yet confirmed against the canonical 2008 article. |
| `BLOCKED-2008` | Requires the canonical 2008 IMA Journal article body, DOI `10.1093/imamat/hxm043`, which is still unavailable locally because noninteractive access hit Cloudflare/403 barriers. |
| `OCR-needs-check` | Visible in the 2001 PDF extraction, but equation typography is degraded enough that any exact symbolic transcription must be checked visually before becoming test data. |

## Sources Inspected

| Source | Status | Notes |
| --- | --- | --- |
| S. P. Decent and A. C. King, "The Recoil of A Broken Liquid Bridge," in *IUTAM Symposium on Free Surface Flows*, 2001. DOI <https://doi.org/10.1007/978-94-010-0796-2_10>. | `C2001` | Local PDF present. Relevant chapter is printed pp. 81-88, PDF pages 93-100. Text extracted with `pdftotext -layout` and `pdftotext -raw`. |
| S. P. Decent and A. C. King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* 73(1), 37-68, 2008. DOI <https://doi.org/10.1093/imamat/hxm043>. | `C2008-meta` / `BLOCKED-2008` | Metadata/abstract verified from institutional records: Birmingham (<https://research.birmingham.ac.uk/en/publications/surface-tension-driven-flow-in-a-slender-cone>) and Dundee (<https://discovery.dundee.ac.uk/en/publications/surface-tension-driven-flow-in-a-slender-cone/>). Article-body equations, page numbers, figures, and numerical constants remain blocked. |
| `docs/papers/README.md` | `C2008-meta` / access record | Records expected filename `DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf` and the current Cloudflare/403 acquisition blocker. |

The 2001 chapter explicitly says it studies a simplified one-dimensional model
and uses that model to show how the full axisymmetric cone calculation can be
approached. The canonical package target must remain the 2008 IMA Journal
article. Do not promote 2001-only constants or current solver outputs to
"Decent-King 2008 benchmarks" without the 2008 extraction.

## Confirmed Metadata And Abstract-Level Facts For 2008

| Fact | Status | Source |
| --- | --- | --- |
| Correct canonical paper is Decent and King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* 73(1), 37-68, 2008. | `C2008-meta` | DOI and Birmingham/Dundee records. |
| DOI is `10.1093/imamat/hxm043`. | `C2008-meta` | DOI URL and institutional metadata. |
| The paper studies recoil after a droplet breaks from a slender thread or jet whose local tip is conical near break-off. | `C2008-meta` | Institutional abstract, paraphrased. |
| The model is an ideal fluid, initially conical, driven by surface tension, with small cone aspect ratio. | `C2008-meta` | Institutional abstract, paraphrased. |
| The paper uses a small-time similarity transformation and identifies a rapidly oscillating nonlinear wave moving away from the tip. | `C2008-meta` | Institutional abstract, paraphrased. |
| Exact nondimensionalization, equation numbers, asymptotic powers, boundary/tip conditions, matching constants, numerical values, and figure data. | `BLOCKED-2008` | Requires the article body. |

## 2001 Precursor Equation And Assumption Ledger

### Physical Problem And Assumptions

| Quantity or assumption | Status | 2001 reference | Current package mapping |
| --- | --- | --- | --- |
| Cylindrical coordinates `(z, r, theta, t)` and axisymmetric velocity potential `Phi(r,z,t)`. | `C2001` | Section 2, printed p. 82 / PDF p. 94. | `src/slender.jl` symbolic variables and comments describe an axisymmetric free surface, but the package evolves primitive variables `R(z,t), u(z,t)` rather than the potential. |
| Free surface `r = mu(z,t)`, tip at `r = 0`, `z = z0(t)`, with `z0(0)=0`. | `C2001` | Section 2, printed pp. 82-83 / PDF pp. 94-95. | Current inner solver uses a constant similarity-tip location `xi0`; PDE uses a truncated left boundary rather than a moving physical tip. |
| Initial cone `mu = epsilon z` at `t = 0`. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | `solve_pde` uses `R0 = epsilon .* z`; inner/outer far-field conditions use `S ~ epsilon xi`. |
| Far-field cone `mu ~ epsilon z` as `z -> infinity`; flow initially stationary and stationary in the far field. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | `solve_inner_bvp` targets far-field slope `epsilon`, velocity `U -> 0`, and curvature decay; `solve_outer` imposes small/zero perturbation at large `xi`. |
| Cone aspect ratio `epsilon` is small. | `C2001` and `C2008-meta` | 2001 abstract and Section 2; 2008 abstract metadata. | `epsilon` appears throughout slender, similarity, outer, composite, PDE, and tests. |
| Ideal fluid, irrotational flow, driven by surface tension alone. | `C2001` and `C2008-meta` | 2001 Section 1, printed p. 82 / PDF p. 94; 2008 abstract metadata. | Package equations are inviscid, no gravity, no viscosity. |
| Gravity neglected because the recoil time scale is short. | `C2001` | Section 1, printed p. 82 / PDF p. 94. | No gravity term in current PDE or similarity equations. |
| Viscosity treated as negligible outside a very small length scale for water. | `C2001` | Section 1, printed pp. 81-82 / PDF pp. 93-94. | No viscous terms implemented. |
| Full curvature must be retained in the 2001 simplified model for robust oscillatory results. | `C2001` | Note 2, printed p. 88 / PDF p. 100. | Current implementation includes `kappa = 1/R - Rzz` and resulting `Rzzz`/`S'''` terms, but the exact 2008 ordering of this term is blocked. |

### Similarity Variables And One-Dimensional Model

| Quantity or equation | Status | 2001 reference | Current package mapping |
| --- | --- | --- | --- |
| Keller-Miksis length scale `(sigma t^2 / rho)^(1/3)` removes time from the problem. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | `src/similarity.jl` uses nondimensional `ell(t)=t^(2/3)`. `src/pde.jl` rescales snapshots with `t^(2/3)`. |
| Similarity coordinates: axial and radial lengths scaled by `(sigma t^2 / rho)^(1/3)`. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | Code variables `xi`, `S(xi)` are equivalent at nondimensional `sigma/rho=1`. |
| Potential scale `(sigma^2 t / rho^2)^(1/3)` for the velocity potential. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | Current package does not expose potential `phi`; it uses axial velocity `u` and similarity velocity `U`. |
| Time-independent tip location in similarity variables, denoted `xi0`/`(0` in OCR. | `C2001` | Section 2, printed p. 83 / PDF p. 95. | `InnerSolution.xi0` and `solve_inner_bvp` shooting parameter. |
| 2001 one-dimensional approximation uses potential amplitude `A(y)`, free-surface amplitude `R(y)`, small-aspect-ratio rescalings in powers of `epsilon^(1/3)`, and constants `y0, y1, y2, ...`. | `C2001`, `OCR-needs-check` | Section 2, printed p. 83 / PDF p. 95. | Current code does not use the 2001 `A(y), R(y)` system directly. It reconstructs a primitive-variable PDE/ODE system. |
| Equation (1): kinematic condition for the 2001 one-dimensional model. | `C2001`, `OCR-needs-check` | Section 2, printed p. 83 / PDF p. 95. | Current `slender_mass_eq` and `similarity_ode_mass` are not a direct transcription of 2001 equation (1); source reconciliation is needed. |
| Equation (2): Bernoulli equation for the 2001 one-dimensional model, with curvature. | `C2001`, `OCR-needs-check` | Section 2, printed p. 83 / PDF p. 95. | Current `slender_momentum_eq`, `similarity_ode_momentum`, and axial-curvature terms are reconstructed. Exact sign/order must be checked against 2008. |
| Near-tip coordinate expansion for equations (1)-(2) gives `R ~ sqrt(2 p0 y)` and `A ~ phi0 + (y0/2)y` as `y -> 0`. | `C2001`, `OCR-needs-check` | Section 2, printed p. 84 / PDF p. 96. | Current inner tip uses `S(xi0)=S0>0`, `S'(xi0)=0`, `U(xi0)=4 xi0/5`; this is not the same 2001 one-dimensional singular-tip condition. |
| Far-field for the 2001 one-dimensional model: `R ~ y + y0`, `A ~ -1/(y+y0)` as `y -> infinity`. | `C2001`, `OCR-needs-check` | Section 2, printed p. 84 / PDF p. 96. | Current far-field uses `S/xi -> epsilon`, `U -> 0`, `S'' -> 0`; mapping to the 2001 `A,R,y` variables is unresolved. |
| Direct shooting of 2001 equations (1)-(2) determines constants `phi0`, `y0`, `p0`; reported limiting values include `y0 -> 1.50` and `p0 -> 2.65` as the truncation parameter tends to zero. | `C2001`, `OCR-needs-check` | Section 2, printed p. 84 / PDF p. 96. | These constants are not represented in current tests. They should not be compared directly to current `xi0,S0,S''0`. |
| Figure 1 shows numerical `R` and `A` for `epsilon=0.01`, including high-frequency oscillations modulated over a longer scale. | `C2001` | Figure 1, printed p. 84 / PDF p. 96. | Current tests check oscillations in `S - epsilon xi`, but against local implementation values only. |

### 2001 Inner/Outer/Matching Asymptotics

| Quantity or equation | Status | 2001 reference | Current package mapping |
| --- | --- | --- | --- |
| One-dimensional model breaks down as `y -> 0` because the potential expansion becomes nonuniform; true leading flow near the tip is three-dimensional axisymmetric. | `C2001` | Section 3, printed pp. 84-85 / PDF pp. 96-97. | Current inner solver is still a one-dimensional similarity ODE in `S,U`; it is not the full 3D axisymmetric potential calculation. |
| Inner and outer regions are introduced close to and far from the moving tip. | `C2001` | Section 3, printed p. 85 / PDF p. 97. | `src/inner.jl`, `src/outer.jl`, `src/composite.jl`. |
| Inner rescaling in 2001 uses `y = epsilon * ybar` and multiple scales. | `C2001`, `OCR-needs-check` | Section 3, printed p. 85 / PDF p. 97. | Current inner coordinate is `xi` directly; no 2001 multiple-scale coordinate is implemented. |
| Inner expansions use `A = phi0 + epsilon A1 + epsilon^2 A2 + ...` and `R = R0 + epsilon R1 + ...`, with both fast and slow variables. | `C2001`, `OCR-needs-check` | Section 3, printed p. 85 / PDF p. 97. | No direct implementation; current `S,U` BVP is an `IMPL-inferred` reconstruction. |
| Equation (3): leading-order 2001 inner equations for `R0` and `A1`, including full nonlinear curvature form. | `C2001`, `OCR-needs-check` | Section 3, printed p. 85 / PDF p. 97. | Current `inner_rhs!` is not equation (3). It solves a third-order primitive similarity ODE with axial curvature. |
| Inner initial behavior: `R0 ~ sqrt(2 p0 y)` and `A1 ~ (y0/2)y` as the fast/inner coordinate tends to zero. | `C2001`, `OCR-needs-check` | Section 3, printed p. 85 / PDF p. 97. | Not used by current inner tip conditions. |
| Inner solution is interpreted as oscillation around a slowly drifting phase-plane centre; over the longer scale it behaves like a stable spiral. | `C2001` | Section 3, printed p. 85 / PDF p. 97. | Current tests check oscillatory capillary-wave shape but not the 2001 phase-plane/secular structure. |
| Inner far-field uses `R0 ~ Rbar0(y) + small oscillatory correction` and `A1 ~ Abar1(y)Y + small oscillatory correction`; it introduces amplitude/phase functions including `b(y)` and `a(y)`. | `C2001`, `OCR-needs-check` | Section 3, printed pp. 85-86 / PDF pp. 97-98. | Current outer matching uses raw interpolated `S,S',S'',U` at `xi_match`, not the 2001 analytic amplitude/phase matching constants. |
| Outer distinguished limit expands `A` in oscillatory harmonics with powers including `epsilon^(3/2)`, `epsilon^2`, and `epsilon^(5/2)`. | `C2001`, `OCR-needs-check` | Equation (4), printed p. 86 / PDF p. 98. | Current `src/outer.jl` uses a single linear perturbation `S = epsilon xi + s1`, `U = u1`; `src/outer_hierarchy.jl` explores powers but has not been reconciled to 2001/2008. |
| Outer surface expansion is equation (5), paired with equation (4). | `C2001`, `OCR-needs-check` | Equation (5), printed p. 86 / PDF p. 98. | Current outer `s1,u1` variables are not source-confirmed against this expansion. |
| Equations (6)-(7): leading outer equations for coefficients such as `B0` and `S0`. | `C2001`, `OCR-needs-check` | Printed pp. 86-87 / PDF pp. 98-99. | Current outer equations differ in variables and must be reconciled. |
| Equations (8)-(9): secularity conditions for outer oscillatory coefficients. | `C2001`, `OCR-needs-check` | Printed p. 87 / PDF p. 99. | No direct implementation. Potential future validation target after 2008 extraction. |
| Matching uses an intermediate coordinate `y = epsilon^M yhat`, `0 < M < 1`, with small-outer/large-inner overlap. | `C2001`, `OCR-needs-check` | Printed pp. 87-88 / PDF pp. 99-100. | Current composite uses additive matching and a linear overlap fit, not the 2001 intermediate-coordinate formula. |
| Matching conditions include `B0(0)=phi0`, `B0'(0)=A1(0)=2y0/3`, small-`y` behavior for `R0`/`S0`, amplitude/phase constants `r0`, `theta0`, and matching-determined `p0,y0,y1`. | `C2001`, `OCR-needs-check` | Printed pp. 87-88 / PDF pp. 99-100. | Not represented in current API or tests. These quantities must be extracted cleanly from 2008 before becoming regression targets. |
| Figure 2 shows inner oscillations and outer coefficients. | `C2001` | Figure 2, printed p. 86 / PDF p. 98. | Current figures are package-generated and are not source-fidelity reproductions yet. |

## Current Implementation Ledger

| Package concept | Current location | Status | Source-fidelity note |
| --- | --- | --- | --- |
| Full curvature `kappa = 1/(R sqrt(1+Rz^2)) - Rzz/(1+Rz^2)^(3/2)` and leading curvature `1/R`. | `src/slender.jl:29-55`; `test/test_bead3.jl:6-12`. | `IMPL-inferred`; curvature retention motivated by `C2001`. | The curvature formula is standard and consistent with 2001's warning to retain full curvature in the simplified model, but the 2008 asymptotic order/sign must still be checked. |
| Slender mass equation `d(R^2)/dt + d(R^2 u)/dz = 0`, equivalent primitive form `Rt = -u Rz - (R/2) uz`. | `src/slender.jl:59-80`; `src/pde.jl:179-181`; tests `test/test_bead3.jl:14-24`. | `IMPL-inferred`. | Not directly extracted from 2008. 2001 equation (1) is a potential/free-surface kinematic equation in different variables. |
| Slender momentum equation `ut + u uz - Rz/R^2 = 0`, optionally with axial curvature term `-Rzzz` in residual form. | `src/slender.jl:82-103`; `src/pde.jl:179-215`; tests `test/test_bead3.jl:26-47`. | `IMPL-inferred`. | Sign and axial-curvature ordering are critical. Current comments/tests are local reconstruction, not source-confirmed. |
| Keller-Miksis scaling `ell(t)=t^(2/3)`, `z=ell xi`, `R=ell S`, `u=ell_dot U`. | `src/similarity.jl:1-7`, `src/similarity.jl:31-77`; `src/pde.jl:275-292`; tests `test/test_bead4.jl:31-33`. | Length scaling `C2001`; velocity convention `IMPL-inferred`. | 2001 confirms length/potential scaling, not this exact primitive velocity normalization. |
| Similarity mass ODE `2S + 2S'(U-xi) + S U' = 0`. | `src/similarity.jl:107-123`; `src/inner.jl:41-54`; tests `test/test_bead4.jl:35-43`, `test/test_numerical_regressions.jl:21-31`. | `IMPL-inferred`. | Must be reconciled against 2008 equations and notation. |
| Similarity momentum ODE `-(2/9)U + (4/9)(U-xi)U' - S'/S^2 - S''' = 0` with axial curvature. | `src/similarity.jl:125-145`; `src/inner.jl:56-58`; tests `test/test_bead4.jl:45-74`, `test/test_bead5.jl:66-82`. | `IMPL-inferred`. | No equation-number source yet. 2001 has different potential equations and says full curvature matters. |
| Inner state vector `[S, S', S'', U]`. | `src/inner.jl:17-31`, `src/inner.jl:37-63`; tests `test/test_bead5.jl`. | `IMPL-inferred`. | Not a direct 2001 variable set; blocked pending 2008. |
| Tip regularity conditions `S'(xi0)=0`, `U(xi0)=4 xi0/5`, free `S0`, free `S''0`, and derived `S'''(xi0)=0`. | `src/inner.jl:65-78`; tests `test/test_bead5.jl:7-14`. | `IMPL-inferred`; `BLOCKED-2008` for source confirmation. | This is one of the highest-priority 2008 extraction items. |
| Inner shooting parameters `(xi0, S0, S''0)` and far-field residuals `S/xi -> epsilon`, `U -> 0`, `S'' -> 0`. | `src/inner.jl:107-126`; tests `test/test_bead5.jl:31-51`. | `IMPL-inferred`. | Current far-field conditions are plausible but not yet source-backed. |
| Outer linearization about `S = epsilon xi`, `U = 0`; state `[s1, s1', s1'', u1]`; driven equation with `1/(epsilon xi^2)` forcing. | `src/outer.jl:1-14`, `src/outer.jl:34-70`; tests `test/test_bead6.jl`. | `IMPL-inferred`; `BLOCKED-2008` for equation/power structure. | 2001 outer expansions use oscillatory coefficients and half-integer powers in different variables. |
| `solve_outer` integrates inward from large `xi`; `solve_outer_matched` extracts perturbation data from the inner solution at `xi_match` and integrates outward. | `src/outer.jl:72-141`; tests `test/test_bead6.jl`, `test/test_numerical_regressions.jl:120-144`. | `IMPL-inferred`. | Needs matching constants and boundary conditions from 2008. |
| Composite solution `S_comp = S_inner + S_outer - S_overlap` with overlap from a linear fit to the inner far field. | `src/composite.jl:1-118`; tests `test/test_bead7.jl`. | `IMPL-inferred`. | 2001 uses an intermediate asymptotic coordinate and analytic matching constants; current composite is a numerical splice. |
| PDE method of lines for primitive `R,u`, nonuniform grid, FBDF time integration, initial cone, truncated-tip and outflow boundary conditions. | `src/pde.jl:1-12`, `src/pde.jl:87-272`; tests `test/test_bead8.jl`. | Governing idea `IMPL-inferred`; IC/far field partly `C2001`; numerical BCs `IMPL-inferred`. | PDE is an internal verification model, not currently a Decent-King source reproduction. |
| Source-status-aware numerical regression values. | `test/test_numerical_regressions.jl:61-172`. | `IMPL-inferred`. | Tests correctly state these are locally blessed, not published Decent-King benchmarks. |
| Fast/slow test partition. | `test/runtests.jl:4-68`. | Implementation infrastructure. | Relevant for future validation gates, not source content. |

## Quantities Needed For Verification

These are the quantities that must be captured before inner, outer, composite,
and PDE verification can claim source fidelity.

Machine-readable reference records currently live in
`test/reference/decent_king_cone_reference.toml`. That file intentionally keeps
2001 precursor values, 2008 metadata, blocked 2008 quantitative placeholders,
and local regression values separate.

| Verification area | Required quantities | Current status |
| --- | --- | --- |
| Slender/PDE model | Dimensional variables, nondimensional scales, exact slender equations, curvature convention, axial-curvature ordering, pressure/Bernoulli sign, initial cone, far-field condition, tip condition, valid time range. | Initial/far-field/scaling partly `C2001`; primitive equations and axial term `IMPL-inferred`; exact 2008 equation IDs `BLOCKED-2008`. |
| Similarity reduction | Exact similarity variables for position, radius, potential/velocity, tip coordinate, transformed equations, regularity conditions, small-time validity statement. | Length scaling `C2001`; current `S,U` ODEs `IMPL-inferred`; 2008 equations `BLOCKED-2008`. |
| Inner BVP | Inner variables, region definition, leading equations, boundary/tip conditions, far-field form, shooting/free constants, reported constants, wave amplitude/phase behavior. | 2001 simplified model supplies `A,R,R0,A1,p0,y0` structure; current `xi0,S0,S''0,U` structure is `IMPL-inferred`; 2008 full-cone details `BLOCKED-2008`. |
| Outer problem | Outer variables, expansion powers, base state, oscillatory harmonics, coefficient equations, boundary conditions at infinity, amplitude/phase variables, numerical coefficient values. | 2001 supplies OCR-limited equations (4)-(9); current `s1,u1` linearization `IMPL-inferred`; 2008 exact ledger `BLOCKED-2008`. |
| Composite/matching | Intermediate coordinate, overlap limits, matching constants, relationships among inner/outer constants, composite formula, error/order estimates. | 2001 supplies OCR-limited matching conditions; current additive composite `IMPL-inferred`; 2008 exact constants `BLOCKED-2008`. |
| Figures and numerical values | Cone aspect ratios used, tabulated constants, plotted profiles, wave envelope, tip/blob measures, comparison tolerances, figure captions. | 2001 reports `epsilon=0.01` and limiting `y0,p0`; current README/test values are local only; 2008 values `BLOCKED-2008`. |

## Unsupported Current Claims To Correct Later

The following claims are allowed as local implementation/provenance notes, but
should not be described as confirmed Decent-King 2008 benchmarks until the 2008
paper is extracted.

| Claim | Location | Status |
| --- | --- | --- |
| Current inner solution values such as `xi0 ~= 2.77`, `S0 ~= 0.21` or current regression values `xi0 ~= 2.7649`, `S0 ~= 0.2383`, `S''0 ~= 1.1581`. | `README.md`, `HANDOFF.md`, `test/test_numerical_regressions.jl:66-85`. | `IMPL-inferred`; not source benchmark. |
| Current 6-10 capillary-wave oscillation count and envelope claims. | `README.md`, `HANDOFF.md`, generated figures, `test/test_numerical_regressions.jl:103-118`. | Qualitatively consistent with 2001/2008 abstract-level wave statements, but quantitative values are local. |
| Current 3D Newton BVP over `(xi0,S0,S''0)` as "the" Decent-King inner BVP. | `README.md`, `HANDOFF.md`, `src/inner.jl`. | `IMPL-inferred`; needs 2008 source equation and boundary-condition extraction. |
| Current outer linearized equation and additive composite as paper-faithful matching. | `README.md`, `HANDOFF.md`, `src/outer.jl`, `src/composite.jl`. | `IMPL-inferred`; 2001 points to a more detailed oscillatory/multiple-scale matching structure. |
| PDE convergence and similarity collapse as Decent-King validation. | `README.md`, `HANDOFF.md`, `src/pde.jl`, `test/test_numerical_regressions.jl:146-172`. | Internal consistency only until source values and problem conditions are extracted. |

## Pending 2008 Extraction Checklist

When `DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf` is acquired locally, extract
the following into this ledger with exact equation numbers and page references:

- Bibliographic confirmation from the PDF front matter: title, authors, journal,
  volume/issue, pages, DOI, received/revised dates if present.
- Nondimensional variables and scaling: length, time, potential/velocity,
  surface-tension coefficient, density, cone aspect ratio, tip coordinate, and
  any notation conversion needed for this package.
- Small aspect-ratio ordering: powers of `epsilon`, distinguished limits,
  region definitions, and whether axial curvature appears at leading or higher
  order in each region.
- Governing equations: full axisymmetric potential problem, free-surface
  kinematic and dynamic/Bernoulli conditions, curvature convention, simplified
  slender equations, and any primitive-variable reduction comparable to
  `R,u`.
- Initial, far-field, boundary, and tip conditions: initially conical surface,
  far-field stationarity, tip regularity/blunting, any contact with the axis,
  and constants/free parameters.
- Similarity equations: exact source ODEs, source variable names, equation
  numbers, algebraic eliminations, and sign conventions.
- Inner problem: inner variables, leading equations, boundary conditions,
  shooting constants, asymptotic far-field form, reported numerical constants,
  wave amplitude/phase definitions.
- Outer problem: expansion powers, coefficient names, equations, source terms,
  boundary conditions at infinity, matching data imported from the inner
  solution.
- Matching constants: intermediate coordinate, overlap assumptions, constants
  analogous to `p0`, `y0`, `y1`, `phi0`, amplitude/phase constants, and any
  normalization choices.
- Numerical values: all tabulated or text-reported constants, aspect ratios,
  solver parameters, tolerances, domain truncations, and comparison values.
- Figures: every profile/wave/composite figure with axis variables, parameter
  values, plotted quantities, and digitizable reference data if no table exists.
- Relation to 2001 precursor: identify which 2001 equations are reused,
  corrected, generalized, or replaced in the 2008 article.

## Immediate Follow-Up For `UP-SF3`

1. Treat `src/slender.jl`, `src/similarity.jl`, `src/inner.jl`, `src/outer.jl`,
   `src/composite.jl`, `src/pde.jl`, and tests as reconstructed until each
   governing equation has a source-ledger row with a 2008 equation ID or an
   explicit "implementation-only verification" label.
2. Do not convert current local regression values into Decent-King benchmark
   tests. Keep them labelled as source-status-aware local regressions.
3. Reconcile the current primitive `S,U` formulation against the 2001
   potential-based `A,R` precursor only as a heuristic bridge. The canonical
   reconciliation target is the 2008 article body.
