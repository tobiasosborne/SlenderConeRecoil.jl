# Ground-Truth Fidelity Review

Review Agent D, 2026-06-01. Scope: mathematical and ground-truth fidelity of the stated slender-cone recoil model, using local docs/source/tests and local papers where available. No Julia package operations, precompilation, manifest updates, or tests were run.

## Findings

### 1. Critical: the cited 2008 primary reference is wrong locally and in the README

Label: source-backed fact.

The README says the project reproduces Decent & King (2008) in QJMAM 61(1), DOI `10.1093/qjmam/hbm028` (`README.md:3`, `README.md:107`). The local file named `docs/papers/DecentKing2008_QJMAM_61_1.pdf` is not that paper: text extraction identifies it as "Buckling of an axisymmetric vesicle under compression" with DOI `10.1093/qjmam/hbm021`. A web lookup of the article title found the actual paper as Decent and King, "Surface-tension-driven flow in a slender cone", IMA Journal of Applied Mathematics 73(1), 37-68, DOI `10.1093/imamat/hxm043`.

Impact: the project currently lacks the main 2008 primary source it claims to reproduce. Any detailed claim about equation fidelity to Decent & King 2008 is blocked until the correct paper is obtained and checked. This also invalidates the README reference block as ground-truth metadata.

### 2. Critical: the implemented governing model is not demonstrably the Decent-King formulation

Label: source-backed fact plus mathematical concern requiring primary-paper verification.

The docs claim a computational reproduction of Decent & King and a symbolic derivation from the full axisymmetric Euler equations (`README.md:3`, `README.md:21`, `docs/method.md:11-16`). The implemented model is instead a hand-coded primitive 1D radius/axial-velocity system,

```text
d(R^2)/dt + d(R^2 u)/dz = 0
u_t + u u_z = -d(1/R)/dz
```

with optional linearized axial curvature (`src/slender.jl:5-8`, `src/slender.jl:83-100`, `src/pde.jl:67-70`). The local Decent-King 2001 proceedings paper formulates the self-similar slender-cone problem through a velocity potential `A(y)` and radius `R(y)`, with kinematic and Bernoulli equations containing epsilon corrections, fourth derivatives of `A`, and the full curvature expression (article p. 83, equations (1)-(2), `docs/papers/DecentKing2001_IUTAM.pdf`). Its later inner/outer asymptotics are written in those variables, not in the code's `[S, S', S'', U]` ODE.

Impact: the code may be a useful Eggers-style 1D slender-jet surrogate, but the repository does not currently show a derivation path from the cited Decent-King equations to the implemented ODE/PDE. The claim that the CAS "derives the 1D slender model from the full axisymmetric Euler equations" is not supported by `verify_slender_derivation`, which only substitutes `R = epsilon*f` into the leading curvature (`src/slender.jl:117-148`).

### 3. High: the tip boundary conditions do not match the local Decent-King 1D source

Label: source-backed fact plus mathematical concern requiring 2008 verification.

The code treats the tip as a finite-radius point with `S'(xi0) = 0`, `S(xi0) = S0 > 0`, and `U(xi0) = 4 xi0 / 5` (`src/inner.jl:56-68`). It then adds a hemispherical cap in figure/composite presentation to connect that finite radius back to the axis (`scripts/figures.jl:96-151`, `scripts/figures.jl:257-262`). The local Decent-King 2001 paper instead states that the cone tip is at `r = 0`, becomes blunted after motion starts, and uses near-tip asymptotics where the one-dimensional radius satisfies `R ~ sqrt(2 rho0 y / epsilon)` as `y -> 0` in the direct 1D equations and `R0 ~ sqrt(2 rho0 yhat)` in the inner equations (article pp. 83-85, `docs/papers/DecentKing2001_IUTAM.pdf`).

Impact: a finite-radius shooting point plus post-hoc cap is not a derived boundary condition for the free-boundary tip. This can change the blob radius, tip location, volume near the nose, and phase of the outgoing capillary waves.

### 4. High: the outer expansion and matching are not the local Decent-King asymptotics

Label: source-backed fact.

The code linearises about `S = epsilon*xi, U = 0` as `S = epsilon*xi + delta*s`, `U = delta*u`, then also documents a hierarchy ansatz `S = epsilon*xi + epsilon^3*sigma1 + epsilon^5*sigma2`, `U = epsilon^2*omega1 + ...` (`src/outer.jl:7-11`, `src/outer_hierarchy.jl:1-4`, `src/outer_hierarchy.jl:51-59`). The local Decent-King 2001 source uses a different outer structure in `y`, with oscillatory exponentials and fractional powers, for example `R = S0(y) + epsilon^(1/2) gamma(y) exp(i epsilon^(-1) alpha(y)) + epsilon(...) + ...` and matching constants obtained from inner far-field data (article pp. 86-88, `docs/papers/DecentKing2001_IUTAM.pdf`).

The code's default outer driver starts from a zero/small artificial seed at large `xi` (`src/outer.jl:67-97`), while the physically matched path exists separately but is not exported (`src/outer.jl:108-122`) and the composite subtracts only the unperturbed cone as overlap (`src/composite.jl:77-84`). The handoff already flags that the "true matching would use inner far-field as the boundary data" (`HANDOFF.md:36-37`).

Impact: the current outer/composite solution is not a matched asymptotic construction in the sense of the local Decent-King source. It is a plausible numerical continuation/linearization exercise, but its phase, amplitude decay, and matching constants are not grounded in the cited asymptotics.

### 5. High: axial curvature is retained only as the small-slope linear term, while the source stresses full curvature retention

Label: source-backed fact plus mathematical concern requiring primary-paper verification.

The code defines the full axisymmetric curvature as

```text
kappa = 1/(R*sqrt(1+Rz^2)) - Rzz/(1+Rz^2)^(3/2)
```

(`src/slender.jl:29-42`), but the actual PDE/ODE solvers use only `kappa ~= 1/R - Rzz`, giving `+Rzzz` in the PDE and `-S'''` in the similarity ODE (`src/slender.jl:90-100`, `src/inner.jl:3-5`, `src/pde.jl:67-70`, `src/pde.jl:88-95`). The local Decent-King 2001 note says their numerical results remain similar "so long as the full expression for the curvature is retained" (article p. 88 note 2, `docs/papers/DecentKing2001_IUTAM.pdf`).

Impact: the `-S'''` term is the right small-slope axial-curvature correction for a slender model, and it is internally sign-consistent in the code. But it is not the full curvature retained in the cited source. This matters especially near the blunted tip/cap where slopes need not be small.

### 6. Medium: sign and scaling comments are internally inconsistent, although the implemented primitive equations are algebraically self-consistent

Label: source-backed fact.

Within the implemented primitive model, the capillary sign is consistent: `u_t + u u_z = -d(1/R)/dz = Rz/R^2`, and with axial curvature `u_t + u u_z = Rz/R^2 + Rzzz` (`src/slender.jl:5-8`, `src/pde.jl:92-95`). The similarity derivation from `u = ldot*U` also leads to `-(2/9)U + (4/9)(U - xi)U' - S'/S^2 - S''' = 0` (`src/similarity.jl:31-76`, `src/similarity.jl:125-144`).

However, the file header in `src/similarity.jl` gives a different mass equation involving `U - (2/3)xi` and `U' - 2/3` (`src/similarity.jl:9-12`), then later corrects itself through several "wait" comments (`src/similarity.jl:57-68`, `src/similarity.jl:80-104`). `test/test_bead3.jl` also describes the "known" model as `u_t + u u_z = d(1/R)/dz` while the implemented/tested expression is the opposite sign convention (`test/test_bead3.jl:63-70`, `test/test_bead3.jl:85-90`). `test/test_bead4.jl` records a corrected far-field velocity interpretation after first writing an inconsistent one (`test/test_bead4.jl:79-103`).

Impact: I did not find a confirmed sign reversal in the executed formulas, but the comments/docs contain enough contradictory sign and scaling statements to make future maintenance risky.

### 7. Medium: PDE "verification" does not verify convergence to the similarity solution

Label: source-backed fact.

The README says the time-dependent PDE snapshots confirm convergence toward the similarity solution (`README.md:27`), and `docs/method.md` says late-time rescaling should collapse onto it (`docs/method.md:37-39`). The PDE tests only check grid monotonicity, finite differences on simple functions, that the solver runs, positivity of `R`, initial condition preservation, and that rescaling returns finite vectors (`test/test_bead8.jl:47-87`). They do not compare against the inner/composite solution or measure convergence rates. The PDE itself uses a truncated left boundary with fixed `R` and `u = 0` (`src/pde.jl:9-11`, `src/pde.jl:97-99`), which is not obviously a faithful moving recoiling tip boundary.

Impact: the PDE layer cannot currently serve as ground-truth validation of the ODE or matching calculations.

### 8. Medium: the far-field boundary condition is treated as a terminal numerical target, not an asymptotic condition with quantified error

Label: source-backed fact plus mathematical concern.

The inner Newton residual enforces only `S(xi_max)/xi_max -> epsilon`, `U(xi_max) -> 0`, and `S''(xi_max) -> 0` at a finite `xi_max` (`src/inner.jl:90-100`). The default Newton tolerance is loose (`1e-4`), and the README states the far-field slope error remains about 2 percent (`README.md:101`). The test allows the far-field slope to differ by `0.05` at `epsilon=0.1`, a 50 percent relative tolerance (`test/test_bead5.jl:39-44`).

Impact: even within the implemented model, the inner BVP is not tightly pinned to the far-field cone. For asymptotic matching, the far-field phase/amplitude data matter, not just a slope ratio at one truncation point.

## Open Questions

1. Obtain the correct Decent & King 2008 IMA Journal paper, DOI `10.1093/imamat/hxm043`, and replace or quarantine the mislabeled local QJMAM vesicle PDF. Then re-run this fidelity check against the actual equations, figures, constants, and notation in that paper.
2. Decide whether the package is intended to reproduce Decent-King's potential/multiple-scales asymptotics, or to implement a simplified Eggers-style 1D slender-jet model inspired by it. The README should not claim the former unless the derivation is present.
3. Verify the exact pressure/curvature sign convention against the correct primary paper. The code is internally consistent for `p ~ kappa` and acceleration `-d(1/R)/dz`, but the local notes and tests contain contradictory sign descriptions.
4. Determine the mathematically correct tip boundary conditions in `R(z)` form for a smooth recoiling axisymmetric nose. The current finite-radius shooting point plus hemispherical cap needs derivation or replacement.
5. Derive the outer problem from the correct paper, including the scaling of oscillatory waves, matching constants, and far-field decay. The current `epsilon^3` hierarchy and seeded linear solver should be treated as provisional.
6. Add quantitative benchmarks from the primary source: reported constants, blob radius/location, far-field phase/amplitude/wavelength, and convergence of PDE rescalings.

## Commands Run

No Julia commands were run.

```bash
bd onboard
bd prime
pwd && rg --files
wc -l README.md docs/method.md src/slender.jl src/similarity.jl src/inner.jl src/outer.jl src/composite.jl src/pde.jl src/outer_hierarchy.jl test/runtests.jl test/test_bead*.jl
rg -n "Decent|King|Keller|Miksis|curvature|similarity|slender|inner|outer|PDE|far|boundary|match|epsilon|kappa|ODE|recoil|cone|axis|tau|rho|eta|xi|R\\(" README.md docs src test
nl -ba README.md
nl -ba docs/method.md
nl -ba src/slender.jl
nl -ba src/similarity.jl
nl -ba src/inner.jl
nl -ba src/outer.jl
nl -ba src/outer_hierarchy.jl
nl -ba src/composite.jl
nl -ba src/pde.jl
nl -ba src/SlenderConeRecoil.jl
nl -ba scripts/figures.jl
nl -ba test/runtests.jl
nl -ba test/test_bead2.jl
nl -ba test/test_bead3.jl
nl -ba test/test_bead4.jl
nl -ba test/test_bead5.jl
nl -ba test/test_bead6.jl
nl -ba test/test_bead7.jl
nl -ba test/test_bead8.jl
nl -ba SlenderConeRecoil_CLAUDE.md
nl -ba HANDOFF.md
nl -ba reviews/2026-06-01-full-audit/README.md
find . -maxdepth 3 -type f \( -iname '*.pdf' -o -iname '*.bib' -o -iname '*paper*' -o -iname '*ref*' \) -print
pdfinfo docs/papers/DecentKing2008_QJMAM_61_1.pdf | sed -n '1,20p'
pdftotext -q docs/papers/DecentKing2008_QJMAM_61_1.pdf - 2>/dev/null | head -80
pdftotext -q docs/papers/DecentKing2008_QJMAM_61_1.pdf - 2>/dev/null | rg -n "BUCKLING|VESICLE|hbm021|Surface-tension-driven flow in a slender cone|Decent|King|Preston|Jensen|Richardson|doi"
pdftotext -q -f 93 -l 100 docs/papers/DecentKing2001_IUTAM.pdf - 2>/dev/null | nl -ba | sed -n '1,430p'
pdftoppm -q -png -r 180 -f 95 -l 100 docs/papers/DecentKing2001_IUTAM.pdf /tmp/dk2001_page
curl -L --fail --silent --show-error -o /tmp/DecentKing2008_hxm043.pdf https://academic.oup.com/imamat/article-pdf/73/1/37/1920099/hxm043.pdf
git status --short
```

External lookup: web search for `10.1093/qjmam/hbm028 Surface-tension-driven flow in a slender cone Decent King` found the Oxford Academic/IMA Journal metadata for DOI `10.1093/imamat/hxm043`; direct PDF download returned HTTP 403.
