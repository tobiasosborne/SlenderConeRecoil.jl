# Capillary-Recoil and Free-Surface Similarity Literature

Date: 2026-06-01

Assigned bead: `scr-x9r`

Scope: source-backed survey of the capillary-recoil/free-surface similarity
literature most relevant to `SlenderConeRecoil.jl`, with local paper inventory
status and candidate package benchmarks.

## Access And Local Inventory

I used local `docs/papers/README.md`, local PDF text where available, publisher
or DOI metadata pages, Crossref metadata, and open-access publisher pages. I did
not download or modify local papers.

| Source | Stable link | Local status | Access notes |
| --- | --- | --- | --- |
| Keller & Miksis (1983), "Surface Tension Driven Flows" | DOI: <https://doi.org/10.1137/0143018> | Present: `KellerMiksis1983_SIAMJAM_43_268.pdf` | Local PDF text and SIAM metadata agree. |
| Peregrine, Shoker & Symon (1990), "The bifurcation of liquid bridges" | DOI: <https://doi.org/10.1017/S0022112090001835> | Not listed in local manifest | Publisher abstract accessible. |
| Keller, King & Ting (1995), "Blob formation" | DOI: <https://doi.org/10.1063/1.868723> | PDF missing; only `KellerKingTing1995_page.png` present | Research found the previous manifest DOI `10.1063/1.868513` was not the blob-formation article; the manifest was corrected in this batch. |
| Billingham (1999), "Surface-tension-driven flow in fat fluid wedges and cones" | DOI: <https://doi.org/10.1017/S0022112099006047> | Missing; the previous local `Billingham1999_JFM_397_45.pdf` was quarantined as `NOT_Billingham_ChenChen1999_JFM_395_327_S0022112099006011.pdf` | Local `pdfinfo` identifies the quarantined file as Chen & Chen, "Effect of gravity modulation on the stability of convection in a vertical slot". Treat Billingham as missing locally until the PDF is reacquired and checksummed. |
| Decent & King (2001), "The Recoil of A Broken Liquid Bridge" | DOI: <https://doi.org/10.1007/978-94-010-0796-2_10> | Present inside `DecentKing2001_IUTAM.pdf` | Local chapter text inspected. |
| Decent & King (2008), "Surface-tension-driven flow in a slender cone" | DOI: <https://doi.org/10.1093/imamat/hxm043> | Missing as of local manifest | OUP and institutional metadata/abstract accessible; primary target PDF should be acquired before source-fidelity benchmark work. |
| Eggers (1997), "Nonlinear dynamics and breakup of free-surface flows" | DOI: <https://doi.org/10.1103/RevModPhys.69.865> | Present: `Eggers1997_RMP_69_865.pdf` | Local PDF text inspected. |

## Source-Backed Literature Facts

### Keller--Miksis Scaling And Potential-Flow Archetypes

Keller & Miksis (1983) is the scaling root for this package. They consider
time-dependent potential flows with surface tension as the only driving force
and analyze two self-similar configurations: a breaking liquid sheet model and a
contact-line/wedge model. Their dimensional result is the familiar inertio-
capillary similarity scale:

- length scale: `(\sigma t^2/\rho)^(1/3)`;
- velocity scale: `(\sigma/(\rho t))^(1/3)`.

They reduce the free-surface problem to an integrodifferential system for the
surface and boundary potential, solve it numerically, and derive far-field
capillary waves analytically. Their local PDF also gives the wave asymptotics
for a linearized planar surface: wave height decays algebraically while the
wavelength shortens with distance. This is the correct conceptual ancestor for
the Decent--King cone scaling and for any outer capillary-wave benchmark.

Peregrine, Shoker & Symon (1990) provide the experimental motivation for the
broken-liquid-bridge model. Their JFM abstract reports detailed liquid-bridge
shapes before, at, and after drop bifurcation; the post-bifurcation recoil is
attributed to unbalanced surface tension, and the observed geometric similarity
supports a Keller--Miksis type inertio-capillary similarity solution.

### Blob Formation

Keller, King & Ting (1995) analyze the blob that grows on the broken ends of a
liquid thread or sheet as surface tension retracts the end. Crossref metadata
for the correct DOI, `10.1063/1.868723`, states that they construct asymptotic
expansions of the blob shape and flow for thread and sheet geometries, extending
earlier mass and velocity calculations. This is closely related to the "blob at
the recoiling tip" observed in the cone problem, but it is not a replacement for
the Decent--King slender-cone matched-asymptotic construction.

### Fat Wedges And Cones

Billingham (1999) treats "fat" viscous wedges and cones with initial
semi-angles close to `pi/2`. Cambridge metadata states that, shortly after
release from rest, there is an inner region where surface tension and viscosity
dominate and an outer region where inertia and viscosity dominate. The initial
tip velocity is singular, `O(log(1/t))`, as `t -> 0`. At long times, the free
surface approaches a similarity form with `O(t^(2/3))` deformations and
capillary waves propagating away from the tip; viscosity damps the capillary
waves a distance `O(t^(3/4))` from the tip. Billingham solves the linearized
problem by double integral transforms and compares the inviscid fat-wedge
asymptotics against nonlinear inviscid numerics for arbitrary wedge semi-angle.

This paper is the main bridge between Keller--Miksis type recoil and a
viscous/non-slender regime. It should not be conflated with Decent & King
(2008), which is the small-aspect-ratio cone calculation.

### Decent--King Slender Cone Recoil

Decent & King (2001) is the IUTAM precursor. The local chapter derives a
simplified one-dimensional model for an initially sharp slender cone of ideal
fluid driven by surface tension. It uses the Keller--Miksis similarity variables
and assumes small cone aspect ratio. The local text records a one-dimensional
shooting solution with a high-frequency oscillation modulated on a longer scale,
then constructs inner and outer asymptotic regions using multiple scales.

Decent & King (2008) is the target source for this package. OUP metadata states
that it studies an ideal fluid initially conical after droplet/thread break-off,
with surface tension as the only force, and finds a small-aspect-ratio
asymptotic solution. Its similarity transformation is valid for small times
after bifurcation and identifies a rapidly oscillating nonlinear wave that
propagates away from the tip, matching the experimental observations.

The package should treat the 2008 IMA paper as the canonical reference and the
2001 IUTAM chapter as a useful but lower-fidelity precursor.

### Slender-Jet Breakup Context

Ting & Keller (1990), "Slender Jets and Thin Sheets with Surface Tension",
DOI <https://doi.org/10.1137/0150090>, derive slender equations for sheets and
jets with surface tension, applicable to breaking jets, merging jets, holes or
slits in sheets, and related slender free-boundary motions. This is a natural
equation-level predecessor for one-dimensional slender-jet tests.

Eggers (1993), DOI <https://doi.org/10.1103/PhysRevLett.71.3458>, identifies
universal pinching in three-dimensional axisymmetric free-surface flow. Eggers
& Dupont (1994), DOI <https://doi.org/10.1017/S0022112094000480>, derive a
one-dimensional approximation of the Navier--Stokes equations for a thin
axisymmetric viscous column, compare with jet and pendant-drop experiments, and
study the singularities formed at neck pinch-off. Eggers (1997), present
locally, reviews the broader breakup literature and emphasizes one-dimensional
models, self-similar singularities, and continuation through topology change.

Papageorgiou (1995), DOI <https://doi.org/10.1063/1.868540>, uses rational
asymptotics under a slender-jet approximation to derive a one-dimensional model
for viscous Newtonian threads driven by capillarity. Its abstract reports
finite-time radius collapse and validity of the slender approximation through
pinch-off. This is the standard viscous-thread benchmark adjacent to, but not
identical with, inviscid cone recoil.

Day, Hinch & Lister (1998), DOI <https://doi.org/10.1103/PhysRevLett.80.704>,
is the key inviscid pinchoff comparison. Public metadata reports a unique
self-similar pinchoff shape, with lengths scaling like `tau^(2/3)` and two
asymptotic cone angles. Lister & Stone (1998), DOI
<https://doi.org/10.1063/1.869799>, show that even small external viscosity can
change the asymptotic balance for capillary breakup of a viscous thread
surrounded by another viscous fluid. Together these papers warn that "surface
tension plus inertia" benchmarks are only asymptotically stable when the
surrounding medium and viscosity regime are controlled.

## Recent Forward Work, 2020--2026

The recent literature does not appear to replace Decent & King as the core
slender-cone target, but it changes what a best-in-class package should
benchmark if it expands beyond the pure inviscid cone.

- Kamat et al. (2020), "Surfactant-driven escape from endpinching during
  contraction of nearly inviscid filaments", DOI
  <https://doi.org/10.1017/jfm.2020.476>, shows that surfactants can prevent
  endpinching in nearly inviscid contracting filaments through Marangoni stress
  and vorticity generation at curved interfaces. This makes surfactant-free
  assumptions a deliberate model boundary, not a harmless omission, for real
  filament recoil.
- Pierson, Magnaudet, Soares & Popinet (2020), "Revisiting the Taylor--Culick
  approximation: Retraction of an axisymmetric filament", DOI
  <https://doi.org/10.1103/PhysRevFluids.5.073602>, numerically revisits
  capillary retraction of an axisymmetric viscous filament in a passive
  surrounding fluid. It is a useful modern retraction benchmark for finite
  Ohnesorge number and finite filament length.
- Gordillo, Onuki & Tagawa (2020), "Impulsive generation of jets by flow
  focusing", DOI <https://doi.org/10.1017/jfm.2020.270>, gives a potential-flow
  model for high-speed jets produced by sudden implosion of locally spherical
  cavities and predicts jet-tip radius, velocity and droplet radii from initial
  interfacial normal velocity data. It is not cone recoil, but it is a modern
  inviscid free-surface focusing benchmark.
- Basak, Farsoiya & Dasgupta (2021), "Jetting in finite-amplitude, free,
  capillary-gravity waves", DOI <https://doi.org/10.1017/jfm.2020.851>, uses
  theory and computation for jet formation from finite-amplitude axisymmetric
  capillary-gravity waves in a cylindrical pool. It demonstrates the utility of
  Bessel-mode initial data and modal comparisons.
- Sen et al. (2021), "The retraction of jetted slender viscoelastic liquid
  filaments", DOI <https://doi.org/10.1017/jfm.2021.855>, extends capillary
  retraction into viscoelastic inkjet-relevant filaments. This is outside the
  current Newtonian/inviscid package scope but relevant to future API design.
- Sanjay et al. (2022), "Taylor--Culick retractions and the influence of the
  surroundings", DOI <https://doi.org/10.1017/jfm.2022.671>, shows that an
  external viscous medium controls the retraction speed of ruptured films in
  two- and three-phase configurations. Surrounding-fluid effects are therefore a
  benchmark axis for any generalized recoil model.
- Constante-Amores et al. (2022), "Role of surfactant-induced Marangoni stresses
  in retracting liquid sheets", DOI <https://doi.org/10.1017/jfm.2022.768>,
  finds in three-dimensional simulations that insoluble surfactants can delay or
  prevent rim breakup and suppress capillary waves ahead of a retracting rim.
- Kayal, Basak & Dasgupta (2022), "Dimples, jets and self-similarity in
  nonlinear capillary waves", DOI <https://doi.org/10.1017/jfm.2022.854>, is
  directly relevant to Keller--Miksis scaling. It presents a minimal inviscid
  capillary model for dimple and jet formation from a collapsing capillary-wave
  trough. For strong nonlinearity it reports a localized space-time window where
  the jet evolves self-similarly with Keller--Miksis inertio-capillary scales.
- Tian, Yang & Thoroddsen (2023), "Conical focusing: mechanism for singular
  jetting from collapsing drop-impact craters", DOI
  <https://doi.org/10.1017/jfm.2022.1085>, reports extremely fast microjets from
  drop-impact crater collapse and identifies a conically converging flow
  mechanism. This suggests conical focusing benchmarks for inviscid free-surface
  codes, but with different physics from a recoiling liquid cone.
- Ismail & Taraki (2024), "Liquid pinching dynamics in an inertial transitioning
  regime", DOI <https://doi.org/10.1103/PhysRevResearch.6.043131>, shows that
  low-viscosity liquid filaments can thin with a capillary-inertial power law
  while the prefactor varies with internal and external viscosity via neck axial
  velocity. This is a modern warning against treating inertial prefactors as
  universal outside the inviscid-air limit.
- Kayal et al. (2025), "Focussing of concentric free-surface waves", DOI
  <https://doi.org/10.1017/jfm.2024.1089>, decomposes initial free-surface
  cavities into Bessel modes and compares theory with DNS for inward-focusing
  gravito-capillary waves. It identifies limits of potential-flow modelling and
  points toward nonlinear viscous theory with boundary-layer effects.
- Kumar et al. (2026), "Slowing of sheet recoil by surface viscosity", DOI
  <https://doi.org/10.1017/jfm.2026.11131>, shows that surface viscosity can
  slow highly viscous sheet recoil. The paper gives a one-dimensional momentum
  equation with capillary, bulk viscous and surface viscous terms, plus an exact
  Lambert-W expression for transient sheet thickness in its asymptotic setting.

## Inferences And Recommendations

These are package-facing conclusions inferred from the sources above, not direct
claims made by every cited paper.

1. The core package benchmark should remain Decent & King (2008). The immediate
   source-fidelity work should reproduce their nondimensionalization,
   small-cone-angle ordering, inner nonlinear similarity problem, outer
   capillary-wave asymptotics, and matching data before advertising numerical
   results as paper-faithful.
2. Keller--Miksis (1983) should be treated as a scaling and far-field-wave
   benchmark layer. Even if the package does not implement their full boundary
   integral formulation, it should expose tests for `t^(2/3)` length scaling,
   `t^(-1/3)` velocity scaling, and capillary-wave phase/amplitude behaviour in
   the relevant asymptotic limit.
3. The current local research inventory needs completion before detailed
   source-fidelity work: Billingham (1999) should be reacquired, and the
   Keller--King--Ting PDF should be acquired under DOI `10.1063/1.868723`.
4. A best-in-class upgrade should separate model families explicitly:
   inviscid potential-flow cone recoil, slender one-dimensional jet/thread
   equations, viscous wedge/cone recoil, retracting finite filaments/sheets, and
   modern extensions with surfactants or surface rheology. Mixing these under
   one "capillary recoil" API would hide regime assumptions.
5. For numerical algorithms, the literature points to three benchmark classes:
   continuation/collocation or shooting for similarity BVPs; one-dimensional
   slender PDE solvers with full curvature for jets/threads/sheets; and
   high-fidelity boundary-integral or VOF references for inviscid/finite-Oh
   validation.
6. The package should store benchmark metadata as dimensionless regimes, not
   only plots: geometry, viscosity ratio, Ohnesorge number, Bond number,
   surface-rheology state, external-fluid model, similarity exponents, fitted
   prefactors, and source DOI.

## Candidate Package Capabilities And Benchmarks

1. **Decent--King cone benchmark**: solve the slender-cone inner similarity BVP
   and outer capillary-wave problem with paper-specified matching quantities.
   Acceptance should be numeric agreement against the 2008 paper once the PDF is
   acquired.
2. **Keller--Miksis scaling checks**: common utilities for converting between
   dimensional and similarity variables, with tests for `L(t) =
   (\sigma t^2/\rho)^(1/3)` and `U(t) = (\sigma/(\rho t))^(1/3)`.
3. **Capillary-wave asymptotic checks**: far-field wave phase and decay tests
   for Keller--Miksis/Decent--King outer solutions, including grid-resolution
   diagnostics for rapidly shortening wavelengths.
4. **Blob/retraction benchmarks**: Keller--King--Ting blob expansion, finite
   filament retraction, and Taylor--Culick sheet/filament speeds as separate
   examples from cone recoil.
5. **Slender-jet breakup suite**: Eggers--Dupont and Papageorgiou one-
   dimensional models, plus Day--Hinch--Lister inviscid pinchoff scaling and
   cone-angle checks, as regression tests for any general slender free-surface
   solver.
6. **Viscosity and ambient-fluid regime tests**: Billingham fat wedge/cone,
   Lister--Stone two-fluid pinch-off, Sanjay et al. external-viscosity sheet
   retraction, and Ismail--Taraki inertial-prefactor transition tests.
7. **Optional physics modules**: surfactant/Marangoni stresses, surface
   viscosity, and viscoelastic filament retraction should be future modules with
   explicit regime flags, not silent corrections to the inviscid cone model.
8. **Reference-data discipline**: each benchmark should cite a DOI, local
   artifact status, extracted numerical values, commands used to reproduce
   figures, and known access gaps.
