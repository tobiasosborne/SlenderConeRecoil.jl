# HANDOFF — SlenderConeRecoil.jl

## What was done

Working implementation aimed at becoming a faithful computational reproduction of Decent and King (2008): slender cone recoil under surface tension using similarity solutions and matched asymptotics. Source fidelity is tracked in `docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md`. Because the canonical 2008 article body is still unavailable, the current primitive-variable equations, BVP constants, and numerical values are reconstructed/local implementation data unless the ledger marks a fact as 2001-backed or 2008-metadata-backed.

### Implemented pipeline from the original bead plan:
1. **Symbolic CAS** (beads 1-2): Expression types mirroring TensorGR.jl, series expansion in ε
2. **1D slender model** (bead 3): Reconstructed primitive mass + momentum with κ = 1/R - Rzz (azimuthal + axial curvature)
3. **Similarity reduction** (bead 4): Source-backed Keller-Miksis t^{2/3} length scaling, plus reconstructed primitive ODE system in S(ξ), U(ξ)
4. **Inner BVP** (bead 5): Local 4-component ODE [S, S', S'', U], 3D damped Newton over (ξ₀, S₀, S''₀), Rodas5P stiff solver. Produces the current implementation's blob + neck + capillary waves.
5. **Outer linearisation** (bead 6): Local reconstructed 4-component [s₁, s₁', s₁'', u₁] with S''' dispersive term
6. **Composite** (bead 7): Local diagnostic additive inner + outer - fitted overlap
7. **PDE verification** (bead 8): Method of lines, FBDF implicit, non-uniform FD, Rzzz term
8. **Figures** (bead 9): 7 figures including matched asymptotic profiles and tip shape
9. **Documentation** (bead 10): docs/method.md

The current Beads tracker now drives review remediation and package upgrades;
the old numbered bead plan is provenance rather than the active task list.

### Bug fix history:
- **15 bugs** found in initial audit (root cause: momentum sign error ∂/∂z(1/R) → -∂/∂z(1/R))
- **7 more** from reviewer round 1 (FD stencils, allocations, canonical ordering)
- **5 more** from reviewer round 2 (left BC, documentation)
- **ξ₀=0 bug**: tip position must be nonzero for physical recoil (2D→3D Newton)
- **Axial curvature**: added -Rzzz dispersive term for capillary waves

### Current implementation outputs:
- Inner solution local regression: ξ₀ ≈ 2.76, S₀ ≈ 0.24, 6-10 capillary wave oscillations
- PDE convergence and similarity collapse are internal consistency checks, not Decent--King benchmarks
- Matched-asymptotic-inspired composite shows the current reconstructed inner/outer structure, pending source-confirmed matching formulae

## What's NOT done / known limitations

1. **3D Newton convergence**: Only reaches ~2% far-field slope error (tolerance 1e-4). A continuation/homotopy approach would improve this.
2. **PDE long-time integration**: Stiffness from capillary pressure (1/(εz)² near tip) limits PDE to t ≈ 0.04. Adaptive mesh refinement or moving mesh needed for longer runs.
3. **PDE blob comparison**: PDE doesn't run long enough to develop full blob structure for direct comparison with similarity solution in the inner region.
4. **Papers/source fidelity**: The correct cone-paper metadata is Decent and King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* 73(1), 37--68, DOI `10.1093/imamat/hxm043`. The source ledger is `docs/research/2026-06-01-similarity-methods/06_decent_king_source_ledger.md`; use its `C2001`, `C2008-meta`, `IMPL-inferred`, and `BLOCKED-2008` labels as authoritative. The previous local file named `DecentKing2008_QJMAM_61_1.pdf` was actually a QJMAM vesicle-compression paper and has been renamed locally to `docs/papers/NOT_DecentKing_PrestonJensenRichardson2008_QJMAM_61_1_hbm021.pdf`. The previous local file named `Billingham1999_JFM_397_45.pdf` was actually a Chen--Chen gravity-modulation paper and has been renamed locally to `docs/papers/NOT_Billingham_ChenChen1999_JFM_395_327_S0022112099006011.pdf`; the correct Billingham 1999 PDF has since been acquired. Licensed PDFs remain ignored by git; see `docs/papers/README.md`. Missing: Decent--King 2008 and Keller--King--Ting 1995 (AIP, DOI `10.1063/1.868723`).
5. **Test gates**: `test/runtests.jl` now defaults to the fast package/API + symbolic/CAS gate. Set `SLENDER_RECOIL_TEST_GROUP=slow` for solver, PDE, composite, and numerical regression coverage, or `all` to run both groups in one Julia process.
6. **Outer solution**: Currently uses zero/small seed at large ξ. True matching would use inner far-field as the boundary data for the outer problem.
7. **Dependency policy**: Core, script, optional-extension, benchmark, and paper-tool dependency boundaries are documented in `docs/dependency_policy.md`.

## How to run

```bash
julia --project test/runtests.jl       # default fast gate
SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl
SLENDER_RECOIL_TEST_GROUP=all julia --project test/runtests.jl
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts scripts/figures.jl    # all figure PDFs/PNGs + figures/metadata.toml
julia --project=scripts scripts/figures.jl --metadata-only  # metadata without plot binaries
julia --project test/test_bead5.jl    # inner solver tests (slowest, ~50s)
julia --project test/test_bead8.jl    # PDE tests
julia --project test/test_bead1.jl    # symbolic tests
```

Do not run multiple Julia package, precompile, or test jobs concurrently in
this project environment. Use the `all` group when both gates are needed.
Figure generation saves both PDF and PNG for every tracked figure stem and
records parameters, package version, git commit, Julia version, and solver
diagnostics in `figures/metadata.toml`. The figure binaries and metadata file
are tracked reproducibility artifacts.

## Key files
| File | LOC | Purpose |
|------|-----|---------|
| src/inner.jl | 130 | Reconstructed inner BVP: 3D Newton, 4-component ODE, Rodas5P |
| src/pde.jl | 130 | PDE: FBDF, non-uniform FD, Rzzz |
| src/similarity.jl | 120 | Keller-Miksis scaling plus reconstructed primitive similarity ODEs with S''' |
| src/outer.jl | 90 | Outer linearisation, 4-component |
| src/expr.jl | 300 | Symbolic expression types |
| src/series.jl | 180 | ε-expansion |
