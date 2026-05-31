# Full Repository Audit - 2026-06-01

Beads issue: `scr-2ml`

This directory stores independent review reports for the full repository audit.
The audit goal is to identify concrete risks and improvement opportunities
before the package upgrade phase.

Review reports should:

- Lead with findings ordered by severity.
- Cite local files and line numbers for code or documentation claims.
- Distinguish confirmed defects from risks, smells, and open questions.
- Include commands run and whether any Julia code was executed.
- Avoid mutating source files.
- Avoid concurrent Julia package operations, precompilation, manifest updates,
  or full test jobs against this shared project environment.

Coverage areas:

- Repository tidiness and packaging hygiene.
- Code smells and maintainability.
- Tests, edge cases, and local quality gates.
- Numerical robustness and solver behavior.
- Ground-truth fidelity to the governing equations and cited literature.
- Architecture and API shape for a future best-in-class similarity-methods
  package.
