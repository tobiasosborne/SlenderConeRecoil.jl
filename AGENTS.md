# Agent Instructions

This project uses **bd** (beads) for issue tracking. Run `bd prime` at the
start of every session, after compaction, or whenever context may have been
lost.

These instructions were adapted from `../BennettVM.jl/AGENTS.md` and the
durable workflow rules in `../Bennett.jl/CLAUDE.md`. `../Bennett.jl` did not
contain an `AGENTS.md` when this file was updated on 2026-06-01.

## Current North Star

Make this repository the place to come to for best-in-class computational
similarity methods: reliable derivations, robust numerical algorithms,
ground-truth fidelity to the literature, reproducible figures, strong tests,
and clear package APIs.

## Durable Master Program

The multi-session program is tracked in Beads. Recover state with `bd ready`,
`bd list --status=open`, and `bd show <id>`.

Umbrella issues created on 2026-06-01:

- `scr-q54` - Persist master review and upgrade programme.
- `scr-2ml` - Run rigorous multi-agent repository review.
- `scr-88v` - Synthesize review findings into remediation plan.
- `scr-6qi` - Execute review remediation serially.
- `scr-x9r` - Research best-in-class similarity methods and literature.
- `scr-ts1` - Create package upgrade plan from research.
- `scr-8l4` - Execute package upgrade plan serially.

Dependency order is intentionally serial: review, synthesis, remediation,
research, upgrade planning, upgrade execution. Do not skip ahead unless Beads
records why the dependency changed.

## Compaction Recovery

When resuming after compaction or a fresh session:

1. Run `bd prime`.
2. Run `bd ready` and inspect the active/blocked issue chain with `bd show`.
3. Read this file, `README.md`, `HANDOFF.md`, and any current
   `reviews/<date>-full-audit/` or `docs/research/<date>-similarity-methods/`
   directory referenced by active Beads issues.
4. Use `bd remember` for durable facts that are not obvious from the code or
   issue descriptions. Do not create `MEMORY.md`.
5. Keep issue notes current when a phase starts, pauses, or changes scope.

## Review And Research Artifacts

Use stable, dated artifact directories:

- Repository review reports: `reviews/YYYY-MM-DD-full-audit/`
- Review synthesis and remediation plan: same review directory.
- Literature and algorithm research: `docs/research/YYYY-MM-DD-similarity-methods/`
- Upgrade plans: `docs/roadmap/` or the research directory if tightly coupled.

Review reports should lead with findings, cite local files and line numbers,
include commands run, and separate confirmed defects from risks or open
questions. Research reports must cite primary sources where possible and
separate source-backed facts from inferences.

## Subagent Orchestration

Use subagents for independent review, synthesis, research, and focused
implementation work. The orchestrating agent owns integration quality:

- Independent review agents should not mutate source files.
- Review agents write reports only in the dated review directory.
- Synthesis agents read every report before creating a plan.
- Implementation proceeds serially through Beads issues.
- Each implementation step is reviewed before closing its Beads issue.

Be careful with Julia races. Do not run concurrent Julia package operations,
manifest updates, precompilation, or test jobs against the same project
environment. This includes the fast, slow, and all test gates. If parallel
exploration is unavoidable, use separate temporary depots and avoid writing to
the project manifest.

## Engineering Rules

- Be skeptical of subagent output, previous work, and your own assumptions.
  Verify claims with code, tests, papers, or reproducible commands.
- Prefer fail-fast behavior over silent fallbacks for mathematical or numerical
  invariants.
- Write or update focused tests for behavior changes. For risky numerical work,
  add regression tests around edge cases and known reference values.
- Run fast focused checks during development and the default fast gate
  (`julia --project test/runtests.jl`) before closing code-changing work.
  When solver, PDE, composite, or numerical regression behavior changes, also
  run `SLENDER_RECOIL_TEST_GROUP=slow julia --project test/runtests.jl`; use
  `SLENDER_RECOIL_TEST_GROUP=all` to run both gates in one Julia process.
- Do not add GitHub Actions or external CI. Quality gates are local.
- Keep generated figures and paper-fetching changes intentional; do not churn
  binary artifacts during unrelated fixes.
- Bundle Beads database/cache changes into the same commit as the source or
  documentation change that caused them. Do not make standalone "bd sync"
  commits.

## Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work atomically
bd close <id>         # Complete work
bd dolt push          # Push beads data to remote
```

## Non-Interactive Shell Commands

**ALWAYS use non-interactive flags** with file operations to avoid hanging on confirmation prompts.

Shell commands like `cp`, `mv`, and `rm` may be aliased to include `-i` (interactive) mode on some systems, causing the agent to hang indefinitely waiting for y/n input.

**Use these forms instead:**
```bash
# Force overwrite without prompting
cp -f source dest           # NOT: cp source dest
mv -f source dest           # NOT: mv source dest
rm -f file                  # NOT: rm file

# For recursive operations
rm -rf directory            # NOT: rm -r directory
cp -rf source dest          # NOT: cp -r source dest
```

**Other commands that may prompt:**
- `scp` - use `-o BatchMode=yes` for non-interactive
- `ssh` - use `-o BatchMode=yes` to fail instead of prompting
- `apt-get` - use `-y` flag
- `brew` - use `HOMEBREW_NO_AUTO_UPDATE=1` env var

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
