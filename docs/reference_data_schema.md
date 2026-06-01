# Reference Data Schema And Artifact Policy

This project keeps source provenance, small reusable reference data, computed
benchmark outputs, and large or licensed artifacts separate. Normal package
loading must not require local papers, generated benchmark output, network
access, or artifact downloads.

## Storage Policy

- `docs/papers/` stores local acquisition notes only. Publisher PDFs, page
  images, browser profiles, and fetch diagnostics remain ignored by git.
- `test/reference/` stores small test fixtures that are safe to redistribute.
  These records may cite licensed papers, but they must not include page images
  or publisher PDF content.
- `data/reference/` is the stable tracked location for small reusable reference
  datasets and benchmark tables that should be consumed outside the test suite.
- Large generated arrays, dense solution dumps, and licensed source artifacts
  must not be committed. If they become required for reproducible reference
  runs, add an `Artifacts.toml` plan that records content hashes and retrieval
  policy before adding package artifact helpers.
- `benchmark/results/` and `benchmarks/results/` are generated output
  locations. They are ignored and are not source-of-truth reference data.

## Source Records

Use this record shape for values extracted from papers, books, source ledgers,
figures, tables, or equation text.

Required fields:

- `id`: stable record identifier.
- `source`: stable source identifier, for example `DecentKing2008-IMAJAM`.
- `status`: source-fidelity status such as `C2001`, `C2008-meta`, or
  `BLOCKED-2008`.
- `fact_type`: `numeric`, `qualitative`, `metadata`, or
  `blocked_quantitative`.
- `quantity`: stable quantity identifier.
- `doi`: canonical source DOI when one exists.
- `page_equation_figure_reference`: page, equation, table, or figure reference;
  use an explicit blocked note when article-body access is unavailable.
- `local_source_artifact_path`: expected local path for the PDF or source
  artifact when one exists.
- `local_source_checksum_algorithm`: currently `sha256`.
- `local_source_checksum`: `sha256:<hex>` for an available local source, or a
  status token such as `pending-acquisition` when unavailable.
- `local_source_checksum_status`: `recorded`, `missing`, `blocked`, or
  `not_applicable`.
- `source_license_access_status`: redistribution and access status, for example
  `licensed-local-ignored` or `institutional-access-blocked`.
- `extraction_method`: manual transcription, text extraction, digitization
  method, or source-ledger handoff.
- `parameters`: table of source parameters or an explicit
  `status = "not_applicable"` / `status = "blocked"` table.
- `solver_settings`: table of solver settings for computed records, or an
  explicit not-applicable reason for pure source extraction.
- `units_scaling`: dimensional or nondimensional interpretation.
- `tolerances`: table with `absolute`, `relative`, or an explicit
  not-applicable/unavailable status plus rationale.
- `tolerance_rationale`: human-readable reason for the tolerance policy.
- `source_status_note`: source-fidelity caveat, especially when a record must
  not be used as a Decent-King 2008 benchmark.
- `data_hash`: `sha256:<hex>` hash of the record payload excluding only the
  `data_hash` field.
- `recorded_date`: ISO date on which the record was added or last verified.

Numeric source records additionally require `value` and `tolerance`.
Qualitative source records require `text_value`. Blocked quantitative records
require `no_numeric_benchmark = true` and must not carry `value` or
`tolerance`.

## Computed Benchmark Records

Computed records must be traceable to both source provenance and the exact
solver configuration that produced them.

Required fields:

- `id`, `quantity`, `status`, `recorded_date`, and `data_hash`.
- `source_doi`, `source_id`, and `source_reference` for the source that defines
  the benchmark target.
- `source_checksum_algorithm`, `source_checksum`, and `source_checksum_status`
  for any local paper or source artifact used during extraction.
- `extraction_method` and `extraction_parameters` for digitization,
  transcription, filtering, or normalization.
- `model_parameters` for physical and nondimensional problem settings.
- `solver_settings` including algorithm, grid/domain, tolerances, stopping
  criteria, package version, package commit, and dirty status when available.
- `comparison_tolerances` separating cheap CI tolerances from slower reference
  tolerances.
- `artifact_path` and `artifact_hash` when the benchmark points to a tracked
  data file or content-addressed external artifact.

## Artifact Policy

Licensed PDFs are local inputs only. They may be named and checksummed in
metadata, but they remain untracked and must never be required by `using
SlenderConeRecoil`, normal tests, or examples.

Small reusable text datasets belong under `data/reference/` and should include
checksums in the metadata record that points to them. Use TOML, CSV, or another
stable text format with deterministic ordering. Large data should stay outside
git until an `Artifacts.toml` entry records content hashes and access policy.

The DOI identifies the source. The local checksum identifies the exact local
copy used for extraction. The data hash identifies the extracted or computed
payload. These three identifiers answer different questions and should not be
collapsed into one field.
