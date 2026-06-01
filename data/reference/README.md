# Reference Data

This directory is the stable tracked location for small reusable reference
datasets and benchmark tables. Files placed here must follow
`docs/reference_data_schema.md` and must be safe to redistribute.

Keep licensed PDFs and source page images out of this directory. Store their
expected local paths, access status, and checksums in metadata instead. Large
generated data should use a future `Artifacts.toml` plan with content hashes
before package code depends on it.
