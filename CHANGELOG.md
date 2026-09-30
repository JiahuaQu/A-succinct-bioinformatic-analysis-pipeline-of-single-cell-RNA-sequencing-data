# Changelog

## v1.0.0 — 2026-09-30

- Reorganized the repository for public reuse and readability.
- Replaced user- and institution-specific filesystem paths with placeholders/configuration variables.
- Replaced hard-coded LSF project/queue values with public placeholders.
- Consolidated sample examples into `config/samples.example.tsv`.
- Parameterized the Cell Ranger LSF job script without changing its core `cellranger count` purpose.
- Renamed the downstream R script for clarity while retaining its analysis logic.
- Expanded README documentation, requirements, workflow scope, privacy notes, and citation guidance.
- Added `.gitignore`, `CITATION.cff`, `VERSION`, and a conservative `COPYRIGHT.md` rights notice.
