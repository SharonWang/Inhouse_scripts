# Project Log

This file records changes to the reusable in-house single-cell RNA-seq
preprocessing pipelines. All project files are maintained under
`D:/Xiaonan/CODEX_projects/Inhouse_scripts/Inhouse_scripts`.

## 2026-09-04

### Added

- Created `Python/scRNAseq_preprocessing.py` with logically separated data
  loading, 10X file preparation, and visualization sections.
- Added `read_one_gsm()` for reading one GEO 10X matrix and preserving barcode
  and GSM metadata.
- Added `convert_genes_to_features()` for converting legacy two-column 10X gene
  files into modern three-column feature files.
- Added `trim_axs()` for removing unused Matplotlib axes from plot grids.
- Added detailed NumPy-style docstrings covering each Python function's input
  parameters, output value, exceptions, side effects, and examples.
- Initialized `R/scRNAseq_preprocessing.R` for future R functions, with Roxygen
  documentation conventions and preprocessing-stage sections. No Python
  functions were translated into R.
- Added the pipeline design and implementation plan.

### Verification

- Automated tests were not added or run, following the user's explicit
  instruction for this project update.
- Reviewed the created files and Git diff before committing.

### Commits

- `7eccf56` — Defined the initial scRNA-seq preprocessing pipeline design.
- `b9f7210` — Added the legacy 10X feature-conversion design.
- `c1584c8` — Added the plotting axes utility design.
- `1e54624` — Added the implementation plan.
- `feat: initialize scRNA-seq preprocessing scripts` — Implementation commit
  containing this log entry.
