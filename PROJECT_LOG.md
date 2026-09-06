# Project Log

This file records changes to the reusable in-house single-cell RNA-seq
preprocessing pipelines. All project files are maintained under
`D:/Xiaonan/CODEX_projects/Inhouse_scripts/Inhouse_scripts`.

## 2026-09-06

### Added

- Added `plot_adata_stacked_bar()` to the Python visualization section for
  stacked cell-composition percentages calculated from AnnData metadata.
- Documented and preserved both exclusion modes: renormalizing after excluded
  cells are removed, or hiding excluded segments after percentages are based on
  all complete observations.
- Added configurable category and bar ordering, colours, internal percentage
  labels, plot styling, output saving, and optional percentage-table return.
- Saved the supplied `MHCIIhi`/`MHCIIlo` colour mapping as the reusable exported
  `MHCII_GROUP_COLORS` palette.
- Kept visualization separate from filtering by operating on a metadata copy;
  caller-supplied colour mappings are also copied before missing colours are
  filled.

### Verification

- Automated tests were not added or run, following the user's standing
  instruction.
- Reviewed function organization, input/output documentation, and the Git diff
  before committing.

### Commit

- `feat: add AnnData composition stacked bar plot` — Added the composition
  plotter and this log entry.

## 2026-09-05

### Added

- Added `plot_anndata_group_umap()` to the Python visualization section for
  publication-sized categorical UMAP plots from AnnData objects.
- Preserved its single-scatter-call behavior so points can be rasterized as one
  layer while axes, labels, title, and legend remain vector objects.
- Documented all layout, palette, styling, legend, coordinate, rasterization,
  saving, and return parameters.

### Changed

- Extended `plot_anndata_group_umap()` with the optional `top_group` parameter.
  A selected category is validated and drawn last so it remains visually above
  all other cells without changing categorical legend order.
- Extended `plot_anndata_group_umap()` with `legend_ncol` for configurable
  multi-column legends and compact inter-column and marker-to-label spacing.
- Expanded `plot_anndata_group_umap()` with `split_by` and `split_categories`.
  It now supports one unsplit UMAP, subset panels split by another annotation,
  or full-data highlight panels when splitting by `group_col` itself.
- Added exact per-panel sizing, configurable panel gaps and titles, shared or
  panel-specific coordinate limits, far-left-only y-axes, and configurable
  legend column spacing. The function now consistently returns a one-dimensional
  NumPy array of axes, including for an unsplit plot.
- Refined the overall UMAP title position so it is horizontally centred over
  the combined panels and vertically centred within the reserved top margin.

### Verification

- Automated tests were not added or run, following the user's standing
  instruction.
- Reviewed the function organization, documentation, and Git diff before
  committing.

### Commit

- `feat: add categorical AnnData UMAP plotter` — Added the UMAP plotting
  utility and this log entry.
- `feat: support top group in AnnData UMAP plots` — Added foreground plotting
  for a selected annotation group.
- `feat: support multi-column UMAP legends` — Added configurable legend columns
  and compact spacing.
- `feat: add split-panel AnnData UMAP plotting` — Added split and highlight
  modes with exact multi-panel layout controls.
- `fix: center UMAP title in top margin` — Refined the overall title placement.

## 2026-09-04

### Added

- Created `Python/scRNAseq_preprocessing.py` with logically separated data
  loading, 10X file preparation, and visualization sections.
- Added `read_one_gsm()` for reading one GEO 10X matrix and preserving barcode
  and GSM metadata.
- Added `convert_genes_to_features()` for converting legacy two-column 10X gene
  files into modern three-column feature files.
- Added `trim_axs()` for removing unused Matplotlib axes from plot grids.
- Added `ordmag_filter()` for approximate Cell Ranger OrdMag Step 1 cell
  calling from dense or SciPy sparse raw UMI matrices.
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
- `feat: add OrdMag cell-calling filter` — Added `ordmag_filter()` and its
  complete input/output documentation.
