# Project Log

This file records changes to the reusable in-house single-cell RNA-seq
preprocessing pipelines. All project files are maintained under
`D:/Xiaonan/CODEX_projects/Inhouse_scripts/Inhouse_scripts`.

## 2026-09-11

### Added

- Added the repository-level `README.md` as the maintained function reference
  for both language-specific single-cell RNA-seq pipelines.
- Catalogued all exported Python functions with their workflow roles,
  inputs, principal return values, practical examples, dependencies, and
  relevant scientific caveats.
- Documented the reusable `MHCII_GROUP_COLORS` and `WT_CKO_COLORS` palettes.
- Documented the R pipeline's structured workflow sections and implemented
  function inventory so users can distinguish available utilities from planned
  sections.
- Added `score_and_assign_two_signatures()` to the Python signature-scoring and
  annotation section. The function scores two supplied gene programs from
  `adata.raw`, a selected layer, or `adata.X`, then assigns each cell to the
  higher score with optional ambiguity handling.
- Added `read_featurecounts_project()` to the R data-loading section for
  importing featureCounts matrices, aligning study metadata, merging optional
  gene annotation and assignment summaries, and optionally constructing an
  edgeR `DGEList`.
- Added complete Roxygen documentation for every input and returned list
  component, including dependencies, side effects, examples, and analytical
  scope.
- Added `check_featurecounts_project()` to the R quality-control section for
  auditing count, metadata, gene, and optional QC ordering; identifier
  uniqueness; missing, negative, and integer-like counts; and assignment-rate
  summaries.
- Preserved the supplied human-readable audit while adding an invisible,
  structured result with overall and individual check outcomes for reuse in
  automated workflows.
- Added `subset_featurecounts_project()` to the R project-manipulation section
  for selecting samples with metadata expressions while synchronizing counts,
  annotations, QC, and featureCounts summary data.
- Rebuilds an optional edgeR `DGEList` for the selected samples while retaining
  every gene for later design-aware expression filtering.

### Changed

- Established the maintenance convention that `README.md` is updated whenever
  a function, palette, parameter interface, or return value is added or
  changed.
- Added validation for signature genes, expression-source conflicts, scoring
  parameters, output-column collisions, and categorical labels while
  preserving first-signature assignment for tied scores.
- Kept scaling isolated from the original expression matrix and documented the
  potential memory cost of zero-centred scaling on sparse matrices.
- Updated the README's R function index and clarified that the featureCounts
  reader supports bulk RNA-seq or sample-level pseudobulk data rather than
  cell-by-gene single-cell matrices.
- Documented the distinction between structural validity and integer-like
  counts, since fractional featureCounts output can be intentional.
- Added an explicit finite-count check so infinite values cannot be omitted
  from aggregate audit results or incorrectly reported as valid.
- Added incoming and outgoing ordering checks plus complete-summary coverage
  validation to prevent silently misaligned featureCounts subsets.
- Rejected duplicate identifiers before name-based subsetting and preserved
  original featureCounts sample-name provenance in selected-sample order.
- Aligned `read_featurecounts_project()` original sample names to retained
  samples so projects imported with `strict=FALSE` remain subset-compatible.
- Prevented stale edgeR `lib.size` and `norm.factors` metadata fields from
  overriding values recalculated for a newly subsetted `DGEList`.

### Verification

- Automated tests were not added or run, following the user's standing
  instruction.
- Reviewed the README against the exported Python interface, the R pipeline,
  and the Git diff before committing.

### Commit

- `docs: add pipeline function reference README` — Added the maintained Python
  and R pipeline reference and this log entry.
- `feat: add two-signature scoring and assignment` — Added reusable two-program
  scoring, assignment, validation, documentation, and the updated README index.
- `feat: add featureCounts project reader` — Added the reusable R importer,
  alignment and QC safeguards, optional annotation and edgeR output, and
  corresponding README documentation.
- `feat: add featureCounts project checker` — Added reusable structural and
  count-value auditing with printed and programmatic results.
- `feat: add featureCounts project subsetting` — Added metadata-driven sample
  selection, synchronized project components, optional edgeR rebuilding, and
  README documentation.

## 2026-09-07

### Added

- Added `plot_seurat_violins()` to the Python visualization section for
  multi-gene, grouped AnnData expression plots combining violins, individual
  cells, central box plots, and optional two-group significance brackets.
- Added two-sided Mann–Whitney comparisons between the first two ordered groups
  and optional Benjamini–Hochberg FDR correction across plotted genes.
- Documented expression-source selection, gene and group ordering, panel layout,
  plot layers, statistical output, saving, and optional expression-data return.
- Documented that cell-level tests are exploratory and do not replace
  replicate-aware pseudobulk or mixed-model differential-expression analysis.
- Saved the supplied `WT`/`cKO` colours as the reusable exported
  `WT_CKO_COLORS` palette.

### Changed

- Updated `plot_seurat_violins()` to accept `mannwhitney`, `mann-whitney`, and
  `mw` as aliases for the two-sided Mann–Whitney U test.
- Added explicit empty-group-order validation, aligned violin hue order with the
  requested group order, and avoided empty-series summary warnings.
- Added the previously missing `legend_fontsize` and `legend_title_fontsize`
  parameters with 12-point defaults, preventing undefined-name errors when a
  shared legend is displayed.
- Removed the automatic `plt.show()` call so callers control interactive display
  and can compose or save the returned figure before showing it.
- Added `gene_fontstyle="italic"` for configurable gene-title styling, rejected
  simultaneous `use_raw=True` and `layer=...`, and hid panels safely when a gene
  has no non-missing expression values in the selected cells.
- Updated the reusable `WT_CKO_COLORS` palette to the newly supplied values:
  `WT="#898C8B"` and `cKO="#9CB4CC"`.

### Verification

- Automated tests were not added or run, following the user's standing
  instruction.
- Reviewed function organization, statistical assumptions, input/output
  documentation, and the Git diff before committing.

### Commit

- `feat: add Seurat-style AnnData violin plots` — Added the violin plotter,
  exploratory statistics, saved palette, and this log entry.
- `fix: refine Seurat-style violin plotting` — Added test aliases, group-order
  safeguards, working legend font controls, and caller-controlled display.
- `feat: refine violin titles and expression source` — Added italic gene titles,
  expression-source validation, empty-panel handling, and the revised WT/cKO
  palette.

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
