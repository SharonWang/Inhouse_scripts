# In-house single-cell RNA-seq utilities

Reusable preprocessing, quality-control, and visualization helpers for
single-cell RNA-seq studies. Python and R code are maintained in separate
language-specific pipelines: a function is added only to the language in which
it was supplied.

This README is the function index for the repository. It must be updated
whenever a function, reusable palette, parameter interface, or return value is
added or changed.

## Repository layout

```text
Inhouse_scripts/
|-- Python/
|   `-- scRNAseq_preprocessing.py
|-- R/
|   `-- scRNAseq_preprocessing.R
|-- PROJECT_LOG.md
`-- README.md
```

## Python pipeline

Source: [`Python/scRNAseq_preprocessing.py`](Python/scRNAseq_preprocessing.py)

### Requirements

The module uses Python 3.10 or newer syntax. Its core dependencies are NumPy,
pandas, SciPy, and Matplotlib. Individual functions load additional packages
only when needed:

- Scanpy and AnnData for 10X import and AnnData expression extraction.
- Seaborn for Seurat-style violin plots.
- statsmodels when Benjamini-Hochberg adjustment is requested for violin-plot
  comparisons.

Install packages in the analysis environment appropriate for the study. This
repository does not pin an environment or alter raw input data automatically.

### Importing the utilities

When Python is started from the repository root:

```python
from Python.scRNAseq_preprocessing import (
    MHCII_GROUP_COLORS,
    WT_CKO_COLORS,
    cluster_expression_summary,
    convert_genes_to_features,
    ordmag_filter,
    plot_adata_stacked_bar,
    plot_anndata_group_umap,
    plot_seurat_violins,
    read_one_gsm,
    score_and_assign_two_signatures,
    trim_axs,
)
```

### Function summary

| Workflow stage | Function | What it does | Main return value(s) |
|---|---|---|---|
| 10X preparation | `convert_genes_to_features()` | Converts legacy compressed 10X gene tables into modern feature tables. | List of output `Path` objects written during the call. |
| Data loading | `read_one_gsm()` | Reads one GEO 10X matrix and records the original barcode and GSM accession. | `AnnData` object. |
| Cell calling | `ordmag_filter()` | Applies an approximate Cell Ranger OrdMag Step 1 UMI threshold. | Cell mask, UMI threshold, and per-barcode UMI totals. |
| Signature scoring | `score_and_assign_two_signatures()` | Scores two gene programs and assigns each cell to the higher-scoring signature. | The annotated input `AnnData` object. |
| Expression summary | `cluster_expression_summary()` | Calculates mean expression and percentage detected for requested genes within each observed cell or spot group. | Long-form pandas `DataFrame`. |
| Figure layout | `trim_axs()` | Removes unused Matplotlib axes from a subplot grid. | Flattened array of retained axes. |
| UMAP visualization | `plot_anndata_group_umap()` | Plots one or more categorical UMAP panels, including split and highlight modes. | Figure and axes array. |
| Composition visualization | `plot_adata_stacked_bar()` | Calculates and plots cell-composition percentages from `adata.obs`. | Figure and axes; optionally the percentage table. |
| Expression visualization | `plot_seurat_violins()` | Draws grouped gene-expression violins with optional points, boxes, and exploratory tests. | Statistics table, figure, and axes; optionally extracted expression data. |

### Data loading and 10X preparation

#### `convert_genes_to_features(data_dir, overwrite=False)`

Converts every `*_genes.tsv.gz` file in a directory from the legacy two-column
10X format into a corresponding `*_features.tsv.gz` file. The function keeps
the gene identifier and symbol, then adds `Gene Expression` as the third
column. Existing feature files are preserved unless `overwrite=True`.

```python
outputs = convert_genes_to_features(
    data_dir="data/10x",
    overwrite=False,
)
```

The returned list contains only files written during the current call; skipped
existing files are not included.

#### `read_one_gsm(data_dir, gsm)`

Finds exactly one `{gsm}_*_matrix.mtx.gz` file, derives its shared 10X filename
prefix, and reads the matrix with `scanpy.read_10x_mtx()`. The matching
`barcodes.tsv.gz` and `features.tsv.gz` files must use the same prefix.

```python
adata = read_one_gsm(
    data_dir="data/geo",
    gsm="GSM7732265",
)
```

The returned object uses unique gene symbols as variable names and includes:

- `adata.obs["barcode"]`: original observation/barcode names.
- `adata.obs["GSM"]`: normalized GSM accession supplied by the caller.

### Quality control and cell calling

#### `ordmag_filter(adata, expect_cells=8000)`

Calculates total UMI counts per barcode, finds the 99th percentile among the
top `expect_cells` barcodes, and calls barcodes at or above one tenth of that
value.

```python
keep, threshold, total_umi = ordmag_filter(
    adata,
    expect_cells=8000,
)

called_adata = adata[keep].copy()
```

`adata.X` must contain raw UMI counts with barcodes in rows and genes in
columns. Sparse and dense matrices are supported. This function returns a mask
and does not subset or otherwise modify the supplied AnnData object.

Important: this is an approximate preliminary cell-calling heuristic. It does
not implement EmptyDrops-style testing, doublet detection, or complete
dataset-specific QC. Inspect count distributions and select QC thresholds for
each study rather than treating the default as universally appropriate.

### Signature scoring and annotation

#### `score_and_assign_two_signatures(...)`

Calculates two Scanpy gene-set scores from `adata.raw`, a named layer, or
`adata.X`, then assigns each cell to the signature with the higher score. Ties
are assigned to the first signature. Optional ambiguity handling assigns cells
to a third label when both scores fall below a specified threshold.

```python
adata = score_and_assign_two_signatures(
    adata,
    gene_list_1=["Cd74", "H2-Ab1"],
    gene_list_2=["S100a8", "S100a9"],
    label_1="MHCIIhi",
    label_2="Inflammatory",
    ambiguous=True,
    ambiguous_threshold=0,
    use_raw=True,
)
```

The function modifies and returns the original AnnData object. It writes the
two score columns and the final assignment column to `adata.obs`; their names
are configurable. Duplicate genes are removed, missing genes are reported, and
each signature must contain at least one gene present in the selected source.

Use an expression representation appropriate for gene-set scoring, commonly
normalized and log-transformed values. Optional scaling is performed on a
temporary copy so the original matrix is protected. However, centred scaling
can densify sparse data and require substantial memory; consider
`zero_center=False` or `scale=False` for large datasets when scientifically
appropriate.

### Cluster-level expression summaries

#### `cluster_expression_summary(adata, genes, groupby, layer=None)`

Returns mean expression, percentage expressing, and cell or spot count for
every available requested gene within each observed group. This is especially
useful for sparse targeted spatial-transcriptomics data, where cluster-level
summaries are generally more stable than gating individual cells on a single
marker.

```python
summary = cluster_expression_summary(
    adata,
    genes=["EPCAM", "KRT8", "KRT18"],
    groupby="leiden",
    layer="log1p",
)
```

The returned long-form table contains the grouping column, `gene`,
`mean_expression`, `pct_expressing`, and `n_cells`. Gene order follows the
request, duplicate requests are removed, unavailable genes are ignored, and
the function raises an error if none are found. Group labels are returned in
alphabetical order. Sparse matrices remain sparse during subsetting, and the
input AnnData object is not modified.

`pct_expressing` is defined as the percentage of selected values greater than
zero. It is therefore most directly interpretable on non-negative count or
normalized-expression data; on centred or scaled layers it instead represents
the percentage above that layer's zero point.

### Visualization

#### `trim_axs(axs, N)`

Flattens a Matplotlib axes collection, removes every axis after the first `N`
from its figure, and returns the retained one-dimensional axes array.

```python
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2)
axs = trim_axs(axs, 3)
```

#### `plot_anndata_group_umap(...)`

Plots categorical annotations stored in `adata.obs` using coordinates from
`adata.obsm["X_umap"]` by default. It supports three modes:

- no `split_by`: one UMAP coloured by `group_col`;
- `split_by` different from `group_col`: one subset panel per split category;
- `split_by == group_col`: all cells remain visible in every panel while one
  category is highlighted per panel.

```python
fig, axes = plot_anndata_group_umap(
    adata,
    group_col="cell_type",
    split_by="condition",
    palette=cell_type_colors,
    title="Cell types by condition",
    save="outputs/cell_type_umap.pdf",
)
```

The function supports explicit split ordering, a group drawn on top, shared or
panel-specific limits, exact panel dimensions, legend layout, rasterized cell
points, high-resolution saving, and transparent output. It does not modify the
AnnData object or a supplied palette.

#### `plot_adata_stacked_bar(...)`

Builds a percentage table from two `adata.obs` columns and plots stacked cell
composition bars. Missing metadata rows are omitted. Excluded categories can
be removed before normalization (`renormalize=True`) or hidden after
percentages are calculated (`renormalize=False`).

```python
fig, ax, percentages = plot_adata_stacked_bar(
    adata,
    x_col="condition",
    stack_col="MHCII_group",
    class_order=["MHCIIhi", "MHCIIlo"],
    colors=MHCII_GROUP_COLORS,
    exclude=["Ambiguous"],
    renormalize=True,
    save="outputs/cell_composition.pdf",
)
```

Set `return_table=False` to return only `(fig, ax)`. Calculations use a copy of
the relevant metadata and do not alter `adata.obs`.

#### `plot_seurat_violins(...)`

Creates one grouped gene-expression panel per requested gene. Each panel can
combine a violin, individual-cell points, a central box plot, and a comparison
between the first two ordered groups. Expression may be read from `adata.X`, a
named layer, or `adata.raw`; `layer` and `use_raw=True` cannot be used together.

```python
stats, fig, axes = plot_seurat_violins(
    adata,
    genes=["Cd74", "H2-Ab1"],
    groupby="genotype",
    order=["WT", "cKO"],
    palette=WT_CKO_COLORS,
    layer="log1p",
    show_stats=True,
    adjust_p=True,
    save="outputs/marker_violins.pdf",
)
```

With statistics enabled, `mannwhitney`, `mann-whitney`, and `mw` all select a
two-sided Mann-Whitney U test. Benjamini-Hochberg FDR adjustment can be applied
across plotted genes. Set `return_data=True` to append the extracted long-form
expression table to the returned tuple.

Important: cell-level tests treat cells as independent observations and are
intended for exploratory visualization. For biological inference with
replicated samples, use a replicate-aware method such as pseudobulk analysis or
an appropriate mixed model.

### Reusable palettes

| Palette | Categories and colours | Intended use |
|---|---|---|
| `MHCII_GROUP_COLORS` | `MHCIIhi: #E7B2B6`, `MHCIIlo: #B7C3E0` | High- and low-MHCII annotation groups. |
| `WT_CKO_COLORS` | `WT: #898C8B`, `cKO: #9CB4CC` | Wild-type and conditional-knockout groups. |

Pass these dictionaries to plotting functions through their `palette` or
`colors` parameter. Supplied reusable palettes are saved in the Python module
and added to this table for future studies.

## R pipeline

Source: [`R/scRNAseq_preprocessing.R`](R/scRNAseq_preprocessing.R)

The R pipeline contains reusable bulk/pseudobulk import, validation,
subsetting, differential-expression, and exploratory visualization utilities.
The file also retains structured sections for future cell-level preprocessing
functions supplied in R.

### Function summary

| Workflow stage | Function | What it does | Main return value(s) |
|---|---|---|---|
| Bulk/pseudobulk count import | `read_featurecounts_project()` | Imports featureCounts counts, aligns metadata, merges optional gene annotation and assignment QC, and optionally creates an edgeR object. | A `featurecounts_project` list containing counts, metadata, genes, QC, summary data, and optional `DGEList`. |
| Bulk/pseudobulk structural QC | `check_featurecounts_project()` | Audits ordering, identifier uniqueness, invalid count values, integer-likeness, and assignment-rate summaries. | An invisible structured audit containing overall flags, individual checks, and assignment summary. |
| Bulk/pseudobulk project selection | `subset_featurecounts_project()` | Selects samples with a metadata expression while synchronizing counts, annotations, QC, and summary data. | A new subsetted `featurecounts_project` with an optional rebuilt `DGEList`. |
| Bulk/pseudobulk expression preparation | `prepare_bulk_expression()` | Selects samples once, applies a documented edgeR library-normalization method, and calculates reusable log2 CPM values. | A `bulk_expression_prepared` list containing aligned counts, metadata, genes, DGEList, log2 CPM, and normalization settings. |
| Pairwise differential expression | `edgeR_pairwise()` | Runs comparison-specific filtering, TMM normalization, dispersion estimation, and an edgeR quasi-likelihood test. | Annotated results plus DGEList, fit, test, design, contrast, and sample details. |
| Pairwise differential expression | `limma_voom_pairwise()` | Runs comparison-specific filtering, TMM normalization, voom weighting, and a limma empirical Bayes test. | Annotated results plus DGEList, voom data, fit, design, contrast, and sample details. |
| Differential-expression collection | `collect_edgeR_pairwise()` | Combines named edgeR pairwise results with analysis, sorting, contrast, provenance, and signed-significance fields. | One row-bound data frame retaining all result and annotation columns. |
| Bulk/pseudobulk visualization | `plot_bulk_pca()` | Filters and normalizes sample-level counts, selects variable genes, calculates PCA, and builds a configurable publication-style sample plot. | Plot, PCA fit and scores, log2 CPM, selected genes, variance summaries, DGEList, and palettes. |
| Bulk/pseudobulk quality-control visualization | `plot_bulk_qc()` | Calculates count-matrix library sizes and plots selected sample-level read, assignment, and library metrics in faceted panels. | Combined patchwork figure, individual metric plots, plotting metadata, and colour mapping. |
| Bulk/pseudobulk expression visualization | `plot_bulk_violin()` | Plots selected genes from a reusable prepared normalized-expression object, with optional boxplots, sample points, and split facets. | Plot, long-format expression data, per-panel summaries, gene mapping, log2 CPM, DGEList, normalization settings, and colours. |
| Gene-program expression heatmap | `plot_gene_set_heatmap()` | Resolves named gene sets from a reusable prepared normalized-expression object, optionally aggregates samples, and displays genes in ordered program slices. | ComplexHeatmap object, displayed and pre-scaled matrices, log2 CPM, gene mapping, missing genes, aligned metadata, row split, colour function, DGEList, and normalization settings. |
| Differential-expression heatmap | `plot_de_heatmap()` | Selects directional DE genes and displays normalized sample-level expression with optional visualization-only batch correction and marked-gene labels. | ComplexHeatmap object, displayed matrix, corrected expression, selected DE rows, ordered metadata, row split, labels, and colour function. |
| Differential-expression visualization | `plot_signed_manhattan()` | Displays signed adjusted-p-value significance for every tested gene across comparisons and sorting groups. | Plot, complete plotting data, selected labels, per-panel summary, colours, and cutoff metadata. |

### `read_featurecounts_project(...)`

Reads a featureCounts gene-level count table and validates the relationship
between its sample columns and a supplied metadata table. Metadata can be a
data frame or a CSV, TSV/TXT, or Excel file. The function can also import a
featureCounts summary, merge an external gene-feature annotation, remove
terminal Ensembl version suffixes, and create an `edgeR::DGEList` when edgeR is
installed.

```r
source("R/scRNAseq_preprocessing.R")

project <- read_featurecounts_project(
  fcounts_file = "counts/featureCounts.txt",
  summary_file = "counts/featureCounts.txt.summary",
  gene_feature_file = "reference/gene_features.txt",
  metadata = "metadata/samples.csv",
  sample_col = "SampleName",
  strict = TRUE,
  make_dge = TRUE
)

project$counts
project$metadata
project$qc
project$dge
```

The returned count matrix, metadata, gene annotation, and QC table are ordered
and checked explicitly. With `strict=TRUE`, the featureCounts and metadata
sample sets must match exactly; `strict=FALSE` retains their intersection in
count-matrix order. Excel metadata requires `readxl`, while DGEList creation
requires edgeR. Original featureCounts sample names are retained in the same
order as the final cleaned sample names.

This function is intended for bulk RNA-seq or sample-level pseudobulk counts,
not a cell-by-gene single-cell matrix. It performs import and structural QC but
does not choose expression filters, normalize libraries, construct a design,
or run differential-expression tests. Those steps must account for biological
replication and the study design.

### `check_featurecounts_project(dat, tolerance=1e-8, verbose=TRUE)`

Checks the internal consistency of an imported featureCounts project without
modifying it. The printed audit retains the original interactive checks, while
the invisible return value makes the result reusable in scripts.

```r
audit <- check_featurecounts_project(project)

audit$valid
audit$all_checks_pass
audit$checks
audit$assignment_summary
```

The function checks that samples and genes have consistent ordering, sample
and gene identifiers are unique, and counts contain no missing or negative
values. It also rejects non-finite values such as positive or negative
infinity. Integer-likeness uses a configurable numerical tolerance. Because
fractional featureCounts output can be intentional, `audit$valid` covers
structural and value integrity without requiring integer-like counts;
`audit$all_checks_pass` includes the integer-likeness result. When
`Assignment_percent` is present in the metadata, its summary is printed and
returned.

### `subset_featurecounts_project(...)`

Creates a new featureCounts project from samples selected by an unquoted
metadata expression. The input project is not modified, and every gene is
retained so comparison-specific expression filtering can be applied later.

```r
lung_project <- subset_featurecounts_project(
  project,
  Tissue == "Lung" & Genotype %in% c("WT", "EPX"),
  make_dge = TRUE,
  drop_levels = TRUE
)

lung_project$sample_names
lung_project$counts
lung_project$dge
```

The function validates the incoming and outgoing count, metadata, gene, and QC
ordering. It also subsets the featureCounts summary when available and creates
a fresh edgeR `DGEList` so library information reflects only the selected
samples. Duplicate sample or gene identifiers are rejected before name-based
indexing, and original featureCounts sample names are retained in selected
sample order. Missing values in the subset expression are treated as `FALSE`,
and the function stops if no samples remain. Existing `lib.size` and
`norm.factors` metadata fields are excluded from DGEList construction so edgeR
recalculates them from the subset counts; the returned metadata itself is left
unchanged.

### `prepare_bulk_expression(...)`

Creates one reusable normalized-expression object for downstream descriptive
bulk or replicate-aware pseudobulk plots. Sample selection happens before
library normalization, and counts, metadata, annotations, the edgeR object,
and log2 CPM values remain in explicit matching order.

```r
bulk <- prepare_bulk_expression(
  project,
  subset = Tissue == "Lung",
  prior_count = 2,
  norm_method = "TMM"
)

bulk$logCPM
bulk$normalization
```

`norm_method` accepts `"TMM"`, `"TMMwsp"`, `"RLE"`, `"upperquartile"`, or
`"none"`, using edgeR's library-size normalization. The return records the
selected method and log-CPM prior count so multiple plots can reuse exactly the
same sample set and transformation. Raw selected counts are retained; the
function does not filter genes, construct a statistical design, or perform
differential-expression testing.

### Pairwise differential expression

Both pairwise functions use the same validated sample-selection and design
helper. `group1` is always the reference and the fitted contrast is always
`group2 - group1`; therefore, a positive `logFC` means higher expression in
`group2`. Optional covariates are included in the model, incomplete model rows
are removed, and rank-deficient designs stop with an error.

#### `edgeR_pairwise(...)`

```r
edger_result <- edgeR_pairwise(
  project,
  group_col = "Genotype",
  group1 = "WT",
  group2 = "cKO",
  subset = Tissue == "Lung",
  covariates = "Batch",
  robust = TRUE,
  fdr_cutoff = 0.05,
  logfc_cutoff = 1
)

head(edger_result$results)
```

This workflow uses `edgeR::filterByExpr()` after selecting the comparison,
followed by TMM normalization, robust dispersion estimation, quasi-likelihood
fitting, and a quasi-likelihood F-test.

#### `limma_voom_pairwise(...)`

```r
voom_result <- limma_voom_pairwise(
  project,
  group_col = "Genotype",
  group1 = "WT",
  group2 = "cKO",
  subset = Tissue == "Lung",
  covariates = "Batch",
  robust = TRUE,
  voom_plot = TRUE
)

head(voom_result$results)
```

This workflow performs the same comparison-specific gene filtering and TMM
normalization, then applies limma-voom precision weights, the requested
contrast, and empirical Bayes moderation. The results retain `adj.P.Val` and
also expose it as `FDR` for consistency with edgeR.

Both functions return every tested gene with annotation, effect size,
significance statistics, and a `Direction` label. Direction thresholds label
results but do not remove rows. Fitted statistics take precedence when an
annotation field has the same name; conflicting annotations are retained with
an `annotation_` prefix. The returned DGEList stores the selected comparison
factor explicitly. These workflows require raw counts and biological
replication. Pairing, donor effects, batch variables, and other covariates must
be specified according to the study design; one group with fewer than two
samples triggers a warning because inference is unreliable.

Dependencies are edgeR for `edgeR_pairwise()` and edgeR plus limma for
`limma_voom_pairwise()`. Robust limma empirical Bayes estimation may also use
statmod through limma.

#### `collect_edgeR_pairwise(de_list, name_sep="__")`

Combines a named list of `edgeR_pairwise()` outputs into one analysis-ready
table. The collector records the list name, reference and comparison groups,
standardized comparison, sorting identifier parsed from the list name, and
original gene row name. It also calculates signed significance as
`-log10(FDR) * sign(logFC)`.

```r
combined_de <- collect_edgeR_pairwise(
  list(
    Myeloid__Treated_vs_Control = myeloid_result,
    Lymphoid__Treated_vs_Control = lymphoid_result
  )
)

head(combined_de[, c(
  "Analysis", "Sorting", "Comparison", "logFC", "FDR",
  "SignedSignificance"
)])
```

Analysis names must be unique and non-empty, and each object must retain its
`group1`, `group2`, and numeric `logFC` and `FDR` fields. FDR values must be in
the interval zero to one. Exact zeros are replaced only when calculating the
finite plotting score; original FDR values remain unchanged. Annotation
columns that differ between results are preserved through a union of fields
with missing entries filled by `NA`.

### `plot_bulk_pca(...)`

Creates an exploratory sample-level PCA from the raw bulk or pseudobulk counts
in a featureCounts project. By default it uses the selected colour variable as
the `edgeR::filterByExpr()` group, performs TMM normalization, calculates log2
CPM, and fits PCA to the 5,000 most variable retained genes.

```r
pca_result <- plot_bulk_pca(
  project,
  color_by = "Genotype",
  shape_by = "Batch",
  subset = Tissue == "Lung",
  color_order = c("WT", "cKO"),
  hull = TRUE,
  label_samples = TRUE,
  save = "outputs/lung_pca.pdf"
)

pca_result$plot
pca_result$variance_explained
pca_result$pca_data
```

The colour and shape columns can be ordered explicitly and supplied with
either named mappings or unnamed vectors in factor order. Convex hulls are
drawn only for groups with at least three samples and three unique PCA
positions. Labels require `ggrepel`; the default `SampleName` label falls back
to count-matrix sample names when that metadata column is unavailable. PDF
output uses Cairo, while other figure formats are saved at 300 dpi.

The returned list retains the plot, complete `prcomp` fit, joined sample-score
table, log2-CPM matrix, variable genes, gene variances, component variance
percentages, filtered DGEList, original-row filter mask, and the exact colour
and shape mappings. The input project is not modified.

This PCA expects raw sample-level bulk or replicate-aware pseudobulk counts and
at least three selected samples. It is an exploratory quality-control view,
not a replacement for design-aware differential-expression analysis. Its
required packages are edgeR and ggplot2, plus ggrepel only when sample labels
are requested.

### `plot_bulk_qc(...)`

Builds a compact overview of sample-level sequencing and assignment QC. The
default panels show total input reads, assigned reads, count-matrix library
size, and assignment percentage. `LibrarySize` is always recalculated from the
aligned raw count matrix, so the panel reflects the project currently being
plotted rather than a potentially stale metadata field.

```r
qc_result <- plot_bulk_qc(
  project,
  color_by = "Condition",
  facet_by = "Sorting",
  color_order = c("Control", "Treated"),
  label_samples = TRUE,
  ncol = 2,
  save = "outputs/bulk_sample_qc.pdf"
)

qc_result$plot
qc_result$plots$Assignment_percent
qc_result$metadata[, c("SampleName", "LibrarySize")]
```

Each requested metric must be numeric and contain at least one observed finite
value; missing individual observations are allowed. Samples are sorted by
facet, colour group, and sample identifier before plotting. The optional
`connect_samples` lines follow this display order and should be enabled only
when connecting sequential samples within a group has a meaningful
interpretation. Sample labels require ggrepel.

The function returns the combined patchwork figure, every individual ggplot,
a sorted metadata copy containing the recalculated library size, and the exact
named colour mapping. It does not modify the input project or apply QC
exclusion thresholds. Required packages are ggplot2, patchwork, and scales;
PDF output uses Cairo and other figure formats are saved at 300 dpi.

### `plot_bulk_violin(...)`

Displays prepared log2-CPM distributions for requested genes across sample
groups. Gene symbols and stable IDs are both accepted. If a requested symbol
resolves to multiple annotation rows, the row with the highest mean expression
is selected and recorded explicitly.

```r
bulk <- prepare_bulk_expression(
  project,
  subset = Tissue == "Lung",
  prior_count = 2,
  norm_method = "TMM"
)

violin_result <- plot_bulk_violin(
  bulk,
  genes = c("GATA1", "SPI1", "CEBPA"),
  group_by = "Condition",
  split_by = "Sorting",
  group_order = c("Control", "Treated"),
  boxplot = TRUE,
  show_points = TRUE,
  save = "outputs/bulk_gene_expression.pdf"
)

violin_result$plot
violin_result$summary
violin_result$gene_mapping
```

The plot inherits the selected samples and normalization recorded in `bulk`.
Requested gene order is preserved, group and split orders can be controlled
explicitly, and point jitter is reproducible. Without
`split_by`, genes use a wrap layout of at most four columns; with `split_by`,
genes form facet rows and split values form columns. The automatic height
adapts to either layout.

`violin_scale` controls how widths are normalized across groups: `"width"`
gives each violin the same maximum width, `"area"` gives each the same total
area, and `"count"` makes width proportional to sample count. The default
`"width"` mode emphasizes distribution shape without encoding unequal group
sizes; individual points and the returned summaries retain sample-count
context.

The returned list retains the plot, its complete long-format data, per-gene
and per-panel descriptive statistics, resolved gene mapping, full normalized
log2-CPM matrix, normalized edgeR object, normalization settings, and exact
colour mapping. Metadata
columns that conflict with generated plotting fields are preserved with a
`metadata_` prefix. The plotting function requires ggplot2; preparation
requires edgeR.

This is a descriptive replicate-level view for raw bulk or replicate-aware
pseudobulk counts. It does not replace study-design-aware differential
expression, effect-size estimation, or biological replication checks.

### `plot_gene_set_heatmap(...)`

Creates a program-oriented expression heatmap from a prepared bulk-expression
object and a named list of gene symbols or stable IDs. Gene-set list order and
gene order are retained by default, and the first gene-set assignment wins
when multiple sets resolve to the same expression row.

```r
programs <- list(
  Stemness = c("GATA2", "KIT", "PROM1"),
  Myeloid = c("SPI1", "CEBPA", "MPO")
)

bulk <- prepare_bulk_expression(
  project,
  subset = Condition != "Excluded",
  prior_count = 2,
  norm_method = "TMM"
)

program_heatmap <- plot_gene_set_heatmap(
  bulk,
  gene_sets = programs,
  aggregate_by = c("Condition", "Sorting"),
  factor_orders = list(
    Condition = c("Control", "Treated")
  ),
  annotation_cols = c("Condition", "Sorting"),
  save = "outputs/gene_set_heatmap.pdf"
)

program_heatmap$gene_mapping
program_heatmap$missing_genes
program_heatmap$matrix
```

The function inherits the selected samples and normalization recorded in
`bulk`. It can show individual samples or aggregate columns using a scalar
summary function such as `mean`. Factor levels control reproducible group
order, and `order_by` can apply a final metadata-based column order. For
aggregated columns, every metadata field used for annotation or ordering must
be constant within its aggregate group; ambiguous annotations stop with an
explanatory error.

Rows can be standardized and capped before plotting. The default five-colour
purple mapping follows the actual z-score cap rather than assuming a fixed
range. `preserve_gene_order = TRUE` disables row clustering; set it to `FALSE`
when row clustering is desired. Column clustering can override the visible
metadata order.

The returned list retains both displayed and pre-scaled expression matrices,
the full normalized log2-CPM matrix, resolved and missing genes, aligned column
metadata, row-slice assignments, colour function, and normalized edgeR object.
The normalization settings are also returned for provenance.
The heatmap is descriptive and does not replace design-aware differential
expression or biological replication checks. The heatmap function requires
ComplexHeatmap and circlize; preparation requires edgeR.

### `plot_de_heatmap(...)`

Builds a direction-split ComplexHeatmap from an existing normalized,
preferably log-transformed gene-by-sample expression matrix and a compatible
differential-expression table. It selects the strongest negative and positive
effects independently, while `force_genes` can retain biologically important
genes that are present and have a non-zero effect direction.

```r
heatmap_result <- plot_de_heatmap(
  expr = voom_result$voom$E,
  meta = project$metadata,
  de = voom_result$results,
  gene_key_col = "gene_id",
  gene_label_col = "gene_name",
  group_col = "Condition",
  group_order = c("Control", "Treated"),
  annotation_cols = c("Condition", "Batch"),
  max_genes_per_side = 25,
  force_genes = c("GATA1", "SPI1"),
  label_force_genes = TRUE,
  save = "outputs/de_heatmap.pdf"
)

heatmap_result$selected_de
heatmap_result$matrix
heatmap_result$metadata
```

Automatic selection requires the adjusted-p-value and absolute log-fold-change
cutoffs and can rank first by adjusted p-value or effect magnitude. Duplicate
gene keys retain the row with the smallest adjusted p-value and then strongest
absolute effect. Gene rows are split into explicitly labelled negative and
positive directions. Row z-scores and symmetric value capping are optional;
zero-variance genes are removed before z-scoring.

`batch_col`, `batch2_col`, and numeric `covariate_cols` invoke
`limma::removeBatchEffect()` only for the displayed expression matrix.
Biological effects listed in `preserve_cols` are retained in that correction
design. Neither the supplied DE statistics nor the original expression matrix
is modified. This visualization-only correction must not be interpreted as a
replacement for modeling batch, donors, pairing, or covariates in the original
differential-expression analysis.

The function returns the assembled heatmap, displayed matrix, complete
corrected expression matrix, selected DE rows, ordered metadata, directional
row split, marked labels, and colour function. Required packages are
ComplexHeatmap and circlize, plus limma only when visualization correction is
requested. Input expression must already be normalized/log-scale; raw counts
must not be passed directly.

### `plot_signed_manhattan(...)`

Creates a signed Manhattan-style overview of complete pairwise DE tables.
Positive values indicate higher expression in the comparison group and
negative values indicate higher expression in the reference group. The
function accepts a list of edgeR results, a named list of compatible result
data frames, or the table returned by `collect_edgeR_pairwise()`.

```r
manhattan_result <- plot_signed_manhattan(
  list(
    Myeloid__Treated_vs_Control = myeloid_result,
    Lymphoid__Treated_vs_Control = lymphoid_result
  ),
  comparison_order = "Treated_vs_Control",
  comparison_labels = c(
    Treated_vs_Control = "Treated vs Control"
  ),
  gene_order = "mean_logFC",
  fdr_cutoff = 0.05,
  logfc_cutoff = 1,
  label_top_up = 5,
  label_top_down = 5,
  cap_y = 25,
  save = "outputs/signed_manhattan.pdf"
)

manhattan_result$plot
manhattan_result$summary
manhattan_result$labels
```

Gene positions can be shared across panels using mean log fold change,
alphabetical gene label, or first input appearance, or ordered independently
within each panel by log fold change. Case-insensitive alphabetical ordering is
the default. Significant and non-significant genes retain the same comparison
colour but use different opacity and point size. Faceting and positive/negative
position provide non-colour distinctions.

Automatic labels are selected separately within every
sorting-by-comparison panel, ranked first by smallest FDR and then by the
stronger directional effect; manually requested genes are always eligible.
Label placement uses a stable internal row identifier and a fixed ggrepel seed,
with expanded unclipped y-axis space for edge labels. Labels require ggrepel,
while optional rasterized points use ggrastr when installed and otherwise fall
back to ordinary ggplot2 points.

The returned data retain uncapped signed scores, capped display scores, cutoff
classes, gene keys, and a `Capped` indicator. When `cap_y` is supplied,
open upward and downward triangles mark clipped values by default. The return
also includes per-panel significance and label summaries plus the applied
cap. This overview does not replace the effect sizes, uncertainty, study
design, replication checks, or complete differential-expression tables
required for interpretation.

### Reusable R palettes

| Palette | Colours | Intended use |
|---|---|---|
| `BULK_PCA_MACARON_COLORS` | 15 muted pastel hexadecimal colours | Default ordered groups in `plot_bulk_pca()` and other bulk/pseudobulk visualizations. |
| `BULK_QC_MACARON_COLORS` | 12 muted colours beginning with neutral grey | Default ordered groups in `plot_bulk_qc()`, especially when the first group is a reference or control. |
| `BULK_VIOLIN_MACARON_COLORS` | 12 muted pastel hexadecimal colours | Default ordered groups in `plot_bulk_violin()`. |
| `DE_HEATMAP_COLORS` | Blue, white, and muted red | Default symmetric low-midpoint-high scale in `plot_de_heatmap()`. |
| `GENE_SET_HEATMAP_COLORS` | Five light-to-dark purple colours | Default expression scale in `plot_gene_set_heatmap()`. |
| `SIGNED_MANHATTAN_MACARON_COLORS` | 12 muted comparison colours | Default comparison labels in `plot_signed_manhattan()`. |

The palette is stored as an unnamed character vector so a plotting function
can assign colours consistently after applying the requested group order. Pass
a named custom vector through `colors` when stable biological labels should
always use the same colours across studies.

When another R function is supplied, it will be added to the appropriate
section in the R pipeline, documented with complete Roxygen comments for inputs
and return values, and listed here with a usage example. Python functions are
not automatically translated into R, and R functions are not automatically
translated into Python.

## Maintenance convention

For every new or updated function or palette:

1. Add it only to its original language's pipeline and place it in the logical
   workflow section.
2. Document Python functions with a complete docstring and R functions with
   complete Roxygen comments, including input and output parameters.
3. Add or revise its README entry, example, dependencies, return values, and
   scientific caveats as applicable.
4. Record the change in [`PROJECT_LOG.md`](PROJECT_LOG.md).
5. Commit the related script and documentation changes together.
