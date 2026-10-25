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
- umap-learn for fitting a reusable reference UMAP model and projecting new
  observations.
- Seaborn for grouped expression and sample-level abundance violin plots.
- statsmodels when Benjamini-Hochberg adjustment is requested for violin-plot
  or donor-level expression comparisons.

Install packages in the analysis environment appropriate for the study. This
repository does not pin an environment or alter raw input data automatically.

### Importing the utilities

When Python is started from the repository root:

```python
from Python.scRNAseq_preprocessing import (
    AGE_EXPRESSION_COLORS,
    AT2_AGE_GROUP_COLORS,
    MHCII_GROUP_COLORS,
    WT_CKO_COLORS,
    cluster_expression_summary,
    convert_genes_to_features,
    get_scanpy_de_table,
    match_genes_to_adata,
    ordmag_filter,
    plot_adata_stacked_bar,
    plot_anndata_group_umap,
    plot_at2_age_cell_style,
    plot_gene_expression_by_age,
    plot_seurat_violins,
    quick_umap,
    quick_umap_proj,
    read_one_gsm,
    score_and_assign_two_signatures,
    trim_axs,
    wilcoxon_am1_vs_am2,
)
```

### Function summary

| Workflow stage | Function | What it does | Main return value(s) |
|---|---|---|---|
| 10X preparation | `convert_genes_to_features()` | Converts legacy compressed 10X gene tables into modern feature tables. | List of output `Path` objects written during the call. |
| Data loading | `read_one_gsm()` | Reads one GEO 10X matrix and records the original barcode and GSM accession. | `AnnData` object. |
| Cell calling | `ordmag_filter()` | Applies an approximate Cell Ranger OrdMag Step 1 UMI threshold. | Cell mask, UMI threshold, and per-barcode UMI totals. |
| Gene matching | `match_genes_to_adata()` | Matches requested genes exactly and then case-insensitively against AnnData variable names. | Ordered unique list of matched AnnData gene names. |
| Signature scoring | `score_and_assign_two_signatures()` | Scores two gene programs and assigns each cell to the higher-scoring signature. | The annotated input `AnnData` object. |
| Expression summary | `cluster_expression_summary()` | Calculates mean expression and percentage detected for requested genes within each observed cell or spot group. | Long-form pandas `DataFrame`. |
| Reference dimensionality reduction | `quick_umap()` | Fits a reusable umap-learn model to reference PCA coordinates and stores its embedding. | Fitted UMAP model; coordinates are also stored in the reference `AnnData`. |
| Reference projection | `quick_umap_proj()` | Projects query PCA coordinates through a fitted reference UMAP model. | NumPy coordinate array. |
| Exploratory marker analysis | `wilcoxon_am1_vs_am2()` | Runs bidirectional Scanpy Wilcoxon comparisons between two annotated AM groups. | Two-group `AnnData`, complete DE tables, and filtered marker tables. |
| Differential-expression reporting | `get_scanpy_de_table()` | Adds group means, detection percentages, average expression, and direction to a stored Scanpy marker table. | Annotated and ordered pandas `DataFrame`. |
| Figure layout | `trim_axs()` | Removes unused Matplotlib axes from a subplot grid. | Flattened array of retained axes. |
| UMAP visualization | `plot_anndata_group_umap()` | Plots one or more categorical UMAP panels, including split and highlight modes. | Figure and axes array. |
| Composition visualization | `plot_adata_stacked_bar()` | Calculates and plots cell-composition percentages from `adata.obs`. | Figure and axes; optionally the percentage table. |
| Expression visualization | `plot_seurat_violins()` | Draws grouped gene-expression violins with optional points, boxes, and exploratory tests. | Statistics table, figure, and axes; optionally extracted expression data. |
| Sample-level abundance visualization | `plot_at2_age_cell_style()` | Compares young and old sample-level AT2 percentages with violins, box plots, replicate points, and a Mann-Whitney test. | Plotting table, figure, axis, and test result. |
| Donor-level expression visualization | `plot_gene_expression_by_age()` | Aggregates cells to donor means and compares two age groups without cell-level pseudoreplication. | Donor table, statistics table, figure, and axes. |

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

### Gene matching

#### `match_genes_to_adata(gene_list, adata)`

Matches marker genes against `adata.var_names`, trying exact spelling first
and then a case-insensitive lookup. Returned names use the spelling stored in
the AnnData object, duplicates are removed without changing request order, and
unavailable genes are omitted.

```python
matched_markers = match_genes_to_adata(
    ["SFTPC", "Krt8", "AGER"],
    adata,
)
```

When multiple variables differ only by capitalization, the first variable in
`adata.var_names` is selected deterministically. The input AnnData object is
not modified.

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

### Reference UMAP fitting and projection

#### `quick_umap(...)`

Fits a `umap-learn` model directly to `adata_ref.obsm["X_pca"]`, stores a copy
of the fitted coordinates in `adata_ref.obsm["X_umap"]`, and returns the fitted
model for later projection. The distance metric and metric keyword arguments
are recovered from `adata_ref.uns["neighbors"]["params"]` when available.

```python
umap_model = quick_umap(
    adata_ref,
    n_neighbors=15,
    min_dist=0.3,
    random_state=0,
)
```

This helper fits from PCA coordinates; it does not optimize from Scanpy's
stored connectivity graph and therefore is not intended to recreate an
existing `sc.tl.umap()` result exactly. Parameters controlling the distance
curve can be derived automatically from `spread` and `min_dist` or supplied as
an explicit `a`, `b` pair.

#### `quick_umap_proj(adata_new, umap_ref)`

Transforms new PCA coordinates through the model returned by `quick_umap()`.
The function returns the projected array without modifying the query object.

```python
query_umap = quick_umap_proj(adata_query, umap_model)
adata_query.obsm["X_umap_ref"] = query_umap
```

The query must be transformed with the same feature set, scaling, PCA loadings,
and component ordering used for the reference. Matching only the number of PCA
columns is insufficient for a scientifically valid reference projection. Save
the fitted UMAP model together with that PCA and preprocessing provenance.

### Exploratory two-group marker analysis

#### `get_scanpy_de_table(...)`

Converts a stored `scanpy.tl.rank_genes_groups()` result into an annotated
table with mean log1p expression and percentage detected in two requested
groups, average expression, detection-frequency difference, and effect
direction.

```python
import scanpy as sc

sc.tl.rank_genes_groups(
    adata,
    groupby="MHCII_group",
    groups=["MHCIIhi"],
    reference="MHCIIlo",
    method="wilcoxon",
    use_raw=True,
    key_added="MHCIIhi_vs_MHCIIlo_explicit",
)

de_table = get_scanpy_de_table(
    adata,
    groupby="MHCII_group",
    group1="MHCIIhi",
    group2="MHCIIlo",
    key="MHCIIhi_vs_MHCIIlo_explicit",
    use_raw=True,
    expression_threshold=0,
)
```

The function does not rerun differential expression or modify the AnnData
object. The selected expression source must contain finite normalized,
log1p-transformed values with unique gene names. It validates the stored
grouping, reference, raw-data setting, and layer against the requested
contrast. Stored comparisons against `rest` are rejected because AnnData
subsetting can preserve results that historically included additional groups;
rerun Scanpy with an explicit `reference=group2` before using this helper.

#### `wilcoxon_am1_vs_am2(...)`

Subsets an AnnData object to two annotation groups, runs Scanpy's Wilcoxon
marker procedure in both directions, and returns complete and filtered tables
for each group. Although the defaults use `MHCIIhi` and `MHCIIlo`, both labels
and the grouping column are configurable.

```python
de = wilcoxon_am1_vs_am2(
    adata_am,
    group_col="MHCII_group",
    am1_label="MHCIIhi",
    am2_label="MHCIIlo",
    use_raw=True,
    min_logfc=0.25,
    max_padj=0.05,
    min_pct=0.10,
)

mhcii_high_markers = de["am1_up"]
mhcii_low_markers = de["am2_up"]
```

The returned dictionary contains the copied two-group subset in `adata`, full
tables in `am1_all` and `am2_all`, filtered positive-marker tables in `am1_up`
and `am2_up`, and the Scanpy result key in `key_added`. Filtering uses positive
log-fold change, adjusted p-value, and—when Scanpy returns it—the fraction of
cells expressing the gene in the focal group. The function also prints the
included cell counts and the number of retained markers for each group.

The selected expression source—`adata.raw`, a named layer, or `adata.X`—must
contain normalized, log-transformed expression rather than raw counts.
Scanpy's approximate log-fold-change estimates assume logarithmized input, and
the function cannot reliably reconstruct that preprocessing provenance from an
AnnData object alone.

This is a cell-level exploratory marker-ranking workflow: cells are treated as
observations, so the Wilcoxon p-values do not model biological replication.
Use a donor-aware model or sample-level pseudobulk method for inference across
replicated samples.

### Visualization

#### `plot_gene_expression_by_age(...)`

Aggregates selected expression values to one mean per biological donor before
comparing two age groups. Donor/group combinations below the configured cell
minimum are excluded, and donors appearing in both groups trigger an error
because the implemented Mann-Whitney test is unpaired.

```python
donor_df, stats_df, fig, axes = plot_gene_expression_by_age(
    adata,
    genes=["PLAUR", "GDF15", "NR4A1"],
    group_col="Age_group",
    donor_col="Sample",
    group_order=["Young", "Old"],
    layer="log1p",
    min_cells_per_donor=20,
    colors=AGE_EXPRESSION_COLORS,
    adjust_p=True,
    save="outputs/age_expression.pdf",
)
```

Tests with fewer than two retained donors in either group are reported as
missing. Benjamini-Hochberg adjustment is applied across plotted genes when
requested. Sparse cell-level matrices remain sparse during donor aggregation.
This donor-level view avoids cell-level pseudoreplication but does not model
batch, covariates, paired donors, or repeated measures.

#### `plot_at2_age_cell_style(...)`

Creates a sample-level AT2 abundance comparison between young and old animals.
Leading numeric ages are parsed from sample names such as `18_Lung_OF1`;
intermediate ages between the two configured thresholds are omitted. Input
abundances must be percentages on the inclusive 0–100 scale.
Unparseable sample ages and nonnumeric or non-finite AT2 percentages stop with
an error naming the affected samples rather than silently changing group size.

```python
plot_df, fig, ax, test = plot_at2_age_cell_style(
    pct,
    at2_row="AT2",
    young_max_age=6,
    old_min_age=18,
    palette=AT2_AGE_GROUP_COLORS,
    save="outputs/AT2_young_old.pdf",
)
```

The figure combines narrow violins, internal box plots, jittered biological-
replicate points, and a two-sided Mann-Whitney U annotation. The statistical
test uses samples—not individual cells—as observations, but it does not model
batch, covariates, pairing, or repeated measures.

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
| `AGE_EXPRESSION_COLORS` | `Young: #78B7A3`, `Old: #C17C3B` | Donor-level Young-versus-Old gene-expression plots. |
| `AT2_AGE_GROUP_COLORS` | `Young: #78B7A3`, `Old: #D98973` | Sample-level young-versus-old AT2 abundance plots. |
| `MHCII_GROUP_COLORS` | `MHCIIhi: #E7B2B6`, `MHCIIlo: #B7C3E0` | High- and low-MHCII annotation groups. |
| `WT_CKO_COLORS` | `WT: #898C8B`, `cKO: #9CB4CC` | Wild-type and conditional-knockout groups. |

Pass these dictionaries to plotting functions through their `palette` or
`colors` parameter. Supplied reusable palettes are saved in the Python module
and added to this table for future studies.

## R pipeline

Source: [`R/scRNAseq_preprocessing.R`](R/scRNAseq_preprocessing.R)

The R pipeline contains reusable single-cell Monocle trajectory utilities and
bulk/pseudobulk import, validation, differential-expression, and visualization
workflows. The Monocle workflow uses monocle3, igraph, Matrix,
SummarizedExperiment, SingleCellExperiment, ggplot2, ComplexHeatmap, circlize,
and viridisLite in the functions that require them.

### Function summary

| Workflow stage | Function | What it does | Main return value(s) |
|---|---|---|---|
| Monocle node mapping | `resolve_principal_node_names()` | Normalizes closest-vertex output returned as principal-node names or one-based indices. | Character vector of principal-node names. |
| Monocle path selection | `get_two_monocle_paths()` | Finds two ordered shortest paths from one root to two principal-graph endpoints. | Node paths, shared nodes, last shared node, and nested-path flag. |
| Monocle path visualization | `plot_two_monocle_paths()` | Draws the shared trajectory trunk and two directional branches over a cell-type UMAP. | ggplot object. |
| Monocle pseudotime visualization | `plot_two_monocle_paths_pseudotime()` | Overlays the selected paths on a continuously pseudotime-coloured UMAP. | ggplot object. |
| Monocle path subsetting | `split_path_by_genotype()` | Retains finite-pseudotime cells and splits the path into configurable WT and cKO labels. | Named list of two cell-data-set subsets. |
| Monocle path QC | `check_path()` | Summarizes path composition and finite pseudotime coverage by genotype. | Per-genotype pseudotime summary data frame. |
| Monocle pseudotime modelling | `run_monocle_Q1()` | Fits full and reduced negative-binomial models to test nonlinear pseudotime associations. | Filtered CDS, model tables, complete/significant results, gene filters, formulas, and settings. |
| Monocle pseudotime heatmap | `make_Q1_pseudotime_heatmap()` | Predicts fitted expression along ordered pseudotime and builds a row-scaled ComplexHeatmap. | Heatmap, fitted matrices, ordered metadata, gene labels, and prediction diagnostics. |
| Monocle persistence | `save_monocle_compatible()` | Saves a Monocle object while preserving legacy Seurat-associated dimensional-reduction entries separately. | Invisible output-directory path. |
| Monocle persistence | `load_monocle_compatible()` | Reloads a saved Monocle object and restores its preserved auxiliary entries. | Restored cell-data-set object. |
| CellChat analysis | `run_cellchat()` | Runs a condition-specific CellChat workflow from a Seurat RNA assay after annotation and minimum-cell filtering. | Aggregated CellChat object. |
| Population expression summary | `calculate_expression_summary()` | Calculates detection fraction and mean normalized expression from aligned count/expression matrices. | Gene-level summary data frame. |
| Custom LRT construction | `build_lrt_from_custom_lr()` | Combines a custom ligand-receptor table, receiver targets, and NicheNet prior into prioritized ligand-receptor-target links. | Link, assignment, ordering, and provenance tables. |
| LRT circos visualization | `plot_lrt_circos()` | Draws the newer order-aware ligand-receptor-target circos implementation with exact sector validation. | Invisible plotting inputs, assignments, colours, and sector order. |
| LR chord visualization | `plot_lr_chord()` | Draws a directional ligand-to-receptor chord diagram with comparable link alpha scaling. | Invisible chord data and resolved colours. |
| NicheNet GO enrichment | `run_nichenet_gobp()` | Maps grouped symbols and runs GO enrichment against an explicit eligible-gene background. | Mappings, enrichment objects/results, and settings. |
| NicheNet GO heatmap | `plot_nichenet_gobp_heatmap()` | Builds a group-by-term adjusted-significance ComplexHeatmap. | Heatmap, matrix, selected terms, and colour limit. |
| NicheNet GO export | `save_nichenet_gobp_heatmap()` | Saves the GO heatmap as PDF, PNG, or TIFF. | Invisible normalized output path. |
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
| TF activity inference | `infer_bulk_tf_activity()` | Infers sample-level transcription-factor activities from prepared expression and a signed regulator-target network. | `bulk_tf_activity` object with activity matrices, coverage, mapping, QC, and provenance. |
| TF differential activity | `limma_tf_pairwise()` | Tests pairwise TF-activity differences with optional covariates and repeated-measure blocking. | Results, limma fit/design/contrast, selected metadata, and model diagnostics. |
| TF activity heatmap | `plot_bulk_tf_heatmap()` | Displays significant TF activities across subtype-condition groups or individual samples. | Plot, raw/z-score matrices, annotations, clustering, selected TFs, and sample counts. |
| Gene DE preparation | `prepare_bulk_gene_de()` | Collects annotated edgeR pairwise results and builds an explicit comparison manifest. | `bulk_gene_de` object with combined table, comparison information, and QC. |
| Significant-gene heatmap | `plot_bulk_gene_expression_heatmap()` | Plots the union of significant DE genes using prepared log2 CPM values. | Heatmap result, selected genes, expression matrices, annotations, and gene mapping. |
| Gene-set enrichment | `run_bulk_gene_enrichment()` | Runs per-contrast GSEA and optional directional ORA with GO or custom gene sets. | Enrichment tables/objects, rankings, universes, QC, settings, and provenance. |
| Comparison heatmap | `plot_bulk_comparison_heatmap()` | Displays gene logFC or pathway NES across an explicit comparison grid while preserving missing cells. | Heatmap, effect/FDR matrices, annotations, selected features, and comparison manifest. |
| Condition-response DE | `edgeR_condition_response()` | Tests an overall categorical condition response or numerical-duration trend in one edgeR model. | Results, DGEList, fit/test, design, coefficients, metadata, and settings. |
| Paired expression plot | `plot_bulk_paired_expression()` | Plots tidy normalized expression with optional within-subject connecting lines and facets. | A ggplot object. |
| Signaling expression preparation | `prepare_bulk_signaling_expression()` | Produces symbol-indexed log2 CPM and unlogged CPM matrices for ligand/receptor availability checks. | Expression matrices, aligned metadata, and gene mapping. |
| NicheNet analysis | `run_bulk_nichenet()` | Combines receiver responses with separate sender, receptor, and target evidence using a species-matched NicheNet prior. | Ligand activities, evidence tables, targets, QC, gene sets/backgrounds, mappings, and settings. |
| NicheNet activity heatmap | `plot_bulk_nichenet_activity()` | Displays corrected AUPR ligand rankings across receiver-condition contexts. | Activity matrix, annotations, and heatmap object. |
| General expression heatmap | `plot_bulk_expression_heatmap()` | Plots requested genes from prepared log2 CPM by sample or metadata-group mean. | Heatmap object and requested row/column clustering trees. |
| Ligand-receptor plot preparation | `prepare_bulk_nichenet_lr_plot_data()` | Joins top ligand rankings with sender ligand and receiver receptor evidence. | Long ligand-receptor plotting data. |
| Ligand-receptor dot plot | `plot_bulk_nichenet_lr_dotplot()` | Shows ligand and receptor expression/effect evidence across senders and conditions. | Plot, plotting data, joined evidence, and pair order. |
| Ligand-receptor connectivity | `plot_bulk_nichenet_lr_connectivity()` | Displays binary prior connectivity with sender-availability annotations. | Heatmap, connectivity matrix, annotations, evidence, and availability rule. |

### Monocle trajectory workflow

Select two rooted principal-graph routes and visualize them over cell type or
pseudotime:

```r
resolved_nodes <- resolve_principal_node_names(
  closest_vertex_values,
  principal_node_names = igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
)

paths <- get_two_monocle_paths(
  cds,
  root_node = "Y_1",
  AT2_1_endpoint = "Y_8",
  AT1_endpoint = "Y_15"
)

celltype_plot <- plot_two_monocle_paths(cds, paths)
pseudotime_plot <- plot_two_monocle_paths_pseudotime(cds, paths)
path_summary <- check_path(cds)
genotype_paths <- split_path_by_genotype(cds)
```

`get_two_monocle_paths()` reports whether the routes diverge or one endpoint is
nested within the other. Plotting functions validate principal-node names and
coordinates and do not mutate the input object. The supplied cell-type and
BuPu pseudotime colours are stored as reusable palettes.
`resolve_principal_node_names()` supports the named and one-based numeric
closest-vertex formats returned by different Monocle versions and rejects
mixed, fractional, out-of-range, or otherwise ambiguous values.
Numeric inputs are always interpreted as indices; character values that match
node names are preserved before digit-string index fallback is considered.

Test genes for nonlinear association with rooted pseudotime, optionally
adjusting for sample identity as a fixed model term, then visualize fitted
profiles:

```r
q1 <- run_monocle_Q1(
  path_cds = cds,
  sample_column = "Sample_ID",
  spline_df = 3,
  min_cells_expressed = 20
)

heatmap_result <- make_Q1_pseudotime_heatmap(
  q1,
  genes = q1$significant$id,
  title = "Genes associated with pseudotime"
)
ComplexHeatmap::draw(heatmap_result$heatmap)
```

The model comparison treats cells as observations. A fixed sample term can
adjust mean sample differences but does not make cells independent biological
replicates; use the results as trajectory-associated gene evidence with
appropriate replicate-aware interpretation. The heatmap shows fitted model
predictions for visualization and is not a separate significance test.

Use the compatible persistence helpers when legacy `gene_loadings` or
`prop_var_expl` entries prevent direct Monocle serialization:

```r
save_monocle_compatible(cds, "outputs/monocle_checkpoint", verify = TRUE)
restored_cds <- load_monocle_compatible("outputs/monocle_checkpoint")
```

Saving requires a new output location and writes a sibling
`_seurat_aux.rds` file. With `verify=TRUE`, the saved object is reloaded and
checked before the function returns.

### Cell-cell communication and NicheNet utilities

Run CellChat for one condition after removing missing annotations and groups
below the requested cell threshold:

```r
cellchat_wt <- run_cellchat(
  seurat_object,
  condition = "WT",
  condition_col = "sample",
  group_col = "cell_type_subtype",
  min_cells = 10
)
```

The helper supports Seurat v5 layered RNA assays, creates normalized data when
needed, aligns metadata to expression columns, and returns the aggregated
CellChat object. Cell-level communication probabilities are exploratory; use
replicated sample-aware analyses for condition-level biological inference.

These utilities use CellChat, Seurat, SeuratObject, Matrix, dplyr, tidyr,
tibble, purrr, scales, circlize, ComplexHeatmap, nichenetr, clusterProfiler,
and an organism annotation package such as org.Mm.eg.db. Install the packages
needed by the functions you call; they are referenced explicitly rather than
attached automatically.

Build and visualize custom ligand-receptor-target evidence:

```r
lrt <- build_lrt_from_custom_lr(
  lr_table,
  AT2_target_genes,
  ligand_target_matrix
)

circos_info <- plot_lrt_circos(
  lrt,
  sender_label = "MHCIIhi",
  sender_colour = "#D98989",
  title_text = "MHCIIhi to AT2 signaling",
  common_lt_range = c(0, 1),
  common_lr_range = c(0, 1),
  lr_order_table = lr_table
)

plot_lr_chord(
  lr_table,
  title = "Ligand-receptor evidence",
  ligand_color = "#D98989"
)
```

The repository retains the newest supplied `plot_lrt_circos()` version. This
version accepts `lr_order_table` as an ordering guide only: fitted pairs and
weights always come from `lrt`, pairs absent from the guide are appended, and
guide-only pairs are ignored. It preserves ligand and fitted-pair order across
panels and validates the exact circos sector set before drawing. The three
identical duplicate attachments and the older implementation were not added.

Run GO enrichment with an explicit eligible-gene universe and create a shared
heatmap:

```r
go_result <- run_nichenet_gobp(
  gene_list = list(MHCIIhi = genes_hi, MHCIIlo = genes_lo),
  background_symbols = tested_background,
  comparison_name = "AT2 response"
)

go_heatmap <- plot_nichenet_gobp_heatmap(go_result$result_table)
save_nichenet_gobp_heatmap(go_heatmap, "outputs/nichenet_go.pdf")
```

The enrichment background must be the genes eligible for selection by the
upstream workflow. GO enrichment and prioritized ligand-target links do not by
themselves establish causal signaling. Foreground identifiers are restricted
to the mapped background universe before enrichment; the conversion summary
reports mapped, eligible, and out-of-universe counts. In the heatmap, zero
means that a group-term combination was absent from the supplied significant
results rather than that it was tested with an exact adjusted p-value of one.

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

### Advanced bulk transcription-factor workflow

`infer_bulk_tf_activity()` uses prepared expression and a signed, species-
matched TF-target network to calculate sample-level activities. Differential
activity is tested separately with `limma_tf_pairwise()`; its
`DeltaActivity` is a score difference, not a log fold change.

```r
tf <- infer_bulk_tf_activity(
  bulk,
  network = collectri_mouse,
  group_cols = c("ATM_subtype", "Condition"),
  min_targets = 10
)

tf_de <- limma_tf_pairwise(
  tf,
  group_col = "Condition",
  group1 = "ND",
  group2 = "2W",
  subset = ATM_subtype == "pATM",
  covariates = "Batch"
)

tf_de_table <- tf_de$results
tf_de_table$ATM_subtype <- "pATM"

tf_heatmap <- plot_bulk_tf_heatmap(
  tf,
  de_results = tf_de_table,
  subtype_order = "pATM",
  condition_order = c("ND", "2W"),
  subtype_colors = subtype_colors,
  condition_colors = condition_colors
)
```

Inference requires edgeR and decoupleR; differential testing requires limma;
the heatmap requires pheatmap. Network organism, identifiers, sign/weight
semantics, sample groups, covariates, and repeated-measure blocks must be
verified for the study design.

### Gene differential-expression and enrichment workflow

`prepare_bulk_gene_de()` consolidates named `edgeR_pairwise()` results into a
comparison-aware object. Its table can drive expression, effect-size, and
enrichment views without treating missing comparisons as zero.

```r
gene_de <- prepare_bulk_gene_de(pairwise_results)

gene_heatmap <- plot_bulk_gene_expression_heatmap(
  bulk,
  gene_de,
  subtype_colors = subtype_colors,
  condition_colors = condition_colors,
  fdr_cutoff = 0.05
)

enrichment <- run_bulk_gene_enrichment(
  gene_de,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  run_ora = TRUE
)

effect_heatmap <- plot_bulk_comparison_heatmap(
  gene_de$table,
  gene_de$comparison_info,
  subtype_colors = subtype_colors,
  condition_colors = condition_colors
)
```

`edgeR_condition_response()` is a separate all-sample model for an omnibus
condition test or a verified numerical-duration trend. Omnibus F statistics
are unsigned and must not be used as signed enrichment rankings. Enrichment
requires clusterProfiler and AnnotationDbi plus either an organism annotation
database or a custom `TERM2GENE` mapping.

### General and paired expression views

`plot_bulk_expression_heatmap()` shows an arbitrary requested gene list from a
prepared object, either by individual sample or metadata-group mean. It does
not rerun normalization or select genes by DE significance.

```r
expression_heatmap <- plot_bulk_expression_heatmap(
  bulk,
  genes = c("Il6", "Ccl2", "Tnf"),
  column_mode = "mean",
  group_by = c("ATM_subtype", "Condition"),
  column_order = list(
    ATM_subtype = c("pATM", "cATM", "sATM"),
    Condition = c("ND", "2W", "1M", "2M")
  )
)

paired_plot <- plot_bulk_paired_expression(
  tidy_expression,
  x = "Sorting",
  y = "Expression",
  pair_by = "MouseID",
  facet_by = "Condition"
)
```

Pairing lines should be used only when `pair_by` is a verified subject ID and
each subject has at most one observation per facet and x-axis group.

### Bulk NicheNet signaling workflow

The signaling workflow keeps receiver response, sender availability, receptor
evidence, and prior connectivity separate. NicheNet activity scores are
predictive rankings rather than p-values, and the unsigned prior does not by
itself establish activation or repression.

```r
sender_signaling <- prepare_bulk_signaling_expression(sender_bulk)
receiver_signaling <- prepare_bulk_signaling_expression(receiver_bulk)

nichenet_result <- run_bulk_nichenet(
  sender_bulk,
  receiver_bulk,
  sender_de,
  receiver_de,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  reference = "ND"
)

activity_plot <- plot_bulk_nichenet_activity(nichenet_result, response = "up")
lr_data <- prepare_bulk_nichenet_lr_plot_data(
  nichenet_result,
  receiver = "ASC",
  response = "up"
)
lr_dotplot <- plot_bulk_nichenet_lr_dotplot(
  nichenet_result,
  receiver = "ASC",
  response = "up"
)
lr_connectivity <- plot_bulk_nichenet_lr_connectivity(
  nichenet_result,
  receiver = "ASC",
  condition = "2W",
  response = "up"
)
```

This workflow requires edgeR, nichenetr, ggplot2, and pheatmap as used by the
individual functions. Sender/receiver DE contrasts must use the declared
condition versus reference orientation, and all identifiers must match the
species-specific NicheNet model.

### Reusable R palettes

| Palette | Colours | Intended use |
|---|---|---|
| `BULK_PCA_MACARON_COLORS` | 15 muted pastel hexadecimal colours | Default ordered groups in `plot_bulk_pca()` and other bulk/pseudobulk visualizations. |
| `BULK_QC_MACARON_COLORS` | 12 muted colours beginning with neutral grey | Default ordered groups in `plot_bulk_qc()`, especially when the first group is a reference or control. |
| `BULK_VIOLIN_MACARON_COLORS` | 12 muted pastel hexadecimal colours | Default ordered groups in `plot_bulk_violin()`. |
| `DE_HEATMAP_COLORS` | Blue, white, and muted red | Default symmetric low-midpoint-high scale in `plot_de_heatmap()`. |
| `GENE_SET_HEATMAP_COLORS` | Five light-to-dark purple colours | Default expression scale in `plot_gene_set_heatmap()`. |
| `SIGNED_MANHATTAN_MACARON_COLORS` | 12 muted comparison colours | Default comparison labels in `plot_signed_manhattan()`. |
| `BULK_PAIRED_EXPRESSION_COLORS` | Six muted macaron colours | Default groups in `plot_bulk_paired_expression()`. |
| `BULK_SIGNALING_DIVERGING_COLORS` | Blue, white, and muted red anchors | Default signed ligand/receptor effect scale. |
| `BULK_ACTIVITY_DIVERGING_COLORS` | 100-colour interpolation of the signaling anchors | Default TF, expression, and comparison heatmaps. |
| `BULK_NICHENET_ACTIVITY_COLORS` | 100-colour white-to-purple gradient | Corrected AUPR ligand-activity heatmaps. |
| `BULK_NICHENET_SENDER_COLORS` | Named pATM, cATM, and sATM colours | Sender annotations in NicheNet connectivity plots. |
| `BULK_NICHENET_CONNECTIVITY_COLORS` | Two colours for absent/present links | Binary NicheNet connectivity matrices. |
| `LRT_RECEPTOR_MACARON_COLORS` | 16 supplied macaron colours | Default receptor assignments in `plot_lrt_circos()`. |
| `MONOCLE_PATH_CELLTYPE_COLORS` | Named AT2-2, AT2-1, AT1, and Club-cell colours | Cell-type backgrounds in selected Monocle path plots. |
| `MONOCLE_HEATMAP_CELLTYPE_COLORS` | Named semi-transparent AT2-2, AT2-1, AT1, and Club-cell colours | Cell-type annotation in Monocle pseudotime heatmaps. |
| `MONOCLE_PSEUDOTIME_COLORS` | 100-colour light-to-dark BuPu gradient | Continuous Monocle pseudotime UMAP scale. |

Unnamed palettes are assigned after applying the requested group order; named
palettes preserve stable biological labels across plots. Pass a documented
custom vector through the relevant colour argument when a study requires a
different mapping.

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
