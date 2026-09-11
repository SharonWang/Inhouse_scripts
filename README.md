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
    convert_genes_to_features,
    ordmag_filter,
    plot_adata_stacked_bar,
    plot_anndata_group_umap,
    plot_seurat_violins,
    read_one_gsm,
    trim_axs,
)
```

### Function summary

| Workflow stage | Function | What it does | Main return value(s) |
|---|---|---|---|
| 10X preparation | `convert_genes_to_features()` | Converts legacy compressed 10X gene tables into modern feature tables. | List of output `Path` objects written during the call. |
| Data loading | `read_one_gsm()` | Reads one GEO 10X matrix and records the original barcode and GSM accession. | `AnnData` object. |
| Cell calling | `ordmag_filter()` | Applies an approximate Cell Ranger OrdMag Step 1 UMI threshold. | Cell mask, UMI threshold, and per-barcode UMI totals. |
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

There are currently no reusable R functions in the pipeline. The file contains
a structured scaffold for data loading, QC and filtering, normalization and
feature selection, dimensionality reduction and clustering, and visualization
and export.

When an R function is supplied, it will be added to the appropriate section in
the R pipeline, documented with complete Roxygen comments for inputs and return
values, and listed here with a usage example. Python functions are not
automatically translated into R, and R functions are not automatically
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
