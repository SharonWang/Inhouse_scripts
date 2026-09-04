# Single-Cell RNA-Seq Preprocessing Pipelines Design

## Purpose

Create two reusable, language-specific in-house scripts for single-cell RNA-seq
preprocessing. Functions supplied in Python will be added only to the Python
pipeline, and functions supplied in R will be added only to the R pipeline.

## Repository Structure

```text
Inhouse_scripts/
|-- Python/
|   `-- scRNAseq_preprocessing.py
|-- R/
|   `-- scRNAseq_preprocessing.R
`-- PROJECT_LOG.md
```

The scripts will be organized into clearly labelled preprocessing stages as
functions are added. Expected stages include data loading, quality control,
normalization, feature selection, dimensionality reduction, and export. Only
stages containing code will be added; unused placeholder implementations will
not be created.

## Initial Python Function

The first Python pipeline function is `read_one_gsm(data_dir, gsm)`. It will:

1. Accept a filesystem directory and a GEO sample accession.
2. Validate that the directory exists and that the accession is not empty.
3. Locate exactly one `{gsm}_*_matrix.mtx.gz` file.
4. Derive the prefix required by `scanpy.read_10x_mtx()`.
5. Load the gene-expression matrix using gene symbols with unique variable
   names and `gex_only=True`.
6. Preserve the original observation names in `adata.obs["barcode"]`.
7. Add the supplied accession to `adata.obs["GSM"]`.
8. Return the resulting `AnnData` object.

The function will raise a clear exception for an invalid directory, an empty
accession, a missing matrix file, or multiple matching matrix files. Scanpy
errors encountered while reading malformed or incomplete 10X input will remain
visible to the caller.

## Legacy 10X Feature Conversion

The Python pipeline will also provide
`convert_genes_to_features(data_dir, overwrite=False)`. It will:

1. Accept a directory containing legacy `*_genes.tsv.gz` files.
2. Validate that the directory exists and contains at least one matching file.
3. Derive each output name by replacing `_genes.tsv.gz` with
   `_features.tsv.gz`.
4. Skip an existing output by default, or replace it when `overwrite=True`.
5. Read the legacy tab-separated file as strings and require at least two
   columns.
6. Retain the first two columns and add `Gene Expression` as the third column
   required by the modern 10X feature format.
7. Write a headerless, index-free, gzip-compressed output file.
8. Return a list of `Path` objects for files converted during that call;
   skipped output files will not be included.

Conversion will process files in sorted order for deterministic results. An
invalid directory, the absence of matching files, or a source file with fewer
than two columns will raise a clear exception. Existing outputs will never be
modified unless the caller explicitly sets `overwrite=True`.

## Documentation and Style

Python functions will use type hints and Google-style docstrings that document
all parameters, return values, raised exceptions, expected input files, and a
short usage example. R functions will use Roxygen comments documenting all
parameters and return values. Both scripts will be divided into logical
preprocessing sections and will be reorganized when necessary as they grow.

## Change Recording and Version Control

`PROJECT_LOG.md` will use dated entries to record files changed, behavior added
or modified, verification performed, and the associated Git commit. Each
completed update will be committed separately. All project edits will remain
under `D:/Xiaonan/CODEX_projects/Inhouse_scripts/Inhouse_scripts`.

## Verification

The GEO reader will be tested with temporary fixture files and a mocked Scanpy
reader so that discovery, prefix derivation, metadata creation, and error paths
can be verified without downloading GEO data. The legacy feature converter will
be tested with temporary gzip-compressed fixtures covering successful
conversion, existing-output skipping, explicit overwriting, missing inputs, and
malformed inputs. The Python module will also be syntax-checked. The R script
will be parse-checked if an R runtime is available; otherwise its initialization
will be recorded as not executable until an R function is added.
