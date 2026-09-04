# Single-Cell RNA-Seq Preprocessing Pipelines Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create reusable Python and R single-cell RNA-seq preprocessing pipeline scripts, initially providing three documented Python utilities and a durable project change log.

**Architecture:** Keep one pipeline script per language and organize functions by preprocessing concern within each script. Python tests use temporary files, mocked optional dependencies, and Matplotlib's non-interactive backend; the R file establishes documentation and organization conventions without translating Python functions.

**Tech Stack:** Python 3, pathlib, pandas, NumPy, Scanpy, Matplotlib, unittest; R with Roxygen comments for future functions; Git and Markdown.

**Spec:** `docs/superpowers/specs/2026-09-04-scrna-preprocessing-pipelines-design.md`

## Global Constraints

- Make all edits under `D:/Xiaonan/CODEX_projects/Inhouse_scripts/Inhouse_scripts`.
- Add supplied functions only to the pipeline for their original language.
- Document every Python function with a detailed docstring and every R function with Roxygen comments.
- Organize functions into clear single-cell RNA-seq preprocessing stages.
- Record completed changes in `PROJECT_LOG.md` and commit each completed update.
- Preserve the public signatures `read_one_gsm(data_dir, gsm)`, `convert_genes_to_features(data_dir, overwrite=False)`, and `trim_axs(axs, N)`.

---

### Task 1: GEO 10X Sample Reader

**Files:**
- Create: `Python/__init__.py`
- Create: `Python/scRNAseq_preprocessing.py`
- Create: `tests/__init__.py`
- Create: `tests/test_scrnaseq_preprocessing.py`

**Interfaces:**
- Consumes: A directory path and GEO sample accession.
- Produces: `read_one_gsm(data_dir: str | Path, gsm: str) -> AnnData`.

- [ ] **Step 1: Write failing tests for valid loading and input discovery errors**

Create a temporary `GSM1_sample_matrix.mtx.gz`, install a fake `scanpy` module in
`sys.modules`, and assert that the reader passes `path`, `prefix`,
`var_names="gene_symbols"`, `make_unique=True`, and `gex_only=True`. Use a fake
AnnData object with a pandas `obs` frame and string observation index; assert
that `barcode` and `GSM` are added. Add cases for a missing directory, blank GSM,
no matching matrix, and multiple matching matrices.

```python
def test_reads_matrix_and_adds_metadata(self):
    (self.data_dir / "GSM1_sample_matrix.mtx.gz").touch()
    adata = SimpleNamespace(obs=pd.DataFrame(index=["AA-1", "BB-1"]))
    adata.obs_names = adata.obs.index
    reader = Mock(return_value=adata)
    fake_scanpy = SimpleNamespace(read_10x_mtx=reader)
    with patch.dict(sys.modules, {"scanpy": fake_scanpy}):
        result = read_one_gsm(self.data_dir, " GSM1 ")
    reader.assert_called_once_with(
        path=self.data_dir, prefix="GSM1_sample_",
        var_names="gene_symbols", make_unique=True, gex_only=True,
    )
    self.assertEqual(result.obs["barcode"].tolist(), ["AA-1", "BB-1"])
    self.assertEqual(result.obs["GSM"].tolist(), ["GSM1", "GSM1"])
```

- [ ] **Step 2: Run the focused tests and verify failure**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.ReadOneGsmTests -v`

Expected: FAIL because `Python/scRNAseq_preprocessing.py` or `read_one_gsm` does
not exist.

- [ ] **Step 3: Implement the reader and its documentation**

Create the module with `from __future__ import annotations`, `Path`,
`TYPE_CHECKING`, and a type-only `AnnData` import. Implement directory and GSM
validation, deterministic matrix
matching, exact prefix removal with `removesuffix("matrix.mtx.gz")`, lazy Scanpy
import, metadata columns, and return. Include a Google-style docstring covering
parameters, return value, expected files, exceptions, and an example.

```python
def _validated_directory(data_dir: str | Path) -> Path:
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_path}")
    if not data_path.is_dir():
        raise NotADirectoryError(f"Data path is not a directory: {data_path}")
    return data_path


def read_one_gsm(data_dir: str | Path, gsm: str) -> AnnData:
    data_path = _validated_directory(data_dir)
    if gsm is None or not str(gsm).strip():
        raise ValueError("gsm must be a non-empty GEO sample accession")
    gsm = str(gsm).strip()
    matrix_files = sorted(data_path.glob(f"{gsm}_*_matrix.mtx.gz"))
    if not matrix_files:
        raise FileNotFoundError(f"No matrix.mtx.gz file found for {gsm}")
    if len(matrix_files) > 1:
        matches = "\n".join(map(str, matrix_files))
        raise RuntimeError(f"Multiple matrix files found for {gsm}:\n{matches}")
    prefix = matrix_files[0].name.removesuffix("matrix.mtx.gz")
    import scanpy as sc
    adata = sc.read_10x_mtx(
        path=data_path, prefix=prefix, var_names="gene_symbols",
        make_unique=True, gex_only=True,
    )
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["GSM"] = gsm
    return adata
```

Add the full Google-style contract described above around this implementation,
including expected `{prefix}barcodes.tsv.gz`, `{prefix}features.tsv.gz`, and
`{prefix}matrix.mtx.gz` files and a concrete call example.

- [ ] **Step 4: Run the focused tests and verify success**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.ReadOneGsmTests -v`

Expected: all reader tests PASS.

- [ ] **Step 5: Commit the reader**

```text
git add Python/__init__.py Python/scRNAseq_preprocessing.py tests/__init__.py tests/test_scrnaseq_preprocessing.py
git commit -m "feat: add GEO 10X sample reader"
```

### Task 2: Legacy 10X Feature Converter

**Files:**
- Modify: `Python/scRNAseq_preprocessing.py`
- Modify: `tests/test_scrnaseq_preprocessing.py`

**Interfaces:**
- Consumes: `convert_genes_to_features(data_dir: str | Path, overwrite: bool = False)`.
- Produces: `list[Path]` containing only outputs written during the call.

- [ ] **Step 1: Write failing conversion tests**

Use temporary gzip-compressed, tab-separated fixtures. Assert that two source
columns become three output columns with `Gene Expression` in the third column,
that output order is deterministic, existing outputs are skipped by default,
`overwrite=True` replaces them, and fewer than two input columns raise
`ValueError`. Add cases for a missing directory and no matching source files.

```python
def test_converts_legacy_genes_file(self):
    source = self.data_dir / "GSM1_sample_genes.tsv.gz"
    pd.DataFrame([["ENSG1", "A"], ["ENSG2", "B"]]).to_csv(
        source, sep="\t", header=False, index=False, compression="gzip"
    )
    expected = self.data_dir / "GSM1_sample_features.tsv.gz"
    self.assertEqual(convert_genes_to_features(self.data_dir), [expected])
    observed = pd.read_csv(expected, sep="\t", header=None, dtype=str)
    self.assertEqual(observed.iloc[:, 2].tolist(), ["Gene Expression"] * 2)
```

- [ ] **Step 2: Run the converter tests and verify failure**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.ConvertGenesTests -v`

Expected: FAIL because `convert_genes_to_features` is not defined.

- [ ] **Step 3: Implement the converter and its documentation**

Add the pandas import and use `Path.glob("*_genes.tsv.gz")`, pandas string input, the first two columns,
a constant third column, and gzip output without headers or index. Validate the
directory and `overwrite` value, preserve existing outputs unless explicitly
overwriting, print a concise summary, and return converted paths. Include a
Google-style docstring covering parameters, return value, file formats,
exceptions, overwrite behavior, and an example.

```python
def convert_genes_to_features(
    data_dir: str | Path,
    overwrite: bool = False,
) -> list[Path]:
    data_path = _validated_directory(data_dir)
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a Boolean")
    gene_files = sorted(data_path.glob("*_genes.tsv.gz"))
    if not gene_files:
        raise FileNotFoundError(f"No *_genes.tsv.gz files found in: {data_path}")
    converted: list[Path] = []
    for gene_file in gene_files:
        output = gene_file.with_name(
            gene_file.name.removesuffix("_genes.tsv.gz") + "_features.tsv.gz"
        )
        if output.exists() and not overwrite:
            continue
        genes = pd.read_csv(gene_file, sep="\t", header=None, dtype=str)
        if genes.shape[1] < 2:
            raise ValueError(
                f"{gene_file.name} has {genes.shape[1]} column(s); expected at least 2"
            )
        features = genes.iloc[:, :2].copy()
        features[2] = "Gene Expression"
        features.to_csv(
            output, sep="\t", header=False, index=False, compression="gzip"
        )
        converted.append(output)
    return converted
```

Wrap this implementation in the complete documented contract described above,
including a concrete call example and the fact that skipped outputs are not
returned.

- [ ] **Step 4: Run all converter tests and verify success**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.ConvertGenesTests -v`

Expected: all converter tests PASS.

- [ ] **Step 5: Commit the converter**

```text
git add Python/scRNAseq_preprocessing.py tests/test_scrnaseq_preprocessing.py
git commit -m "feat: add legacy 10X feature converter"
```

### Task 3: Matplotlib Axes Trimmer

**Files:**
- Modify: `Python/scRNAseq_preprocessing.py`
- Modify: `tests/test_scrnaseq_preprocessing.py`

**Interfaces:**
- Consumes: `trim_axs(axs, N: int)` where `axs` is an Axes or array-like axes collection.
- Produces: A one-dimensional NumPy array containing exactly the first `N` axes.

- [ ] **Step 1: Write failing axes tests**

Set the Matplotlib backend to `Agg`. Create a 2-by-2 subplot grid, retain two
axes, and assert that the returned array has shape `(2,)` and the figure now has
two axes. Add tests for a single Axes input, zero axes retained, negative and
oversized counts, Boolean/non-integer counts, and elements without `remove()`.

```python
def test_flattens_and_removes_surplus_axes(self):
    figure, axs = plt.subplots(2, 2)
    retained = trim_axs(axs, 2)
    self.assertEqual(retained.shape, (2,))
    self.assertEqual(len(figure.axes), 2)
    self.assertIs(retained[0], axs.flat[0])
    plt.close(figure)
```

- [ ] **Step 2: Run the plotting tests and verify failure**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.TrimAxesTests -v`

Expected: FAIL because `trim_axs` is not defined.

- [ ] **Step 3: Implement the axes utility and its documentation**

Add the NumPy import and convert input with
`numpy.asarray(axs, dtype=object).reshape(-1)`. Reject Boolean
or non-integral `N` with `TypeError`, reject values outside `0..len(axes)` with
`ValueError`, validate that every element has a callable `remove`, remove
surplus axes, and return `axes[:N]`. Include a Google-style docstring covering
accepted inputs, exact output shape/type, mutation of the figure, exceptions,
and an example.

```python
def trim_axs(axs: object, N: int) -> np.ndarray:
    if isinstance(N, bool) or not isinstance(N, Integral):
        raise TypeError("N must be an integer")
    axes = np.asarray(axs, dtype=object).reshape(-1)
    if N < 0 or N > len(axes):
        raise ValueError(f"N must be between 0 and {len(axes)}, inclusive")
    if any(not callable(getattr(ax, "remove", None)) for ax in axes):
        raise TypeError("axs must contain only Axes-like objects")
    for ax in axes[N:]:
        ax.remove()
    return axes[:N]
```

Wrap this implementation in the complete documented contract described above,
including a concrete subplot example and an explicit note that the source
figure is mutated.

- [ ] **Step 4: Run the plotting tests and verify success**

Run: `python -m unittest tests.test_scrnaseq_preprocessing.TrimAxesTests -v`

Expected: all plotting tests PASS.

- [ ] **Step 5: Commit the plotting utility**

```text
git add Python/scRNAseq_preprocessing.py tests/test_scrnaseq_preprocessing.py
git commit -m "feat: add plotting axes utility"
```

### Task 4: R Pipeline Scaffold, Project Log, and Final Verification

**Files:**
- Create: `R/scRNAseq_preprocessing.R`
- Create: `PROJECT_LOG.md`
- Modify: `tests/test_scrnaseq_preprocessing.py` only if final verification exposes a test defect.

**Interfaces:**
- Consumes: The committed Python pipeline utilities and design history.
- Produces: A documented R pipeline entry point and dated repository change record.

- [ ] **Step 1: Create the R pipeline scaffold**

Add a file-level description explaining that R functions supplied later will be
organized by preprocessing stage and documented with Roxygen. Do not translate
or pre-create versions of the Python functions.

```r
#' In-house single-cell RNA-seq preprocessing utilities
#'
#' R functions supplied to this project are organized here by preprocessing
#' stage. Each function must document its inputs and outputs with Roxygen.
#'
#' @keywords internal
NULL
```

- [ ] **Step 2: Run complete Python verification**

Run: `python -m unittest discover -s tests -v`

Expected: all tests PASS.

Run: `python -m py_compile Python/scRNAseq_preprocessing.py tests/test_scrnaseq_preprocessing.py`

Expected: exit status 0 with no syntax errors.

- [ ] **Step 3: Check the R scaffold when R is available**

Run: `Rscript -e "parse(file='R/scRNAseq_preprocessing.R')"`

Expected: exit status 0. If `Rscript` is unavailable, record that fact in the
project log without installing software.

- [ ] **Step 4: Write the project log**

Create a `2026-09-04` entry listing the specification and plan commits, each
implemented function, the R scaffold, test commands and outcomes, and the
individual implementation commit identifiers. State that all changes were made
inside the D: repository.

```markdown
# Project Log

## 2026-09-04

### Added

- Initialized language-specific single-cell RNA-seq preprocessing pipelines.
- Added the documented Python functions and automated tests.
- Initialized the R pipeline without translating Python functions.

### Verification

- Record each exact verification command and its outcome.

### Commits

- Record the design, plan, and implementation commit identifiers.
```

- [ ] **Step 5: Check the final diff and repository state**

Run: `git diff --check`

Expected: no whitespace errors.

Run: `git status --short`

Expected: only `R/scRNAseq_preprocessing.R` and `PROJECT_LOG.md` are uncommitted.

- [ ] **Step 6: Commit the scaffold and project log**

```text
git add R/scRNAseq_preprocessing.R PROJECT_LOG.md
git commit -m "docs: initialize R pipeline and project log"
```

- [ ] **Step 7: Confirm clean completion**

Run: `git status --short --branch`

Expected: branch `main` with a clean working tree.
