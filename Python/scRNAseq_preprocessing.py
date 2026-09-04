"""Reusable Python utilities for single-cell RNA-seq preprocessing.

Functions are grouped by their role in a typical preprocessing workflow. New
Python functions should be added to the relevant section and documented with a
complete docstring describing their inputs, outputs, side effects, and errors.
"""

from __future__ import annotations

from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData


__all__ = [
    "convert_genes_to_features",
    "read_one_gsm",
    "trim_axs",
]


# =============================================================================
# Shared filesystem helpers
# =============================================================================


def _validate_data_directory(data_dir: str | Path) -> Path:
    """Return a validated directory path used by file-based utilities."""
    data_path = Path(data_dir)

    if not data_path.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_path}")

    if not data_path.is_dir():
        raise NotADirectoryError(f"Data path is not a directory: {data_path}")

    return data_path


# =============================================================================
# Data loading and 10X file preparation
# =============================================================================


def read_one_gsm(data_dir: str | Path, gsm: str) -> AnnData:
    """Read one GEO 10X sample and attach its barcode and GSM metadata.

    The input directory must contain one matrix file whose name matches
    ``{gsm}_*_matrix.mtx.gz``. Scanpy uses the portion before
    ``matrix.mtx.gz`` as the filename prefix and expects the corresponding
    barcode and feature files to use that same prefix::

        {prefix}barcodes.tsv.gz
        {prefix}features.tsv.gz
        {prefix}matrix.mtx.gz

    Parameters
    ----------
    data_dir : str or pathlib.Path
        Directory containing the gzip-compressed 10X matrix, barcode, and
        feature files. The directory must already exist.
    gsm : str
        GEO Sample accession, such as ``"GSM7732265"``. Leading and trailing
        whitespace is removed before matching files and recording metadata.

    Returns
    -------
    anndata.AnnData
        Gene-expression matrix returned by :func:`scanpy.read_10x_mtx`.
        ``adata.obs["barcode"]`` contains the original observation names and
        ``adata.obs["GSM"]`` contains the normalized ``gsm`` value.

    Raises
    ------
    FileNotFoundError
        If ``data_dir`` does not exist or no matching matrix file is found.
    NotADirectoryError
        If ``data_dir`` exists but is not a directory.
    ValueError
        If ``gsm`` is empty or contains only whitespace.
    RuntimeError
        If more than one matrix file matches the supplied ``gsm``.
    ImportError
        If Scanpy is not installed in the active Python environment.

    Examples
    --------
    >>> adata = read_one_gsm("data/geo", "GSM7732265")
    >>> adata.obs[["barcode", "GSM"]].head()

    Notes
    -----
    Gene symbols are used as variable names and made unique. Only
    gene-expression features are retained by Scanpy.
    """
    data_path = _validate_data_directory(data_dir)

    # Normalize the accession before using it in filenames and metadata.
    if gsm is None or not str(gsm).strip():
        raise ValueError("gsm must be a non-empty GEO sample accession")
    gsm = str(gsm).strip()

    # Exactly one matrix is required so that sample selection is unambiguous.
    matrix_files = sorted(data_path.glob(f"{gsm}_*_matrix.mtx.gz"))

    if not matrix_files:
        raise FileNotFoundError(
            f"No matrix.mtx.gz file found for {gsm} in: {data_path}"
        )

    if len(matrix_files) > 1:
        matches = "\n".join(map(str, matrix_files))
        raise RuntimeError(f"Multiple matrix files found for {gsm}:\n{matches}")

    matrix_file = matrix_files[0]
    prefix = matrix_file.name.removesuffix("matrix.mtx.gz")

    print(f"Reading {gsm}")
    print(f"  prefix: {prefix}")

    # Import Scanpy only when this reader is used. This keeps file-conversion
    # and plotting helpers importable in environments that do not need Scanpy.
    import scanpy as sc

    sample_adata = sc.read_10x_mtx(
        path=data_path,
        prefix=prefix,
        var_names="gene_symbols",
        make_unique=True,
        gex_only=True,
    )

    # Preserve original barcodes before later processing can rename cells.
    sample_adata.obs["barcode"] = sample_adata.obs_names.astype(str)
    sample_adata.obs["GSM"] = gsm

    return sample_adata


def convert_genes_to_features(
    data_dir: str | Path,
    overwrite: bool = False,
) -> list[Path]:
    """Convert legacy 10X gene files to modern 10X feature files.

    Each ``*_genes.tsv.gz`` input is expected to contain at least two
    tab-separated columns: gene identifier and gene symbol. The first two
    columns are retained and a third column containing ``"Gene Expression"``
    is added to create the corresponding ``*_features.tsv.gz`` file.

    Parameters
    ----------
    data_dir : str or pathlib.Path
        Directory containing one or more ``*_genes.tsv.gz`` files. The
        directory must already exist.
    overwrite : bool, default=False
        Whether an existing ``*_features.tsv.gz`` output may be replaced. When
        ``False``, existing outputs are preserved and reported as skipped.

    Returns
    -------
    list[pathlib.Path]
        Output paths written during this call, sorted by their source
        filenames. Existing files skipped because ``overwrite=False`` are not
        included.

    Raises
    ------
    FileNotFoundError
        If ``data_dir`` does not exist or contains no ``*_genes.tsv.gz`` files.
    NotADirectoryError
        If ``data_dir`` exists but is not a directory.
    TypeError
        If ``overwrite`` is not a Boolean value.
    ValueError
        If an input file contains fewer than two columns.
    pandas.errors.ParserError
        If pandas cannot parse an input file as tab-separated data.

    Examples
    --------
    >>> outputs = convert_genes_to_features("data/10x")
    >>> for output in outputs:
    ...     print(output.name)

    Notes
    -----
    Outputs are written as headerless, index-free, gzip-compressed TSV files.
    """
    data_path = _validate_data_directory(data_dir)

    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a Boolean value")

    # Sort inputs to make conversion and return order deterministic.
    gene_files = sorted(data_path.glob("*_genes.tsv.gz"))

    if not gene_files:
        raise FileNotFoundError(
            f"No *_genes.tsv.gz files found in: {data_path}"
        )

    print(f"Found {len(gene_files)} genes.tsv.gz files\n")

    converted: list[Path] = []
    skipped: list[Path] = []

    for gene_file in gene_files:
        # Replace only the terminal legacy suffix when deriving the output.
        feature_file = gene_file.with_name(
            gene_file.name.removesuffix("_genes.tsv.gz") + "_features.tsv.gz"
        )

        if feature_file.exists() and not overwrite:
            print(f"Skipping existing: {feature_file.name}")
            skipped.append(feature_file)
            continue

        # Read every field as text so identifiers retain their exact spelling.
        genes = pd.read_csv(gene_file, sep="\t", header=None, dtype=str)

        if genes.shape[1] < 2:
            raise ValueError(
                f"{gene_file.name} has only {genes.shape[1]} column(s). "
                "Expected at least 2 columns."
            )

        # Modern 10X feature files use ID, symbol, and feature-type columns.
        features = genes.iloc[:, :2].copy()
        features[2] = "Gene Expression"

        features.to_csv(
            feature_file,
            sep="\t",
            header=False,
            index=False,
            compression="gzip",
        )

        converted.append(feature_file)
        print(
            "Converted:\n"
            f"  {gene_file.name}\n"
            f"  -> {feature_file.name}\n"
            f"  {len(features):,} features"
        )

    print(f"\nDone.\nConverted: {len(converted)}\nSkipped:   {len(skipped)}")

    return converted


# =============================================================================
# Visualization utilities
# =============================================================================


def trim_axs(axs: Any, N: int) -> np.ndarray:
    """Flatten an axes collection and remove axes beyond the first ``N``.

    This helper is useful when a rectangular subplot grid contains more panels
    than are needed for single-cell RNA-seq quality-control or exploratory
    plots. Removing the unused axes prevents empty panels from appearing.

    Parameters
    ----------
    axs : matplotlib.axes.Axes or array-like of matplotlib.axes.Axes
        One Matplotlib Axes object or a collection returned by functions such
        as :func:`matplotlib.pyplot.subplots`. Multidimensional collections are
        flattened in row-major order.
    N : int
        Number of axes to retain. Must be between zero and the total number of
        supplied axes, inclusive. Boolean values are not accepted as integers.

    Returns
    -------
    numpy.ndarray
        One-dimensional object array containing exactly the first ``N`` axes.

    Raises
    ------
    TypeError
        If ``N`` is not an integer, or if an element of ``axs`` does not expose
        a callable ``remove`` method.
    ValueError
        If ``N`` is negative or exceeds the number of supplied axes.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, axs = plt.subplots(2, 2)
    >>> axs = trim_axs(axs, 3)
    >>> len(axs)
    3

    Notes
    -----
    The source figure is modified in place because every surplus axis is
    removed from its figure.
    """
    if isinstance(N, bool) or not isinstance(N, Integral):
        raise TypeError("N must be an integer")

    # Convert scalar, list-like, and multidimensional inputs to one shape.
    axes = np.asarray(axs, dtype=object).reshape(-1)

    if N < 0 or N > len(axes):
        raise ValueError(f"N must be between 0 and {len(axes)}, inclusive")

    if any(not callable(getattr(ax, "remove", None)) for ax in axes):
        raise TypeError("axs must contain only Axes-like objects")

    for ax in axes[N:]:
        ax.remove()

    return axes[:N]
