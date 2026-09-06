"""Reusable Python utilities for single-cell RNA-seq preprocessing.

Functions are grouped by their role in a typical preprocessing workflow. New
Python functions should be added to the relevant section and documented with a
complete docstring describing their inputs, outputs, side effects, and errors.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, PercentFormatter
from scipy import sparse

if TYPE_CHECKING:
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


__all__ = [
    "MHCII_GROUP_COLORS",
    "WT_CKO_COLORS",
    "convert_genes_to_features",
    "ordmag_filter",
    "plot_adata_stacked_bar",
    "plot_anndata_group_umap",
    "plot_seurat_violins",
    "read_one_gsm",
    "trim_axs",
]


# =============================================================================
# Reusable colour palettes
# =============================================================================


MHCII_GROUP_COLORS = {
    "MHCIIhi": "#E7B2B6",
    "MHCIIlo": "#B7C3E0",
}
"""Reusable colour mapping for high- and low-MHCII annotation groups."""


WT_CKO_COLORS = {
    "WT": "#DCE8F2",
    "cKO": "#F3C77F",
}
"""Reusable colour mapping for wild-type and conditional-knockout groups."""


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
# Quality control and cell calling
# =============================================================================


def ordmag_filter(
    adata: AnnData,
    expect_cells: int = 8000,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Call cells using an approximate Cell Ranger OrdMag Step 1 filter.

    The function ranks barcodes by their total UMI counts, calculates the 99th
    percentile among the highest-ranked ``expect_cells`` barcodes, and retains
    barcodes with total counts of at least one tenth of that value.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing the raw UMI count matrix in ``adata.X``.
        Rows must represent barcodes or cells and columns must represent genes.
        Both dense arrays and SciPy sparse matrices are supported.
    expect_cells : int, default=8000
        Approximate number of expected cells. The calculation uses the smaller
        of this value and the number of available barcodes. Must be a positive
        integer; Boolean values are not accepted.

    Returns
    -------
    keep : numpy.ndarray
        One-dimensional Boolean array with one value per barcode. ``True``
        identifies a barcode whose total UMI count meets the OrdMag cutoff.
    threshold : float
        OrdMag UMI cutoff, calculated as one tenth of the 99th percentile of
        the top-ranked barcodes.
    total_umi : numpy.ndarray
        One-dimensional array containing the total UMI count for every barcode
        in the original row order.

    Raises
    ------
    TypeError
        If ``expect_cells`` is not an integer.
    ValueError
        If ``expect_cells`` is not positive or ``adata.X`` has no barcode rows.
    AttributeError
        If ``adata`` does not provide an ``X`` count matrix.

    Examples
    --------
    >>> keep, threshold, total_umi = ordmag_filter(adata, expect_cells=8000)
    >>> filtered_adata = adata[keep].copy()

    Notes
    -----
    This is an approximation of Cell Ranger's OrdMag Step 1 procedure. It does
    not perform EmptyDrops-style testing or other later cell-calling steps.
    """
    if isinstance(expect_cells, bool) or not isinstance(expect_cells, Integral):
        raise TypeError("expect_cells must be an integer")

    if expect_cells <= 0:
        raise ValueError("expect_cells must be greater than zero")

    # Sum raw UMI counts across genes while supporting both sparse and dense X.
    if sparse.issparse(adata.X):
        total_umi = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        total_umi = np.asarray(adata.X.sum(axis=1)).ravel()

    if total_umi.size == 0:
        raise ValueError("adata.X must contain at least one barcode row")

    # Rank totals without changing the order of values returned to the caller.
    ranked = np.sort(total_umi)[::-1]
    n_barcodes = min(expect_cells, len(ranked))
    top_n = ranked[:n_barcodes]

    # Cell Ranger's OrdMag heuristic uses one tenth of a robust upper count.
    robust_maximum = float(np.percentile(top_n, 99))
    threshold = robust_maximum / 10
    keep = total_umi >= threshold

    print(f"Expected cells: {expect_cells:,}")
    print(f"99th percentile m: {robust_maximum:,.1f} UMIs")
    print(f"OrdMag cutoff m/10: {threshold:,.1f} UMIs")
    print(f"Called cells: {keep.sum():,}")

    return keep, threshold, total_umi


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


def plot_anndata_group_umap(
    adata: AnnData,
    group_col: str,
    split_by: str | None = None,
    split_categories: Iterable[Any] | None = None,
    palette: dict[Any, Any] | None = None,
    point_size: float = 18,
    point_alpha: float = 0.5,
    na_color: Any = "#D3D3D3",
    umap_key: str = "X_umap",
    top_group: Any | None = None,
    width_cm: float = 4,
    height_cm: float = 4,
    left_cm: float = 1.20,
    right_cm: float = 0.25,
    bottom_cm: float = 1.05,
    top_cm: float = 0.45,
    panel_gap_cm: float = 0.50,
    show_legend: bool = True,
    legend_title: str | None = None,
    legend_marker_size: float = 6,
    legend_fontsize: float = 8,
    legend_title_size: float = 9,
    legend_gap_cm: float = 0.35,
    legend_width_cm: float = 3.5,
    legend_ncol: int = 2,
    legend_columnspacing: float = 1.2,
    legend_handletextpad: float = 0.5,
    title: str | None = None,
    xlabel: str = "UMAP1",
    ylabel: str = "UMAP2",
    axis_label_size: float = 10,
    title_size: float = 11,
    split_title_size: float = 10,
    tick_label_size: float = 8,
    spine_width: float = 1.0,
    tick_width: float = 1.0,
    tick_length: float = 3,
    tick_nbins: int = 4,
    padding_fraction: float = 0.03,
    shared_limits: bool = True,
    rasterized: bool = False,
    transparent: bool = True,
    save: str | Path | None = None,
    dpi: float = 600,
) -> tuple[Figure, np.ndarray]:
    """Plot categorical AnnData annotations on one or more UMAP panels.

    With no split, one UMAP is coloured by ``group_col``. Splitting by another
    annotation creates one cell-subset panel per split category. Splitting by
    ``group_col`` instead creates highlight panels: every panel contains all
    cells, but only its selected category is coloured and drawn on top. Each
    panel uses one scatter call, allowing points to be rasterized together while
    axes, labels, titles, ticks, and the legend remain vector objects.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object whose ``obs`` table contains ``group_col`` and whose
        ``obsm`` mapping contains two-dimensional UMAP coordinates.
    group_col : str
        Name of the categorical column in ``adata.obs`` used to colour cells.
        Non-categorical values are converted to a pandas categorical series.
    split_by : str, optional
        Annotation column used to create multiple panels. If different from
        ``group_col``, each panel contains only cells from one split category.
        If equal to ``group_col``, each panel contains all cells and highlights
        one group. If omitted, a single unsplit panel is drawn.
    split_categories : iterable, optional
        Ordered subset of observed ``split_by`` categories to plot. By default,
        every observed category is used in categorical order. Ignored when
        ``split_by`` is omitted.
    palette : dict, optional
        Mapping from every observed non-missing group to a Matplotlib-compatible
        colour. If omitted, colours are read from
        ``adata.uns[f"{group_col}_colors"]`` when available; otherwise the
        Matplotlib ``tab20`` colour map is used. A supplied mapping is not
        modified or stored in ``adata``.
    point_size : float, default=18
        Marker area passed to ``Axes.scatter`` in points squared.
    point_alpha : float, default=0.5
        Opacity of plotted cell markers.
    na_color : color, default="#D3D3D3"
        Matplotlib-compatible colour used for cells with missing group values.
    umap_key : str, default="X_umap"
        Key in ``adata.obsm`` containing at least two UMAP coordinate columns.
    top_group : optional
        Category to draw after all other groups so its cells appear visually on
        top. The value must exactly match one of the observed, non-missing
        categories in ``group_col``. This affects unsplit and normal split
        panels; highlight panels always draw their selected group last.
    width_cm : float, default=4
        Physical width of each UMAP plotting box in centimetres.
    height_cm : float, default=4
        Physical height of each UMAP plotting box in centimetres.
    left_cm : float, default=1.20
        Space to the left of the first plotting box in centimetres.
    right_cm : float, default=0.25
        Space after the last plotting box in centimetres.
    bottom_cm : float, default=1.05
        Space below the plotting box in centimetres.
    top_cm : float, default=0.45
        Space above the plotting box in centimetres.
    panel_gap_cm : float, default=0.50
        Horizontal gap between adjacent UMAP panels in centimetres.
    show_legend : bool, default=True
        Whether to create a figure-level legend for groups and missing values.
    legend_title : str, optional
        Legend heading. Defaults to ``group_col`` when the legend is displayed.
    legend_marker_size : float, default=6
        Diameter of legend markers in points.
    legend_fontsize : float, default=8
        Font size of legend labels in points.
    legend_title_size : float, default=9
        Font size of the legend title in points.
    legend_gap_cm : float, default=0.35
        Horizontal gap before the legend area in centimetres.
    legend_width_cm : float, default=3.5
        Width reserved for the legend in centimetres.
    legend_ncol : int, default=2
        Number of columns used to arrange legend entries.
    legend_columnspacing : float, default=1.2
        Horizontal spacing between legend columns in font-size units.
    legend_handletextpad : float, default=0.5
        Gap between each legend marker and its label in font-size units.
    title : str, optional
        Overall plot title, horizontally centred over the combined panel area
        and vertically centred within the reserved top margin. No title is
        drawn when omitted.
    xlabel : str, default="UMAP1"
        Label displayed on the horizontal axis.
    ylabel : str, default="UMAP2"
        Label displayed on the vertical axis.
    axis_label_size : float, default=10
        Font size of both axis labels in points.
    title_size : float, default=11
        Font size of the overall figure title in points.
    split_title_size : float, default=10
        Font size of category titles shown above split panels.
    tick_label_size : float, default=8
        Font size of axis tick labels in points.
    spine_width : float, default=1.0
        Line width of the visible left and bottom axes spines.
    tick_width : float, default=1.0
        Line width of major tick marks.
    tick_length : float, default=3
        Length of major tick marks in points.
    tick_nbins : int, default=4
        Maximum approximate number of integer major tick intervals per axis.
    padding_fraction : float, default=0.03
        Fraction of each coordinate range added to both sides of its axis. A
        fixed padding of one is used when a coordinate range is zero.
    shared_limits : bool, default=True
        Whether every split panel uses limits calculated from all cells. When
        ``False``, each panel derives limits from the cells it displays.
    rasterized : bool, default=False
        Whether to rasterize each panel's scatter collection in vector output.
    transparent : bool, default=True
        Whether saved figures use a transparent background.
    save : str or pathlib.Path, optional
        Output filename passed to ``Figure.savefig``. The figure is not written
        when omitted.
    dpi : float, default=600
        Resolution passed to ``Figure.savefig`` when ``save`` is provided.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Newly created figure, sized from the requested centimetre dimensions.
    axes : numpy.ndarray
        One-dimensional object array containing the UMAP axes in split-category
        order. An unsplit plot returns an array containing one Axes object.

    Raises
    ------
    KeyError
        If ``group_col`` or ``split_by`` is absent from ``adata.obs``,
        ``umap_key`` is absent from ``adata.obsm``, or a supplied palette lacks
        an observed group.
    ValueError
        If the UMAP array is not two-dimensional with at least two columns, its
        row count differs from ``adata.n_obs``, no non-missing groups exist, or
        stored AnnData colours do not match all categorical levels. Also raised
        when ``top_group`` is not an observed group, a requested split category
        is unavailable, or no split panels can be created.

    Examples
    --------
    >>> fig, axes = plot_anndata_group_umap(
    ...     adata,
    ...     group_col="cell_type",
    ...     split_by="sample",
    ...     rasterized=True,
    ...     save="cell_type_umap.pdf",
    ... )

    Notes
    -----
    In split mode only the far-left panel retains its y-axis. Category and split
    order follow pandas categorical order. ``top_group`` changes drawing order
    without changing the legend order.
    """
    # -------------------------------------------------------------------------
    # Validate the AnnData inputs and UMAP coordinates.
    # -------------------------------------------------------------------------
    if group_col not in adata.obs.columns:
        raise KeyError(f"'{group_col}' was not found in adata.obs.")

    if split_by is not None and split_by not in adata.obs.columns:
        raise KeyError(f"'{split_by}' was not found in adata.obs.")

    if umap_key not in adata.obsm:
        raise KeyError(f"'{umap_key}' was not found in adata.obsm.")

    umap_xy = np.asarray(adata.obsm[umap_key])

    if umap_xy.ndim != 2 or umap_xy.shape[1] < 2:
        raise ValueError(
            f"adata.obsm['{umap_key}'] must contain at least two columns."
        )

    if umap_xy.shape[0] != adata.n_obs:
        raise ValueError(
            f"adata.obsm['{umap_key}'] has {umap_xy.shape[0]} rows, "
            f"but adata has {adata.n_obs} cells."
        )

    # -------------------------------------------------------------------------
    # Normalize group information and omit unused categorical levels.
    # -------------------------------------------------------------------------
    group_values = adata.obs[group_col].copy()

    if not isinstance(group_values.dtype, pd.CategoricalDtype):
        group_values = group_values.astype("category")

    all_categories = list(group_values.cat.categories)
    groups = [group for group in all_categories if (group_values == group).any()]

    if not groups:
        raise ValueError(f"No non-missing groups found in '{group_col}'.")

    if top_group is not None and top_group not in groups:
        raise ValueError(
            f"top_group={top_group!r} was not found in "
            f"adata.obs['{group_col}']. Available groups: {groups}"
        )

    # -------------------------------------------------------------------------
    # Resolve a complete palette for the groups that are present.
    # -------------------------------------------------------------------------
    color_key = f"{group_col}_colors"

    if palette is None:
        if color_key in adata.uns:
            stored_colors = list(adata.uns[color_key])

            if len(stored_colors) != len(all_categories):
                raise ValueError(
                    f"adata.obs['{group_col}'] has {len(all_categories)} "
                    f"categories, but adata.uns['{color_key}'] has "
                    f"{len(stored_colors)} colours."
                )

            full_palette = dict(zip(all_categories, stored_colors))
            palette = {group: full_palette[group] for group in groups}
        else:
            cmap = plt.get_cmap("tab20")
            palette = {group: cmap(i % 20) for i, group in enumerate(groups)}
    else:
        missing = [group for group in groups if group not in palette]
        if missing:
            raise KeyError(
                "No colour supplied for: " + ", ".join(map(str, missing))
            )

    # -------------------------------------------------------------------------
    # Resolve the ordered set of panels.
    # -------------------------------------------------------------------------
    if split_by is None:
        split_levels = [None]
    else:
        split_values = adata.obs[split_by].copy()
        if not isinstance(split_values.dtype, pd.CategoricalDtype):
            split_values = split_values.astype("category")

        available_split_levels = [
            level
            for level in split_values.cat.categories
            if (split_values == level).any()
        ]

        if split_categories is None:
            split_levels = available_split_levels
        else:
            requested_split_levels = list(split_categories)
            missing_split = [
                level
                for level in requested_split_levels
                if level not in available_split_levels
            ]
            if missing_split:
                raise ValueError(
                    "These split categories were not found: "
                    + ", ".join(map(str, missing_split))
                )
            split_levels = requested_split_levels

    if not split_levels:
        raise ValueError("No split categories available to plot.")

    # -------------------------------------------------------------------------
    # Calculate shared limits and exact multi-panel figure dimensions.
    # -------------------------------------------------------------------------
    x_all = umap_xy[:, 0]
    y_all = umap_xy[:, 1]
    x_range = np.ptp(x_all)
    y_range = np.ptp(y_all)
    x_pad = padding_fraction * x_range if x_range > 0 else 1
    y_pad = padding_fraction * y_range if y_range > 0 else 1
    global_xlim = (np.nanmin(x_all) - x_pad, np.nanmax(x_all) + x_pad)
    global_ylim = (np.nanmin(y_all) - y_pad, np.nanmax(y_all) + y_pad)

    n_panels = len(split_levels)
    plot_area_width_cm = n_panels * width_cm + (n_panels - 1) * panel_gap_cm
    legend_extra_cm = legend_gap_cm + legend_width_cm if show_legend else 0
    figure_width_cm = (
        left_cm + plot_area_width_cm + right_cm + legend_extra_cm
    )
    figure_height_cm = bottom_cm + height_cm + top_cm
    cm_to_inch = 1 / 2.54
    fig = plt.figure(
        figsize=(figure_width_cm * cm_to_inch, figure_height_cm * cm_to_inch)
    )

    axes: list[Axes] = []

    # -------------------------------------------------------------------------
    # Build and format each UMAP panel.
    # -------------------------------------------------------------------------
    for panel_i, split_level in enumerate(split_levels):
        panel_left_cm = left_cm + panel_i * (width_cm + panel_gap_cm)
        ax = fig.add_axes(
            [
                panel_left_cm / figure_width_cm,
                bottom_cm / figure_height_cm,
                width_cm / figure_width_cm,
                height_cm / figure_height_cm,
            ]
        )
        axes.append(ax)

        # Highlight mode keeps every cell and colours only the selected group.
        if split_by == group_col:
            x = x_all
            y = y_all
            panel_groups = group_values.copy()
            colors = np.full(adata.n_obs, na_color, dtype=object)
            selected_mask = (panel_groups == split_level).to_numpy()
            colors[selected_mask] = palette[split_level]
            order_rank = np.zeros(adata.n_obs, dtype=int)
            order_rank[selected_mask] = 1
        else:
            # Normal mode uses every cell when unsplit or a split-level subset.
            if split_by is None:
                panel_mask = np.ones(adata.n_obs, dtype=bool)
            else:
                panel_mask = (adata.obs[split_by] == split_level).to_numpy()

            panel_idx = np.where(panel_mask)[0]
            x = x_all[panel_idx]
            y = y_all[panel_idx]
            panel_groups = group_values.iloc[panel_idx]
            colors = np.full(len(panel_idx), na_color, dtype=object)
            order_rank = np.zeros(len(panel_idx), dtype=int)

            for group_i, group in enumerate(groups, start=1):
                mask = (panel_groups == group).to_numpy()
                colors[mask] = palette[group]
                order_rank[mask] = group_i

            if top_group is not None:
                top_mask = (panel_groups == top_group).to_numpy()
                order_rank[top_mask] = len(groups) + 1

        # Stable sorting preserves original cell order within each category.
        order = np.argsort(order_rank, kind="stable")
        points = ax.scatter(
            x[order],
            y[order],
            s=point_size,
            c=colors[order],
            alpha=point_alpha,
            linewidths=0,
            edgecolors="none",
            rasterized=rasterized,
            zorder=2,
        )
        points.set_gid(
            "all_umap_dots"
            if split_by is None
            else f"umap_dots_{split_level}"
        )

        if shared_limits:
            ax.set_xlim(global_xlim)
            ax.set_ylim(global_ylim)
        else:
            local_x_range = np.ptp(x)
            local_y_range = np.ptp(y)
            local_x_pad = (
                padding_fraction * local_x_range if local_x_range > 0 else 1
            )
            local_y_pad = (
                padding_fraction * local_y_range if local_y_range > 0 else 1
            )
            ax.set_xlim(np.nanmin(x) - local_x_pad, np.nanmax(x) + local_x_pad)
            ax.set_ylim(np.nanmin(y) - local_y_pad, np.nanmax(y) + local_y_pad)

        ax.set_box_aspect(height_cm / width_cm)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=tick_nbins, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=tick_nbins, integer=True))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(spine_width)
        ax.spines["bottom"].set_linewidth(spine_width)
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=tick_label_size,
            width=tick_width,
            length=tick_length,
            direction="out",
        )
        ax.grid(False)
        ax.set_xlabel(xlabel, fontsize=axis_label_size)

        # Only the far-left split panel retains the y-axis and its label.
        if panel_i == 0:
            ax.set_ylabel(ylabel, fontsize=axis_label_size)
        else:
            ax.set_ylabel("")
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)

        if split_by is not None:
            ax.set_title(
                str(split_level),
                fontsize=split_title_size,
                fontweight="normal",
                pad=7,
            )

    axes_array = np.asarray(axes, dtype=object)

    # Place an overall title above the centre of the combined plotting area.
    if title is not None:
        plot_center_cm = left_cm + plot_area_width_cm / 2
        title_y_cm = bottom_cm + height_cm + top_cm * 0.55
        fig.text(
            plot_center_cm / figure_width_cm,
            title_y_cm / figure_height_cm,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
        )

    # Create one vector legend shared by all panels.
    if show_legend:
        if legend_title is None:
            legend_title = group_col

        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor=palette[group],
                markeredgecolor="none",
                alpha=1.0,
                markersize=legend_marker_size,
                label=str(group),
            )
            for group in groups
        ]

        if group_values.isna().any():
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=na_color,
                    markeredgecolor="none",
                    alpha=1.0,
                    markersize=legend_marker_size,
                    label="NA",
                )
            )

        legend_left_cm = (
            left_cm + plot_area_width_cm + right_cm + legend_gap_cm
        )
        fig.legend(
            handles=legend_handles,
            title=legend_title,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(legend_left_cm / figure_width_cm, 0.5),
            bbox_transform=fig.transFigure,
            borderaxespad=0,
            fontsize=legend_fontsize,
            title_fontsize=legend_title_size,
            ncol=legend_ncol,
            columnspacing=legend_columnspacing,
            handletextpad=legend_handletextpad,
        )

    if save is not None:
        fig.savefig(save, dpi=dpi, transparent=transparent)

    return fig, axes_array


def plot_adata_stacked_bar(
    adata: AnnData,
    x_col: str,
    stack_col: str,
    class_order: Iterable[Any] | None = None,
    colors: dict[Any, Any] | None = None,
    x_order: Iterable[Any] | None = None,
    exclude: Iterable[Any] | None = None,
    renormalize: bool = True,
    width: float = 0.68,
    figsize: tuple[float, float] = (4.8, 5.2),
    ylabel: str = "Cell proportion (%)",
    xlabel: str | None = None,
    title: str | None = None,
    legend_title: str | None = None,
    show_labels: bool = True,
    label_min_pct: float = 8,
    label_fontsize: float = 10,
    xtick_rotation: float = 0,
    xtick_fontsize: float = 12,
    ytick_fontsize: float = 11,
    legend_fontsize: float = 11,
    legend_title_fontsize: float = 11,
    edgecolor: Any = "white",
    linewidth: float = 1.2,
    save: str | Path | None = None,
    dpi: float = 600,
    transparent: bool = True,
    return_table: bool = True,
) -> tuple[Figure, Axes, pd.DataFrame] | tuple[Figure, Axes]:
    """Plot cell-composition percentages from AnnData metadata as stacked bars.

    The function tabulates two columns from ``adata.obs`` without modifying the
    AnnData object. Missing metadata rows are removed before percentages are
    calculated. Excluded stack categories may be removed either before or after
    normalization, depending on ``renormalize``.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing the grouping metadata in ``adata.obs``.
    x_col : str
        Observation column defining individual bars, for example ``"sample"``,
        ``"Condition_FatType"``, or ``"Study"``.
    stack_col : str
        Observation column defining the categories stacked within each bar, for
        example ``"MHCII_group"`` or ``"CellType_Broad"``.
    class_order : iterable, optional
        Stack-category order from the bottom to the top of each bar. Categories
        absent from the data are retained as zero-height segments. By default,
        the crosstab column order is used.
    colors : dict, optional
        Mapping from stack categories to Matplotlib-compatible colours. Missing
        mappings are filled from ``tab20``. The supplied mapping is copied and
        is therefore not modified. If omitted, all colours come from ``tab20``.
    x_order : iterable, optional
        Desired left-to-right order of bars. Requested values absent from the
        percentage table are ignored. By default, crosstab index order is used.
    exclude : iterable, optional
        Stack categories to omit, such as ``["Ambiguous"]``.
    renormalize : bool, default=True
        If ``True``, excluded cells are removed before calculating percentages,
        so each non-empty bar sums to 100%. If ``False``, percentages are first
        calculated from all cells and excluded segments are then hidden, so a
        displayed bar may sum to less than 100%.
    width : float, default=0.68
        Width of each bar in Matplotlib x-axis units.
    figsize : tuple of float, default=(4.8, 5.2)
        Figure width and height in inches.
    ylabel : str, default="Cell proportion (%)"
        Label displayed on the vertical axis.
    xlabel : str, optional
        Label displayed on the horizontal axis. No label is added when omitted.
    title : str, optional
        Plot title. No title is added when omitted.
    legend_title : str, optional
        Legend heading. Defaults to ``stack_col`` with underscores replaced by
        spaces.
    show_labels : bool, default=True
        Whether sufficiently large bar segments display internal percentages.
    label_min_pct : float, default=8
        Minimum segment percentage required for an internal text label.
    label_fontsize : float, default=10
        Font size of internal percentage labels in points.
    xtick_rotation : float, default=0
        Rotation of horizontal-axis labels in degrees. Rotated labels are
        right-aligned; unrotated labels are centred.
    xtick_fontsize : float, default=12
        Font size of horizontal-axis tick labels in points.
    ytick_fontsize : float, default=11
        Font size of vertical-axis tick labels in points.
    legend_fontsize : float, default=11
        Font size of legend entries in points.
    legend_title_fontsize : float, default=11
        Font size of the legend title in points.
    edgecolor : color, default="white"
        Matplotlib-compatible colour used for segment borders.
    linewidth : float, default=1.2
        Width of segment borders in points.
    save : str or pathlib.Path, optional
        Output filename passed to ``Figure.savefig``. Nothing is written when
        omitted.
    dpi : float, default=600
        Resolution used when saving the figure.
    transparent : bool, default=True
        Whether a saved figure uses a transparent background.
    return_table : bool, default=True
        Whether to return the percentage table together with the figure and
        axes.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the stacked percentage bar plot.
    ax : matplotlib.axes.Axes
        Axes containing the stacked bars and legend.
    pct_df : pandas.DataFrame
        Percentage table with bars as rows and stack categories as columns.
        Returned only when ``return_table=True``.

    Raises
    ------
    ValueError
        If ``x_col`` or ``stack_col`` is absent from ``adata.obs``, or if no
        complete observations remain after removing missing and excluded data.

    Examples
    --------
    >>> fig, ax, percentages = plot_adata_stacked_bar(
    ...     adata,
    ...     x_col="sample",
    ...     stack_col="MHCII_group",
    ...     colors=MHCII_GROUP_COLORS,
    ...     exclude=["Ambiguous"],
    ... )

    Notes
    -----
    ``plt.show()`` is called before returning. The source AnnData object and a
    caller-supplied ``colors`` mapping are not modified.
    """
    # -------------------------------------------------------------------------
    # Validate and copy the required observation metadata.
    # -------------------------------------------------------------------------
    for column in (x_col, stack_col):
        if column not in adata.obs.columns:
            raise ValueError(
                f"'{column}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    metadata = adata.obs[[x_col, stack_col]].copy()
    metadata = metadata.dropna(subset=[x_col, stack_col])
    metadata[x_col] = metadata[x_col].astype(str)
    metadata[stack_col] = metadata[stack_col].astype(str)

    if metadata.empty:
        raise ValueError(
            f"No complete observations remain for '{x_col}' and '{stack_col}'."
        )

    excluded_classes = [] if exclude is None else [str(x) for x in exclude]

    # -------------------------------------------------------------------------
    # Calculate percentages before or after category exclusion as requested.
    # -------------------------------------------------------------------------
    if renormalize:
        if excluded_classes:
            metadata = metadata[~metadata[stack_col].isin(excluded_classes)]
        if metadata.empty:
            raise ValueError("No observations remain after excluding categories.")
        pct_df = pd.crosstab(
            metadata[x_col], metadata[stack_col], normalize="index"
        ) * 100
    else:
        pct_df = pd.crosstab(
            metadata[x_col], metadata[stack_col], normalize="index"
        ) * 100
        if excluded_classes:
            pct_df = pct_df.drop(columns=excluded_classes, errors="ignore")

    # Use explicit stack and x-axis orderings when the caller supplies them.
    if class_order is None:
        ordered_classes = list(pct_df.columns)
    else:
        ordered_classes = [str(category) for category in class_order]
    pct_df = pct_df.reindex(columns=ordered_classes, fill_value=0)

    if x_order is not None:
        requested_x_order = [str(value) for value in x_order]
        existing_x_order = [
            value for value in requested_x_order if value in pct_df.index
        ]
        pct_df = pct_df.reindex(existing_x_order)

    if pct_df.empty:
        raise ValueError("No bars remain after applying the requested ordering.")

    # -------------------------------------------------------------------------
    # Resolve a complete colour mapping without mutating caller-owned input.
    # -------------------------------------------------------------------------
    if colors is None:
        cmap = plt.get_cmap("tab20")
        resolved_colors = {
            category: cmap(i % 20)
            for i, category in enumerate(ordered_classes)
        }
    else:
        resolved_colors = dict(colors)
        missing_colors = [
            category
            for category in ordered_classes
            if category not in resolved_colors
        ]
        cmap = plt.get_cmap("tab20")
        for i, category in enumerate(missing_colors):
            resolved_colors[category] = cmap(i % 20)

    # -------------------------------------------------------------------------
    # Draw each category from the bottom to the top of every stacked bar.
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    x_positions = np.arange(len(pct_df))
    bottom = np.zeros(len(pct_df))

    for state in ordered_classes:
        values = pct_df[state].to_numpy(dtype=float)
        ax.bar(
            x_positions,
            values,
            bottom=bottom,
            width=width,
            color=resolved_colors[state],
            edgecolor=edgecolor,
            linewidth=linewidth,
            label=str(state).replace("_", "-"),
        )

        if show_labels:
            for bar_i, value in enumerate(values):
                if value >= label_min_pct:
                    ax.text(
                        x_positions[bar_i],
                        bottom[bar_i] + value / 2,
                        f"{value:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=label_fontsize,
                        fontweight="bold",
                        color="black",
                    )

        bottom += values

    # -------------------------------------------------------------------------
    # Format axes, title, and the external legend.
    # -------------------------------------------------------------------------
    ax.set_xticks(x_positions)
    tick_alignment = "right" if xtick_rotation != 0 else "center"
    ax.set_xticklabels(
        pct_df.index.astype(str),
        fontsize=xtick_fontsize,
        fontweight="bold",
        rotation=xtick_rotation,
        ha=tick_alignment,
    )
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=13, fontweight="bold")

    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax.tick_params(axis="y", labelsize=ytick_fontsize, width=1, length=4)

    if title is not None:
        ax.set_title(title, fontsize=14, fontweight="bold", pad=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.grid(False)

    if legend_title is None:
        legend_title = stack_col.replace("_", " ")

    legend = ax.legend(
        title=legend_title,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        handlelength=1.2,
        labelspacing=0.7,
    )
    legend.get_title().set_fontweight("bold")
    plt.tight_layout()

    if save is not None:
        fig.savefig(
            save,
            dpi=dpi,
            bbox_inches="tight",
            transparent=transparent,
        )

    plt.show()

    if return_table:
        return fig, ax, pct_df
    return fig, ax


def plot_seurat_violins(
    adata: AnnData,
    genes: str | Iterable[str],
    groupby: str,
    order: Iterable[Any] | None = None,
    palette: dict[Any, Any] | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    ncols: int = 3,
    show_points: bool = True,
    point_size: float = 1.6,
    point_alpha: float = 0.65,
    point_marker: str = "D",
    point_color: Any = "black",
    jitter: float = 0.28,
    show_box: bool = True,
    box_width: float = 0.12,
    show_stats: bool = True,
    test: str = "mannwhitney",
    adjust_p: bool = True,
    show_pvalue: bool = False,
    violin_linewidth: float = 1.0,
    title_fontsize: float = 17,
    axis_label_fontsize: float = 14,
    tick_fontsize: float = 14,
    ylabel: str = "Expression Level",
    xlabel: str = "Identity",
    show_legend: bool = True,
    legend_title: str | None = None,
    legend_fontsize: float = 12,
    legend_title_fontsize: float = 12,
    panel_width: float = 4.2,
    panel_height: float = 4.5,
    save: str | Path | None = None,
    dpi: float = 300,
    transparent: bool = False,
    return_data: bool = False,
) -> (
    tuple[pd.DataFrame, Figure, np.ndarray]
    | tuple[pd.DataFrame, Figure, np.ndarray, pd.DataFrame]
):
    """Draw Seurat-style gene-expression violin plots from an AnnData object.

    Each gene panel combines a violin, optional individual-cell points, an
    optional central box plot, and an optional comparison between the first two
    groups in ``order``. Mann–Whitney p-values can be adjusted across all genes
    with the Benjamini–Hochberg false-discovery-rate procedure.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing expression values and grouping metadata.
    genes : str or iterable of str
        One gene or an ordered collection of genes to plot. Genes unavailable
        from the selected expression source are reported and skipped.
    groupby : str
        Column in ``adata.obs`` used to group cells along the horizontal axis.
    order : iterable, optional
        Group order for plotting. When statistics are enabled, only the first
        two groups are compared. By default, categorical order is preserved;
        otherwise first-observed order is used.
    palette : dict, optional
        Mapping from group labels to Matplotlib-compatible colours. Missing
        mappings are filled from Seaborn ``Set2``. The supplied mapping is
        copied and is not modified.
    layer : str, optional
        AnnData layer passed to :func:`scanpy.get.obs_df`. By default,
        expression is read from ``adata.X`` unless ``use_raw=True``.
    use_raw : bool, default=False
        Whether to obtain expression from ``adata.raw``.
    ncols : int, default=3
        Maximum number of gene panels per figure row. Must be positive.
    show_points : bool, default=True
        Whether to overlay individual-cell points.
    point_size : float, default=1.6
        Size of individual-cell markers in Seaborn units.
    point_alpha : float, default=0.65
        Opacity of individual-cell markers.
    point_marker : str, default="D"
        Matplotlib marker code used for individual cells.
    point_color : color, default="black"
        Matplotlib-compatible colour used for individual-cell markers.
    jitter : float, default=0.28
        Horizontal jitter applied to individual-cell markers.
    show_box : bool, default=True
        Whether to overlay a white box plot with whiskers and median.
    box_width : float, default=0.12
        Width of the central box plot in categorical-axis units.
    show_stats : bool, default=True
        Whether to compare ``order[0]`` and ``order[1]`` for every plotted gene
        and draw a significance bracket.
    test : str, default="mannwhitney"
        Statistical test name. ``"mannwhitney"``, ``"mann-whitney"``, and
        ``"mw"`` select the same two-sided Mann–Whitney U test.
    adjust_p : bool, default=True
        Whether to apply Benjamini–Hochberg FDR adjustment across plotted genes.
        When ``False``, the unadjusted p-value is copied to the ``FDR`` column.
    show_pvalue : bool, default=False
        Whether brackets show numeric ``P``/``FDR`` values instead of
        significance symbols.
    violin_linewidth : float, default=1.0
        Width of violin boundary lines.
    title_fontsize : float, default=17
        Font size of gene titles in points.
    axis_label_fontsize : float, default=14
        Font size of horizontal and vertical axis labels in points.
    tick_fontsize : float, default=14
        Font size of axis tick labels in points.
    ylabel : str, default="Expression Level"
        Vertical-axis label shown on the left-most panel of each row.
    xlabel : str, default="Identity"
        Horizontal-axis label shown on every used panel.
    show_legend : bool, default=True
        Whether to add a shared figure-level group legend.
    legend_title : str, optional
        Shared legend title. No title is displayed when omitted.
    legend_fontsize : float, default=12
        Font size of shared legend entries in points.
    legend_title_fontsize : float, default=12
        Font size of the shared legend title in points.
    panel_width : float, default=4.2
        Width of each gene panel in inches.
    panel_height : float, default=4.5
        Height of each gene panel in inches.
    save : str or pathlib.Path, optional
        Output filename passed to ``Figure.savefig``. Nothing is written when
        omitted.
    dpi : float, default=300
        Resolution used when saving the figure.
    transparent : bool, default=False
        Whether a saved figure uses a transparent background.
    return_data : bool, default=False
        Whether to append the extracted long-form expression table to the
        returned values.

    Returns
    -------
    stats_df : pandas.DataFrame
        One row per plotted gene with compared groups, cell counts, means,
        medians, Mann–Whitney statistic, raw p-value, and FDR. Empty when
        ``show_stats=False``.
    fig : matplotlib.figure.Figure
        Figure containing the gene-expression panels and optional shared legend.
    axes : numpy.ndarray
        Flattened object array containing all allocated axes. Unused axes are
        hidden when the panel grid is larger than the number of genes.
    expression_df : pandas.DataFrame
        Extracted expression values and ordered group labels. Returned only when
        ``return_data=True``.

    Raises
    ------
    ValueError
        If ``groupby`` is unavailable, no genes are supplied or found, raw data
        are requested but absent, fewer than two groups are available for
        enabled statistics, ``ncols`` is invalid, no selected-group cells remain,
        or ``test`` is unsupported.
    ImportError
        If Scanpy, Seaborn, SciPy statistics, or statsmodels is unavailable in
        the active Python environment.

    Examples
    --------
    >>> stats, fig, axes = plot_seurat_violins(
    ...     adata,
    ...     genes=["Gata3", "Il5"],
    ...     groupby="genotype",
    ...     order=["WT", "cKO"],
    ...     palette=WT_CKO_COLORS,
    ... )

    Notes
    -----
    The statistical comparison treats cells as independent observations. It is
    suitable for exploratory visualization but does not account for biological
    replication or donor/sample structure and should not replace replicate-aware
    pseudobulk or mixed-model differential-expression analysis. The function
    does not call ``plt.show()``; display remains under caller control.
    """
    # -------------------------------------------------------------------------
    # Validate grouping, genes, expression source, and panel layout.
    # -------------------------------------------------------------------------
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"'{groupby}' is not present in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    requested_genes = [genes] if isinstance(genes, str) else list(genes)
    if not requested_genes:
        raise ValueError("No genes were provided.")

    if isinstance(ncols, bool) or not isinstance(ncols, Integral) or ncols <= 0:
        raise ValueError("ncols must be a positive integer")

    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        available_genes = set(adata.raw.var_names.astype(str))
    else:
        available_genes = set(adata.var_names.astype(str))

    missing_genes = [
        gene for gene in requested_genes if gene not in available_genes
    ]
    if missing_genes:
        print("Warning: the following genes were not found and will be skipped:")
        print(missing_genes)

    selected_genes = [
        gene for gene in requested_genes if gene in available_genes
    ]
    if not selected_genes:
        raise ValueError("None of the requested genes were found.")

    # Preserve categorical order unless the caller provides an explicit order.
    if order is None:
        if isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
            group_order = list(adata.obs[groupby].cat.categories)
        else:
            group_order = list(pd.unique(adata.obs[groupby].dropna()))
    else:
        group_order = list(order)
    group_order = [str(group) for group in group_order]

    if not group_order:
        raise ValueError("No groups were found.")

    if show_stats and len(group_order) < 2:
        raise ValueError("At least two groups are required for statistical testing.")

    supported_tests = {"mannwhitney", "mann-whitney", "mw"}
    if show_stats and test.lower() not in supported_tests:
        raise ValueError(
            f"Unsupported test: {test}. Currently supported: 'mannwhitney'."
        )

    # Keep optional visualization/statistics dependencies local to this helper.
    import scanpy as sc
    import seaborn as sns
    from matplotlib.patches import Patch
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Resolve colours without mutating a palette owned by the caller.
    if palette is None:
        default_colors = sns.color_palette("Set2", n_colors=len(group_order))
        resolved_palette = dict(zip(group_order, default_colors))
    else:
        resolved_palette = {str(key): value for key, value in palette.items()}
        missing_colors = [
            group for group in group_order if group not in resolved_palette
        ]
        extra_colors = sns.color_palette("Set2", n_colors=len(missing_colors))
        resolved_palette.update(zip(missing_colors, extra_colors))

    expression_df = sc.get.obs_df(
        adata,
        keys=[groupby] + selected_genes,
        layer=layer,
        use_raw=use_raw,
    )
    expression_df[groupby] = expression_df[groupby].astype(str)
    expression_df = expression_df[
        expression_df[groupby].isin(group_order)
    ].copy()

    if expression_df.empty:
        raise ValueError("No cells remain after applying the requested group order.")

    expression_df[groupby] = pd.Categorical(
        expression_df[groupby], categories=group_order, ordered=True
    )

    # -------------------------------------------------------------------------
    # Compare the first two ordered groups for every selected gene.
    # -------------------------------------------------------------------------
    stat_results: list[dict[str, Any]] = []
    if show_stats:
        group1, group2 = group_order[:2]

        for gene in selected_genes:
            values1 = expression_df.loc[
                expression_df[groupby] == group1, gene
            ].dropna()
            values2 = expression_df.loc[
                expression_df[groupby] == group2, gene
            ].dropna()

            if len(values1) > 0 and len(values2) > 0:
                result = mannwhitneyu(values1, values2, alternative="two-sided")
                pvalue = result.pvalue
                statistic = result.statistic
            else:
                pvalue = np.nan
                statistic = np.nan

            stat_results.append(
                {
                    "Gene": gene,
                    "Group1": group1,
                    "Group2": group2,
                    "n_Group1": len(values1),
                    "n_Group2": len(values2),
                    "mean_Group1": values1.mean() if len(values1) > 0 else np.nan,
                    "mean_Group2": values2.mean() if len(values2) > 0 else np.nan,
                    "median_Group1": (
                        values1.median() if len(values1) > 0 else np.nan
                    ),
                    "median_Group2": (
                        values2.median() if len(values2) > 0 else np.nan
                    ),
                    "statistic": statistic,
                    "pvalue": pvalue,
                }
            )

        stats_df = pd.DataFrame(stat_results)
        stats_df["FDR"] = np.nan
        valid_pvalues = stats_df["pvalue"].notna()
        if valid_pvalues.any():
            if adjust_p:
                stats_df.loc[valid_pvalues, "FDR"] = multipletests(
                    stats_df.loc[valid_pvalues, "pvalue"], method="fdr_bh"
                )[1]
            else:
                stats_df.loc[valid_pvalues, "FDR"] = stats_df.loc[
                    valid_pvalues, "pvalue"
                ]
    else:
        stats_df = pd.DataFrame()

    def _pvalue_label(pvalue: float) -> str:
        """Format one raw or adjusted p-value for a plot annotation."""
        if pd.isna(pvalue):
            return "NA"
        if show_pvalue:
            label = "FDR" if adjust_p else "P"
            return f"{label} = {pvalue:.2g}"
        if pvalue < 0.0001:
            return "****"
        if pvalue < 0.001:
            return "***"
        if pvalue < 0.01:
            return "**"
        if pvalue < 0.05:
            return "*"
        return "ns"

    # -------------------------------------------------------------------------
    # Create the panel grid and render each selected gene.
    # -------------------------------------------------------------------------
    ncols_use = min(ncols, len(selected_genes))
    nrows = math.ceil(len(selected_genes) / ncols_use)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols_use,
        figsize=(panel_width * ncols_use, panel_height * nrows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    for gene_i, gene in enumerate(selected_genes):
        ax = axes_flat[gene_i]
        gene_df = expression_df[[groupby, gene]].copy()
        gene_df.columns = [groupby, "Expression"]
        gene_df = gene_df.dropna()

        sns.violinplot(
            data=gene_df,
            x=groupby,
            y="Expression",
            order=group_order,
            hue=groupby,
            hue_order=group_order,
            palette=resolved_palette,
            legend=False,
            inner=None,
            cut=0,
            linewidth=violin_linewidth,
            saturation=1,
            ax=ax,
        )

        if show_points:
            sns.stripplot(
                data=gene_df,
                x=groupby,
                y="Expression",
                order=group_order,
                color=point_color,
                marker=point_marker,
                size=point_size,
                alpha=point_alpha,
                jitter=jitter,
                ax=ax,
                zorder=2,
            )

        if show_box:
            sns.boxplot(
                data=gene_df,
                x=groupby,
                y="Expression",
                order=group_order,
                width=box_width,
                showfliers=False,
                showcaps=True,
                boxprops={
                    "facecolor": "white",
                    "edgecolor": "black",
                    "linewidth": 1.2,
                    "zorder": 3,
                },
                whiskerprops={"color": "black", "linewidth": 1.2},
                capprops={"color": "black", "linewidth": 1.2},
                medianprops={"color": "black", "linewidth": 1.8},
                ax=ax,
                zorder=3,
            )

        ymin = gene_df["Expression"].min()
        ymax = gene_df["Expression"].max()
        yrange = ymax - ymin
        if not np.isfinite(yrange) or yrange == 0:
            yrange = 1

        if show_stats:
            adjusted_pvalue = stats_df.loc[
                stats_df["Gene"] == gene, "FDR"
            ].iloc[0]
            bracket_y = ymax + 0.07 * yrange
            bracket_height = 0.035 * yrange
            ax.plot(
                [0, 0, 1, 1],
                [
                    bracket_y,
                    bracket_y + bracket_height,
                    bracket_y + bracket_height,
                    bracket_y,
                ],
                color="black",
                linewidth=1.2,
                clip_on=False,
            )
            ax.text(
                0.5,
                bracket_y + bracket_height,
                _pvalue_label(adjusted_pvalue),
                ha="center",
                va="bottom",
                fontsize=11,
            )
            upper_margin = 0.22
        else:
            upper_margin = 0.08

        ax.set_ylim(ymin - 0.03 * yrange, ymax + upper_margin * yrange)
        ax.set_title(gene, fontsize=title_fontsize, pad=14)
        ax.set_xlabel(xlabel, fontsize=axis_label_fontsize)
        if gene_i % ncols_use == 0:
            ax.set_ylabel(ylabel, fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=tick_fontsize, width=1)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.1)
        ax.spines["bottom"].set_linewidth(1.1)

        if ax.get_legend() is not None:
            ax.get_legend().remove()

    # Hide surplus axes while retaining them in the returned flattened array.
    for unused_i in range(len(selected_genes), len(axes_flat)):
        axes_flat[unused_i].set_visible(False)

    if show_legend:
        handles = [
            Patch(
                facecolor=resolved_palette[group],
                edgecolor="black",
                label=group,
            )
            for group in group_order
        ]
        legend = fig.legend(
            handles=handles,
            title=legend_title,
            loc="center right",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=legend_fontsize,
            title_fontsize=legend_title_fontsize,
        )
        if legend_title is not None:
            legend.get_title().set_fontweight("bold")
        fig.tight_layout(rect=[0, 0, 0.94, 1])
    else:
        fig.tight_layout()

    if save is not None:
        fig.savefig(
            save,
            dpi=dpi,
            bbox_inches="tight",
            transparent=transparent,
            facecolor="none" if transparent else "white",
        )

    if return_data:
        return stats_df, fig, axes_flat, expression_df
    return stats_df, fig, axes_flat
