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
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy import sparse

if TYPE_CHECKING:
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


__all__ = [
    "convert_genes_to_features",
    "ordmag_filter",
    "plot_anndata_group_umap",
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
    palette: dict[Any, Any] | None = None,
    point_size: float = 18,
    point_alpha: float = 0.5,
    na_color: Any = "#D3D3D3",
    umap_key: str = "X_umap",
    top_group: Any | None = None,
    width_cm: float = 4,
    height_cm: float = 4,
    left_cm: float = 1.30,
    right_cm: float = 0.30,
    bottom_cm: float = 1.10,
    top_cm: float = 0.35,
    show_legend: bool = True,
    legend_title: str | None = None,
    legend_marker_size: float = 6,
    legend_fontsize: float = 8,
    legend_title_size: float = 9,
    legend_gap_cm: float = 0.25,
    legend_width_cm: float = 2.5,
    title: str | None = None,
    xlabel: str = "UMAP1",
    ylabel: str = "UMAP2",
    axis_label_size: float = 10,
    title_size: float = 11,
    tick_label_size: float = 8,
    spine_width: float = 1.0,
    tick_width: float = 1.0,
    tick_length: float = 3,
    tick_nbins: int = 4,
    padding_fraction: float = 0.03,
    rasterized: bool = False,
    transparent: bool = True,
    save: str | Path | None = None,
    dpi: float = 600,
) -> tuple[Figure, Axes]:
    """Plot one categorical AnnData annotation on precomputed UMAP coordinates.

    All cells are drawn with one scatter call. When ``rasterized=True``, this
    normally produces one raster layer for the points in PDF or SVG output,
    while axes, labels, title, ticks, and legend remain editable vector objects.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object whose ``obs`` table contains ``group_col`` and whose
        ``obsm`` mapping contains two-dimensional UMAP coordinates.
    group_col : str
        Name of the categorical column in ``adata.obs`` used to colour cells.
        Non-categorical values are converted to a pandas categorical series.
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
        categories in ``group_col``. By default, normal categorical plotting
        order is used without prioritizing one group.
    width_cm : float, default=4
        Physical width of the UMAP plotting box in centimetres.
    height_cm : float, default=4
        Physical height of the UMAP plotting box in centimetres.
    left_cm : float, default=1.30
        Space to the left of the plotting box in centimetres.
    right_cm : float, default=0.30
        Space between the plotting box and optional legend area in centimetres.
    bottom_cm : float, default=1.10
        Space below the plotting box in centimetres.
    top_cm : float, default=0.35
        Space above the plotting box in centimetres.
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
    legend_gap_cm : float, default=0.25
        Horizontal gap before the legend area in centimetres.
    legend_width_cm : float, default=2.5
        Width reserved for the legend in centimetres.
    title : str, optional
        Plot title. No title is drawn when omitted.
    xlabel : str, default="UMAP1"
        Label displayed on the horizontal axis.
    ylabel : str, default="UMAP2"
        Label displayed on the vertical axis.
    axis_label_size : float, default=10
        Font size of both axis labels in points.
    title_size : float, default=11
        Font size of the title in points.
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
    rasterized : bool, default=False
        Whether to rasterize the single scatter collection in vector output.
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
    ax : matplotlib.axes.Axes
        Axes containing the UMAP scatter plot.

    Raises
    ------
    KeyError
        If ``group_col`` is absent from ``adata.obs``, ``umap_key`` is absent
        from ``adata.obsm``, or a supplied palette lacks an observed group.
    ValueError
        If the UMAP array is not two-dimensional with at least two columns, its
        row count differs from ``adata.n_obs``, no non-missing groups exist, or
        stored AnnData colours do not match all categorical levels. Also raised
        when ``top_group`` is not an observed group.

    Examples
    --------
    >>> fig, ax = plot_anndata_group_umap(
    ...     adata,
    ...     group_col="cell_type",
    ...     rasterized=True,
    ...     save="cell_type_umap.pdf",
    ... )

    Notes
    -----
    Cells with missing group values are drawn first. Remaining cells are drawn
    in categorical order so non-missing annotations remain visible above them.
    When supplied, ``top_group`` is drawn last without changing legend order.
    """
    # -------------------------------------------------------------------------
    # Validate the AnnData inputs and UMAP coordinates.
    # -------------------------------------------------------------------------
    if group_col not in adata.obs.columns:
        raise KeyError(f"'{group_col}' was not found in adata.obs.")

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
    # Calculate coordinate limits and the exact physical figure dimensions.
    # -------------------------------------------------------------------------
    x = umap_xy[:, 0]
    y = umap_xy[:, 1]
    x_range = np.ptp(x)
    y_range = np.ptp(y)
    x_pad = padding_fraction * x_range if x_range > 0 else 1
    y_pad = padding_fraction * y_range if y_range > 0 else 1
    xlim = (np.nanmin(x) - x_pad, np.nanmax(x) + x_pad)
    ylim = (np.nanmin(y) - y_pad, np.nanmax(y) + y_pad)

    extra_legend_cm = legend_gap_cm + legend_width_cm if show_legend else 0
    figure_width_cm = left_cm + width_cm + right_cm + extra_legend_cm
    figure_height_cm = bottom_cm + height_cm + top_cm
    cm_to_inch = 1 / 2.54

    fig = plt.figure(
        figsize=(figure_width_cm * cm_to_inch, figure_height_cm * cm_to_inch)
    )
    ax = fig.add_axes(
        [
            left_cm / figure_width_cm,
            bottom_cm / figure_height_cm,
            width_cm / figure_width_cm,
            height_cm / figure_height_cm,
        ]
    )

    # -------------------------------------------------------------------------
    # Build per-cell colours and draw missing values before annotated groups.
    # -------------------------------------------------------------------------
    colors = np.array([na_color] * adata.n_obs, dtype=object)
    order_rank = np.full(adata.n_obs, fill_value=-1, dtype=int)
    order_rank[group_values.isna().to_numpy()] = 0

    for i, group in enumerate(groups, start=1):
        mask = (group_values == group).to_numpy()
        colors[mask] = palette[group]
        order_rank[mask] = i

    # Promote the selected group above every other category while retaining
    # stable cell order within that group.
    if top_group is not None:
        top_mask = (group_values == top_group).to_numpy()
        order_rank[top_mask] = len(groups) + 1

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
    points.set_gid("all_umap_dots")

    # -------------------------------------------------------------------------
    # Format axes, labels, and optional title.
    # -------------------------------------------------------------------------
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
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
    ax.set_ylabel(ylabel, fontsize=axis_label_size)

    if title is not None:
        ax.set_title(title, fontsize=title_size, fontweight="normal", pad=8)

    # -------------------------------------------------------------------------
    # Create a vector legend independently of the scatter collection.
    # -------------------------------------------------------------------------
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

        legend_left_cm = left_cm + width_cm + right_cm + legend_gap_cm
        legend_x = legend_left_cm / figure_width_cm
        fig.legend(
            handles=legend_handles,
            title=legend_title,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(legend_x, 0.5),
            bbox_transform=fig.transFigure,
            borderaxespad=0,
            fontsize=legend_fontsize,
            title_fontsize=legend_title_size,
        )

    # Save only when the caller explicitly supplies an output path.
    if save is not None:
        fig.savefig(save, dpi=dpi, transparent=transparent)

    return fig, ax
