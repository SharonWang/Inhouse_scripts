#' In-house single-cell RNA-seq preprocessing utilities
#'
#' This script contains reusable R functions for single-cell RNA-seq
#' preprocessing. Functions supplied in R will be placed in the appropriate
#' workflow section and documented with Roxygen comments describing all input
#' parameters and return values.
#'
#' Python functions are maintained separately in
#' `Python/scRNAseq_preprocessing.py` and are not translated automatically.
#'
#' @keywords internal
NULL


# =============================================================================
# Reusable palettes
# =============================================================================


#' Macaron colour palette for bulk RNA-seq visualizations
#'
#' A reusable vector of 15 muted pastel colours originally supplied with
#' [plot_bulk_pca()]. Unnamed colours are intentional so plotting functions can
#' assign them in the requested factor-level order.
#'
#' @return A character vector containing 15 hexadecimal colour values.
#'
#' @examples
#' BULK_PCA_MACARON_COLORS
#' BULK_PCA_MACARON_COLORS[seq_len(3)]
BULK_PCA_MACARON_COLORS <- c(
    "#E6A4A8", # strawberry pink
    "#9DB7D5", # blueberry blue
    "#A8C8A0", # pistachio green
    "#E5C07B", # vanilla yellow
    "#B8A1C8", # lavender
    "#8FC7C3", # mint
    "#D7A6C2", # raspberry
    "#C7B299", # caramel
    "#A9B8C6", # blue-grey
    "#D8B4A0", # peach
    "#A5C3D5", # sky blue
    "#C3B6D8", # violet
    "#D8C79F", # biscuit
    "#9FBDB0", # sage
    "#D2A7A0"  # dusty rose
)


#' Macaron colour palette for bulk RNA-seq QC plots
#'
#' A reusable vector of 12 muted colours supplied with [plot_bulk_qc()]. The
#' leading neutral grey is useful for reference or control groups. Colours are
#' intentionally unnamed so functions can assign them in a requested group
#' order.
#'
#' @return A character vector containing 12 hexadecimal colour values.
#'
#' @examples
#' BULK_QC_MACARON_COLORS
#' BULK_QC_MACARON_COLORS[seq_len(3)]
BULK_QC_MACARON_COLORS <- c(
    "#BDBDBD", # grey
    "#9DB7D5", # blueberry
    "#A8C8A0", # pistachio
    "#E6A4A8", # strawberry
    "#E5C07B", # vanilla
    "#B8A1C8", # lavender
    "#8FC7C3", # mint
    "#D7A6C2", # raspberry
    "#C7B299", # caramel
    "#A5C3D5", # sky
    "#D8B4A0", # peach
    "#9FBDB0"  # sage
)


#' Macaron colour palette for bulk RNA-seq violin plots
#'
#' A reusable vector of 12 muted pastel colours supplied with
#' [plot_bulk_violin()]. Colours are intentionally unnamed so the function can
#' assign them after resolving the requested group order.
#'
#' @return A character vector containing 12 hexadecimal colour values.
#'
#' @examples
#' BULK_VIOLIN_MACARON_COLORS
#' BULK_VIOLIN_MACARON_COLORS[seq_len(3)]
BULK_VIOLIN_MACARON_COLORS <- c(
    "#E6A4A8", # pink
    "#9DB7D5", # blue
    "#A8C8A0", # green
    "#E5C07B", # yellow
    "#B8A1C8", # lavender
    "#8FC7C3", # mint
    "#D7A6C2", # raspberry
    "#D8B4A0", # peach
    "#A5C3D5", # sky
    "#9FBDB0", # sage
    "#C7B299", # caramel
    "#C3B6D8"  # violet
)


#' Diverging palette for differential-expression heatmaps
#'
#' A reusable blue-white-red palette supplied with [plot_de_heatmap()]. The
#' low, midpoint, and high colours are intentionally unnamed because they map
#' directly to the symmetric heatmap scale.
#'
#' @return A character vector containing three hexadecimal colour values.
#'
#' @examples
#' DE_HEATMAP_COLORS
DE_HEATMAP_COLORS <- c(
    "#5B78A5", # lower expression
    "#FFFFFF", # midpoint
    "#C76B6B"  # higher expression
)


#' Sequential purple palette for gene-set expression heatmaps
#'
#' A reusable five-colour light-to-dark purple palette supplied with
#' [plot_gene_set_heatmap()]. For row z-scores it is mapped symmetrically from
#' the negative display limit through zero to the positive display limit.
#'
#' @return A character vector containing five hexadecimal colour values.
#'
#' @examples
#' GENE_SET_HEATMAP_COLORS
GENE_SET_HEATMAP_COLORS <- c(
    "#F2F2F2", # near-white
    "#D6D6E7", # pale lavender
    "#A89BD0", # purple midpoint
    "#6A51A3", # deep purple
    "#301934"  # darkest purple
)


#' Macaron comparison palette for signed Manhattan plots
#'
#' A reusable vector of 12 muted colours supplied with
#' [plot_signed_manhattan()]. Colours are intentionally unnamed so they can be
#' assigned after comparison labels and their display order have been resolved.
#'
#' @return A character vector containing 12 hexadecimal colour values.
#'
#' @examples
#' SIGNED_MANHATTAN_MACARON_COLORS
#' SIGNED_MANHATTAN_MACARON_COLORS[seq_len(3)]
SIGNED_MANHATTAN_MACARON_COLORS <- c(
    "#9DB7D5", # blueberry
    "#A8C8A0", # pistachio
    "#E6A4A8", # strawberry
    "#B8A1C8", # lavender
    "#E5C07B", # vanilla
    "#8FC7C3", # mint
    "#D7A6C2", # raspberry
    "#D8B4A0", # peach
    "#A5C3D5", # sky
    "#9FBDB0", # sage
    "#C7B299", # caramel
    "#C3B6D8"  # violet
)


#' Macaron palette for paired bulk-expression plots
#'
#' Six muted colours supplied with [plot_bulk_paired_expression()].
#'
#' @return A character vector containing six hexadecimal colour values.
BULK_PAIRED_EXPRESSION_COLORS <- c(
    "#E8B4B8", "#B8C7E0", "#B9D8C2",
    "#F3C77F", "#D6B8D8", "#AFCFD3"
)


#' Diverging palettes for bulk activity and expression heatmaps
#'
#' `BULK_SIGNALING_DIVERGING_COLORS` stores the supplied blue, white, and
#' muted-red anchors. `BULK_ACTIVITY_DIVERGING_COLORS` expands those anchors to
#' 100 colours for TF, gene-expression, comparison, and signaling heatmaps.
#'
#' @return Character vectors containing three and 100 hexadecimal colours,
#'   respectively.
BULK_SIGNALING_DIVERGING_COLORS <- c("#8199C1", "#FFFFFF", "#D68188")

BULK_ACTIVITY_DIVERGING_COLORS <- grDevices::colorRampPalette(
    BULK_SIGNALING_DIVERGING_COLORS
)(100)


#' Sequential palette for NicheNet ligand-activity heatmaps
#'
#' A 100-colour white-to-purple gradient for corrected AUPR scores.
#'
#' @return A character vector containing 100 hexadecimal colour values.
BULK_NICHENET_ACTIVITY_COLORS <- grDevices::colorRampPalette(
    c("#FFFFFF", "#B4B9D9", "#59348A")
)(100)


#' Sender and connectivity palettes for NicheNet plots
#'
#' Named sender-population colours and binary absent/present link colours used
#' by [plot_bulk_nichenet_lr_connectivity()].
#'
#' @return `BULK_NICHENET_SENDER_COLORS` is a named three-colour vector;
#'   `BULK_NICHENET_CONNECTIVITY_COLORS` is a two-colour vector.
BULK_NICHENET_SENDER_COLORS <- c(
    pATM = "#D9A6A6",
    cATM = "#AFC8B8",
    sATM = "#B8B3D8"
)

BULK_NICHENET_CONNECTIVITY_COLORS <- c("#F5F5F5", "#8199C1")


# =============================================================================
# Data loading and import
# =============================================================================


#' Import a featureCounts project with metadata and assignment QC
#'
#' Read a featureCounts gene-level count table, reconcile its sample columns
#' with study metadata, optionally merge external gene annotations and the
#' featureCounts summary file, and optionally construct an
#' [edgeR::DGEList()] object. Sample order is validated explicitly across the
#' returned count matrix, metadata, and QC table.
#'
#' Although this utility is stored in the shared R pipeline, featureCounts
#' gene-level matrices are normally used for bulk or pseudobulk RNA-seq rather
#' than cell-level single-cell matrices.
#'
#' @param fcounts_file Character scalar. Path to a tab-delimited featureCounts
#'   count file containing a `Geneid` column, standard featureCounts annotation
#'   columns, and one count column per sample.
#' @param summary_file `NULL` or a character scalar. Optional path to the
#'   featureCounts `.summary` file. When supplied, assignment metrics and all
#'   featureCounts status categories are parsed and aligned to the retained
#'   samples.
#' @param gene_feature_file `NULL` or a character scalar. Optional path to a
#'   whitespace-delimited, headerless gene-annotation file. The first three
#'   columns must be gene ID, gene type, and gene name; optional columns four
#'   through six are chromosome, start, and end. Further columns are retained
#'   as `annotation_extra_*`.
#' @param metadata A data frame or a character scalar giving a metadata file.
#'   CSV, tab-delimited text, and Excel (`.xlsx` or `.xls`) files are supported.
#'   Excel input requires the `readxl` package.
#' @param sample_col Character scalar. Name of the metadata column containing
#'   sample identifiers after featureCounts sample names have been cleaned.
#' @param strip_ensembl_version Logical scalar. If `TRUE`, remove a terminal
#'   numeric Ensembl version suffix such as `.5` from count-table and external
#'   annotation gene IDs.
#' @param sample_name_fun `NULL` or a function accepting and returning a
#'   character vector. It is applied to count-table and summary-file sample
#'   names. The default removes directory paths, STAR
#'   `_Aligned.sortedByCoord.out.bam` suffixes, and `.bam`, `.sam`, or `.cram`
#'   extensions.
#' @param strict Logical scalar. If `TRUE`, stop unless count-table and metadata
#'   sample sets match exactly. If `FALSE`, retain their intersection in count
#'   matrix order and report unmatched samples.
#' @param make_dge Logical scalar. Whether to construct an [edgeR::DGEList()]
#'   in the returned `dge` element. If edgeR is unavailable, a warning is
#'   issued and `dge` remains `NULL`.
#' @param add_all_fc_qc_to_metadata Logical scalar. When summary data are
#'   available, whether every `Unassigned_*` featureCounts status column should
#'   also be copied into the aligned metadata.
#' @param verbose Logical scalar. Whether to print progress, matching
#'   diagnostics, and an import summary.
#'
#' @return An object of class `featurecounts_project`, implemented as a named
#'   list with the following elements:
#'   \describe{
#'     \item{counts}{Numeric gene-by-sample count matrix.}
#'     \item{metadata}{Metadata aligned exactly to the count-matrix columns,
#'       with summary QC columns added when available.}
#'     \item{genes}{Gene annotation aligned exactly to the count-matrix rows.}
#'     \item{qc}{Per-sample featureCounts assignment QC table, or `NULL` when a
#'       usable summary was not supplied.}
#'     \item{featurecounts_summary}{Imported featureCounts summary data frame,
#'       or `NULL`.}
#'     \item{dge}{Optional [edgeR::DGEList()] object, or `NULL`.}
#'     \item{original_sample_names}{Original featureCounts sample column names
#'       for retained samples, aligned to `sample_names`.}
#'     \item{sample_names}{Cleaned sample names retained in the final object.}
#'   }
#'
#' @details
#' The function checks for duplicate gene IDs after optional Ensembl-version
#' removal, duplicate or empty cleaned sample names, missing or negative count
#' values, and ordering inconsistencies. Fractional counts are retained but
#' trigger a warning because they may indicate fractional featureCounts
#' assignment. Existing metadata columns with generated QC names may be
#' replaced.
#'
#' This function imports and validates counts but does not perform library-size
#' filtering, normalization, exploratory QC plots, design construction, or
#' differential-expression testing. Those decisions should reflect the study
#' design, biological replicates, strandedness, reference annotation, and
#' featureCounts configuration.
#'
#' @examples
#' \dontrun{
#' project <- read_featurecounts_project(
#'   fcounts_file = "counts/featureCounts.txt",
#'   summary_file = "counts/featureCounts.txt.summary",
#'   metadata = "metadata/samples.csv",
#'   gene_feature_file = "reference/gene_features.txt",
#'   sample_col = "SampleName"
#' )
#'
#' project$counts[, 1:3]
#' project$qc
#' }
#'
#' @export
read_featurecounts_project <- function(
    fcounts_file,
    summary_file = NULL,
    gene_feature_file = NULL,
    metadata,
    sample_col = "SampleName",
    strip_ensembl_version = TRUE,
    sample_name_fun = NULL,
    strict = TRUE,
    make_dge = TRUE,
    add_all_fc_qc_to_metadata = TRUE,
    verbose = TRUE
) {
    # -------------------------------------------------------------------------
    # Validate scalar arguments and input paths
    # -------------------------------------------------------------------------
    logical_arguments <- list(
        strip_ensembl_version = strip_ensembl_version,
        strict = strict,
        make_dge = make_dge,
        add_all_fc_qc_to_metadata = add_all_fc_qc_to_metadata,
        verbose = verbose
    )

    invalid_logical <- vapply(
        logical_arguments,
        function(x) length(x) != 1L || is.na(x) || !is.logical(x),
        logical(1)
    )
    if (any(invalid_logical)) {
        stop(
            "The following arguments must be non-missing logical scalars: ",
            paste(names(logical_arguments)[invalid_logical], collapse = ", ")
        )
    }

    if (!is.character(sample_col) ||
        length(sample_col) != 1L ||
        is.na(sample_col) ||
        sample_col == "") {
        stop("sample_col must be a non-empty character scalar.")
    }

    validate_input_file <- function(path, argument_name) {
        if (!is.character(path) ||
            length(path) != 1L ||
            is.na(path) ||
            path == "") {
            stop(argument_name, " must be a non-empty file path.")
        }
        if (!file.exists(path)) {
            stop(argument_name, " does not exist: ", path)
        }
        if (dir.exists(path)) {
            stop(argument_name, " must refer to a file, not a directory: ", path)
        }
        invisible(path)
    }

    validate_input_file(fcounts_file, "fcounts_file")
    if (!is.null(summary_file)) {
        validate_input_file(summary_file, "summary_file")
    }
    if (!is.null(gene_feature_file)) {
        validate_input_file(gene_feature_file, "gene_feature_file")
    }

    # -------------------------------------------------------------------------
    # Local helpers
    # -------------------------------------------------------------------------
    default_sample_name_fun <- function(x) {
        x <- basename(x)
        x <- sub(
            "_Aligned\\.sortedByCoord\\.out\\.bam$",
            "",
            x,
            ignore.case = TRUE
        )
        x <- sub(
            "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
            "",
            x,
            ignore.case = TRUE
        )
        sub("\\.(bam|sam|cram)$", "", x, ignore.case = TRUE)
    }

    if (is.null(sample_name_fun)) {
        sample_name_fun <- default_sample_name_fun
    }
    if (!is.function(sample_name_fun)) {
        stop("sample_name_fun must be NULL or a function.")
    }

    clean_gene_id <- function(x) {
        x <- as.character(x)
        if (strip_ensembl_version) {
            x <- sub("\\.[0-9]+$", "", x)
        }
        x
    }

    clean_sample_names <- function(x, source_name) {
        cleaned <- sample_name_fun(x)
        if (length(cleaned) != length(x)) {
            stop("sample_name_fun changed the number of ", source_name, " names.")
        }
        cleaned <- as.character(cleaned)
        if (anyNA(cleaned) || any(cleaned == "")) {
            stop("Cleaning produced empty or NA ", source_name, " sample names.")
        }
        if (anyDuplicated(cleaned)) {
            duplicates <- unique(cleaned[duplicated(cleaned)])
            stop(
                "Cleaning created duplicate ", source_name, " sample names:\n",
                paste(duplicates, collapse = "\n")
            )
        }
        cleaned
    }

    read_metadata <- function(x) {
        if (is.data.frame(x)) {
            return(as.data.frame(
                x,
                stringsAsFactors = FALSE,
                check.names = FALSE
            ))
        }

        validate_input_file(x, "metadata")
        extension <- tolower(tools::file_ext(x))
        if (extension == "csv") {
            return(read.csv(
                x,
                stringsAsFactors = FALSE,
                check.names = FALSE
            ))
        }
        if (extension %in% c("xlsx", "xls")) {
            if (!requireNamespace("readxl", quietly = TRUE)) {
                stop(
                    "Package 'readxl' is required for Excel metadata.\n",
                    "Install with: install.packages('readxl')"
                )
            }
            return(as.data.frame(
                readxl::read_excel(x),
                stringsAsFactors = FALSE,
                check.names = FALSE
            ))
        }
        read.delim(
            x,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    # -------------------------------------------------------------------------
    # Read featureCounts output and construct the count matrix
    # -------------------------------------------------------------------------
    if (verbose) {
        message("1. Reading featureCounts count file...")
    }
    featurecounts <- read.delim(
        fcounts_file,
        header = TRUE,
        sep = "\t",
        comment.char = "#",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (!"Geneid" %in% colnames(featurecounts)) {
        stop("'Geneid' column was not found in:\n", fcounts_file)
    }

    standard_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    annotation_columns <- intersect(standard_columns, colnames(featurecounts))
    count_columns <- setdiff(colnames(featurecounts), annotation_columns)
    if (length(count_columns) == 0L) {
        stop("No sample count columns were found.")
    }

    gene_ids <- clean_gene_id(featurecounts$Geneid)
    if (anyNA(gene_ids) || any(gene_ids == "")) {
        stop("Empty or NA gene IDs were found after cleaning.")
    }
    if (anyDuplicated(gene_ids)) {
        duplicates <- unique(gene_ids[duplicated(gene_ids)])
        stop(
            "Duplicated gene IDs after cleaning. Examples:\n",
            paste(head(duplicates, 10L), collapse = "\n")
        )
    }

    counts <- as.matrix(featurecounts[, count_columns, drop = FALSE])
    suppressWarnings(storage.mode(counts) <- "numeric")
    rownames(counts) <- gene_ids
    if (anyNA(counts)) {
        stop(
            "NA values were detected after converting sample columns to numeric ",
            "counts. Check for non-count columns or malformed values."
        )
    }
    if (any(counts < 0)) {
        stop("Negative values were detected in the count matrix.")
    }

    original_sample_names <- colnames(counts)
    cleaned_sample_names <- clean_sample_names(
        original_sample_names,
        "featureCounts"
    )
    original_names_by_cleaned <- stats::setNames(
        original_sample_names,
        cleaned_sample_names
    )
    colnames(counts) <- cleaned_sample_names

    # -------------------------------------------------------------------------
    # Read metadata and align shared samples in count-matrix order
    # -------------------------------------------------------------------------
    if (verbose) {
        message("2. Reading metadata...")
    }
    metadata_table <- read_metadata(metadata)
    if (!sample_col %in% colnames(metadata_table)) {
        stop("Metadata does not contain column: ", sample_col)
    }

    metadata_table[[sample_col]] <- as.character(metadata_table[[sample_col]])
    if (anyNA(metadata_table[[sample_col]]) ||
        any(metadata_table[[sample_col]] == "")) {
        stop("Metadata contains empty or NA sample identifiers.")
    }
    if (anyDuplicated(metadata_table[[sample_col]])) {
        duplicates <- unique(
            metadata_table[[sample_col]][duplicated(metadata_table[[sample_col]])]
        )
        stop(
            "Duplicated sample names were found in metadata:\n",
            paste(duplicates, collapse = "\n")
        )
    }

    counts_not_metadata <- setdiff(
        colnames(counts),
        metadata_table[[sample_col]]
    )
    metadata_not_counts <- setdiff(
        metadata_table[[sample_col]],
        colnames(counts)
    )
    if (verbose && length(counts_not_metadata) > 0L) {
        message(
            "\nSamples in featureCounts but not metadata:\n",
            paste(counts_not_metadata, collapse = "\n")
        )
    }
    if (verbose && length(metadata_not_counts) > 0L) {
        message(
            "\nSamples in metadata but not featureCounts:\n",
            paste(metadata_not_counts, collapse = "\n")
        )
    }
    if (strict &&
        (length(counts_not_metadata) > 0L || length(metadata_not_counts) > 0L)) {
        stop(
            "\nfeatureCounts and metadata sample names do not match.\n",
            "Use strict = FALSE to retain only samples present in both."
        )
    }

    common_samples <- colnames(counts)[
        colnames(counts) %in% metadata_table[[sample_col]]
    ]
    if (length(common_samples) == 0L) {
        stop("No matching samples between featureCounts and metadata.")
    }
    counts <- counts[, common_samples, drop = FALSE]
    retained_original_sample_names <- unname(
        original_names_by_cleaned[colnames(counts)]
    )
    metadata_table <- metadata_table[
        match(colnames(counts), metadata_table[[sample_col]]),
        ,
        drop = FALSE
    ]
    rownames(metadata_table) <- metadata_table[[sample_col]]

    # -------------------------------------------------------------------------
    # Build featureCounts and optional external gene annotation
    # -------------------------------------------------------------------------
    if (verbose) {
        message("3. Building gene annotation...")
    }
    genes <- data.frame(
        gene_id = gene_ids,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    extra_annotation <- setdiff(annotation_columns, "Geneid")
    for (column in extra_annotation) {
        genes[[paste0("featureCounts_", column)]] <- featurecounts[[column]]
    }

    if (!is.null(gene_feature_file)) {
        if (verbose) {
            message("4. Reading external gene annotation...")
        }
        gene_features <- read.table(
            gene_feature_file,
            header = FALSE,
            sep = "",
            quote = "\"",
            comment.char = "",
            fill = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        if (ncol(gene_features) < 3L) {
            stop(
                "Gene feature file must contain at least gene_id, gene_type, ",
                "and gene_name."
            )
        }

        standard_names <- c("gene_id", "gene_type", "gene_name", "chr", "start", "end")
        named_columns <- min(ncol(gene_features), length(standard_names))
        colnames(gene_features)[seq_len(named_columns)] <-
            standard_names[seq_len(named_columns)]
        if (ncol(gene_features) > length(standard_names)) {
            extra_count <- ncol(gene_features) - length(standard_names)
            colnames(gene_features)[
                (length(standard_names) + 1L):ncol(gene_features)
            ] <- paste0("annotation_extra_", seq_len(extra_count))
        }

        gene_features$gene_id <- clean_gene_id(gene_features$gene_id)
        if (anyDuplicated(gene_features$gene_id)) {
            warning(
                "Duplicated gene IDs found in external annotation; ",
                "keeping the first occurrence."
            )
            gene_features <- gene_features[
                !duplicated(gene_features$gene_id),
                ,
                drop = FALSE
            ]
        }

        annotation_match <- match(genes$gene_id, gene_features$gene_id)
        annotation_fields <- setdiff(colnames(gene_features), "gene_id")
        for (column in annotation_fields) {
            genes[[column]] <- gene_features[[column]][annotation_match]
        }
    }
    rownames(genes) <- genes$gene_id

    # -------------------------------------------------------------------------
    # Import the optional featureCounts summary and construct sample QC
    # -------------------------------------------------------------------------
    featurecounts_summary <- NULL
    qc <- NULL
    if (!is.null(summary_file)) {
        if (verbose) {
            message("5. Reading featureCounts summary...")
        }
        featurecounts_summary <- read.delim(
            summary_file,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        if (ncol(featurecounts_summary) < 2L) {
            stop("featureCounts summary must contain status and sample columns.")
        }
        colnames(featurecounts_summary)[1L] <- "Status"
        if (anyDuplicated(featurecounts_summary$Status)) {
            stop("Duplicated status rows were found in the featureCounts summary.")
        }

        colnames(featurecounts_summary)[-1L] <- clean_sample_names(
            colnames(featurecounts_summary)[-1L],
            "featureCounts summary"
        )
        summary_matrix <- as.matrix(
            featurecounts_summary[, -1L, drop = FALSE]
        )
        suppressWarnings(storage.mode(summary_matrix) <- "numeric")
        rownames(summary_matrix) <- featurecounts_summary$Status
        if (anyNA(summary_matrix)) {
            stop("Non-numeric or missing values were found in the summary counts.")
        }

        missing_summary <- setdiff(colnames(counts), colnames(summary_matrix))
        if (length(missing_summary) > 0L) {
            stop(
                "Samples in counts but missing from featureCounts summary:\n",
                paste(missing_summary, collapse = "\n")
            )
        }
        summary_matrix <- summary_matrix[, colnames(counts), drop = FALSE]

        if (!"Assigned" %in% rownames(summary_matrix)) {
            warning("'Assigned' was not found in the featureCounts summary.")
        } else {
            assigned_reads <- summary_matrix["Assigned", , drop = TRUE]
            total_input_reads <- colSums(summary_matrix, na.rm = TRUE)
            assignment_percent <- assigned_reads / total_input_reads * 100
            qc <- data.frame(
                SampleName = colnames(counts),
                Total_input_reads = as.numeric(
                    total_input_reads[colnames(counts)]
                ),
                Assigned_reads = as.numeric(assigned_reads[colnames(counts)]),
                Assignment_percent = as.numeric(
                    assignment_percent[colnames(counts)]
                ),
                Count_matrix_total = as.numeric(colSums(counts)),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            rownames(qc) <- qc$SampleName

            status_categories <- setdiff(rownames(summary_matrix), "Assigned")
            for (status in status_categories) {
                qc[[status]] <- as.numeric(
                    summary_matrix[status, rownames(qc)]
                )
            }
            qc <- qc[colnames(counts), , drop = FALSE]

            core_qc_columns <- c(
                "Total_input_reads",
                "Assigned_reads",
                "Assignment_percent",
                "Count_matrix_total"
            )
            for (column in core_qc_columns) {
                metadata_table[[column]] <- qc[rownames(metadata_table), column]
            }
            if (add_all_fc_qc_to_metadata) {
                unassigned_columns <- grep(
                    "^Unassigned_",
                    colnames(qc),
                    value = TRUE
                )
                for (column in unassigned_columns) {
                    metadata_table[[column]] <- qc[
                        rownames(metadata_table),
                        column
                    ]
                }
            }
        }
    }

    # -------------------------------------------------------------------------
    # Validate final ordering and count characteristics
    # -------------------------------------------------------------------------
    if (verbose) {
        message("6. Validating object ordering...")
    }
    if (!identical(colnames(counts), rownames(metadata_table))) {
        stop("Internal error: counts and metadata are not in the same order.")
    }
    if (!is.null(qc) && !identical(colnames(counts), rownames(qc))) {
        stop("Internal error: counts and QC are not in the same order.")
    }
    if (!identical(rownames(counts), rownames(genes))) {
        stop("Internal error: counts and gene annotation are not in the same order.")
    }

    integer_like <- all(abs(counts - round(counts)) < 1e-8)
    if (!integer_like) {
        warning(
            "Count matrix contains non-integer values; check whether ",
            "featureCounts used fractional counting."
        )
    }

    # -------------------------------------------------------------------------
    # Optionally construct an edgeR object and return all aligned components
    # -------------------------------------------------------------------------
    dge <- NULL
    if (make_dge) {
        if (!requireNamespace("edgeR", quietly = TRUE)) {
            warning("Package 'edgeR' is not installed; skipping DGEList creation.")
        } else {
            dge <- edgeR::DGEList(
                counts = counts,
                samples = metadata_table,
                genes = genes
            )
        }
    }

    result <- list(
        counts = counts,
        metadata = metadata_table,
        genes = genes,
        qc = qc,
        featurecounts_summary = featurecounts_summary,
        dge = dge,
        original_sample_names = retained_original_sample_names,
        sample_names = colnames(counts)
    )
    class(result) <- c("featurecounts_project", "list")

    if (verbose) {
        message("")
        message("==============================================")
        message(" featureCounts import complete")
        message("==============================================")
        message("Genes:                   ", format(nrow(counts), big.mark = ","))
        message("Samples:                 ", ncol(counts))
        message("Integer count matrix:    ", integer_like)
        message(
            "Counts <-> metadata:     ",
            identical(colnames(counts), rownames(metadata_table))
        )
        message(
            "Counts <-> genes:        ",
            identical(rownames(counts), rownames(genes))
        )
        if (!is.null(qc)) {
            message(
                "Counts <-> QC:           ",
                identical(colnames(counts), rownames(qc))
            )
            message(
                "Mean assignment rate:    ",
                round(mean(metadata_table$Assignment_percent, na.rm = TRUE), 1L),
                "%"
            )
            message(
                "Min assignment rate:     ",
                round(min(metadata_table$Assignment_percent, na.rm = TRUE), 1L),
                "%"
            )
            message(
                "Max assignment rate:     ",
                round(max(metadata_table$Assignment_percent, na.rm = TRUE), 1L),
                "%"
            )
        }
        if (!is.null(gene_feature_file) && "gene_name" %in% colnames(genes)) {
            message(
                "Annotated gene symbols:  ",
                sum(!is.na(genes$gene_name) & genes$gene_name != ""),
                " / ",
                nrow(genes)
            )
        }
        message("==============================================")
        message("")
    }

    result
}

# =============================================================================
# Quality control and filtering
# =============================================================================


#' Check the internal consistency of a featureCounts project
#'
#' Audit the count matrix, metadata, gene annotation, and optional QC table
#' returned by [read_featurecounts_project()]. The function checks ordering,
#' identifier uniqueness, missing and non-finite values, negative counts, and
#' whether count values are integer-like. Assignment-rate summary statistics
#' are included when the metadata contains `Assignment_percent`.
#'
#' @param dat A `featurecounts_project` object, or a list with `counts`,
#'   `metadata`, and `genes` elements and an optional `qc` element. `counts`
#'   must be a numeric matrix with gene row names and sample column names;
#'   metadata and gene annotation must have row names.
#' @param tolerance Non-negative numeric scalar. Maximum absolute difference
#'   between a count and its rounded value for that count to be considered
#'   integer-like. The default allows negligible floating-point error.
#' @param verbose Logical scalar. Whether to print the formatted audit and
#'   assignment-rate summary. The structured result is returned regardless.
#'
#' @return Invisibly returns a named list with these elements:
#'   \describe{
#'     \item{valid}{Logical scalar indicating whether all structural and value
#'       checks pass. Integer-likeness is reported separately and does not make
#'       fractional featureCounts output structurally invalid.}
#'     \item{all_checks_pass}{Logical scalar indicating whether every
#'       applicable check, including integer-likeness, passes.}
#'     \item{checks}{Named logical vector of positive checks. The QC-order
#'       element is `NA` when `dat$qc` is `NULL`.}
#'     \item{assignment_summary}{The result of [summary()] for
#'       `metadata$Assignment_percent`, or `NULL` when that column is absent.}
#'   }
#'
#' @details
#' The audit is read-only and does not modify `dat`. A `FALSE` result identifies
#' a consistency issue but does not repair it. Fractional counts can be valid
#' when featureCounts was deliberately run with fractional assignment, so the
#' `valid` result excludes the integer-likeness check while
#' `all_checks_pass` includes it.
#'
#' @examples
#' \dontrun{
#' project <- read_featurecounts_project(
#'   fcounts_file = "counts/featureCounts.txt",
#'   metadata = "metadata/samples.csv"
#' )
#'
#' audit <- check_featurecounts_project(project)
#' audit$valid
#' audit$checks
#' }
#'
#' @export
check_featurecounts_project <- function(
    dat,
    tolerance = 1e-8,
    verbose = TRUE
) {
    # -------------------------------------------------------------------------
    # Validate the supplied project structure
    # -------------------------------------------------------------------------
    if (!is.list(dat)) {
        stop("dat must be a featurecounts_project object or compatible list.")
    }

    required_elements <- c("counts", "metadata", "genes")
    missing_elements <- setdiff(required_elements, names(dat))
    if (length(missing_elements) > 0L) {
        stop(
            "dat is missing required elements: ",
            paste(missing_elements, collapse = ", ")
        )
    }

    if (!is.matrix(dat$counts) || !is.numeric(dat$counts)) {
        stop("dat$counts must be a numeric matrix.")
    }
    if (is.null(rownames(dat$counts)) || is.null(colnames(dat$counts))) {
        stop("dat$counts must have gene row names and sample column names.")
    }
    if (!is.data.frame(dat$metadata) || is.null(rownames(dat$metadata))) {
        stop("dat$metadata must be a data frame with sample row names.")
    }
    if (!is.data.frame(dat$genes) || is.null(rownames(dat$genes))) {
        stop("dat$genes must be a data frame with gene row names.")
    }
    if (!is.null(dat$qc) &&
        (!is.data.frame(dat$qc) || is.null(rownames(dat$qc)))) {
        stop("dat$qc must be NULL or a data frame with sample row names.")
    }
    if (!is.numeric(tolerance) ||
        length(tolerance) != 1L ||
        is.na(tolerance) ||
        !is.finite(tolerance) ||
        tolerance < 0) {
        stop("tolerance must be a finite, non-negative numeric scalar.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be a non-missing logical scalar.")
    }

    # -------------------------------------------------------------------------
    # Calculate positive consistency checks
    # -------------------------------------------------------------------------
    counts_metadata_order <- identical(
        colnames(dat$counts),
        rownames(dat$metadata)
    )
    counts_qc_order <- if (is.null(dat$qc)) {
        NA
    } else {
        identical(colnames(dat$counts), rownames(dat$qc))
    }
    counts_genes_order <- identical(
        rownames(dat$counts),
        rownames(dat$genes)
    )
    sample_names_unique <- anyDuplicated(colnames(dat$counts)) == 0L
    gene_ids_unique <- anyDuplicated(rownames(dat$counts)) == 0L
    counts_complete <- !anyNA(dat$counts)
    counts_finite <- all(is.finite(dat$counts))
    counts_nonnegative <- if (counts_complete && counts_finite) {
        !any(dat$counts < 0)
    } else {
        NA
    }
    integer_like_counts <- if (counts_complete && counts_finite) {
        all(abs(dat$counts - round(dat$counts)) <= tolerance)
    } else {
        NA
    }

    checks <- c(
        counts_metadata_order = counts_metadata_order,
        counts_qc_order = counts_qc_order,
        counts_genes_order = counts_genes_order,
        sample_names_unique = sample_names_unique,
        gene_ids_unique = gene_ids_unique,
        counts_complete = counts_complete,
        counts_finite = counts_finite,
        counts_nonnegative = counts_nonnegative,
        integer_like_counts = integer_like_counts
    )

    structural_checks <- checks[names(checks) != "integer_like_counts"]
    valid <- all(structural_checks[!is.na(structural_checks)])
    all_checks_pass <- all(checks[!is.na(checks)])

    assignment_summary <- NULL
    if ("Assignment_percent" %in% colnames(dat$metadata)) {
        assignment_summary <- summary(dat$metadata$Assignment_percent)
    }

    # -------------------------------------------------------------------------
    # Print a human-readable audit matching the original interactive workflow
    # -------------------------------------------------------------------------
    if (verbose) {
        cat("====================================\n")
        cat("Checking featureCounts project\n")
        cat("====================================\n\n")
        cat(
            "counts <-> metadata order: ",
            counts_metadata_order,
            "\n",
            sep = ""
        )
        if (!is.null(dat$qc)) {
            cat(
                "counts <-> QC order:       ",
                counts_qc_order,
                "\n",
                sep = ""
            )
        }
        cat(
            "counts <-> genes order:    ",
            counts_genes_order,
            "\n",
            sep = ""
        )
        cat(
            "Duplicated sample names:   ",
            !sample_names_unique,
            "\n",
            sep = ""
        )
        cat(
            "Duplicated gene IDs:       ",
            !gene_ids_unique,
            "\n",
            sep = ""
        )
        cat(
            "NA values in counts:       ",
            !counts_complete,
            "\n",
            sep = ""
        )
        cat(
            "Non-finite counts:         ",
            !counts_finite,
            "\n",
            sep = ""
        )
        cat(
            "Negative counts:           ",
            !counts_nonnegative,
            "\n",
            sep = ""
        )
        cat(
            "Integer counts:            ",
            integer_like_counts,
            "\n",
            sep = ""
        )

        if (!is.null(assignment_summary)) {
            cat("\nAssignment rate:\n")
            print(assignment_summary)
        }
        cat("\n====================================\n")
    }

    result <- list(
        valid = valid,
        all_checks_pass = all_checks_pass,
        checks = checks,
        assignment_summary = assignment_summary
    )
    invisible(result)
}


# =============================================================================
# Project subsetting and manipulation
# =============================================================================


#' Subset a featureCounts project using sample metadata
#'
#' Select samples from a `featurecounts_project` with an expression evaluated
#' against its metadata, then subset the count matrix, metadata, optional QC,
#' and optional featureCounts summary in a synchronized order. All genes are
#' retained, and a new [edgeR::DGEList()] can be constructed for the selected
#' comparison.
#'
#' @param dat A `featurecounts_project` object, or a compatible list containing
#'   `counts`, `metadata`, and `genes`, with optional `qc` and
#'   `featurecounts_summary` elements. Count-matrix columns and metadata row
#'   names must be identical; count-matrix and gene-annotation row names must
#'   also be identical.
#' @param subset An unquoted expression evaluated in `dat$metadata` that returns
#'   one logical value per sample. For example,
#'   `Tissue == "Lung" & Genotype %in% c("WT", "EPX")`. Values from the
#'   calling environment may also be referenced. Missing results are treated as
#'   `FALSE`.
#' @param make_dge Logical scalar. Whether to build a new [edgeR::DGEList()]
#'   from the subset. Rebuilding rather than subsetting an existing DGEList
#'   resets its library information for the selected comparison. If edgeR is
#'   unavailable, a warning is issued and the returned `dge` is `NULL`.
#' @param drop_levels Logical scalar. Whether to remove unused levels from every
#'   factor column in the subset metadata.
#' @param verbose Logical scalar. Whether to print retained sample and gene
#'   counts together with the selected sample identifiers.
#'
#' @return A new object of class `featurecounts_project`, implemented as a
#'   named list with the following elements:
#'   \describe{
#'     \item{counts}{Gene-by-sample count matrix containing only selected
#'       samples and all input genes.}
#'     \item{metadata}{Selected sample metadata in count-matrix column order.}
#'     \item{genes}{Unfiltered gene annotation in count-matrix row order.}
#'     \item{qc}{Selected sample QC in count-matrix column order, or `NULL`.}
#'     \item{featurecounts_summary}{FeatureCounts status table restricted to
#'       selected sample columns, or `NULL`.}
#'     \item{dge}{New optional [edgeR::DGEList()] for the subset, or `NULL`.}
#'     \item{original_sample_names}{Original featureCounts sample column names
#'       for selected samples in `sample_names` order, or `NULL` when the input
#'       compatible list does not provide this provenance.}
#'     \item{sample_names}{Character vector of retained sample identifiers in
#'       count-matrix order.}
#'   }
#'
#' @details
#' This function does not filter genes. Retaining the complete gene set allows
#' expression filtering, such as [edgeR::filterByExpr()], to be performed later
#' for the specific design and comparison. The input object is not modified.
#'
#' @examples
#' \dontrun{
#' lung_project <- subset_featurecounts_project(
#'   project,
#'   Tissue == "Lung" & Genotype %in% c("WT", "EPX")
#' )
#'
#' lung_project$sample_names
#' lung_project$dge
#' }
#'
#' @export
subset_featurecounts_project <- function(
    dat,
    subset,
    make_dge = TRUE,
    drop_levels = TRUE,
    verbose = TRUE
) {
    # -------------------------------------------------------------------------
    # Validate the project components and their incoming order
    # -------------------------------------------------------------------------
    if (!is.list(dat)) {
        stop("dat must be a featurecounts_project object or compatible list.")
    }

    required_elements <- c("counts", "metadata", "genes")
    missing_elements <- setdiff(required_elements, names(dat))
    if (length(missing_elements) > 0L) {
        stop(
            "Missing components in dat: ",
            paste(missing_elements, collapse = ", ")
        )
    }

    logical_arguments <- list(
        make_dge = make_dge,
        drop_levels = drop_levels,
        verbose = verbose
    )
    invalid_logical <- vapply(
        logical_arguments,
        function(x) length(x) != 1L || is.na(x) || !is.logical(x),
        logical(1)
    )
    if (any(invalid_logical)) {
        stop(
            "The following arguments must be non-missing logical scalars: ",
            paste(names(logical_arguments)[invalid_logical], collapse = ", ")
        )
    }

    if (!is.matrix(dat$counts) || !is.numeric(dat$counts) ||
        is.null(colnames(dat$counts)) ||
        is.null(rownames(dat$counts))) {
        stop("dat$counts must be a numeric matrix with sample and gene names.")
    }
    if (!is.data.frame(dat$metadata) || is.null(rownames(dat$metadata))) {
        stop("dat$metadata must be a data frame with sample row names.")
    }
    if (!is.data.frame(dat$genes) || is.null(rownames(dat$genes))) {
        stop("dat$genes must be a data frame with gene row names.")
    }
    if (!identical(colnames(dat$counts), rownames(dat$metadata))) {
        stop("Input counts and metadata are not in the same sample order.")
    }
    if (!identical(rownames(dat$counts), rownames(dat$genes))) {
        stop("Input counts and gene annotation are not in the same gene order.")
    }
    if (anyDuplicated(colnames(dat$counts))) {
        stop("Input count matrix contains duplicated sample names.")
    }
    if (anyDuplicated(rownames(dat$counts))) {
        stop("Input count matrix contains duplicated gene IDs.")
    }
    if (!is.null(dat$qc)) {
        if (!is.data.frame(dat$qc) || is.null(rownames(dat$qc))) {
            stop("dat$qc must be NULL or a data frame with sample row names.")
        }
        if (!identical(colnames(dat$counts), rownames(dat$qc))) {
            stop("Input counts and QC are not in the same sample order.")
        }
    }
    if (!is.null(dat$sample_names) &&
        !identical(colnames(dat$counts), as.character(dat$sample_names))) {
        stop("dat$sample_names is not aligned with the count matrix.")
    }
    if (!is.null(dat$original_sample_names) &&
        length(dat$original_sample_names) != ncol(dat$counts)) {
        stop(
            "dat$original_sample_names must contain one value per count-matrix ",
            "column. Re-import legacy objects before subsetting."
        )
    }

    metadata <- dat$metadata

    # -------------------------------------------------------------------------
    # Evaluate the caller's expression against sample metadata
    # -------------------------------------------------------------------------
    if (missing(subset)) {
        stop("subset must be supplied as a metadata expression.")
    }
    subset_expression <- substitute(subset)
    keep <- eval(
        subset_expression,
        envir = metadata,
        enclos = parent.frame()
    )
    if (!is.logical(keep)) {
        stop("subset must return a logical vector.")
    }
    if (length(keep) != nrow(metadata)) {
        stop(
            "subset returned ",
            length(keep),
            " values but metadata contains ",
            nrow(metadata),
            " samples."
        )
    }

    keep[is.na(keep)] <- FALSE
    if (!any(keep)) {
        stop("No samples remain after subsetting.")
    }
    selected_samples <- rownames(metadata)[keep]

    # -------------------------------------------------------------------------
    # Subset metadata, counts, genes, and optional QC in a shared order
    # -------------------------------------------------------------------------
    metadata_subset <- metadata[selected_samples, , drop = FALSE]
    if (drop_levels) {
        metadata_subset[] <- lapply(metadata_subset, function(column) {
            if (is.factor(column)) {
                droplevels(column)
            } else {
                column
            }
        })
    }

    counts_subset <- dat$counts[, selected_samples, drop = FALSE]

    # Keep every gene until a design-aware filter is applied downstream.
    genes_subset <- dat$genes[rownames(counts_subset), , drop = FALSE]

    qc_subset <- NULL
    if (!is.null(dat$qc)) {
        qc_subset <- dat$qc[selected_samples, , drop = FALSE]
    }

    featurecounts_summary_subset <- NULL
    if (!is.null(dat$featurecounts_summary)) {
        featurecounts_summary <- dat$featurecounts_summary
        if (!is.data.frame(featurecounts_summary) ||
            !"Status" %in% colnames(featurecounts_summary)) {
            stop(
                "dat$featurecounts_summary must be NULL or a data frame ",
                "containing a Status column."
            )
        }
        if (anyDuplicated(colnames(featurecounts_summary))) {
            stop("dat$featurecounts_summary contains duplicated column names.")
        }
        missing_summary_samples <- setdiff(
            selected_samples,
            colnames(featurecounts_summary)
        )
        if (length(missing_summary_samples) > 0L) {
            stop(
                "Selected samples missing from featureCounts summary: ",
                paste(missing_summary_samples, collapse = ", ")
            )
        }
        featurecounts_summary_subset <- featurecounts_summary[
            ,
            c("Status", selected_samples),
            drop = FALSE
        ]
    }

    # -------------------------------------------------------------------------
    # Validate the subset before constructing downstream objects
    # -------------------------------------------------------------------------
    if (!identical(colnames(counts_subset), rownames(metadata_subset))) {
        stop("Internal error: subset counts and metadata order do not match.")
    }
    if (!is.null(qc_subset) &&
        !identical(colnames(counts_subset), rownames(qc_subset))) {
        stop("Internal error: subset counts and QC order do not match.")
    }
    if (!identical(rownames(counts_subset), rownames(genes_subset))) {
        stop("Internal error: subset counts and gene annotation order do not match.")
    }

    # -------------------------------------------------------------------------
    # Build a fresh DGEList so library information reflects selected samples
    # -------------------------------------------------------------------------
    dge_subset <- NULL
    if (make_dge) {
        if (!requireNamespace("edgeR", quietly = TRUE)) {
            warning("edgeR is not installed; DGEList was not generated.")
        } else {
            # Do not let library fields from a previous edgeR object override
            # values that must be recalculated from this subset's counts.
            dge_samples <- metadata_subset[
                ,
                setdiff(
                    colnames(metadata_subset),
                    c("lib.size", "norm.factors")
                ),
                drop = FALSE
            ]
            dge_subset <- edgeR::DGEList(
                counts = counts_subset,
                samples = dge_samples,
                genes = genes_subset
            )
        }
    }

    original_sample_names_subset <- NULL
    if (!is.null(dat$original_sample_names)) {
        selected_positions <- match(selected_samples, colnames(dat$counts))
        original_sample_names_subset <- as.character(
            dat$original_sample_names[selected_positions]
        )
    }

    result <- list(
        counts = counts_subset,
        metadata = metadata_subset,
        genes = genes_subset,
        qc = qc_subset,
        featurecounts_summary = featurecounts_summary_subset,
        dge = dge_subset,
        original_sample_names = original_sample_names_subset,
        sample_names = selected_samples
    )
    class(result) <- c("featurecounts_project", "list")

    if (verbose) {
        message("")
        message("==========================================")
        message(" featureCounts project subset")
        message("==========================================")
        message(
            "Samples retained: ",
            ncol(counts_subset),
            " / ",
            ncol(dat$counts)
        )
        message("Genes retained:   ", nrow(counts_subset))
        message("")
        message("Samples:")
        message(paste(selected_samples, collapse = "\n"))
        message("==========================================")
        message("")
    }

    result
}


# =============================================================================
# Shared bulk/pseudobulk expression preparation
# =============================================================================


#' Prepare normalized bulk or pseudobulk expression for reusable plots
#'
#' Validate a featureCounts-style project, optionally select samples, calculate
#' library-composition normalization factors, and produce one aligned log2-CPM
#' matrix that can be reused by multiple exploratory visualization functions.
#'
#' @param dat A list-like featureCounts project containing `counts`, `metadata`,
#'   and `genes`. `counts` must be a non-empty numeric gene-by-sample matrix or
#'   Matrix object. Count columns must be identical to metadata row names, and
#'   count rows must be identical to gene-annotation row names.
#' @param subset Optional unquoted logical expression evaluated within
#'   `dat$metadata`, with access to objects in the calling environment. It must
#'   return one value per sample; missing results are treated as `FALSE`.
#' @param prior_count Non-negative finite numeric scalar passed to
#'   [edgeR::cpm()] when calculating log2 counts per million.
#' @param norm_method One of `"TMM"`, `"TMMwsp"`, `"RLE"`,
#'   `"upperquartile"`, or `"none"`, passed to
#'   [edgeR::calcNormFactors()]. `"none"` retains library-size normalization
#'   while disabling an additional composition-normalization adjustment.
#' @param verbose Logical scalar. Whether to print the number of genes and
#'   samples, normalization label, prior count, and edgeR library summary.
#'
#' @return An object of class `bulk_expression_prepared`, implemented as a
#'   named list containing:
#'   \describe{
#'     \item{counts}{The original or subsetted raw count matrix.}
#'     \item{metadata}{Sample metadata aligned to count-matrix columns.}
#'     \item{genes}{Gene annotation aligned to count-matrix rows.}
#'     \item{dge}{The normalized [edgeR::DGEList()] object.}
#'     \item{logCPM}{The normalized gene-by-sample log2-CPM matrix.}
#'     \item{normalization}{A list containing the normalization `method`,
#'       `prior_count`, and human-readable `label`.}
#'   }
#'
#' @details
#' This function centralizes normalization for descriptive bulk RNA-seq and
#' replicate-aware pseudobulk visualizations. Prepare the subset once and pass
#' the returned object to [plot_bulk_violin()] or
#' [plot_gene_set_heatmap()] so both figures use exactly the same samples,
#' library factors, and prior count. Raw counts are retained unchanged in the
#' return object. This preparation does not filter genes, construct a design,
#' or replace design-aware differential-expression analysis. The `edgeR`
#' package is required.
#'
#' @examples
#' \dontrun{
#' bulk <- prepare_bulk_expression(
#'     project,
#'     subset = Tissue == "Lung",
#'     prior_count = 2,
#'     norm_method = "TMM"
#' )
#' bulk$logCPM
#' bulk$normalization
#' }
#'
#' @export
prepare_bulk_expression <- function(
    dat,
    subset = NULL,
    prior_count = 2,
    norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
    verbose = TRUE
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package `edgeR` is required.", call. = FALSE)
    }
    norm_method <- match.arg(norm_method)

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
    }
    if (!is.numeric(prior_count) || length(prior_count) != 1L ||
        is.na(prior_count) || !is.finite(prior_count) || prior_count < 0) {
        stop(
            "`prior_count` must be one finite non-negative number.",
            call. = FALSE
        )
    }
    if (!is.list(dat)) {
        stop("`dat` must be a list-like featureCounts project.", call. = FALSE)
    }
    required_objects <- c("counts", "metadata", "genes")
    missing_objects <- setdiff(required_objects, names(dat))
    if (length(missing_objects) > 0L) {
        stop(
            "`dat` is missing: ", paste(missing_objects, collapse = ", "),
            call. = FALSE
        )
    }

    counts <- dat$counts
    meta <- as.data.frame(dat$metadata)
    anno <- as.data.frame(dat$genes)
    if ((!is.matrix(counts) && !inherits(counts, "Matrix")) ||
        !is.numeric(counts) || nrow(counts) == 0L || ncol(counts) == 0L) {
        stop(
            "`dat$counts` must be a non-empty numeric matrix-like object.",
            call. = FALSE
        )
    }
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
        is.null(rownames(meta)) || is.null(rownames(anno))) {
        stop(
            "Counts, metadata, and gene annotations require row/column names.",
            call. = FALSE
        )
    }
    if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
        anyDuplicated(rownames(meta)) || anyDuplicated(rownames(anno))) {
        stop("Sample and gene identifiers must be unique.", call. = FALSE)
    }
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
        stop(
            "`dat$counts` must contain finite, non-negative values.",
            call. = FALSE
        )
    }
    if (!identical(colnames(counts), rownames(meta))) {
        stop(
            "colnames(dat$counts) and rownames(dat$metadata) must be ",
            "identical and in the same order.", call. = FALSE
        )
    }
    if (!identical(rownames(counts), rownames(anno))) {
        stop(
            "rownames(dat$counts) and rownames(dat$genes) must be identical ",
            "and in the same order.", call. = FALSE
        )
    }

    subset_expr <- substitute(subset)
    if (!identical(subset_expr, quote(NULL))) {
        keep_samples <- eval(
            subset_expr,
            envir = meta,
            enclos = parent.frame()
        )
        if (!is.logical(keep_samples) || length(keep_samples) != nrow(meta)) {
            stop(
                "`subset` must return one logical value per metadata row.",
                call. = FALSE
            )
        }
        keep_samples[is.na(keep_samples)] <- FALSE
        if (!any(keep_samples)) {
            stop("No samples remain after subsetting.", call. = FALSE)
        }
        meta <- meta[keep_samples, , drop = FALSE]
        counts <- counts[, rownames(meta), drop = FALSE]
    }

    library_sizes <- colSums(counts)
    if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
        stop(
            "At least one selected sample has library size <= 0.",
            call. = FALSE
        )
    }

    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = norm_method)
    logcpm <- edgeR::cpm(
        dge,
        log = TRUE,
        prior.count = prior_count,
        normalized.lib.sizes = TRUE
    )
    if (!identical(colnames(logcpm), rownames(meta)) ||
        !identical(rownames(logcpm), rownames(anno))) {
        stop(
            "Internal error: normalized expression is not aligned with ",
            "metadata and gene annotations.", call. = FALSE
        )
    }
    if (anyNA(logcpm) || any(!is.finite(logcpm))) {
        stop(
            "Internal error: normalized expression contains non-finite values.",
            call. = FALSE
        )
    }

    normalization_label <- if (norm_method == "TMM") {
        "TMM-normalized log2 CPM"
    } else if (norm_method == "none") {
        "library-size-normalized log2 CPM"
    } else {
        paste0(norm_method, "-normalized log2 CPM")
    }
    result <- list(
        counts = counts,
        metadata = meta,
        genes = anno,
        dge = dge,
        logCPM = logcpm,
        normalization = list(
            method = norm_method,
            prior_count = prior_count,
            label = normalization_label
        )
    )
    class(result) <- c("bulk_expression_prepared", "list")

    if (verbose) {
        message("")
        message("========================================")
        message(" Bulk RNA-seq expression preparation")
        message("========================================")
        message("Genes:          ", nrow(logcpm))
        message("Samples:        ", ncol(logcpm))
        message("Normalization:  ", normalization_label)
        message("Prior count:    ", prior_count)
        message("")
        message("Library sizes:")
        print(dge$samples)
        message("========================================")
    }

    result
}


# =============================================================================
# Pairwise differential expression
# =============================================================================


#' Prepare a two-group RNA-seq differential-expression comparison
#'
#' Validate and select samples, remove incomplete model rows, construct a
#' two-group design with optional covariates, and define a contrast for
#' `group2 - group1`. This internal helper centralizes the design semantics used
#' by [edgeR_pairwise()] and [limma_voom_pairwise()].
#'
#' @param dat A `featurecounts_project` or compatible list containing a numeric
#'   gene-by-sample `counts` matrix, sample `metadata`, and gene annotation
#'   `genes` data frame in matching orders.
#' @param group_col Character scalar naming the metadata column that defines the
#'   comparison groups.
#' @param group1 Character scalar naming the reference group.
#' @param group2 Character scalar naming the comparison group. Positive log-fold
#'   changes represent higher expression in `group2` than `group1`.
#' @param subset_expr `NULL` or a quoted expression that returns one logical
#'   value per metadata row. Missing values are treated as `FALSE`.
#' @param subset_env Environment used to resolve names in `subset_expr` that are
#'   not metadata columns.
#' @param covariates `NULL` or a character vector naming additional metadata
#'   variables to include in the model.
#' @param verbose Logical scalar. Whether to print the comparison, replicate
#'   counts, removed incomplete samples, and covariates.
#'
#' @return A list containing selected `counts`, `metadata`, and `genes`; the
#'   ordered group factor; design matrix; numeric `group2 - group1` contrast;
#'   group names and group-column name; and sample counts for both groups.
#'
#' @keywords internal
.prepare_pairwise_dat <- function(
    dat,
    group_col,
    group1,
    group2,
    subset_expr = NULL,
    subset_env = parent.frame(),
    covariates = NULL,
    verbose = TRUE
) {
    # -------------------------------------------------------------------------
    # Validate the project and comparison specification
    # -------------------------------------------------------------------------
    if (!is.list(dat)) {
        stop("dat must be a featurecounts_project object or compatible list.")
    }
    required_elements <- c("counts", "metadata", "genes")
    missing_elements <- setdiff(required_elements, names(dat))
    if (length(missing_elements) > 0L) {
        stop("dat is missing: ", paste(missing_elements, collapse = ", "))
    }
    if (!is.matrix(dat$counts) || !is.numeric(dat$counts) ||
        is.null(rownames(dat$counts)) || is.null(colnames(dat$counts))) {
        stop("dat$counts must be a numeric matrix with gene and sample names.")
    }
    if (anyNA(dat$counts) || any(!is.finite(dat$counts)) || any(dat$counts < 0)) {
        stop("dat$counts must contain finite, non-missing, non-negative values.")
    }
    if (anyDuplicated(rownames(dat$counts)) ||
        anyDuplicated(colnames(dat$counts))) {
        stop("dat$counts must have unique gene and sample names.")
    }
    if (!is.data.frame(dat$metadata) || is.null(rownames(dat$metadata))) {
        stop("dat$metadata must be a data frame with sample row names.")
    }
    if (!is.data.frame(dat$genes) || is.null(rownames(dat$genes))) {
        stop("dat$genes must be a data frame with gene row names.")
    }
    if (!identical(colnames(dat$counts), rownames(dat$metadata))) {
        stop("dat$counts and dat$metadata are not in the same sample order.")
    }
    if (!identical(rownames(dat$counts), rownames(dat$genes))) {
        stop("dat$counts and dat$genes are not in the same gene order.")
    }

    comparison_values <- list(
        group_col = group_col,
        group1 = group1,
        group2 = group2
    )
    invalid_comparison <- vapply(
        comparison_values,
        function(x) !is.character(x) || length(x) != 1L || is.na(x) || x == "",
        logical(1)
    )
    if (any(invalid_comparison)) {
        stop(
            "The following arguments must be non-empty character scalars: ",
            paste(names(comparison_values)[invalid_comparison], collapse = ", ")
        )
    }
    if (identical(group1, group2)) {
        stop("group1 and group2 must be different.")
    }
    if (!group_col %in% colnames(dat$metadata)) {
        stop("group_col '", group_col, "' was not found in metadata.")
    }
    if (!is.null(covariates)) {
        if (!is.character(covariates) || anyNA(covariates) ||
            any(covariates == "") || anyDuplicated(covariates)) {
            stop("covariates must be NULL or unique, non-empty column names.")
        }
        if (group_col %in% covariates) {
            stop("group_col must not also be listed in covariates.")
        }
    }
    if (!is.environment(subset_env)) {
        stop("subset_env must be an environment.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be a non-missing logical scalar.")
    }

    counts <- dat$counts
    metadata <- dat$metadata
    genes <- dat$genes

    # -------------------------------------------------------------------------
    # Apply the optional metadata expression, then retain the requested pair
    # -------------------------------------------------------------------------
    if (is.null(subset_expr)) {
        keep_subset <- rep(TRUE, nrow(metadata))
    } else {
        keep_subset <- eval(
            subset_expr,
            envir = metadata,
            enclos = subset_env
        )
        if (!is.logical(keep_subset)) {
            stop("subset must evaluate to a logical vector.")
        }
        if (length(keep_subset) != nrow(metadata)) {
            stop(
                "subset returned ",
                length(keep_subset),
                " values but metadata contains ",
                nrow(metadata),
                " samples."
            )
        }
        keep_subset[is.na(keep_subset)] <- FALSE
    }

    keep_group <- as.character(metadata[[group_col]]) %in% c(group1, group2)
    keep <- keep_subset & keep_group
    if (!any(keep)) {
        stop("No samples remain after metadata and group selection.")
    }

    selected_samples <- rownames(metadata)[keep]
    metadata <- metadata[selected_samples, , drop = FALSE]
    counts <- counts[, selected_samples, drop = FALSE]

    # -------------------------------------------------------------------------
    # Remove samples missing any model variable
    # -------------------------------------------------------------------------
    variables_needed <- c(group_col, covariates)
    missing_variables <- setdiff(variables_needed, colnames(metadata))
    if (length(missing_variables) > 0L) {
        stop(
            "Variables missing from metadata: ",
            paste(missing_variables, collapse = ", ")
        )
    }

    complete_samples <- complete.cases(
        metadata[, variables_needed, drop = FALSE]
    )
    if (any(!complete_samples)) {
        if (verbose) {
            message(
                "Removing ",
                sum(!complete_samples),
                " sample(s) with missing model metadata."
            )
        }
        metadata <- metadata[complete_samples, , drop = FALSE]
        counts <- counts[, rownames(metadata), drop = FALSE]
    }

    group_values <- as.character(metadata[[group_col]])
    n_group1 <- sum(group_values == group1)
    n_group2 <- sum(group_values == group2)
    if (n_group1 == 0L) {
        stop("No complete samples were found for group1: ", group1)
    }
    if (n_group2 == 0L) {
        stop("No complete samples were found for group2: ", group2)
    }
    if (n_group1 < 2L || n_group2 < 2L) {
        warning(
            "One group contains fewer than two biological samples; ",
            "differential-expression inference may not be reliable."
        )
    }

    # -------------------------------------------------------------------------
    # Build a no-intercept group design and optional covariate columns
    # -------------------------------------------------------------------------
    group <- factor(group_values, levels = c(group1, group2))
    design_group <- cbind(
        G1 = as.numeric(group == group1),
        G2 = as.numeric(group == group2)
    )
    design <- design_group

    if (!is.null(covariates) && length(covariates) > 0L) {
        covariate_data <- droplevels(
            metadata[, covariates, drop = FALSE]
        )
        no_variation <- vapply(
            covariate_data,
            function(x) length(unique(x[!is.na(x)])) <= 1L,
            logical(1)
        )
        if (any(no_variation)) {
            stop(
                "These covariates have no variation within selected samples: ",
                paste(names(no_variation)[no_variation], collapse = ", ")
            )
        }

        covariate_design <- model.matrix(~ ., data = covariate_data)
        if ("(Intercept)" %in% colnames(covariate_design)) {
            covariate_design <- covariate_design[
                ,
                colnames(covariate_design) != "(Intercept)",
                drop = FALSE
            ]
        }
        reserved_columns <- intersect(
            colnames(covariate_design),
            colnames(design_group)
        )
        if (length(reserved_columns) > 0L) {
            stop(
                "Covariate design uses reserved columns: ",
                paste(reserved_columns, collapse = ", ")
            )
        }
        design <- cbind(design_group, covariate_design)
    }

    rownames(design) <- rownames(metadata)
    if (qr(design)$rank < ncol(design)) {
        stop(
            "Design matrix is not full rank. Check whether group and ",
            "covariates are confounded."
        )
    }

    contrast <- stats::setNames(rep(0, ncol(design)), colnames(design))
    contrast["G1"] <- -1
    contrast["G2"] <- 1
    genes <- genes[rownames(counts), , drop = FALSE]

    if (verbose) {
        message("")
        message("Comparison: ", group2, " vs ", group1)
        message(group1, ": n = ", n_group1)
        message(group2, ": n = ", n_group2)
        if (!is.null(covariates) && length(covariates) > 0L) {
            message("Covariates: ", paste(covariates, collapse = ", "))
        }
    }

    list(
        counts = counts,
        metadata = metadata,
        genes = genes,
        group = group,
        design = design,
        contrast = contrast,
        n_group1 = n_group1,
        n_group2 = n_group2,
        group1 = group1,
        group2 = group2,
        group_col = group_col
    )
}


#' Run a pairwise edgeR quasi-likelihood analysis
#'
#' Select two biological groups, optionally restrict samples and adjust for
#' covariates, filter genes for that comparison, apply TMM normalization, fit
#' an edgeR quasi-likelihood model, and return an annotated ranked table. The
#' contrast is always `group2 - group1`, so positive `logFC` means higher
#' expression in `group2`.
#'
#' @param dat A `featurecounts_project` or compatible list containing raw
#'   gene-by-sample counts, aligned sample metadata, and aligned gene
#'   annotation. Counts must be finite, non-missing, and non-negative.
#' @param group_col Character scalar naming the metadata grouping column.
#' @param group1 Character scalar naming the reference group.
#' @param group2 Character scalar naming the comparison group.
#' @param subset `NULL` or an unquoted metadata expression used to restrict
#'   samples before retaining `group1` and `group2`. Missing expression results
#'   are treated as `FALSE`.
#' @param covariates `NULL` or a character vector naming additional metadata
#'   variables in the design. Incomplete samples are removed and the final
#'   design must be full rank.
#' @param robust Logical scalar passed to [edgeR::estimateDisp()] and
#'   [edgeR::glmQLFit()] for robust empirical Bayes estimation.
#' @param fdr_cutoff Numeric scalar in `[0, 1]` used with `logfc_cutoff` to label
#'   genes in the returned `Direction` column. It does not remove result rows.
#' @param logfc_cutoff Non-negative numeric scalar giving the absolute log2-fold
#'   change required for a directional label. When zero, only strictly positive
#'   or negative significant effects receive directional labels.
#' @param verbose Logical scalar. Whether to print sample counts, filtering
#'   totals, covariates, and significant-gene summaries.
#'
#' @return A named list containing `method`, comparison label, annotated
#'   `results`, filtered and normalized `dge`, fitted `fit`, quasi-likelihood
#'   `test`, `design`, numeric `contrast`, selected `metadata`, the original
#'   gene-length `keep_genes` filter, group names, and group sample counts.
#'
#' @details
#' Expression filtering is performed with [edgeR::filterByExpr()] only after
#' the requested samples and model have been selected. Existing `group`,
#' `lib.size`, and `norm.factors` metadata columns are excluded when constructing
#' the DGEList, and the prepared comparison factor is stored as its group.
#' Annotation fields that conflict with fitted statistic names receive an
#' `annotation_` prefix so fitted statistics remain authoritative. Genes are
#' labelled `Up_<group2>`, `Up_<group1>`, or `NS`; the complete ranked table is
#' retained. Replication, batch adjustment, pairing variables, and other
#' covariates must reflect the study design. A warning is issued when a group
#' contains fewer than two samples.
#'
#' @examples
#' \dontrun{
#' result <- edgeR_pairwise(
#'   project,
#'   group_col = "Genotype",
#'   group1 = "WT",
#'   group2 = "cKO",
#'   subset = Tissue == "Lung",
#'   covariates = "Batch"
#' )
#'
#' head(result$results)
#' }
#'
#' @export
edgeR_pairwise <- function(
    dat,
    group_col,
    group1,
    group2,
    subset = NULL,
    covariates = NULL,
    robust = TRUE,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0,
    verbose = TRUE
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required.")
    }

    logical_arguments <- list(robust = robust, verbose = verbose)
    invalid_logical <- vapply(
        logical_arguments,
        function(x) length(x) != 1L || is.na(x) || !is.logical(x),
        logical(1)
    )
    if (any(invalid_logical)) {
        stop(
            "The following arguments must be non-missing logical scalars: ",
            paste(names(logical_arguments)[invalid_logical], collapse = ", ")
        )
    }
    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        is.na(fdr_cutoff) || !is.finite(fdr_cutoff) ||
        fdr_cutoff < 0 || fdr_cutoff > 1) {
        stop("fdr_cutoff must be a finite numeric scalar between zero and one.")
    }
    if (!is.numeric(logfc_cutoff) || length(logfc_cutoff) != 1L ||
        is.na(logfc_cutoff) || !is.finite(logfc_cutoff) || logfc_cutoff < 0) {
        stop("logfc_cutoff must be a finite, non-negative numeric scalar.")
    }

    subset_expression <- substitute(subset)
    if (identical(subset_expression, quote(NULL))) {
        subset_expression <- NULL
    }
    prep <- .prepare_pairwise_dat(
        dat = dat,
        group_col = group_col,
        group1 = group1,
        group2 = group2,
        subset_expr = subset_expression,
        subset_env = parent.frame(),
        covariates = covariates,
        verbose = verbose
    )

    dge_samples <- prep$metadata[
        ,
        setdiff(
            colnames(prep$metadata),
            c("group", "lib.size", "norm.factors")
        ),
        drop = FALSE
    ]
    dge <- edgeR::DGEList(
        counts = prep$counts,
        samples = dge_samples,
        group = prep$group,
        genes = prep$genes
    )
    n_genes_before <- nrow(dge)

    keep <- edgeR::filterByExpr(dge, design = prep$design)
    if (!any(keep)) {
        stop("No genes passed edgeR::filterByExpr for this comparison.")
    }
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    n_genes_after <- nrow(dge)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge <- edgeR::estimateDisp(
        dge,
        design = prep$design,
        robust = robust
    )
    fit <- edgeR::glmQLFit(
        dge,
        design = prep$design,
        robust = robust
    )
    test <- edgeR::glmQLFTest(fit, contrast = prep$contrast)
    table <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table

    statistic_columns <- c("logFC", "logCPM", "F", "PValue", "FDR")
    missing_statistics <- setdiff(statistic_columns, colnames(table))
    if (length(missing_statistics) > 0L) {
        stop(
            "edgeR result is missing expected statistics: ",
            paste(missing_statistics, collapse = ", ")
        )
    }
    statistics <- table[, statistic_columns, drop = FALSE]
    gene_annotation <- prep$genes[rownames(statistics), , drop = FALSE]
    annotation_names <- colnames(gene_annotation)
    conflicting_annotation <- annotation_names %in%
        c(statistic_columns, "Direction")
    annotation_names[conflicting_annotation] <- paste0(
        "annotation_",
        annotation_names[conflicting_annotation]
    )
    colnames(gene_annotation) <- make.unique(annotation_names)
    results <- cbind(gene_annotation, statistics)

    results$Direction <- "NS"
    significant <- !is.na(results$FDR) & results$FDR < fdr_cutoff
    if (logfc_cutoff == 0) {
        up_group2 <- significant & !is.na(results$logFC) & results$logFC > 0
        up_group1 <- significant & !is.na(results$logFC) & results$logFC < 0
    } else {
        up_group2 <- significant &
            !is.na(results$logFC) & results$logFC >= logfc_cutoff
        up_group1 <- significant &
            !is.na(results$logFC) & results$logFC <= -logfc_cutoff
    }
    results$Direction[which(up_group2)] <- paste0("Up_", prep$group2)
    results$Direction[which(up_group1)] <- paste0("Up_", prep$group1)

    if (verbose) {
        message("")
        message("edgeR QL analysis complete")
        message("Genes before filtering: ", n_genes_before)
        message("Genes after filtering:  ", n_genes_after)
        message("FDR < ", fdr_cutoff, ": ", sum(significant))
        message(
            prep$group2,
            " up: ",
            sum(results$Direction == paste0("Up_", prep$group2))
        )
        message(
            prep$group1,
            " up: ",
            sum(results$Direction == paste0("Up_", prep$group1))
        )
    }

    list(
        method = "edgeR_QL",
        comparison = paste0(prep$group2, "_vs_", prep$group1),
        results = results,
        dge = dge,
        fit = fit,
        test = test,
        design = prep$design,
        contrast = prep$contrast,
        metadata = prep$metadata,
        keep_genes = keep,
        group1 = prep$group1,
        group2 = prep$group2,
        n_group1 = prep$n_group1,
        n_group2 = prep$n_group2
    )
}


#' Run a pairwise limma-voom differential-expression analysis
#'
#' Select two biological groups, optionally restrict samples and adjust for
#' covariates, filter genes, apply TMM normalization and voom precision weights,
#' fit a limma model, and return an annotated ranked table. The contrast is
#' always `group2 - group1`, so positive `logFC` means higher expression in
#' `group2`.
#'
#' @inheritParams edgeR_pairwise
#' @param robust Logical scalar passed to [limma::eBayes()] for robust
#'   empirical Bayes estimation.
#' @param trend Logical scalar passed to [limma::eBayes()] to enable an
#'   intensity-dependent prior trend after voom.
#' @param voom_plot Logical scalar. Whether [limma::voom()] draws its
#'   mean-variance trend plot on the active graphics device.
#'
#' @return A named list containing `method`, comparison label, annotated
#'   `results`, filtered and normalized `dge`, voom object `voom`, empirical
#'   Bayes model `fit`, `design`, named contrast matrix, selected `metadata`, the
#'   original gene-length `keep_genes` filter, group names, and sample counts.
#'
#' @details
#' Genes are filtered with [edgeR::filterByExpr()] after sample selection, then
#' normalized with TMM before [limma::voom()]. Existing `group`, `lib.size`, and
#' `norm.factors` metadata fields are excluded from DGEList construction, and
#' the prepared comparison factor is stored as its group. Annotation fields
#' that conflict with fitted statistic names receive an `annotation_` prefix.
#' `adj.P.Val` is retained and copied to the standardized `FDR` result column.
#' Direction labels use the same thresholds and interpretation as
#' [edgeR_pairwise()]. Robust empirical Bayes estimation may require the
#' `statmod` package through limma.
#'
#' @examples
#' \dontrun{
#' result <- limma_voom_pairwise(
#'   project,
#'   group_col = "Genotype",
#'   group1 = "WT",
#'   group2 = "cKO",
#'   subset = Tissue == "Lung",
#'   covariates = "Batch"
#' )
#'
#' head(result$results)
#' }
#'
#' @export
limma_voom_pairwise <- function(
    dat,
    group_col,
    group1,
    group2,
    subset = NULL,
    covariates = NULL,
    robust = TRUE,
    trend = FALSE,
    voom_plot = FALSE,
    fdr_cutoff = 0.05,
    logfc_cutoff = 0,
    verbose = TRUE
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for DGEList construction, filtering, and TMM.")
    }
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma is required.")
    }

    logical_arguments <- list(
        robust = robust,
        trend = trend,
        voom_plot = voom_plot,
        verbose = verbose
    )
    invalid_logical <- vapply(
        logical_arguments,
        function(x) length(x) != 1L || is.na(x) || !is.logical(x),
        logical(1)
    )
    if (any(invalid_logical)) {
        stop(
            "The following arguments must be non-missing logical scalars: ",
            paste(names(logical_arguments)[invalid_logical], collapse = ", ")
        )
    }
    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        is.na(fdr_cutoff) || !is.finite(fdr_cutoff) ||
        fdr_cutoff < 0 || fdr_cutoff > 1) {
        stop("fdr_cutoff must be a finite numeric scalar between zero and one.")
    }
    if (!is.numeric(logfc_cutoff) || length(logfc_cutoff) != 1L ||
        is.na(logfc_cutoff) || !is.finite(logfc_cutoff) || logfc_cutoff < 0) {
        stop("logfc_cutoff must be a finite, non-negative numeric scalar.")
    }

    subset_expression <- substitute(subset)
    if (identical(subset_expression, quote(NULL))) {
        subset_expression <- NULL
    }
    prep <- .prepare_pairwise_dat(
        dat = dat,
        group_col = group_col,
        group1 = group1,
        group2 = group2,
        subset_expr = subset_expression,
        subset_env = parent.frame(),
        covariates = covariates,
        verbose = verbose
    )

    dge_samples <- prep$metadata[
        ,
        setdiff(
            colnames(prep$metadata),
            c("group", "lib.size", "norm.factors")
        ),
        drop = FALSE
    ]
    dge <- edgeR::DGEList(
        counts = prep$counts,
        samples = dge_samples,
        group = prep$group,
        genes = prep$genes
    )
    n_genes_before <- nrow(dge)

    keep <- edgeR::filterByExpr(dge, design = prep$design)
    if (!any(keep)) {
        stop("No genes passed edgeR::filterByExpr for this comparison.")
    }
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    n_genes_after <- nrow(dge)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    voom <- limma::voom(dge, design = prep$design, plot = voom_plot)
    fit <- limma::lmFit(voom, prep$design)

    comparison <- paste0(prep$group2, "_vs_", prep$group1)
    contrast_matrix <- matrix(
        prep$contrast,
        ncol = 1L,
        dimnames = list(colnames(prep$design), comparison)
    )
    fit <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
    fit <- limma::eBayes(fit, robust = robust, trend = trend)
    table <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "P")

    statistic_columns <- c(
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
    )
    missing_statistics <- setdiff(statistic_columns, colnames(table))
    if (length(missing_statistics) > 0L) {
        stop(
            "limma result is missing expected statistics: ",
            paste(missing_statistics, collapse = ", ")
        )
    }
    statistics <- table[, statistic_columns, drop = FALSE]
    gene_annotation <- prep$genes[rownames(statistics), , drop = FALSE]
    annotation_names <- colnames(gene_annotation)
    conflicting_annotation <- annotation_names %in%
        c(statistic_columns, "FDR", "Direction")
    annotation_names[conflicting_annotation] <- paste0(
        "annotation_",
        annotation_names[conflicting_annotation]
    )
    colnames(gene_annotation) <- make.unique(annotation_names)
    results <- cbind(gene_annotation, statistics)
    results$FDR <- results$adj.P.Val
    results$Direction <- "NS"

    significant <- !is.na(results$FDR) & results$FDR < fdr_cutoff
    if (logfc_cutoff == 0) {
        up_group2 <- significant & !is.na(results$logFC) & results$logFC > 0
        up_group1 <- significant & !is.na(results$logFC) & results$logFC < 0
    } else {
        up_group2 <- significant &
            !is.na(results$logFC) & results$logFC >= logfc_cutoff
        up_group1 <- significant &
            !is.na(results$logFC) & results$logFC <= -logfc_cutoff
    }
    results$Direction[which(up_group2)] <- paste0("Up_", prep$group2)
    results$Direction[which(up_group1)] <- paste0("Up_", prep$group1)

    if (verbose) {
        message("")
        message("limma-voom analysis complete")
        message("Genes before filtering: ", n_genes_before)
        message("Genes after filtering:  ", n_genes_after)
        message("FDR < ", fdr_cutoff, ": ", sum(significant))
        message(
            prep$group2,
            " up: ",
            sum(results$Direction == paste0("Up_", prep$group2))
        )
        message(
            prep$group1,
            " up: ",
            sum(results$Direction == paste0("Up_", prep$group1))
        )
    }

    list(
        method = "limma_voom",
        comparison = comparison,
        results = results,
        dge = dge,
        voom = voom,
        fit = fit,
        design = prep$design,
        contrast = contrast_matrix,
        metadata = prep$metadata,
        keep_genes = keep,
        group1 = prep$group1,
        group2 = prep$group2,
        n_group1 = prep$n_group1,
        n_group2 = prep$n_group2
    )
}


# Internal collector shared by the public edgeR collector and plotter.
.collect_pairwise_results <- function(
    de_list,
    name_sep,
    logfc_col,
    padj_col,
    allow_data_frames = FALSE,
    require_names = TRUE
) {
    if (!is.list(de_list) || length(de_list) == 0L) {
        stop("de_list must be a non-empty list of pairwise result objects.")
    }
    list_names <- names(de_list)
    if (is.null(list_names)) {
        if (require_names) {
            stop("de_list must be named so analyses can be identified.")
        }
        list_names <- paste0("Analysis_", seq_along(de_list))
    }
    if (length(list_names) != length(de_list) || anyNA(list_names) ||
        any(!nzchar(list_names)) || anyDuplicated(list_names)) {
        stop("de_list names must be unique, non-missing, and non-empty.")
    }

    collected <- lapply(seq_along(de_list), function(index) {
        analysis_name <- list_names[index]
        object <- de_list[[index]]
        is_pairwise_object <- is.list(object) &&
            "results" %in% names(object) && is.data.frame(object$results)
        if (is_pairwise_object) {
            result <- object$results
        } else if (allow_data_frames && is.data.frame(object)) {
            result <- object
        } else {
            stop(
                "'", analysis_name,
                "' does not look like an edgeR_pairwise result object."
            )
        }
        missing_statistics <- setdiff(c(logfc_col, padj_col), colnames(result))
        if (length(missing_statistics) > 0L) {
            stop(
                "Analysis '", analysis_name, "' is missing columns: ",
                paste(missing_statistics, collapse = ", ")
            )
        }
        if (!is.numeric(result[[logfc_col]]) ||
            !is.numeric(result[[padj_col]])) {
            stop(
                "Analysis '", analysis_name, "' must contain numeric ",
                logfc_col, " and ", padj_col, " columns."
            )
        }
        observed_padj <- result[[padj_col]][!is.na(result[[padj_col]])]
        if (any(!is.finite(observed_padj)) ||
            any(observed_padj < 0 | observed_padj > 1)) {
            stop(
                "Analysis '", analysis_name, "' contains invalid ",
                padj_col, " values; observed values must be finite and in [0, 1]."
            )
        }

        scalar_from_object <- function(field, fallback = NA_character_) {
            if (is_pairwise_object && field %in% names(object)) {
                value <- as.character(object[[field]])
                if (length(value) == 1L && !is.na(value) && nzchar(value)) {
                    return(value)
                }
            }
            if (field %in% colnames(result)) {
                value <- unique(as.character(result[[field]]))
                value <- value[!is.na(value) & nzchar(value)]
                if (length(value) == 1L) {
                    return(value)
                }
            }
            fallback
        }
        reference <- scalar_from_object("group1")
        if (is.na(reference)) {
            reference <- scalar_from_object("Reference")
        }
        comparison_group <- scalar_from_object("group2")
        if (is.na(comparison_group)) {
            comparison_group <- scalar_from_object("ComparisonGroup")
        }
        if (!allow_data_frames &&
            (is.na(reference) || is.na(comparison_group))) {
            stop(
                "Analysis '", analysis_name,
                "' must contain scalar group1 and group2 identifiers."
            )
        }
        comparison <- scalar_from_object("comparison")
        if (is.na(comparison)) {
            comparison <- scalar_from_object("Comparison")
        }
        if (is.na(comparison) && !is.na(reference) &&
            !is.na(comparison_group)) {
            comparison <- paste0(comparison_group, "_vs_", reference)
        }
        if (is.na(comparison)) {
            comparison <- analysis_name
        }

        parts <- strsplit(analysis_name, name_sep, fixed = TRUE)[[1]]
        sorting <- parts[1]
        result$Analysis <- analysis_name
        result$Reference <- reference
        result$ComparisonGroup <- comparison_group
        result$Comparison <- comparison
        result$Sorting <- sorting

        plotting_p <- result[[padj_col]]
        positive_p <- plotting_p[
            !is.na(plotting_p) & is.finite(plotting_p) & plotting_p > 0
        ]
        replacement <- if (length(positive_p) > 0L) {
            max(min(positive_p) / 10, .Machine$double.xmin)
        } else {
            .Machine$double.xmin
        }
        plotting_p[is.na(plotting_p)] <- 1
        plotting_p[plotting_p <= 0] <- replacement
        result$SignedSignificance <-
            -log10(plotting_p) * sign(result[[logfc_col]])
        result$GeneRowName <- rownames(result)
        rownames(result) <- NULL
        result
    })

    all_columns <- unique(unlist(lapply(collected, colnames), use.names = FALSE))
    collected <- lapply(collected, function(result) {
        missing_columns <- setdiff(all_columns, colnames(result))
        for (column in missing_columns) {
            result[[column]] <- NA
        }
        result[, all_columns, drop = FALSE]
    })
    output <- do.call(rbind, collected)
    rownames(output) <- NULL
    output
}


#' Collect named edgeR pairwise analyses into one result table
#'
#' Combine multiple [edgeR_pairwise()] return objects while recording analysis,
#' reference, comparison, and sorting identifiers. The function also calculates
#' a signed significance score, `-log10(FDR) * sign(logFC)`, for downstream
#' visualization. Positive values indicate higher expression in the comparison
#' group and negative values indicate higher expression in the reference group.
#'
#' @param de_list Non-empty named list of [edgeR_pairwise()] return objects.
#'   Every object must contain a data-frame `results` element with numeric `FDR`
#'   and `logFC` columns plus scalar `group1` and `group2` identifiers.
#' @param name_sep Non-empty character scalar separating fields in each analysis
#'   name. Text before the first separator is stored as `Sorting`; when the
#'   separator is absent, the complete analysis name is used.
#'
#' @return A data frame formed from all input result tables. Existing result
#'   columns are retained and the following columns are added or refreshed:
#'   `Analysis`, `Reference`, `ComparisonGroup`, `Comparison`, `Sorting`,
#'   `SignedSignificance`, and `GeneRowName`. Rows retain list and within-table
#'   order, and row names are reset.
#'
#' @details
#' FDR values must be missing or finite values in `[0, 1]`. Missing FDR values
#' are assigned a plotting value of one, while exact zeros use one tenth of the
#' smallest positive FDR across their own analysis, bounded below by machine
#' precision. These substitutions affect only `SignedSignificance`; the
#' original FDR column is unchanged. Result tables with different annotation
#' columns are combined by taking their union and filling absent fields with
#' `NA`.
#'
#' @examples
#' \dontrun{
#' combined <- collect_edgeR_pairwise(
#'     list(
#'         Myeloid__Treated_vs_Control = myeloid_result,
#'         Lymphoid__Treated_vs_Control = lymphoid_result
#'     )
#' )
#' head(combined)
#' }
#'
#' @export
collect_edgeR_pairwise <- function(de_list, name_sep = "__") {
    if (!is.character(name_sep) || length(name_sep) != 1L ||
        is.na(name_sep) || !nzchar(name_sep)) {
        stop("name_sep must be one non-empty character string.")
    }
    output <- .collect_pairwise_results(
        de_list = de_list,
        name_sep = name_sep,
        logfc_col = "logFC",
        padj_col = "FDR",
        allow_data_frames = FALSE,
        require_names = TRUE
    )
    output
}


# =============================================================================
# Normalization and feature selection
# =============================================================================


# =============================================================================
# Dimensionality reduction and clustering
# =============================================================================


# =============================================================================
# Visualization and export
# =============================================================================


#' Plot a principal-component analysis of bulk or pseudobulk RNA-seq samples
#'
#' Prepare a validated featureCounts project for exploratory sample-level PCA,
#' optionally filter low-expression genes with [edgeR::filterByExpr()], apply
#' TMM library normalization, select the most variable log2-CPM genes, and
#' create a publication-style `ggplot2` scatter plot. Optional shapes, convex
#' hulls, sample labels, custom factor orders, palettes, and figure export are
#' supported.
#'
#' @param dat A `featurecounts_project`-like named list containing `counts`, a
#'   numeric gene-by-sample raw-count matrix, and `metadata`, a data frame whose
#'   row names are identical to the count-matrix column names and order.
#' @param color_by Character scalar. Metadata column used to colour samples and,
#'   when `filter_genes = TRUE`, define biological groups for
#'   [edgeR::filterByExpr()].
#' @param shape_by `NULL` or a character scalar. Optional metadata column used
#'   to map plotting symbols.
#' @param subset An optional unquoted logical expression evaluated within
#'   `dat$metadata`, for example `Tissue == "Lung" & Condition != "Control"`.
#'   Missing results are treated as `FALSE`. The default retains every sample.
#' @param n_variable_genes Positive integer. Maximum number of genes with the
#'   largest variance in log2 CPM to use as PCA features.
#' @param filter_genes Logical scalar. Whether to remove low-expression genes
#'   using [edgeR::filterByExpr()] with `color_by` as its grouping factor.
#' @param normalize Logical scalar. Whether to calculate TMM normalization
#'   factors and use normalized library sizes for log2-CPM calculation.
#' @param pc_x,pc_y Distinct positive integers identifying the principal
#'   components to plot on the horizontal and vertical axes.
#' @param point_size Positive numeric scalar controlling point size.
#' @param point_alpha Numeric scalar between zero and one controlling point
#'   transparency.
#' @param hull Logical scalar. Whether to draw filled convex hulls around colour
#'   groups having at least three samples and three unique coordinate pairs.
#' @param hull_alpha Numeric scalar between zero and one controlling convex-hull
#'   fill transparency.
#' @param hull_linewidth Non-negative numeric scalar controlling convex-hull
#'   border width.
#' @param label_samples Logical scalar. Whether to label individual samples with
#'   `ggrepel`; this option requires the `ggrepel` package.
#' @param label_col Character scalar. Column in the assembled PCA data used for
#'   labels. The default `"SampleName"` falls back to count-matrix sample names
#'   when that metadata column is absent.
#' @param color_order `NULL` or a unique character vector giving the desired
#'   order of colour groups. Every observed group must be included.
#' @param shape_order `NULL` or a unique character vector giving the desired
#'   order of shape groups. Every observed group must be included.
#' @param colors `NULL` or a character vector of colours. An unnamed vector is
#'   assigned in `color_order`; a named vector must contain every observed
#'   colour group. By default, `BULK_PCA_MACARON_COLORS` is used for up to 15
#'   groups and `grDevices::hcl.colors(..., palette = "Pastel 1")` thereafter.
#' @param shapes `NULL` or a vector of plotting symbols. An unnamed vector is
#'   assigned in `shape_order`; a named vector must contain every observed shape
#'   group. The default supports up to eight groups.
#' @param title `NULL` or a character scalar used as the plot title.
#' @param width,height Positive numeric scalars giving saved figure dimensions
#'   in inches.
#' @param save `NULL` or a character scalar giving an output figure path. PDF
#'   files use `grDevices::cairo_pdf`; other formats are written at 300 dpi.
#' @param verbose Logical scalar. Whether to print a concise PCA summary.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{plot}{The assembled `ggplot` object.}
#'     \item{pca}{The fitted `prcomp` object.}
#'     \item{pca_data}{Sample scores joined to aligned sample metadata.}
#'     \item{logCPM}{The filtered gene-by-sample log2-CPM matrix.}
#'     \item{variable_genes}{Gene identifiers used as PCA features.}
#'     \item{gene_variance}{Named, decreasingly sorted log2-CPM variances for
#'       all retained genes.}
#'     \item{variance_explained}{Percentage of variance explained by each
#'       available principal component.}
#'     \item{dge}{The filtered and optionally TMM-normalized `DGEList`.}
#'     \item{keep_genes}{Logical vector aligned to the original count rows,
#'       identifying genes retained by expression filtering.}
#'     \item{colors}{Named colour vector used in the plot.}
#'     \item{shapes}{Named shape vector used in the plot, or `NULL`.}
#'   }
#'
#' @details
#' PCA is exploratory and sample-level: supplied counts must represent bulk
#' libraries or replicate-aware pseudobulk profiles, not individual cells.
#' The function does not modify `dat`. If `save` is supplied, the figure is
#' written as a side effect. At least three selected samples and two retained
#' genes are required.
#'
#' @examples
#' \dontrun{
#' pca_result <- plot_bulk_pca(
#'     project,
#'     color_by = "Genotype",
#'     shape_by = "Batch",
#'     subset = Tissue == "Lung",
#'     color_order = c("WT", "cKO"),
#'     hull = TRUE,
#'     label_samples = TRUE,
#'     save = "outputs/lung_pca.pdf"
#' )
#' pca_result$plot
#' }
#'
#' @export
plot_bulk_pca <- function(
    dat,
    color_by = "Sorting",
    shape_by = NULL,
    subset = NULL,
    n_variable_genes = 5000,
    filter_genes = TRUE,
    normalize = TRUE,
    pc_x = 1,
    pc_y = 2,
    point_size = 3.2,
    point_alpha = 0.85,
    hull = FALSE,
    hull_alpha = 0.10,
    hull_linewidth = 0.6,
    label_samples = FALSE,
    label_col = "SampleName",
    color_order = NULL,
    shape_order = NULL,
    colors = NULL,
    shapes = NULL,
    title = NULL,
    width = 6,
    height = 5,
    save = NULL,
    verbose = TRUE
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }
    if (isTRUE(label_samples) &&
        !requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required when label_samples = TRUE.")
    }

    required <- c("counts", "metadata")
    missing_required <- setdiff(required, names(dat))
    if (length(missing_required) > 0L) {
        stop("dat is missing: ", paste(missing_required, collapse = ", "))
    }

    counts <- dat$counts
    meta <- dat$metadata
    if (!is.matrix(counts) || !is.numeric(counts)) {
        stop("dat$counts must be a numeric matrix.")
    }
    if (!is.data.frame(meta)) {
        stop("dat$metadata must be a data frame.")
    }
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
        anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
        any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
        anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
        stop("dat$counts must have unique, non-empty gene and sample names.")
    }
    if (!identical(colnames(counts), rownames(meta))) {
        stop(
            "dat$counts columns and dat$metadata rows are not in identical ",
            "sample order."
        )
    }
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
        stop("dat$counts must contain finite, non-negative values without NA.")
    }

    scalar_character <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    scalar_logical <- function(x) {
        is.logical(x) && length(x) == 1L && !is.na(x)
    }
    if (!scalar_character(color_by) || !color_by %in% colnames(meta)) {
        stop("color_by must name one metadata column.")
    }
    if (!is.null(shape_by) &&
        (!scalar_character(shape_by) || !shape_by %in% colnames(meta))) {
        stop("shape_by must be NULL or name one metadata column.")
    }
    if (!scalar_character(label_col)) {
        stop("label_col must be one non-empty character string.")
    }
    if (!is.null(title) && !scalar_character(title)) {
        stop("title must be NULL or one non-empty character string.")
    }
    if (!is.null(save) && !scalar_character(save)) {
        stop("save must be NULL or one non-empty character path.")
    }
    logical_arguments <- list(
        filter_genes = filter_genes,
        normalize = normalize,
        hull = hull,
        label_samples = label_samples,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[
        !vapply(logical_arguments, scalar_logical, logical(1))
    ]
    if (length(invalid_logical) > 0L) {
        stop(
            paste(invalid_logical, collapse = ", "),
            " must be non-missing logical scalars."
        )
    }
    if (length(n_variable_genes) != 1L || is.na(n_variable_genes) ||
        !is.numeric(n_variable_genes) || !is.finite(n_variable_genes) ||
        n_variable_genes < 1 ||
        n_variable_genes != round(n_variable_genes)) {
        stop("n_variable_genes must be a positive integer.")
    }
    for (pc_name in c("pc_x", "pc_y")) {
        pc_value <- get(pc_name)
        if (length(pc_value) != 1L || is.na(pc_value) ||
            !is.numeric(pc_value) || !is.finite(pc_value) || pc_value < 1 ||
            pc_value != round(pc_value)) {
            stop(pc_name, " must be a positive integer.")
        }
    }
    if (pc_x == pc_y) {
        stop("pc_x and pc_y must identify different principal components.")
    }
    numeric_values <- list(
        point_size = point_size,
        point_alpha = point_alpha,
        hull_alpha = hull_alpha,
        hull_linewidth = hull_linewidth,
        width = width,
        height = height
    )
    numeric_ranges <- list(
        point_size = c(0, Inf),
        point_alpha = c(0, 1),
        hull_alpha = c(0, 1),
        hull_linewidth = c(0, Inf),
        width = c(0, Inf),
        height = c(0, Inf)
    )
    for (argument_name in names(numeric_values)) {
        value <- numeric_values[[argument_name]]
        supported_range <- numeric_ranges[[argument_name]]
        lower <- supported_range[1]
        upper <- supported_range[2]
        if (length(value) != 1L || is.na(value) || !is.finite(value) ||
            value < lower || value > upper ||
            (argument_name %in% c("point_size", "width", "height") &&
                value == 0)) {
            stop(argument_name, " is outside its supported numeric range.")
        }
    }

    subset_expr <- substitute(subset)
    if (!identical(subset_expr, quote(NULL))) {
        keep_samples <- eval(subset_expr, envir = meta, enclos = parent.frame())
        if (!is.logical(keep_samples) || length(keep_samples) != nrow(meta)) {
            stop("subset must evaluate to one logical value per metadata row.")
        }
        keep_samples[is.na(keep_samples)] <- FALSE
        if (!any(keep_samples)) {
            stop("No samples remain after applying subset.")
        }
        meta <- meta[keep_samples, , drop = FALSE]
        counts <- counts[, rownames(meta), drop = FALSE]
    }
    if (ncol(counts) < 3L) {
        stop("At least 3 samples are required for PCA.")
    }
    if (anyNA(meta[[color_by]])) {
        stop("color_by contains missing values in the selected samples.")
    }
    if (!is.null(shape_by) && anyNA(meta[[shape_by]])) {
        stop("shape_by contains missing values in the selected samples.")
    }

    make_ordered_factor <- function(values, requested_order, argument_name) {
        observed <- unique(as.character(values))
        if (any(!nzchar(observed))) {
            stop(argument_name, " contains empty group labels.")
        }
        if (is.null(requested_order)) {
            return(factor(as.character(values), levels = observed))
        }
        if (!is.character(requested_order) || anyNA(requested_order) ||
            any(!nzchar(requested_order)) || anyDuplicated(requested_order)) {
            stop(argument_name, "_order must contain unique, non-empty labels.")
        }
        missing_levels <- setdiff(observed, requested_order)
        if (length(missing_levels) > 0L) {
            stop(
                argument_name,
                "_order is missing observed groups: ",
                paste(missing_levels, collapse = ", ")
            )
        }
        factor(as.character(values), levels = requested_order)
    }
    if (!is.null(shape_by) && identical(shape_by, color_by)) {
        if (!is.null(color_order) && !is.null(shape_order) &&
            !identical(color_order, shape_order)) {
            stop(
                "color_order and shape_order must be identical when color_by ",
                "and shape_by name the same metadata column."
            )
        }
        shared_order <- if (!is.null(color_order)) color_order else shape_order
        meta[[color_by]] <- make_ordered_factor(
            meta[[color_by]], shared_order, "color/shape"
        )
    } else {
        meta[[color_by]] <- make_ordered_factor(
            meta[[color_by]], color_order, "color"
        )
    }
    if (!is.null(shape_by) && !identical(shape_by, color_by)) {
        meta[[shape_by]] <- make_ordered_factor(
            meta[[shape_by]], shape_order, "shape"
        )
    }

    dge <- edgeR::DGEList(counts = counts)
    n_genes_before <- nrow(dge)
    if (filter_genes) {
        keep_genes <- edgeR::filterByExpr(dge, group = meta[[color_by]])
        dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    } else {
        keep_genes <- rep(TRUE, nrow(dge))
        names(keep_genes) <- rownames(dge)
    }
    n_genes_after <- nrow(dge)
    if (n_genes_after < 2L) {
        stop("Too few genes remain after filtering.")
    }
    if (normalize) {
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
    }
    logcpm <- edgeR::cpm(
        dge,
        log = TRUE,
        prior.count = 2,
        normalized.lib.sizes = normalize
    )

    gene_variance <- apply(logcpm, 1L, stats::var, na.rm = TRUE)
    gene_variance[!is.finite(gene_variance)] <- 0
    gene_variance <- sort(gene_variance, decreasing = TRUE)
    n_variable_genes <- min(as.integer(n_variable_genes), length(gene_variance))
    variable_genes <- names(gene_variance)[seq_len(n_variable_genes)]
    pca_input <- logcpm[variable_genes, , drop = FALSE]
    pca <- stats::prcomp(t(pca_input), center = TRUE, scale. = FALSE)
    total_component_variance <- sum(pca$sdev^2)
    if (!is.finite(total_component_variance) || total_component_variance <= 0) {
        stop("PCA has no positive finite variance to display.")
    }
    variance_explained <- 100 * pca$sdev^2 / total_component_variance
    if (max(pc_x, pc_y) > ncol(pca$x)) {
        stop(
            "Requested PC does not exist. Available PCs: 1-",
            ncol(pca$x)
        )
    }

    pca_df <- as.data.frame(pca$x, check.names = FALSE)
    pca_df$PCA_sample_name <- rownames(pca_df)
    meta_tmp <- meta[pca_df$PCA_sample_name, , drop = FALSE]
    score_name_collisions <- intersect(colnames(meta_tmp), colnames(pca$x))
    if (length(score_name_collisions) > 0L) {
        stop(
            "Metadata column names conflict with PCA score columns: ",
            paste(score_name_collisions, collapse = ", ")
        )
    }
    meta_tmp$PCA_sample_name <- NULL
    pca_df <- cbind(pca_df, meta_tmp)

    color_levels <- levels(droplevels(meta[[color_by]]))
    if (is.null(colors)) {
        if (length(color_levels) <= length(BULK_PCA_MACARON_COLORS)) {
            colors <- BULK_PCA_MACARON_COLORS[seq_along(color_levels)]
        } else {
            colors <- grDevices::hcl.colors(
                length(color_levels), palette = "Pastel 1"
            )
        }
        names(colors) <- color_levels
    } else {
        if (!is.character(colors) || length(colors) == 0L || anyNA(colors)) {
            stop("colors must be a non-empty vector without missing values.")
        }
        if (is.null(names(colors))) {
            if (length(colors) < length(color_levels)) {
                stop("Not enough colors supplied. Need ", length(color_levels), ".")
            }
            colors <- colors[seq_along(color_levels)]
            names(colors) <- color_levels
        } else {
            if (anyDuplicated(names(colors))) {
                stop("Named colors must have unique names.")
            }
            missing_colors <- setdiff(color_levels, names(colors))
            if (length(missing_colors) > 0L) {
                stop("Missing colors for: ", paste(missing_colors, collapse = ", "))
            }
            colors <- colors[color_levels]
        }
    }

    shape_levels <- NULL
    if (!is.null(shape_by)) {
        shape_levels <- levels(droplevels(meta[[shape_by]]))
        if (is.null(shapes)) {
            default_shapes <- c(16, 17, 15, 18, 8, 3, 7, 4)
            if (length(shape_levels) > length(default_shapes)) {
                stop("Too many shape groups; supply shapes manually.")
            }
            shapes <- default_shapes[seq_along(shape_levels)]
            names(shapes) <- shape_levels
        } else {
            if (!is.atomic(shapes) || length(shapes) == 0L || anyNA(shapes)) {
                stop("shapes must be a non-empty vector without missing values.")
            }
            if (is.null(names(shapes))) {
                if (length(shapes) < length(shape_levels)) {
                    stop("Not enough shapes supplied. Need ", length(shape_levels), ".")
                }
                shapes <- shapes[seq_along(shape_levels)]
                names(shapes) <- shape_levels
            } else {
                if (anyDuplicated(names(shapes))) {
                    stop("Named shapes must have unique names.")
                }
                missing_shapes <- setdiff(shape_levels, names(shapes))
                if (length(missing_shapes) > 0L) {
                    stop("Missing shapes for: ", paste(missing_shapes, collapse = ", "))
                }
                shapes <- shapes[shape_levels]
            }
        }
    }

    x_col <- paste0("PC", pc_x)
    y_col <- paste0("PC", pc_y)
    x_lab <- paste0(x_col, " (", round(variance_explained[pc_x], 1), "%)")
    y_lab <- paste0(y_col, " (", round(variance_explained[pc_y], 1), "%)")
    p <- ggplot2::ggplot(
        pca_df,
        ggplot2::aes(
            x = .data[[x_col]],
            y = .data[[y_col]],
            colour = .data[[color_by]]
        )
    )

    if (hull) {
        hull_list <- lapply(
            split(pca_df, pca_df[[color_by]], drop = TRUE),
            function(df) {
                if (nrow(df) < 3L) {
                    return(NULL)
                }
                coordinates <- unique(df[, c(x_col, y_col), drop = FALSE])
                if (nrow(coordinates) < 3L) {
                    return(NULL)
                }
                df[grDevices::chull(df[[x_col]], df[[y_col]]), , drop = FALSE]
            }
        )
        hull_list <- hull_list[
            !vapply(hull_list, is.null, logical(1))
        ]
        if (length(hull_list) > 0L) {
            hull_df <- do.call(rbind, hull_list)
            rownames(hull_df) <- NULL
            p <- p +
                ggplot2::geom_polygon(
                    data = hull_df,
                    mapping = ggplot2::aes(
                        x = .data[[x_col]],
                        y = .data[[y_col]],
                        group = .data[[color_by]],
                        fill = .data[[color_by]]
                    ),
                    inherit.aes = FALSE,
                    alpha = hull_alpha,
                    colour = NA,
                    show.legend = FALSE
                ) +
                ggplot2::geom_polygon(
                    data = hull_df,
                    mapping = ggplot2::aes(
                        x = .data[[x_col]],
                        y = .data[[y_col]],
                        group = .data[[color_by]],
                        colour = .data[[color_by]]
                    ),
                    inherit.aes = FALSE,
                    fill = NA,
                    linewidth = hull_linewidth,
                    show.legend = FALSE
                )
        }
    }

    if (is.null(shape_by)) {
        p <- p + ggplot2::geom_point(
            size = point_size, alpha = point_alpha, stroke = 0.4
        )
    } else {
        p <- p + ggplot2::geom_point(
            mapping = ggplot2::aes(shape = .data[[shape_by]]),
            size = point_size,
            alpha = point_alpha,
            stroke = 0.6
        )
    }

    if (label_samples) {
        if (!label_col %in% colnames(pca_df)) {
            if (identical(label_col, "SampleName")) {
                pca_df$SampleName <- pca_df$PCA_sample_name
            } else {
                stop("label_col '", label_col, "' was not found in metadata.")
            }
        }
        p <- p + ggrepel::geom_text_repel(
            data = pca_df,
            mapping = ggplot2::aes(
                x = .data[[x_col]],
                y = .data[[y_col]],
                label = .data[[label_col]]
            ),
            inherit.aes = FALSE,
            size = 2.6,
            box.padding = 0.35,
            point.padding = 0.25,
            segment.linewidth = 0.3,
            max.overlaps = Inf,
            show.legend = FALSE
        )
    }

    p <- p + ggplot2::scale_colour_manual(
        values = colors, breaks = color_levels, drop = FALSE
    )
    if (hull) {
        p <- p + ggplot2::scale_fill_manual(
            values = colors, breaks = color_levels, drop = FALSE
        )
    }
    if (!is.null(shape_by)) {
        p <- p + ggplot2::scale_shape_manual(
            values = shapes, breaks = shape_levels, drop = FALSE
        )
    }

    p <- p +
        ggplot2::labs(
            x = x_lab,
            y = y_lab,
            colour = color_by,
            shape = shape_by,
            title = title
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(linewidth = 0.55, colour = "black"),
            axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "black"),
            axis.ticks.length = grid::unit(2, "mm"),
            axis.text = ggplot2::element_text(colour = "black", size = 9),
            axis.title = ggplot2::element_text(colour = "black", size = 10),
            plot.title = ggplot2::element_text(
                hjust = 0.5, size = 11, face = "bold", colour = "black"
            ),
            legend.title = ggplot2::element_text(size = 9, colour = "black"),
            legend.text = ggplot2::element_text(size = 8, colour = "black"),
            legend.key = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
        )

    if (!is.null(save)) {
        if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
            ggplot2::ggsave(
                filename = save,
                plot = p,
                width = width,
                height = height,
                units = "in",
                device = grDevices::cairo_pdf
            )
        } else {
            ggplot2::ggsave(
                filename = save,
                plot = p,
                width = width,
                height = height,
                units = "in",
                dpi = 300
            )
        }
    }

    if (verbose) {
        message("")
        message("========================================")
        message(" Bulk RNA-seq PCA")
        message("========================================")
        message("Samples:                ", ncol(counts))
        message("Genes before filtering: ", n_genes_before)
        message("Genes after filtering:  ", n_genes_after)
        message("Genes used for PCA:     ", nrow(pca_input))
        message(x_col, ":                    ", round(variance_explained[pc_x], 1), "%")
        message(y_col, ":                    ", round(variance_explained[pc_y], 1), "%")
        message("Color:                  ", color_by)
        if (!is.null(shape_by)) {
            message("Shape:                  ", shape_by)
        }
        message("========================================")
        message("")
    }

    list(
        plot = p,
        pca = pca,
        pca_data = pca_df,
        logCPM = logcpm,
        variable_genes = variable_genes,
        gene_variance = gene_variance,
        variance_explained = variance_explained,
        dge = dge,
        keep_genes = keep_genes,
        colors = colors,
        shapes = if (!is.null(shape_by)) shapes else NULL
    )
}


#' Plot sample-level quality metrics for a bulk RNA-seq project
#'
#' Calculate count-matrix library sizes, validate requested sample-level QC
#' metrics, and create one faceted point plot per metric. Plots are combined
#' with a shared colour legend and can optionally include within-facet group
#' lines, sample labels, custom group ordering, custom colours, and figure
#' export.
#'
#' @param dat A `featurecounts_project`-like named list containing `counts`, a
#'   numeric gene-by-sample count matrix, and `metadata`, a data frame whose row
#'   names are identical to the count-matrix column names and order.
#' @param metrics Non-empty character vector naming numeric metadata columns to
#'   plot. `"LibrarySize"` is calculated from `colSums(dat$counts)` and replaces
#'   any existing metadata column with that name. Duplicate metric names are
#'   plotted once, in order of first appearance.
#' @param color_by Character scalar naming a categorical metadata column used
#'   to colour points and optional connecting lines.
#' @param facet_by Character scalar naming a categorical metadata column used
#'   to create sample panels within every metric plot.
#' @param sample_col Character scalar naming the sample-label metadata column.
#'   If absent, it is created from metadata row names. Values must be unique and
#'   non-missing. The reserved name `"LibrarySize"` is not supported here.
#' @param color_order `NULL` or a unique character vector giving the desired
#'   colour-group order. Every observed group must be included.
#' @param colors `NULL` or a character vector of colours. An unnamed vector is
#'   assigned in `color_order`; a named vector must contain every observed
#'   colour group. By default, `BULK_QC_MACARON_COLORS` is used for up to 12
#'   groups and `grDevices::hcl.colors(..., palette = "Pastel 1")` thereafter.
#' @param point_size Positive numeric scalar controlling point size.
#' @param point_alpha Numeric scalar between zero and one controlling point
#'   transparency.
#' @param label_samples Logical scalar. Whether to label points using
#'   `ggrepel`; this option requires the `ggrepel` package.
#' @param label_size Positive numeric scalar controlling sample-label text size.
#' @param connect_samples Logical scalar. Whether to connect samples belonging
#'   to the same colour group within each facet, following the displayed sample
#'   order. Use only when that ordering has a meaningful interpretation.
#' @param ncol Positive integer giving the number of columns in the combined
#'   patchwork layout.
#' @param width,height Positive numeric scalars giving saved figure dimensions
#'   in inches.
#' @param save `NULL` or a character scalar giving an output figure path. PDF
#'   files use `grDevices::cairo_pdf`; other formats are written at 300 dpi.
#' @param verbose Logical scalar. Whether to print the sample count and the
#'   observed assignment-rate range when available.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{plot}{The combined `patchwork` figure.}
#'     \item{plots}{A named list of individual `ggplot` objects, one per unique
#'       requested metric.}
#'     \item{metadata}{A sorted copy of the aligned metadata containing the
#'       calculated `LibrarySize` column and plotting factor levels.}
#'     \item{colors}{The named colour vector used in every panel.}
#'   }
#'
#' @details
#' This function summarizes available featureCounts and count-matrix QC fields;
#' it does not alter `dat` or apply exclusion thresholds. Missing metric values
#' are permitted, but each requested metric must contain at least one finite
#' value. If `save` is supplied, the combined figure is written as a side
#' effect. Required packages are `ggplot2`, `patchwork`, and `scales`, with
#' `ggrepel` additionally required for labels.
#'
#' @examples
#' \dontrun{
#' qc_result <- plot_bulk_qc(
#'     project,
#'     color_by = "Condition",
#'     facet_by = "Sorting",
#'     color_order = c("Control", "Treated"),
#'     label_samples = TRUE,
#'     save = "outputs/bulk_sample_qc.pdf"
#' )
#' qc_result$plot
#' qc_result$metadata
#' }
#'
#' @export
plot_bulk_qc <- function(
    dat,
    metrics = c(
        "Total_input_reads",
        "Assigned_reads",
        "LibrarySize",
        "Assignment_percent"
    ),
    color_by = "Condition",
    facet_by = "Sorting",
    sample_col = "SampleName",
    color_order = NULL,
    colors = NULL,
    point_size = 3.0,
    point_alpha = 0.9,
    label_samples = FALSE,
    label_size = 2.5,
    connect_samples = FALSE,
    ncol = 2,
    width = 10,
    height = 7,
    save = NULL,
    verbose = TRUE
) {
    required_packages <- c("ggplot2", "patchwork", "scales")
    missing_packages <- required_packages[
        !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_packages) > 0L) {
        stop(
            "Required package(s) not installed: ",
            paste(missing_packages, collapse = ", ")
        )
    }
    if (isTRUE(label_samples) &&
        !requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required when label_samples = TRUE.")
    }

    required_elements <- c("counts", "metadata")
    missing_elements <- setdiff(required_elements, names(dat))
    if (length(missing_elements) > 0L) {
        stop("dat is missing: ", paste(missing_elements, collapse = ", "))
    }
    counts <- dat$counts
    meta <- dat$metadata
    if (!is.matrix(counts) || !is.numeric(counts)) {
        stop("dat$counts must be a numeric matrix.")
    }
    if (!is.data.frame(meta)) {
        stop("dat$metadata must be a data frame.")
    }
    if (is.null(colnames(counts)) || anyNA(colnames(counts)) ||
        any(!nzchar(colnames(counts))) || anyDuplicated(colnames(counts))) {
        stop("dat$counts must have unique, non-empty sample column names.")
    }
    if (!identical(colnames(counts), rownames(meta))) {
        stop(
            "dat$counts columns and dat$metadata rows are not in identical ",
            "sample order."
        )
    }
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
        stop("dat$counts must contain finite, non-negative values without NA.")
    }
    if (base::ncol(counts) == 0L) {
        stop("dat must contain at least one sample.")
    }

    scalar_character <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    scalar_logical <- function(x) {
        is.logical(x) && length(x) == 1L && !is.na(x)
    }
    if (!is.character(metrics) || length(metrics) == 0L || anyNA(metrics) ||
        any(!nzchar(metrics))) {
        stop("metrics must be a non-empty character vector without missing values.")
    }
    metrics <- unique(metrics)
    column_arguments <- list(
        color_by = color_by,
        facet_by = facet_by,
        sample_col = sample_col
    )
    invalid_columns <- names(column_arguments)[
        !vapply(column_arguments, scalar_character, logical(1))
    ]
    if (length(invalid_columns) > 0L) {
        stop(
            paste(invalid_columns, collapse = ", "),
            " must be non-empty character scalars."
        )
    }
    if (identical(sample_col, "LibrarySize")) {
        stop("sample_col cannot use the reserved calculated column 'LibrarySize'.")
    }
    if (!is.null(save) && !scalar_character(save)) {
        stop("save must be NULL or one non-empty character path.")
    }
    logical_arguments <- list(
        label_samples = label_samples,
        connect_samples = connect_samples,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[
        !vapply(logical_arguments, scalar_logical, logical(1))
    ]
    if (length(invalid_logical) > 0L) {
        stop(
            paste(invalid_logical, collapse = ", "),
            " must be non-missing logical scalars."
        )
    }
    numeric_values <- list(
        point_size = point_size,
        point_alpha = point_alpha,
        label_size = label_size,
        width = width,
        height = height
    )
    numeric_ranges <- list(
        point_size = c(0, Inf),
        point_alpha = c(0, 1),
        label_size = c(0, Inf),
        width = c(0, Inf),
        height = c(0, Inf)
    )
    for (argument_name in names(numeric_values)) {
        value <- numeric_values[[argument_name]]
        supported_range <- numeric_ranges[[argument_name]]
        if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
            !is.finite(value) || value < supported_range[1] ||
            value > supported_range[2] ||
            (argument_name != "point_alpha" && value == 0)) {
            stop(argument_name, " is outside its supported numeric range.")
        }
    }
    if (length(ncol) != 1L || !is.numeric(ncol) || is.na(ncol) ||
        !is.finite(ncol) || ncol < 1 || ncol != round(ncol)) {
        stop("ncol must be a positive integer.")
    }
    ncol <- as.integer(ncol)

    if (!sample_col %in% colnames(meta)) {
        meta[[sample_col]] <- rownames(meta)
    }
    meta[[sample_col]] <- as.character(meta[[sample_col]])
    if (anyNA(meta[[sample_col]]) || any(!nzchar(meta[[sample_col]])) ||
        anyDuplicated(meta[[sample_col]])) {
        stop("sample_col must contain unique, non-empty sample identifiers.")
    }
    meta$LibrarySize <- as.numeric(colSums(counts))

    required_columns <- unique(c(metrics, color_by, facet_by, sample_col))
    missing_columns <- setdiff(required_columns, colnames(meta))
    if (length(missing_columns) > 0L) {
        stop(
            "These variables are missing from metadata: ",
            paste(missing_columns, collapse = ", ")
        )
    }
    if (identical(color_by, "LibrarySize") ||
        identical(facet_by, "LibrarySize")) {
        stop("color_by and facet_by must name categorical metadata columns.")
    }
    metric_group_overlap <- intersect(metrics, c(color_by, facet_by))
    if (length(metric_group_overlap) > 0L) {
        stop(
            "QC metrics cannot also be used as categorical grouping columns: ",
            paste(metric_group_overlap, collapse = ", ")
        )
    }
    if (anyNA(meta[[color_by]]) || anyNA(meta[[facet_by]])) {
        stop("color_by and facet_by cannot contain missing values.")
    }
    if (any(!nzchar(as.character(meta[[color_by]]))) ||
        any(!nzchar(as.character(meta[[facet_by]])))) {
        stop("color_by and facet_by cannot contain empty group labels.")
    }
    invalid_metrics <- metrics[
        !vapply(meta[metrics], is.numeric, logical(1))
    ]
    if (length(invalid_metrics) > 0L) {
        stop(
            "Requested metrics must be numeric: ",
            paste(invalid_metrics, collapse = ", ")
        )
    }
    nonfinite_metrics <- metrics[
        vapply(
            meta[metrics],
            function(x) any(!is.finite(x[!is.na(x)])),
            logical(1)
        )
    ]
    if (length(nonfinite_metrics) > 0L) {
        stop(
            "Requested metrics contain non-finite values: ",
            paste(nonfinite_metrics, collapse = ", ")
        )
    }
    empty_metrics <- metrics[
        vapply(meta[metrics], function(x) all(is.na(x)), logical(1))
    ]
    if (length(empty_metrics) > 0L) {
        stop(
            "Requested metrics contain no observed values: ",
            paste(empty_metrics, collapse = ", ")
        )
    }

    observed_colors <- unique(as.character(meta[[color_by]]))
    if (is.null(color_order)) {
        color_order <- observed_colors
    } else {
        if (!is.character(color_order) || anyNA(color_order) ||
            any(!nzchar(color_order)) || anyDuplicated(color_order)) {
            stop("color_order must contain unique, non-empty labels.")
        }
        missing_levels <- setdiff(observed_colors, color_order)
        if (length(missing_levels) > 0L) {
            stop(
                "color_order is missing observed groups: ",
                paste(missing_levels, collapse = ", ")
            )
        }
    }
    meta[[color_by]] <- factor(
        as.character(meta[[color_by]]), levels = color_order
    )
    color_levels <- levels(droplevels(meta[[color_by]]))

    if (is.null(colors)) {
        if (length(color_levels) <= length(BULK_QC_MACARON_COLORS)) {
            colors <- BULK_QC_MACARON_COLORS[seq_along(color_levels)]
        } else {
            colors <- grDevices::hcl.colors(
                length(color_levels), palette = "Pastel 1"
            )
        }
        names(colors) <- color_levels
    } else {
        if (!is.character(colors) || length(colors) == 0L || anyNA(colors)) {
            stop("colors must be a non-empty character vector without NA.")
        }
        if (is.null(names(colors))) {
            if (length(colors) < length(color_levels)) {
                stop("Not enough colors supplied. Need ", length(color_levels), ".")
            }
            colors <- colors[seq_along(color_levels)]
            names(colors) <- color_levels
        } else {
            if (anyDuplicated(names(colors))) {
                stop("Named colors must have unique names.")
            }
            missing_colors <- setdiff(color_levels, names(colors))
            if (length(missing_colors) > 0L) {
                stop("Missing colors for: ", paste(missing_colors, collapse = ", "))
            }
            colors <- colors[color_levels]
        }
    }

    meta <- meta[
        order(meta[[facet_by]], meta[[color_by]], meta[[sample_col]]),
        ,
        drop = FALSE
    ]
    cell_theme <- ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
            axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
            axis.text = ggplot2::element_text(colour = "black", size = 8),
            axis.title = ggplot2::element_text(colour = "black", size = 9),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(
                face = "bold", colour = "black", size = 9
            ),
            legend.title = ggplot2::element_text(size = 9),
            legend.text = ggplot2::element_text(size = 8),
            legend.key = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(
                face = "bold", size = 10, hjust = 0.5
            ),
            plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
        )
    metric_titles <- c(
        Total_input_reads = "Total input reads",
        Assigned_reads = "Assigned reads",
        LibrarySize = "Count matrix library size",
        Assignment_percent = "Assignment rate (%)"
    )

    plot_list <- list()
    for (metric in metrics) {
        plot_data <- meta
        plot_data$SamplePlot <- factor(
            plot_data[[sample_col]], levels = unique(plot_data[[sample_col]])
        )
        plot <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(
                x = SamplePlot,
                y = .data[[metric]],
                colour = .data[[color_by]]
            )
        )
        if (connect_samples) {
            plot <- plot + ggplot2::geom_line(
                ggplot2::aes(group = .data[[color_by]]),
                linewidth = 0.4,
                alpha = 0.45
            )
        }
        plot <- plot + ggplot2::geom_point(
            size = point_size,
            alpha = point_alpha,
            stroke = 0.3
        )
        if (label_samples) {
            plot <- plot + ggrepel::geom_text_repel(
                ggplot2::aes(label = .data[[sample_col]]),
                size = label_size,
                box.padding = 0.25,
                point.padding = 0.2,
                max.overlaps = Inf,
                show.legend = FALSE
            )
        }
        plot <- plot +
            ggplot2::facet_wrap(
                stats::reformulate(facet_by), scales = "free_x"
            ) +
            ggplot2::scale_colour_manual(
                values = colors, breaks = color_levels, drop = FALSE
            )
        if (metric %in% c(
            "Total_input_reads", "Assigned_reads", "LibrarySize"
        )) {
            plot <- plot + ggplot2::scale_y_continuous(
                labels = scales::label_number(
                    scale_cut = scales::cut_short_scale()
                ),
                expand = ggplot2::expansion(mult = c(0.05, 0.12))
            )
        } else {
            plot <- plot + ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.05, 0.12))
            )
        }
        pretty_title <- if (metric %in% names(metric_titles)) {
            metric_titles[[metric]]
        } else {
            metric
        }
        plot <- plot +
            ggplot2::labs(
                x = NULL,
                y = pretty_title,
                colour = color_by,
                title = pretty_title
            ) +
            cell_theme +
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )
        plot_list[[metric]] <- plot
    }

    combined <- patchwork::wrap_plots(
        plot_list, ncol = ncol, guides = "collect"
    ) & ggplot2::theme(legend.position = "right")
    if (!is.null(save)) {
        if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
            ggplot2::ggsave(
                filename = save,
                plot = combined,
                width = width,
                height = height,
                units = "in",
                device = grDevices::cairo_pdf
            )
        } else {
            ggplot2::ggsave(
                filename = save,
                plot = combined,
                width = width,
                height = height,
                units = "in",
                dpi = 300
            )
        }
    }

    if (verbose) {
        message("")
        message("======================================")
        message(" Bulk RNA-seq sample QC")
        message("======================================")
        message("Samples: ", nrow(meta))
        assignment_values <- if (
            "Assignment_percent" %in% colnames(meta) &&
            is.numeric(meta$Assignment_percent)
        ) {
            meta$Assignment_percent[is.finite(meta$Assignment_percent)]
        } else {
            numeric(0)
        }
        if (length(assignment_values) > 0L) {
            message(
                "Assignment rate: ",
                round(min(assignment_values), 1),
                "% - ",
                round(max(assignment_values), 1),
                "%"
            )
        }
        message("======================================")
    }

    list(
        plot = combined,
        plots = plot_list,
        metadata = meta,
        colors = colors
    )
}


#' Plot normalized bulk RNA-seq expression distributions
#'
#' Display selected genes from a prepared normalized-expression object as
#' faceted sample-level violin plots. Optional boxplots and jittered sample
#' points retain distribution summaries and individual biological-replicate
#' context without recalculating normalization for each figure.
#'
#' @param bulk A `bulk_expression_prepared` object returned by
#'   [prepare_bulk_expression()]. Its `logCPM`, `metadata`, `genes`, and `dge`
#'   components must remain aligned.
#' @param genes Non-empty character vector of gene symbols or stable gene IDs.
#'   Requested order is preserved. Missing genes produce warnings; when one
#'   request matches multiple annotation rows, the row with the highest mean
#'   log2 CPM is used.
#' @param group_by Character scalar naming the metadata column displayed on the
#'   x-axis and used for violin colours.
#' @param split_by `NULL` or a character scalar naming a metadata column used
#'   for facet columns. Genes form facet rows when this is supplied.
#' @param gene_col Character scalar naming the preferred gene-symbol column in
#'   `bulk$genes`. Stable IDs are used when this column is absent or blank.
#' @param gene_id_col Character scalar naming the preferred stable-ID column in
#'   `bulk$genes`. Annotation row names are used when this column is absent.
#' @param group_order `NULL` or a unique character vector defining x-axis group
#'   order. Observed groups omitted from the vector are appended in first-seen
#'   order; unobserved supplied values are dropped from the plotted levels.
#' @param split_order `NULL` or a unique character vector defining split-panel
#'   order. Observed values omitted from the vector are appended.
#' @param colors `NULL` or a character colour vector. Unnamed colours are
#'   assigned in plotted group order. Named colours must cover every observed
#'   group. The default uses `BULK_VIOLIN_MACARON_COLORS` for up to 12 groups
#'   and a pastel HCL palette for larger sets.
#' @param violin_width Positive numeric scalar controlling violin width.
#' @param violin_alpha Numeric scalar in `[0, 1]` controlling violin opacity.
#' @param violin_linewidth Non-negative numeric scalar controlling violin
#'   outline width.
#' @param violin_scale One of `"width"`, `"area"`, or `"count"`, passed to
#'   [ggplot2::geom_violin()]. `"width"` gives every violin the same maximum
#'   width, `"area"` gives every violin the same total area, and `"count"`
#'   scales width according to the number of observations.
#' @param trim Logical scalar. Whether violin densities are trimmed to the
#'   observed expression range.
#' @param boxplot Logical scalar. Whether to overlay a white boxplot.
#' @param box_width Positive numeric scalar controlling boxplot width.
#' @param box_linewidth Non-negative numeric scalar controlling box outlines.
#' @param show_points Logical scalar. Whether to overlay individual samples.
#' @param point_size Positive numeric scalar controlling sample-point size.
#' @param point_alpha Numeric scalar in `[0, 1]` controlling point opacity.
#' @param jitter_width Non-negative numeric scalar controlling horizontal
#'   jitter. A fixed seed makes point placement reproducible.
#' @param facet_ncol `NULL` or a positive integer giving the number of gene
#'   columns when `split_by` is `NULL`. The default uses at most four columns.
#' @param free_y Logical scalar. Whether each facet row may use an independent
#'   y-axis scale.
#' @param italic_gene_title Logical scalar. Whether gene facet labels use
#'   italic type.
#' @param rotate_x Numeric scalar giving the x-axis label angle in degrees.
#' @param y_title Character, expression, or other ggplot2-compatible y-axis
#'   title. The default displays log2 CPM.
#' @param title `NULL` or a character scalar used as the plot title.
#' @param width Positive numeric scalar giving saved figure width in inches.
#' @param height `NULL` or a positive numeric scalar giving saved figure height
#'   in inches. When omitted, height is derived from the facet layout.
#' @param save `NULL` or a character scalar giving an output figure path. PDF
#'   files use `grDevices::cairo_pdf`; other formats are written at 300 dpi.
#' @param verbose Logical scalar. Whether to print sample, gene, grouping, and
#'   gene-resolution summaries.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{plot}{The assembled `ggplot` object.}
#'     \item{data}{Long-format plotting data with expression and sample
#'       metadata. Metadata fields that collide with generated columns are
#'       retained with a `metadata_` prefix.}
#'     \item{summary}{Per-gene, per-group, and, when requested, per-split
#'       sample count, mean, median, and standard deviation of log2 CPM.}
#'     \item{gene_mapping}{Resolved requested gene, symbol, stable ID, and
#'       original count-row mapping.}
#'     \item{logCPM}{The complete prepared log2-CPM matrix before gene
#'       selection.}
#'     \item{dge}{The normalized [edgeR::DGEList()] object.}
#'     \item{colors}{The named group-colour mapping used by the plot.}
#'     \item{normalization}{The normalization settings inherited from `bulk`.}
#'   }
#'
#' @details
#' Prepare raw sample-level bulk RNA-seq or replicate-aware pseudobulk counts
#' with [prepare_bulk_expression()] before plotting. The plot is descriptive
#' and does not replace a design-aware differential-expression model, effect
#' estimates, or biological replication checks. The input object is not
#' modified. The `ggplot2` package is required.
#'
#' @examples
#' \dontrun{
#' bulk <- prepare_bulk_expression(project, subset = Tissue == "Lung")
#' violin_result <- plot_bulk_violin(
#'     bulk,
#'     genes = c("GATA1", "SPI1", "CEBPA"),
#'     group_by = "Condition",
#'     split_by = "Sorting",
#'     group_order = c("Control", "Treated"),
#'     save = "outputs/bulk_gene_expression.pdf"
#' )
#' violin_result$plot
#' violin_result$summary
#' violin_result$gene_mapping
#' }
#'
#' @export
plot_bulk_violin <- function(
    bulk,
    genes,
    group_by = "Sorting",
    split_by = NULL,
    gene_col = "gene_name",
    gene_id_col = "gene_id",
    group_order = NULL,
    split_order = NULL,
    colors = NULL,
    violin_width = 0.85,
    violin_alpha = 0.65,
    violin_linewidth = 0.45,
    violin_scale = c("width", "area", "count"),
    trim = FALSE,
    boxplot = TRUE,
    box_width = 0.16,
    box_linewidth = 0.4,
    show_points = TRUE,
    point_size = 2.2,
    point_alpha = 0.9,
    jitter_width = 0.08,
    facet_ncol = NULL,
    free_y = FALSE,
    italic_gene_title = TRUE,
    rotate_x = 45,
    y_title = expression(log[2] ~ CPM),
    title = NULL,
    width = 7,
    height = NULL,
    save = NULL,
    verbose = TRUE
) {
    violin_scale <- match.arg(violin_scale)
    required_packages <- "ggplot2"
    missing_packages <- required_packages[!vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )]
    if (length(missing_packages) > 0L) {
        stop(
            "Please install required package(s): ",
            paste(missing_packages, collapse = ", "),
            call. = FALSE
        )
    }

    required_objects <- c("logCPM", "metadata", "genes", "dge")
    if (!is.list(bulk)) {
        stop(
            "`bulk` must be the output of prepare_bulk_expression().",
            call. = FALSE
        )
    }
    missing_objects <- setdiff(required_objects, names(bulk))
    if (length(missing_objects) > 0L) {
        stop(
            "`bulk` is missing: ", paste(missing_objects, collapse = ", "),
            ". Run prepare_bulk_expression() first.",
            call. = FALSE
        )
    }

    logcpm <- bulk$logCPM
    meta <- as.data.frame(bulk$metadata)
    anno <- as.data.frame(bulk$genes)
    dge <- bulk$dge
    normalization_label <- if (!is.null(bulk$normalization$label)) {
        bulk$normalization$label
    } else {
        "precomputed log2 CPM"
    }
    if ((!is.matrix(logcpm) && !inherits(logcpm, "Matrix")) ||
        !is.numeric(logcpm) || nrow(logcpm) == 0L || ncol(logcpm) == 0L) {
        stop("`bulk$logCPM` must be a non-empty numeric matrix-like object.",
             call. = FALSE)
    }
    if (is.null(colnames(logcpm)) || is.null(rownames(logcpm)) ||
        is.null(rownames(meta)) || is.null(rownames(anno))) {
        stop("Expression, metadata, and annotations require row/column names.",
             call. = FALSE)
    }
    if (anyDuplicated(colnames(logcpm)) || anyDuplicated(rownames(logcpm)) ||
        anyDuplicated(rownames(meta)) || anyDuplicated(rownames(anno))) {
        stop("Sample and gene row identifiers must be unique.", call. = FALSE)
    }
    if (!identical(colnames(logcpm), rownames(meta))) {
        stop(
            "colnames(bulk$logCPM) and rownames(bulk$metadata) must be ",
            "identical and in the same order.", call. = FALSE
        )
    }
    if (!identical(rownames(logcpm), rownames(anno))) {
        stop(
            "rownames(bulk$logCPM) and rownames(bulk$genes) must be identical ",
            "and in the same order.", call. = FALSE
        )
    }
    if (anyNA(logcpm) || any(!is.finite(logcpm))) {
        stop("`bulk$logCPM` must contain finite values.",
             call. = FALSE)
    }

    scalar_text <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    if (!scalar_text(group_by) ||
        (!is.null(split_by) && !scalar_text(split_by)) ||
        !scalar_text(gene_col) || !scalar_text(gene_id_col)) {
        stop("Grouping and annotation column names must be non-empty strings.",
             call. = FALSE)
    }
    reserved_columns <- c(".SampleID", "Gene", "GeneName", "GeneID", "Expression")
    if (group_by %in% reserved_columns ||
        (!is.null(split_by) && split_by %in% reserved_columns)) {
        stop("`group_by` and `split_by` cannot use generated output names.",
             call. = FALSE)
    }
    if (!group_by %in% colnames(meta)) {
        stop("`group_by = \"", group_by, "\"` was not found in bulk$metadata.",
             call. = FALSE)
    }
    if (!is.null(split_by) && !split_by %in% colnames(meta)) {
        stop("`split_by = \"", split_by, "\"` was not found in bulk$metadata.",
             call. = FALSE)
    }
    if (!is.null(split_by) && identical(group_by, split_by)) {
        stop("`split_by` must differ from `group_by`.", call. = FALSE)
    }

    logical_arguments <- list(
        trim = trim,
        boxplot = boxplot,
        show_points = show_points,
        free_y = free_y,
        italic_gene_title = italic_gene_title,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[!vapply(
        logical_arguments,
        function(x) is.logical(x) && length(x) == 1L && !is.na(x),
        logical(1)
    )]
    if (length(invalid_logical) > 0L) {
        stop(
            "These arguments must be single TRUE/FALSE values: ",
            paste(invalid_logical, collapse = ", "), call. = FALSE
        )
    }

    numeric_arguments <- list(
        violin_width = violin_width,
        violin_alpha = violin_alpha,
        violin_linewidth = violin_linewidth,
        box_width = box_width,
        box_linewidth = box_linewidth,
        point_size = point_size,
        point_alpha = point_alpha,
        jitter_width = jitter_width,
        width = width,
        rotate_x = rotate_x
    )
    invalid_numeric <- names(numeric_arguments)[!vapply(
        numeric_arguments,
        function(x) is.numeric(x) && length(x) == 1L &&
            !is.na(x) && is.finite(x),
        logical(1)
    )]
    if (length(invalid_numeric) > 0L) {
        stop(
            "These arguments must be finite numeric scalars: ",
            paste(invalid_numeric, collapse = ", "), call. = FALSE
        )
    }
    if (violin_width <= 0 || violin_linewidth < 0 ||
        box_width <= 0 || box_linewidth < 0 || point_size <= 0 ||
        jitter_width < 0 || width <= 0 ||
        violin_alpha < 0 || violin_alpha > 1 ||
        point_alpha < 0 || point_alpha > 1) {
        stop("Plot dimensions, widths, sizes, and alpha values are out of range.",
             call. = FALSE)
    }
    if (!is.null(height) && (!is.numeric(height) || length(height) != 1L ||
        is.na(height) || !is.finite(height) || height <= 0)) {
        stop("`height` must be NULL or a positive numeric scalar.", call. = FALSE)
    }
    if (!is.null(facet_ncol) && (!is.numeric(facet_ncol) ||
        length(facet_ncol) != 1L || is.na(facet_ncol) ||
        facet_ncol < 1 || facet_ncol != as.integer(facet_ncol))) {
        stop("`facet_ncol` must be NULL or a positive integer.", call. = FALSE)
    }
    if (!is.null(save) && !scalar_text(save)) {
        stop("`save` must be NULL or a non-empty file path.", call. = FALSE)
    }

    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0L) {
        stop("`genes` must contain at least one non-empty gene label or ID.",
             call. = FALSE)
    }
    for (order_name in c("group_order", "split_order")) {
        order_value <- get(order_name)
        if (!is.null(order_value) &&
            (!is.character(order_value) || anyNA(order_value) ||
             any(!nzchar(order_value)) || anyDuplicated(order_value))) {
            stop("`", order_name, "` must be NULL or unique non-empty strings.",
                 call. = FALSE)
        }
    }

    group_values <- as.character(meta[[group_by]])
    if (anyNA(group_values) || any(!nzchar(group_values))) {
        stop("The selected `group_by` column contains missing or empty values.",
             call. = FALSE)
    }
    if (!is.null(split_by)) {
        split_values <- as.character(meta[[split_by]])
        if (anyNA(split_values) || any(!nzchar(split_values))) {
            stop("The selected `split_by` column contains missing or empty values.",
                 call. = FALSE)
        }
    }

    anno$.GeneID_internal <- if (gene_id_col %in% colnames(anno)) {
        as.character(anno[[gene_id_col]])
    } else {
        rownames(anno)
    }
    bad_gene_ids <- is.na(anno$.GeneID_internal) |
        !nzchar(anno$.GeneID_internal)
    anno$.GeneID_internal[bad_gene_ids] <- rownames(anno)[bad_gene_ids]
    anno$.GeneName_internal <- if (gene_col %in% colnames(anno)) {
        as.character(anno[[gene_col]])
    } else {
        anno$.GeneID_internal
    }
    bad_gene_names <- is.na(anno$.GeneName_internal) |
        !nzchar(anno$.GeneName_internal)
    anno$.GeneName_internal[bad_gene_names] <-
        anno$.GeneID_internal[bad_gene_names]

    selected_rows <- list()
    gene_match_rows <- list()
    for (gene in genes) {
        candidates <- unique(c(
            which(anno$.GeneName_internal == gene),
            which(anno$.GeneID_internal == gene)
        ))
        if (length(candidates) == 0L) {
            warning("Gene not found: ", gene, call. = FALSE)
            next
        }
        if (length(candidates) > 1L) {
            mean_expression <- rowMeans(
                logcpm[candidates, , drop = FALSE], na.rm = TRUE
            )
            chosen <- candidates[which.max(mean_expression)]
            if (verbose) {
                message(
                    "Gene '", gene, "' matched ", length(candidates),
                    " rows; using ", anno$.GeneID_internal[chosen],
                    " (highest mean logCPM)."
                )
            }
            candidates <- chosen
        }
        selected_rows[[gene]] <- candidates
        gene_match_rows[[gene]] <- data.frame(
            RequestedGene = gene,
            GeneName = anno$.GeneName_internal[candidates],
            GeneID = anno$.GeneID_internal[candidates],
            CountRow = rownames(anno)[candidates],
            stringsAsFactors = FALSE
        )
    }
    if (length(selected_rows) == 0L) {
        stop("None of the requested genes were found.", call. = FALSE)
    }
    gene_mapping <- do.call(rbind, gene_match_rows)
    rownames(gene_mapping) <- NULL

    expression_rows <- lapply(names(selected_rows), function(gene) {
        index <- selected_rows[[gene]]
        data.frame(
            .SampleID = colnames(logcpm),
            Gene = gene,
            GeneName = anno$.GeneName_internal[index],
            GeneID = anno$.GeneID_internal[index],
            Expression = as.numeric(logcpm[index, , drop = TRUE]),
            stringsAsFactors = FALSE
        )
    })
    expression_data <- do.call(rbind, expression_rows)
    rownames(expression_data) <- NULL

    metadata_copy <- meta
    metadata_copy$.SampleID <- NULL
    collisions <- intersect(colnames(metadata_copy), colnames(expression_data))
    if (length(collisions) > 0L) {
        replacement_names <- paste0("metadata_", collisions)
        while (any(replacement_names %in% colnames(metadata_copy))) {
            replacement_names <- paste0("metadata_", replacement_names)
        }
        colnames(metadata_copy)[match(collisions, colnames(metadata_copy))] <-
            replacement_names
    }
    metadata_index <- match(expression_data$.SampleID, rownames(metadata_copy))
    if (anyNA(metadata_index)) {
        stop(
            "Metadata join failed for at least one sample. Check project order.",
            call. = FALSE
        )
    }
    expression_data <- cbind(
        expression_data,
        metadata_copy[metadata_index, , drop = FALSE]
    )
    rownames(expression_data) <- NULL

    observed_groups <- unique(as.character(expression_data[[group_by]]))
    if (is.null(group_order)) {
        group_order <- observed_groups
    } else {
        group_order <- c(group_order, setdiff(observed_groups, group_order))
    }
    group_order <- group_order[group_order %in% observed_groups]
    expression_data[[group_by]] <- factor(
        as.character(expression_data[[group_by]]), levels = group_order
    )

    if (!is.null(split_by)) {
        observed_splits <- unique(as.character(expression_data[[split_by]]))
        if (is.null(split_order)) {
            split_order <- observed_splits
        } else {
            split_order <- c(split_order, setdiff(observed_splits, split_order))
        }
        split_order <- split_order[split_order %in% observed_splits]
        expression_data[[split_by]] <- factor(
            as.character(expression_data[[split_by]]), levels = split_order
        )
    }

    genes_found <- names(selected_rows)
    expression_data$Gene <- factor(expression_data$Gene, levels = genes_found)
    group_levels <- levels(droplevels(expression_data[[group_by]]))
    if (is.null(colors)) {
        colors <- if (length(group_levels) <=
            length(BULK_VIOLIN_MACARON_COLORS)) {
            BULK_VIOLIN_MACARON_COLORS[seq_along(group_levels)]
        } else {
            grDevices::hcl.colors(length(group_levels), palette = "Pastel 1")
        }
        names(colors) <- group_levels
    } else {
        if (!is.character(colors) || anyNA(colors) || any(!nzchar(colors))) {
            stop("`colors` must contain valid non-missing colour strings.",
                 call. = FALSE)
        }
        if (is.null(names(colors))) {
            if (length(colors) < length(group_levels)) {
                stop("Not enough colors supplied.", call. = FALSE)
            }
            colors <- colors[seq_along(group_levels)]
            names(colors) <- group_levels
        } else {
            missing_colors <- setdiff(group_levels, names(colors))
            if (length(missing_colors) > 0L) {
                stop(
                    "Missing colors for group(s): ",
                    paste(missing_colors, collapse = ", "), call. = FALSE
                )
            }
            colors <- colors[group_levels]
        }
    }

    plot <- ggplot2::ggplot(
        expression_data,
        ggplot2::aes(
            x = .data[[group_by]],
            y = Expression,
            fill = .data[[group_by]],
            colour = .data[[group_by]]
        )
    ) +
        ggplot2::geom_violin(
            width = violin_width,
            alpha = violin_alpha,
            trim = trim,
            linewidth = violin_linewidth,
            scale = violin_scale
        )
    if (boxplot) {
        plot <- plot + ggplot2::geom_boxplot(
            width = box_width,
            outlier.shape = NA,
            fill = "white",
            colour = "black",
            linewidth = box_linewidth
        )
    }
    if (show_points) {
        plot <- plot + ggplot2::geom_point(
            position = ggplot2::position_jitter(
                width = jitter_width, height = 0, seed = 123
            ),
            size = point_size,
            alpha = point_alpha,
            stroke = 0.3
        )
    }

    scale_setting <- if (free_y) "free_y" else "fixed"
    n_genes <- length(genes_found)
    if (is.null(split_by)) {
        facet_columns <- if (is.null(facet_ncol)) {
            min(n_genes, 4L)
        } else {
            as.integer(facet_ncol)
        }
        plot <- plot + ggplot2::facet_wrap(
            ~Gene, ncol = facet_columns, scales = scale_setting
        )
    } else {
        split_formula <- stats::reformulate(split_by, response = "Gene")
        plot <- plot + ggplot2::facet_grid(
            split_formula, scales = scale_setting
        )
    }

    gene_face_x <- if (is.null(split_by) && italic_gene_title) "italic" else "bold"
    gene_face_y <- if (!is.null(split_by) && italic_gene_title) "italic" else "bold"
    plot <- plot +
        ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
        ggplot2::scale_colour_manual(values = colors, drop = FALSE) +
        ggplot2::labs(
            x = NULL, y = y_title, fill = group_by,
            colour = group_by, title = title
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
            axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
            axis.ticks.length = grid::unit(1.5, "mm"),
            axis.text = ggplot2::element_text(colour = "black", size = 8.5),
            axis.text.x = ggplot2::element_text(
                angle = rotate_x,
                hjust = if (rotate_x == 0) 0.5 else 1,
                vjust = if (rotate_x == 0) 0.5 else 1
            ),
            axis.title = ggplot2::element_text(colour = "black", size = 9.5),
            strip.background = ggplot2::element_blank(),
            strip.text.y = ggplot2::element_text(
                colour = "black", face = gene_face_y, size = 9
            ),
            strip.text.x = ggplot2::element_text(
                colour = "black", face = gene_face_x, size = 9.5
            ),
            plot.title = ggplot2::element_text(
                colour = "black", face = "bold", size = 10.5, hjust = 0.5
            ),
            legend.position = "none",
            plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
        )

    if (is.null(height)) {
        height <- if (is.null(split_by)) {
            max(3.5, ceiling(n_genes / facet_columns) * 3)
        } else {
            max(3.5, n_genes * 2.5)
        }
    }
    if (!is.null(save)) {
        save_arguments <- list(
            filename = save,
            plot = plot,
            width = width,
            height = height,
            units = "in"
        )
        if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
            save_arguments$device <- grDevices::cairo_pdf
        } else {
            save_arguments$dpi <- 300
        }
        do.call(ggplot2::ggsave, save_arguments)
    }

    summary_fields <- c("Gene", group_by)
    if (!is.null(split_by)) {
        summary_fields <- c(summary_fields, split_by)
    }
    summary_indices <- split(
        seq_len(nrow(expression_data)),
        do.call(
            interaction,
            c(expression_data[summary_fields], list(drop = TRUE, lex.order = TRUE))
        )
    )
    summary_data <- do.call(rbind, lapply(summary_indices, function(index) {
        group_values <- expression_data[index[1], summary_fields, drop = FALSE]
        values <- expression_data$Expression[index]
        cbind(
            group_values,
            data.frame(
                n = length(values),
                mean = mean(values, na.rm = TRUE),
                median = stats::median(values, na.rm = TRUE),
                sd = stats::sd(values, na.rm = TRUE),
                row.names = NULL
            )
        )
    }))
    rownames(summary_data) <- NULL

    if (verbose) {
        message("")
        message("========================================")
        message(" Bulk RNA-seq violin plot")
        message("========================================")
        message("Samples:  ", length(unique(expression_data$.SampleID)))
        message("Genes:    ", paste(genes_found, collapse = ", "))
        message("Group by: ", group_by)
        if (!is.null(split_by)) {
            message("Split by: ", split_by)
        }
        message("Violin scale: ", violin_scale)
        message("Normalization: ", normalization_label)
        message("")
        message("Gene mapping:")
        print(gene_mapping)
        message("")
        message("Group counts:")
        print(table(expression_data[[group_by]]))
        message("========================================")
    }

    list(
        plot = plot,
        data = expression_data,
        summary = summary_data,
        gene_mapping = gene_mapping,
        logCPM = logcpm,
        dge = dge,
        normalization = bulk$normalization,
        colors = colors
    )
}


#' Plot expression programs as a gene-set heatmap
#'
#' Resolve ordered gene sets against a prepared bulk-expression object,
#' optionally aggregate samples by metadata combinations, and display the
#' selected genes in named row slices.
#'
#' @param bulk A `bulk_expression_prepared` object returned by
#'   [prepare_bulk_expression()]. Its `logCPM`, `metadata`, `genes`, and `dge`
#'   components must retain identical sample and gene ordering.
#' @param gene_sets Named non-empty list of character vectors containing gene
#'   symbols or stable IDs. List order defines row-slice order. When the same
#'   expression row appears in multiple sets, its first set assignment wins.
#' @param gene_col Character scalar naming the preferred gene-symbol column in
#'   `bulk$genes`. Stable IDs are used when this column is absent or blank.
#' @param gene_id_col Character scalar naming the preferred stable-ID column in
#'   `bulk$genes`. Annotation row names are used when this column is absent.
#' @param aggregate_by `NULL` or a unique character vector naming metadata
#'   columns whose combinations define aggregated heatmap columns.
#' @param aggregate_fun Function used to combine replicate expression values
#'   within each aggregate group. It must accept a numeric vector and
#'   `na.rm = TRUE`, and return one finite numeric value.
#' @param factor_orders `NULL` or a named list mapping metadata columns to
#'   unique character vectors defining their level order. Observed values not
#'   listed are appended in first-seen order.
#' @param order_by `NULL` or a unique character vector naming metadata columns
#'   used for explicit heatmap-column ordering after optional aggregation.
#' @param annotation_cols `NULL` or a unique character vector naming metadata
#'   columns shown above the heatmap.
#' @param annotation_colors `NULL` or a named list of annotation colour maps
#'   passed to [ComplexHeatmap::HeatmapAnnotation()].
#' @param show_annotation_name Logical scalar. Whether top-annotation names are
#'   displayed.
#' @param show_annotation_legend Logical scalar. Whether top-annotation legends
#'   are displayed.
#' @param row_zscore Logical scalar. Whether each selected gene is standardized
#'   across displayed columns. Zero-variance genes are removed before scaling.
#' @param z_cap `NULL` or a positive numeric scalar used to cap row z-scores
#'   symmetrically. It is ignored when `row_zscore = FALSE`.
#' @param col_fun `NULL` or a colour-mapping function accepted by
#'   [ComplexHeatmap::Heatmap()]. The default uses
#'   `GENE_SET_HEATMAP_COLORS` across the observed or capped range.
#' @param preserve_gene_order Logical scalar. Whether to force the requested
#'   gene order and disable row clustering. Set to `FALSE` for `cluster_rows`
#'   to take effect.
#' @param cluster_rows Logical scalar. Whether rows are clustered when
#'   `preserve_gene_order = FALSE`.
#' @param cluster_row_slices Logical scalar. Whether named gene-set slices are
#'   clustered relative to one another.
#' @param cluster_columns Logical scalar. Whether displayed columns are
#'   clustered. Enabling this can override the visible `order_by` arrangement.
#' @param show_row_dend,show_column_dend Logical scalars controlling row and
#'   column dendrogram display.
#' @param clustering_distance_rows Character scalar or distance function passed
#'   to [ComplexHeatmap::Heatmap()] for row clustering.
#' @param clustering_method_rows Character scalar naming the row-clustering
#'   linkage method.
#' @param show_row_names Logical scalar. Whether gene labels are displayed.
#' @param row_names_italic Logical scalar. Whether gene labels use italic text.
#' @param row_names_fontsize Positive numeric scalar controlling gene-label
#'   text size.
#' @param show_row_titles Logical scalar. Whether gene-set names are displayed
#'   as row-slice titles.
#' @param row_title_fontsize Positive numeric scalar controlling gene-set title
#'   text size.
#' @param row_gap_mm Non-negative numeric scalar controlling space between
#'   gene-set slices in millimetres.
#' @param show_column_names Logical scalar. Whether sample or aggregate column
#'   labels are displayed.
#' @param column_names_fontsize Positive numeric scalar controlling column-label
#'   text size.
#' @param column_names_rot Numeric scalar giving the column-label rotation in
#'   degrees.
#' @param border Logical scalar. Whether to draw a heatmap border.
#' @param cell_border Character scalar specifying cell-border colour.
#' @param cell_border_lwd Non-negative numeric scalar controlling cell-border
#'   line width.
#' @param use_raster Logical scalar passed to [ComplexHeatmap::Heatmap()].
#' @param heatmap_name Non-empty character scalar used as the heatmap name.
#' @param heatmap_legend_title Non-empty character scalar used as the heatmap
#'   legend title.
#' @param width Positive numeric scalar giving saved figure width in inches.
#' @param height `NULL` or a positive numeric scalar giving saved figure height
#'   in inches. The default scales with the number of displayed genes.
#' @param save `NULL` or a character scalar giving an output path. PDF files use
#'   `grDevices::cairo_pdf`; other paths are written as 300-dpi PNG files.
#' @param draw_plot Logical scalar. Whether to draw the heatmap on the current
#'   graphics device in addition to any saved output.
#' @param verbose Logical scalar. Whether to report gene, column, program,
#'   normalization, aggregation, scaling, and missing-gene summaries.
#'
#' @return Invisibly returns a named list containing:
#'   \describe{
#'     \item{heatmap}{The assembled `ComplexHeatmap::Heatmap` object.}
#'     \item{matrix}{The displayed matrix after optional aggregation, row
#'       standardization, and capping.}
#'     \item{expression}{The selected prepared log2-CPM matrix after optional
#'       aggregation and before row standardization.}
#'     \item{logCPM}{The complete prepared log2-CPM matrix before gene
#'       selection.}
#'     \item{gene_mapping}{Requested set, requested gene, resolved symbol,
#'       stable ID, source count row, and unique display name for every plotted
#'       gene.}
#'     \item{missing_genes}{Requested gene-set members that could not be
#'       resolved, with their originating gene-set names.}
#'     \item{metadata}{Metadata aligned to displayed heatmap columns.}
#'     \item{row_split}{Gene-set factor aligned to displayed matrix rows.}
#'     \item{col_fun}{The colour-mapping function used by the heatmap.}
#'     \item{dge}{The prepared [edgeR::DGEList()] object.}
#'     \item{normalization}{The normalization settings inherited from `bulk`.}
#'   }
#'
#' @details
#' Prepare raw sample-level bulk RNA-seq counts or replicate-aware pseudobulk
#' counts with [prepare_bulk_expression()] before calling this function. This
#' keeps normalization and sample selection consistent across reusable plots.
#' When aggregation is requested, any
#' metadata fields used for ordering or annotation must be constant within
#' every aggregate group; otherwise the function stops rather than displaying
#' a misleading annotation. This heatmap is descriptive and does not replace
#' design-aware differential-expression analysis or biological replication
#' checks. The supplied object is not modified. Required packages are
#' `ComplexHeatmap` and `circlize`.
#'
#' @examples
#' \dontrun{
#' programs <- list(
#'     Stemness = c("GATA2", "KIT", "PROM1"),
#'     Myeloid = c("SPI1", "CEBPA", "MPO")
#' )
#' bulk <- prepare_bulk_expression(
#'     project,
#'     subset = Condition != "Excluded"
#' )
#' program_heatmap <- plot_gene_set_heatmap(
#'     bulk,
#'     gene_sets = programs,
#'     aggregate_by = c("Condition", "Sorting"),
#'     factor_orders = list(Condition = c("Control", "Treated")),
#'     annotation_cols = c("Condition", "Sorting"),
#'     save = "outputs/gene_set_heatmap.pdf"
#' )
#' program_heatmap$gene_mapping
#' program_heatmap$matrix
#' }
#'
#' @export
plot_gene_set_heatmap <- function(
    bulk,
    gene_sets,
    gene_col = "gene_name",
    gene_id_col = "gene_id",
    aggregate_by = NULL,
    aggregate_fun = mean,
    factor_orders = NULL,
    order_by = NULL,
    annotation_cols = NULL,
    annotation_colors = NULL,
    show_annotation_name = FALSE,
    show_annotation_legend = TRUE,
    row_zscore = TRUE,
    z_cap = 2.5,
    col_fun = NULL,
    preserve_gene_order = TRUE,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    show_row_names = TRUE,
    row_names_italic = TRUE,
    row_names_fontsize = 9,
    show_row_titles = TRUE,
    row_title_fontsize = 9,
    row_gap_mm = 2,
    show_column_names = TRUE,
    column_names_fontsize = 8,
    column_names_rot = 45,
    border = TRUE,
    cell_border = "white",
    cell_border_lwd = 0.4,
    use_raster = FALSE,
    heatmap_name = "Row z-score",
    heatmap_legend_title = "Expression\nrow z-score",
    width = 8,
    height = NULL,
    save = NULL,
    draw_plot = TRUE,
    verbose = TRUE
) {
    required_packages <- c("ComplexHeatmap", "circlize")
    missing_packages <- required_packages[!vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )]
    if (length(missing_packages) > 0L) {
        stop(
            "Please install required package(s): ",
            paste(missing_packages, collapse = ", "), call. = FALSE
        )
    }

    scalar_text <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    if (!scalar_text(gene_col) || !scalar_text(gene_id_col) ||
        !scalar_text(heatmap_name) || !scalar_text(heatmap_legend_title) ||
        !scalar_text(cell_border) ||
        (!is.null(save) && !scalar_text(save))) {
        stop("Column, display, and path arguments must be non-empty strings.",
             call. = FALSE)
    }
    if (!is.function(aggregate_fun)) {
        stop("`aggregate_fun` must be a function.", call. = FALSE)
    }
    if (!is.null(col_fun) && !is.function(col_fun)) {
        stop("`col_fun` must be NULL or a colour-mapping function.",
             call. = FALSE)
    }
    if (!(scalar_text(clustering_distance_rows) ||
        is.function(clustering_distance_rows)) ||
        !scalar_text(clustering_method_rows)) {
        stop("Row clustering distance and method arguments are invalid.",
             call. = FALSE)
    }

    logical_arguments <- list(
        show_annotation_name = show_annotation_name,
        show_annotation_legend = show_annotation_legend,
        row_zscore = row_zscore,
        preserve_gene_order = preserve_gene_order,
        cluster_rows = cluster_rows,
        cluster_row_slices = cluster_row_slices,
        cluster_columns = cluster_columns,
        show_row_dend = show_row_dend,
        show_column_dend = show_column_dend,
        show_row_names = show_row_names,
        row_names_italic = row_names_italic,
        show_row_titles = show_row_titles,
        show_column_names = show_column_names,
        border = border,
        use_raster = use_raster,
        draw_plot = draw_plot,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[!vapply(
        logical_arguments,
        function(x) is.logical(x) && length(x) == 1L && !is.na(x),
        logical(1)
    )]
    if (length(invalid_logical) > 0L) {
        stop(
            "These arguments must be single TRUE/FALSE values: ",
            paste(invalid_logical, collapse = ", "), call. = FALSE
        )
    }
    numeric_arguments <- list(
        row_names_fontsize = row_names_fontsize,
        row_title_fontsize = row_title_fontsize,
        row_gap_mm = row_gap_mm,
        column_names_fontsize = column_names_fontsize,
        column_names_rot = column_names_rot,
        cell_border_lwd = cell_border_lwd,
        width = width
    )
    invalid_numeric <- names(numeric_arguments)[!vapply(
        numeric_arguments,
        function(x) is.numeric(x) && length(x) == 1L &&
            !is.na(x) && is.finite(x),
        logical(1)
    )]
    if (length(invalid_numeric) > 0L) {
        stop(
            "These arguments must be finite numeric scalars: ",
            paste(invalid_numeric, collapse = ", "), call. = FALSE
        )
    }
    if (row_names_fontsize <= 0 ||
        row_title_fontsize <= 0 || row_gap_mm < 0 ||
        column_names_fontsize <= 0 || cell_border_lwd < 0 || width <= 0) {
        stop("Display values are outside allowed ranges.",
             call. = FALSE)
    }
    if (!is.null(height) && (!is.numeric(height) || length(height) != 1L ||
        is.na(height) || !is.finite(height) || height <= 0)) {
        stop("`height` must be NULL or a positive numeric scalar.", call. = FALSE)
    }
    if (!is.null(z_cap) && (!is.numeric(z_cap) || length(z_cap) != 1L ||
        is.na(z_cap) || !is.finite(z_cap) || z_cap <= 0)) {
        stop("`z_cap` must be NULL or a positive numeric scalar.", call. = FALSE)
    }

    if (!is.list(bulk)) {
        stop(
            "`bulk` must be the output of prepare_bulk_expression().",
            call. = FALSE
        )
    }
    required_objects <- c("logCPM", "metadata", "genes", "dge")
    missing_objects <- setdiff(required_objects, names(bulk))
    if (length(missing_objects) > 0L) {
        stop(
            "`bulk` is missing: ", paste(missing_objects, collapse = ", "),
            ". Run prepare_bulk_expression() first.",
            call. = FALSE
        )
    }
    logcpm <- bulk$logCPM
    meta <- as.data.frame(bulk$metadata)
    anno <- as.data.frame(bulk$genes)
    dge <- bulk$dge
    normalization_label <- if (!is.null(bulk$normalization$label)) {
        bulk$normalization$label
    } else {
        "precomputed log2 CPM"
    }
    if ((!is.matrix(logcpm) && !inherits(logcpm, "Matrix")) ||
        !is.numeric(logcpm) || nrow(logcpm) == 0L || ncol(logcpm) == 0L) {
        stop("`bulk$logCPM` must be a non-empty numeric matrix-like object.",
             call. = FALSE)
    }
    if (is.null(rownames(logcpm)) || is.null(colnames(logcpm)) ||
        is.null(rownames(meta)) || is.null(rownames(anno))) {
        stop("Expression, metadata, and annotations require row/column names.",
             call. = FALSE)
    }
    if (anyDuplicated(rownames(logcpm)) || anyDuplicated(colnames(logcpm)) ||
        anyDuplicated(rownames(meta)) || anyDuplicated(rownames(anno))) {
        stop("Sample and gene identifiers must be unique.", call. = FALSE)
    }
    if (anyNA(logcpm) || any(!is.finite(logcpm))) {
        stop("`bulk$logCPM` must contain finite values.",
             call. = FALSE)
    }
    if (!identical(colnames(logcpm), rownames(meta))) {
        stop("Expression columns and metadata rows are not aligned.",
             call. = FALSE)
    }
    if (!identical(rownames(logcpm), rownames(anno))) {
        stop("Expression rows and gene annotation rows are not aligned.",
             call. = FALSE)
    }

    if (!is.list(gene_sets) || length(gene_sets) == 0L ||
        is.null(names(gene_sets)) || anyNA(names(gene_sets)) ||
        any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets))) {
        stop("`gene_sets` must be a non-empty list with unique names.",
             call. = FALSE)
    }
    gene_sets <- lapply(gene_sets, function(genes) {
        genes <- unique(as.character(genes))
        genes[!is.na(genes) & nzchar(genes)]
    })
    empty_sets <- names(gene_sets)[lengths(gene_sets) == 0L]
    if (length(empty_sets) > 0L) {
        stop(
            "Gene sets cannot be empty after cleaning: ",
            paste(empty_sets, collapse = ", "), call. = FALSE
        )
    }

    if (!is.null(factor_orders)) {
        if (!is.list(factor_orders) || is.null(names(factor_orders)) ||
            anyNA(names(factor_orders)) || any(!nzchar(names(factor_orders))) ||
            anyDuplicated(names(factor_orders))) {
            stop("`factor_orders` must be NULL or a uniquely named list.",
                 call. = FALSE)
        }
        missing_factor_columns <- setdiff(names(factor_orders), colnames(meta))
        if (length(missing_factor_columns) > 0L) {
            stop(
                "factor_orders column(s) not found: ",
                paste(missing_factor_columns, collapse = ", "), call. = FALSE
            )
        }
        for (column in names(factor_orders)) {
            requested_levels <- factor_orders[[column]]
            if (!is.character(requested_levels) || anyNA(requested_levels) ||
                any(!nzchar(requested_levels)) || anyDuplicated(requested_levels)) {
                stop("Every factor order must contain unique non-empty strings.",
                     call. = FALSE)
            }
            observed_levels <- unique(as.character(meta[[column]]))
            if (anyNA(observed_levels) || any(!nzchar(observed_levels))) {
                stop("Ordered metadata fields cannot contain missing values.",
                     call. = FALSE)
            }
            levels_use <- c(
                requested_levels, setdiff(observed_levels, requested_levels)
            )
            levels_use <- levels_use[levels_use %in% observed_levels]
            meta[[column]] <- factor(
                as.character(meta[[column]]), levels = levels_use
            )
        }
    }

    anno$.GeneID <- if (gene_id_col %in% colnames(anno)) {
        as.character(anno[[gene_id_col]])
    } else {
        rownames(anno)
    }
    bad_ids <- is.na(anno$.GeneID) | !nzchar(anno$.GeneID)
    anno$.GeneID[bad_ids] <- rownames(anno)[bad_ids]
    anno$.GeneName <- if (gene_col %in% colnames(anno)) {
        as.character(anno[[gene_col]])
    } else {
        anno$.GeneID
    }
    bad_names <- is.na(anno$.GeneName) | !nzchar(anno$.GeneName)
    anno$.GeneName[bad_names] <- anno$.GeneID[bad_names]

    mapping_rows <- list()
    mapping_index <- 0L
    missing_gene_rows <- list()
    for (set_name in names(gene_sets)) {
        for (gene in gene_sets[[set_name]]) {
            candidates <- which(anno$.GeneName == gene | anno$.GeneID == gene)
            if (length(candidates) == 0L) {
                missing_gene_rows[[length(missing_gene_rows) + 1L]] <- data.frame(
                    GeneSet = set_name,
                    RequestedGene = gene,
                    stringsAsFactors = FALSE
                )
                warning("Gene not found: ", gene, " [", set_name, "]",
                        call. = FALSE)
                next
            }
            if (length(candidates) > 1L) {
                mean_expression <- rowMeans(
                    logcpm[candidates, , drop = FALSE], na.rm = TRUE
                )
                candidates <- candidates[which.max(mean_expression)]
                if (verbose) {
                    message(
                        "Gene ", gene, " matched multiple rows; using ",
                        anno$.GeneID[candidates]
                    )
                }
            }
            mapping_index <- mapping_index + 1L
            mapping_rows[[mapping_index]] <- data.frame(
                GeneSet = set_name,
                RequestedGene = gene,
                GeneName = anno$.GeneName[candidates],
                GeneID = anno$.GeneID[candidates],
                MatrixRow = rownames(anno)[candidates],
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(mapping_rows) == 0L) {
        stop("None of the requested genes were found.", call. = FALSE)
    }
    gene_mapping <- do.call(rbind, mapping_rows)
    rownames(gene_mapping) <- NULL
    duplicated_rows <- duplicated(gene_mapping$MatrixRow)
    if (any(duplicated_rows) && verbose) {
        message(sum(duplicated_rows), " duplicate gene-set assignment(s) removed.")
    }
    gene_mapping <- gene_mapping[!duplicated_rows, , drop = FALSE]

    expression_matrix <- logcpm[gene_mapping$MatrixRow, , drop = FALSE]
    display_names <- make.unique(gene_mapping$GeneName)
    rownames(expression_matrix) <- display_names
    gene_mapping$DisplayName <- display_names

    character_columns <- function(x, argument) {
        if (is.null(x)) {
            return(NULL)
        }
        if (!is.character(x) || length(x) == 0L || anyNA(x) ||
            any(!nzchar(x)) || anyDuplicated(x)) {
            stop("`", argument, "` must contain unique non-empty names.",
                 call. = FALSE)
        }
        x
    }
    aggregate_by <- character_columns(aggregate_by, "aggregate_by")
    order_by <- character_columns(order_by, "order_by")
    annotation_cols <- character_columns(annotation_cols, "annotation_cols")
    if (!is.null(annotation_colors) && is.null(annotation_cols)) {
        stop("`annotation_colors` requires `annotation_cols`.", call. = FALSE)
    }
    requested_metadata <- unique(c(aggregate_by, order_by, annotation_cols))
    missing_metadata <- setdiff(requested_metadata, colnames(meta))
    if (length(missing_metadata) > 0L) {
        stop(
            "Requested metadata column(s) not found: ",
            paste(missing_metadata, collapse = ", "), call. = FALSE
        )
    }

    if (!is.null(aggregate_by)) {
        label_parts <- lapply(meta[aggregate_by], as.character)
        invalid_aggregate_values <- vapply(
            label_parts,
            function(x) anyNA(x) || any(!nzchar(x)),
            logical(1)
        )
        if (any(invalid_aggregate_values)) {
            stop("Aggregation columns cannot contain missing or empty values.",
                 call. = FALSE)
        }
        group_labels <- do.call(paste, c(label_parts, sep = " | "))
        first_group_rows <- which(!duplicated(group_labels))
        ordering_data <- meta[first_group_rows, aggregate_by, drop = FALSE]
        group_order_index <- do.call(
            order, c(ordering_data, list(na.last = TRUE))
        )
        group_levels <- group_labels[first_group_rows][group_order_index]

        aggregate_one_group <- function(group_label) {
            sample_indices <- which(group_labels == group_label)
            if (length(sample_indices) == 1L) {
                return(as.numeric(expression_matrix[, sample_indices]))
            }
            vapply(seq_len(nrow(expression_matrix)), function(row_index) {
                value <- aggregate_fun(
                    expression_matrix[row_index, sample_indices], na.rm = TRUE
                )
                if (!is.numeric(value) || length(value) != 1L ||
                    is.na(value) || !is.finite(value)) {
                    stop(
                        "`aggregate_fun` must return one finite numeric value ",
                        "per gene and group.", call. = FALSE
                    )
                }
                as.numeric(value)
            }, numeric(1))
        }
        aggregated_columns <- lapply(group_levels, aggregate_one_group)
        matrix_plot <- do.call(cbind, aggregated_columns)
        rownames(matrix_plot) <- rownames(expression_matrix)
        colnames(matrix_plot) <- group_levels

        metadata_fields <- unique(c(aggregate_by, order_by, annotation_cols))
        aggregate_metadata_rows <- lapply(group_levels, function(group_label) {
            sample_indices <- which(group_labels == group_label)
            for (column in metadata_fields) {
                values <- meta[[column]][sample_indices]
                observed <- unique(as.character(values[!is.na(values)]))
                if (length(observed) != 1L || anyNA(values)) {
                    stop(
                        "Metadata column '", column,
                        "' is not constant within aggregate group '",
                        group_label, "'.", call. = FALSE
                    )
                }
            }
            output <- meta[sample_indices[1L], metadata_fields, drop = FALSE]
            output$.HeatmapID <- group_label
            output
        })
        column_metadata <- do.call(rbind, aggregate_metadata_rows)
        rownames(column_metadata) <- group_levels
    } else {
        matrix_plot <- expression_matrix
        column_metadata <- meta
        column_metadata$.HeatmapID <- rownames(column_metadata)
    }

    if (!is.null(order_by)) {
        if (anyNA(column_metadata[order_by])) {
            stop("Ordering columns cannot contain missing values.", call. = FALSE)
        }
        order_index <- do.call(
            order, c(column_metadata[order_by], list(na.last = TRUE))
        )
        column_metadata <- column_metadata[order_index, , drop = FALSE]
        matrix_plot <- matrix_plot[
            , column_metadata$.HeatmapID, drop = FALSE
        ]
    }

    if (row_zscore) {
        if (ncol(matrix_plot) < 2L) {
            stop("At least two displayed columns are required for row z-scores.",
                 call. = FALSE)
        }
        row_sd <- apply(matrix_plot, 1L, stats::sd)
        keep_gene <- is.finite(row_sd) & row_sd > 0
        if (any(!keep_gene)) {
            warning(
                sum(!keep_gene),
                " gene(s) had zero or undefined variance and were removed.",
                call. = FALSE
            )
            matrix_plot <- matrix_plot[keep_gene, , drop = FALSE]
            gene_mapping <- gene_mapping[keep_gene, , drop = FALSE]
        }
        if (nrow(matrix_plot) == 0L) {
            stop("No variable genes remain for row z-scoring.", call. = FALSE)
        }
        matrix_scaled <- t(scale(t(matrix_plot)))
    } else {
        matrix_scaled <- matrix_plot
    }
    if (row_zscore && !is.null(z_cap)) {
        matrix_scaled <- pmax(pmin(matrix_scaled, z_cap), -z_cap)
    }

    gene_set_order <- names(gene_sets)
    row_split <- factor(gene_mapping$GeneSet, levels = gene_set_order)
    cluster_rows_use <- if (preserve_gene_order) FALSE else cluster_rows

    top_annotation <- NULL
    if (!is.null(annotation_cols)) {
        if (!is.null(annotation_colors) &&
            (!is.list(annotation_colors) || is.null(names(annotation_colors)) ||
             anyNA(names(annotation_colors)) || any(!nzchar(names(annotation_colors))))) {
            stop("`annotation_colors` must be NULL or a named list.",
                 call. = FALSE)
        }
        annotation_data <- column_metadata[, annotation_cols, drop = FALSE]
        annotation_arguments <- list(
            df = annotation_data,
            show_annotation_name = show_annotation_name,
            show_legend = show_annotation_legend
        )
        if (!is.null(annotation_colors)) {
            annotation_arguments$col <- annotation_colors
        }
        top_annotation <- do.call(
            ComplexHeatmap::HeatmapAnnotation, annotation_arguments
        )
    }

    legend_breaks <- NULL
    if (is.null(col_fun)) {
        if (row_zscore) {
            color_limit <- if (is.null(z_cap)) max(abs(matrix_scaled)) else z_cap
            if (!is.finite(color_limit) || color_limit == 0) {
                color_limit <- 1
            }
            legend_breaks <- c(
                -color_limit, -0.4 * color_limit, 0,
                0.4 * color_limit, color_limit
            )
            col_fun <- circlize::colorRamp2(
                legend_breaks, GENE_SET_HEATMAP_COLORS
            )
        } else {
            expression_range <- range(matrix_scaled)
            if (diff(expression_range) == 0) {
                expression_range <- expression_range + c(-1, 1)
            }
            color_breaks <- seq(
                expression_range[1L], expression_range[2L], length.out = 5L
            )
            col_fun <- circlize::colorRamp2(
                color_breaks, GENE_SET_HEATMAP_COLORS
            )
        }
    }

    heatmap_legend_parameters <- list(
        title = heatmap_legend_title,
        title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9)
    )
    if (!is.null(legend_breaks)) {
        heatmap_legend_parameters$at <- legend_breaks
    }
    row_title_value <- if (show_row_titles) levels(row_split) else NULL
    heatmap <- ComplexHeatmap::Heatmap(
        matrix_scaled,
        name = heatmap_name,
        col = col_fun,
        top_annotation = top_annotation,
        split = row_split,
        row_title = row_title_value,
        row_title_gp = grid::gpar(
            fontsize = row_title_fontsize, fontface = "bold"
        ),
        row_gap = grid::unit(row_gap_mm, "mm"),
        cluster_rows = cluster_rows_use,
        cluster_row_slices = cluster_row_slices,
        clustering_distance_rows = clustering_distance_rows,
        clustering_method_rows = clustering_method_rows,
        show_row_dend = show_row_dend,
        show_row_names = show_row_names,
        row_names_gp = grid::gpar(
            fontsize = row_names_fontsize,
            fontface = if (row_names_italic) "italic" else "plain"
        ),
        cluster_columns = cluster_columns,
        show_column_dend = show_column_dend,
        show_column_names = show_column_names,
        column_names_rot = column_names_rot,
        column_names_gp = grid::gpar(fontsize = column_names_fontsize),
        border = border,
        rect_gp = grid::gpar(col = cell_border, lwd = cell_border_lwd),
        use_raster = use_raster,
        heatmap_legend_param = heatmap_legend_parameters
    )

    if (is.null(height)) {
        height <- max(5, nrow(matrix_scaled) * 0.23 + 2)
    }
    draw_heatmap <- function() {
        ComplexHeatmap::draw(
            heatmap,
            heatmap_legend_side = "right",
            annotation_legend_side = "right"
        )
    }
    if (!is.null(save)) {
        save_heatmap <- function() {
            if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
                grDevices::cairo_pdf(filename = if (is.null(save)) NA else save, width = width, height = height)
            } else {
                grDevices::png(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300
                )
            }
            on.exit(grDevices::dev.off(), add = TRUE)
            draw_heatmap()
        }
        save_heatmap()
    }
    if (draw_plot) {
        draw_heatmap()
    }

    missing_gene_data <- if (length(missing_gene_rows) == 0L) {
        data.frame(GeneSet = character(0), RequestedGene = character(0))
    } else {
        do.call(rbind, missing_gene_rows)
    }
    rownames(missing_gene_data) <- NULL
    if (verbose) {
        message("")
        message("==========================================")
        message(" Gene-set expression heatmap")
        message("==========================================")
        message("Genes plotted:       ", nrow(matrix_scaled))
        message("Columns plotted:     ", ncol(matrix_scaled))
        message("Gene programs:       ", length(unique(gene_mapping$GeneSet)))
        if (is.null(aggregate_by)) {
            message("Expression level:    individual samples")
        } else {
            message("Aggregated by:       ", paste(aggregate_by, collapse = " + "))
        }
        message("Normalization:        ", normalization_label)
        if (row_zscore) {
            message("Scaling:             row z-score")
            message(
                "Z-score cap:         ",
                if (is.null(z_cap)) "none" else paste0("+/-", z_cap)
            )
        }
        message("")
        message("Genes per program:")
        print(table(gene_mapping$GeneSet))
        if (nrow(missing_gene_data) > 0L) {
            message("")
            message("Missing genes:")
            print(missing_gene_data)
        }
        message("==========================================")
    }

    invisible(list(
        heatmap = heatmap,
        matrix = matrix_scaled,
        expression = matrix_plot,
        logCPM = logcpm,
        gene_mapping = gene_mapping,
        missing_genes = missing_gene_data,
        metadata = column_metadata,
        row_split = row_split,
        col_fun = col_fun,
        dge = dge,
        normalization = bulk$normalization
    ))
}


#' Plot a differential-expression-selected expression heatmap
#'
#' Select the strongest negatively and positively changing genes from a
#' differential-expression table, optionally force specified genes into the
#' selection, and display their normalized sample-level expression in a
#' direction-split ComplexHeatmap. Optional batch correction affects the
#' visualization matrix only and never changes the supplied DE statistics.
#'
#' @param expr Numeric gene-by-sample matrix of normalized, preferably
#'   log-transformed expression values. Row names must match `gene_key_col` in
#'   `de`, and column names must match metadata sample identifiers. Do not pass
#'   raw counts directly.
#' @param meta Data frame containing one row per sample and the grouping,
#'   ordering, annotation, batch, covariate, and preservation fields requested
#'   below.
#' @param de Data frame containing one row per tested gene and at least the
#'   columns named by `gene_key_col`, `logfc_col`, and `padj_col`.
#' @param gene_key_col Character scalar naming the stable gene key in `de` that
#'   matches `rownames(expr)`.
#' @param gene_label_col `NULL` or a character scalar naming the preferred
#'   display-label column in `de`. Missing or blank labels fall back to the
#'   stable gene key.
#' @param logfc_col Character scalar naming the numeric log2-fold-change column.
#'   Positive values are assigned the positive side label and negative values
#'   the negative side label.
#' @param padj_col Character scalar naming the numeric adjusted-p-value column.
#'   Values must be finite and between zero and one.
#' @param sample_id_col `NULL` or a character scalar naming sample identifiers
#'   in `meta`. Metadata row names are used when this is `NULL`.
#' @param group_col Character scalar naming the biological group column in
#'   `meta`.
#' @param group_order `NULL` or a unique character vector defining biological
#'   group order. Observed groups omitted from the vector are appended in
#'   first-seen order; unobserved supplied values are dropped.
#' @param side_labels Named character vector containing unique `negative` and
#'   `positive` labels for the two directional row slices.
#' @param fdr_cutoff Numeric scalar in `(0, 1]` defining adjusted-p-value
#'   significance for automatic gene selection.
#' @param logfc_cutoff Non-negative numeric scalar defining the minimum absolute
#'   log2 fold change for automatic gene selection.
#' @param max_genes_per_side Non-negative integer giving the target maximum
#'   number of selected genes in each direction. Forced genes reserve slots but
#'   can exceed this maximum when more are requested.
#' @param force_genes `NULL` or a character vector of gene keys or display
#'   labels to include when present and directionally non-zero, regardless of
#'   the DE cutoffs.
#' @param rank_by Character scalar. `"padj"` ranks by adjusted p-value then
#'   absolute effect; `"abs_logFC"` reverses those priorities.
#' @param batch_col,batch2_col `NULL` or character scalars naming categorical
#'   nuisance variables removed from the visualization matrix with
#'   [limma::removeBatchEffect()].
#' @param covariate_cols `NULL` or a character vector naming numeric nuisance
#'   covariates removed from the visualization matrix.
#' @param preserve_cols `NULL` or a character vector naming biological effects
#'   represented in the design supplied to [limma::removeBatchEffect()]. These
#'   variables are preserved during visualization-only correction.
#' @param row_zscore Logical scalar. Whether to standardize each selected gene
#'   across samples. Zero-variance genes are removed before scaling.
#' @param z_cap `NULL` or a positive numeric scalar used to symmetrically cap
#'   displayed matrix values after optional row standardization.
#' @param order_by `NULL` or a character vector naming metadata columns used to
#'   order heatmap columns. The default uses `group_col`.
#' @param annotation_cols `NULL` or a character vector naming metadata columns
#'   shown above the heatmap. The default uses `group_col`.
#' @param annotation_colors `NULL` or a named list of annotation colour maps
#'   passed to [ComplexHeatmap::HeatmapAnnotation()].
#' @param show_annotation_name Logical scalar. Whether top-annotation names are
#'   displayed.
#' @param label_genes `NULL` or a character vector of selected gene keys or
#'   display labels to mark beside the heatmap.
#' @param label_top_per_side Non-negative integer giving the maximum number of
#'   automatically marked genes in each direction.
#' @param label_force_genes Logical scalar. Whether successfully selected forced
#'   genes are marked automatically.
#' @param label_side One of `"left"` or `"right"`, controlling the side used by
#'   the marked-gene row annotation.
#' @param label_fontsize Positive numeric scalar controlling marked-gene text.
#' @param label_italic Logical scalar. Whether marked-gene labels use italic
#'   text.
#' @param label_link_width_mm,label_padding_mm Non-negative numeric scalars
#'   controlling marked-gene connector width and label padding in millimetres.
#' @param cluster_rows Logical scalar. Whether rows are clustered within the
#'   heatmap.
#' @param cluster_row_slices Logical scalar. Whether the negative and positive
#'   directional slices are clustered relative to each other.
#' @param show_row_dend Logical scalar. Whether to display the row dendrogram.
#' @param cluster_columns Logical scalar. Whether sample columns are clustered.
#'   Enabling this can override the visible `order_by` arrangement.
#' @param show_column_dend Logical scalar. Whether to display the column
#'   dendrogram.
#' @param clustering_distance_rows Character scalar or distance function passed
#'   to [ComplexHeatmap::Heatmap()] for row clustering.
#' @param clustering_method_rows Character scalar naming the row-clustering
#'   linkage method passed to [ComplexHeatmap::Heatmap()].
#' @param heatmap_colors Character vector of exactly three valid colours for
#'   low, midpoint, and high values. The default is `DE_HEATMAP_COLORS`.
#' @param heatmap_name Non-empty character scalar used as the heatmap name.
#' @param show_column_names,show_row_names Logical scalars controlling sample
#'   and gene names drawn by the main heatmap body.
#' @param border Logical scalar. Whether to draw a heatmap border.
#' @param cell_border Character scalar specifying cell-border colour.
#' @param cell_border_lwd Non-negative numeric scalar controlling cell-border
#'   line width.
#' @param row_gap_mm Non-negative numeric scalar controlling directional-slice
#'   spacing in millimetres.
#' @param use_raster Logical scalar passed to [ComplexHeatmap::Heatmap()].
#' @param heatmap_legend_title Character scalar used as the heatmap legend title.
#' @param width,height Positive numeric scalars giving saved figure dimensions
#'   in inches.
#' @param save `NULL` or a character scalar giving an output path. PDF files use
#'   `grDevices::cairo_pdf`; other paths are written as 300-dpi PNG files.
#' @param draw_plot Logical scalar. Whether to draw the heatmap on the current
#'   graphics device in addition to any saved output.
#' @param verbose Logical scalar. Whether to report directional gene counts and
#'   visualization-correction details.
#'
#' @return Invisibly returns a named list containing:
#'   \describe{
#'     \item{heatmap}{The assembled `Heatmap` or `HeatmapList` object.}
#'     \item{matrix}{The selected, optionally corrected, row-standardized, and
#'       capped matrix displayed by the heatmap.}
#'     \item{expression_corrected}{The complete normalized input matrix after
#'       optional visualization-only batch correction and before selection or
#'       row standardization.}
#'     \item{selected_de}{The resolved DE rows used in the heatmap, including
#'       internal gene keys, labels, effects, adjusted p-values, directions,
#'       and unique display labels.}
#'     \item{metadata}{Expression-aligned metadata in displayed column order.}
#'     \item{row_split}{Directional factor aligned to heatmap rows.}
#'     \item{labels}{Display labels marked by `anno_mark`, or `character(0)`.}
#'     \item{col_fun}{The colour-mapping function used by the heatmap.}
#'   }
#'
#' @details
#' DE-based gene selection is a visualization choice and must not be reused as
#' an independent statistical test. Batch correction is applied only to the
#' expression matrix displayed here; the original DE table remains unchanged.
#' Required packages are `ComplexHeatmap` and `circlize`; `limma` is additionally
#' required when any visualization-correction variable is supplied.
#'
#' @examples
#' \dontrun{
#' heatmap_result <- plot_de_heatmap(
#'     expr = voom_result$voom$E,
#'     meta = project$metadata,
#'     de = voom_result$results,
#'     gene_key_col = "gene_id",
#'     gene_label_col = "gene_name",
#'     group_col = "Condition",
#'     group_order = c("Control", "Treated"),
#'     annotation_cols = c("Condition", "Batch"),
#'     max_genes_per_side = 25,
#'     draw_plot = TRUE,
#'     save = "outputs/de_heatmap.pdf"
#' )
#' heatmap_result$selected_de
#' heatmap_result$matrix
#' }
#'
#' @export
plot_de_heatmap <- function(
    expr,
    meta,
    de,
    gene_key_col = "gene",
    gene_label_col = NULL,
    logfc_col = "logFC",
    padj_col = "adj.P.Val",
    sample_id_col = NULL,
    group_col,
    group_order = NULL,
    side_labels = c(negative = "Lower", positive = "Higher"),
    fdr_cutoff = 0.05,
    logfc_cutoff = 0,
    max_genes_per_side = 25,
    force_genes = NULL,
    rank_by = c("padj", "abs_logFC"),
    batch_col = NULL,
    batch2_col = NULL,
    covariate_cols = NULL,
    preserve_cols = NULL,
    row_zscore = TRUE,
    z_cap = 2.5,
    order_by = NULL,
    annotation_cols = NULL,
    annotation_colors = NULL,
    show_annotation_name = FALSE,
    label_genes = NULL,
    label_top_per_side = 0,
    label_force_genes = TRUE,
    label_side = c("left", "right"),
    label_fontsize = 9,
    label_italic = TRUE,
    label_link_width_mm = 4,
    label_padding_mm = 1,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    heatmap_colors = DE_HEATMAP_COLORS,
    heatmap_name = "Row z-score",
    show_column_names = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    cell_border = "white",
    cell_border_lwd = 0.4,
    row_gap_mm = 1.5,
    use_raster = FALSE,
    heatmap_legend_title = "Expression\nrow z-score",
    width = 7,
    height = 8,
    save = NULL,
    draw_plot = TRUE,
    verbose = TRUE
) {
    required_packages <- c("ComplexHeatmap", "circlize")
    missing_packages <- required_packages[!vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )]
    if (length(missing_packages) > 0L) {
        stop(
            "Please install required package(s): ",
            paste(missing_packages, collapse = ", "), call. = FALSE
        )
    }
    correction_requested <- !is.null(batch_col) || !is.null(batch2_col) ||
        !is.null(covariate_cols)
    if (correction_requested && !requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required for removeBatchEffect().",
             call. = FALSE)
    }

    rank_by <- match.arg(rank_by)
    label_side <- match.arg(label_side)
    scalar_text <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    column_arguments <- list(
        gene_key_col = gene_key_col,
        logfc_col = logfc_col,
        padj_col = padj_col,
        group_col = group_col
    )
    invalid_columns <- names(column_arguments)[!vapply(
        column_arguments, scalar_text, logical(1)
    )]
    if (length(invalid_columns) > 0L ||
        (!is.null(gene_label_col) && !scalar_text(gene_label_col)) ||
        (!is.null(sample_id_col) && !scalar_text(sample_id_col)) ||
        (!is.null(batch_col) && !scalar_text(batch_col)) ||
        (!is.null(batch2_col) && !scalar_text(batch2_col))) {
        stop("Column-name arguments must be NULL or non-empty strings as allowed.",
             call. = FALSE)
    }
    if (!is.null(preserve_cols) && !correction_requested) {
        stop(
            "`preserve_cols` is only used when a batch or covariate correction ",
            "is requested.", call. = FALSE
        )
    }
    if (!(scalar_text(clustering_distance_rows) ||
        is.function(clustering_distance_rows)) ||
        !scalar_text(clustering_method_rows)) {
        stop("Row clustering distance and method arguments are invalid.",
             call. = FALSE)
    }

    logical_arguments <- list(
        row_zscore = row_zscore,
        show_annotation_name = show_annotation_name,
        label_force_genes = label_force_genes,
        label_italic = label_italic,
        cluster_rows = cluster_rows,
        cluster_row_slices = cluster_row_slices,
        show_row_dend = show_row_dend,
        cluster_columns = cluster_columns,
        show_column_dend = show_column_dend,
        show_column_names = show_column_names,
        show_row_names = show_row_names,
        border = border,
        use_raster = use_raster,
        draw_plot = draw_plot,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[!vapply(
        logical_arguments,
        function(x) is.logical(x) && length(x) == 1L && !is.na(x),
        logical(1)
    )]
    if (length(invalid_logical) > 0L) {
        stop(
            "These arguments must be single TRUE/FALSE values: ",
            paste(invalid_logical, collapse = ", "), call. = FALSE
        )
    }

    numeric_arguments <- list(
        fdr_cutoff = fdr_cutoff,
        logfc_cutoff = logfc_cutoff,
        label_fontsize = label_fontsize,
        label_link_width_mm = label_link_width_mm,
        label_padding_mm = label_padding_mm,
        cell_border_lwd = cell_border_lwd,
        row_gap_mm = row_gap_mm,
        width = width,
        height = height
    )
    invalid_numeric <- names(numeric_arguments)[!vapply(
        numeric_arguments,
        function(x) is.numeric(x) && length(x) == 1L &&
            !is.na(x) && is.finite(x),
        logical(1)
    )]
    if (length(invalid_numeric) > 0L) {
        stop(
            "These arguments must be finite numeric scalars: ",
            paste(invalid_numeric, collapse = ", "), call. = FALSE
        )
    }
    if (fdr_cutoff <= 0 || fdr_cutoff > 1 || logfc_cutoff < 0 ||
        label_fontsize <= 0 || label_link_width_mm < 0 ||
        label_padding_mm < 0 || cell_border_lwd < 0 || row_gap_mm < 0 ||
        width <= 0 || height <= 0) {
        stop("Cutoffs and display dimensions are outside their allowed ranges.",
             call. = FALSE)
    }
    integer_arguments <- list(
        max_genes_per_side = max_genes_per_side,
        label_top_per_side = label_top_per_side
    )
    invalid_integer <- names(integer_arguments)[!vapply(
        integer_arguments,
        function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
            is.finite(x) && x >= 0 && x == as.integer(x),
        logical(1)
    )]
    if (length(invalid_integer) > 0L) {
        stop(
            "These arguments must be non-negative integers: ",
            paste(invalid_integer, collapse = ", "), call. = FALSE
        )
    }
    if (!is.null(z_cap) && (!is.numeric(z_cap) || length(z_cap) != 1L ||
        is.na(z_cap) || !is.finite(z_cap) || z_cap <= 0)) {
        stop("`z_cap` must be NULL or a positive numeric scalar.", call. = FALSE)
    }
    if (!is.character(heatmap_colors) || length(heatmap_colors) != 3L ||
        anyNA(heatmap_colors) || any(!nzchar(heatmap_colors))) {
        stop("`heatmap_colors` must contain exactly three colour strings.",
             call. = FALSE)
    }
    if (!scalar_text(heatmap_name) || !scalar_text(heatmap_legend_title) ||
        !scalar_text(cell_border) ||
        (!is.null(save) && !scalar_text(save))) {
        stop("Heatmap text and output-path arguments must be non-empty strings.",
             call. = FALSE)
    }

    expr <- as.matrix(expr)
    meta <- as.data.frame(meta)
    de_original <- as.data.frame(de)
    if (!is.numeric(expr) || nrow(expr) == 0L || ncol(expr) == 0L) {
        stop("`expr` must be a non-empty numeric gene-by-sample matrix.",
             call. = FALSE)
    }
    if (is.null(rownames(expr)) || is.null(colnames(expr)) ||
        anyNA(rownames(expr)) || anyNA(colnames(expr)) ||
        any(!nzchar(rownames(expr))) || any(!nzchar(colnames(expr)))) {
        stop("`expr` requires non-empty gene row names and sample column names.",
             call. = FALSE)
    }
    if (anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) {
        stop("Expression gene and sample identifiers must be unique.",
             call. = FALSE)
    }
    if (anyNA(expr) || any(!is.finite(expr))) {
        stop("`expr` must contain only finite, non-missing values.",
             call. = FALSE)
    }
    if (row_zscore && ncol(expr) < 2L) {
        stop("At least two expression samples are required for row z-scores.",
             call. = FALSE)
    }

    required_de_columns <- c(gene_key_col, logfc_col, padj_col)
    missing_de_columns <- setdiff(required_de_columns, colnames(de_original))
    if (length(missing_de_columns) > 0L) {
        stop(
            "DE table is missing: ",
            paste(missing_de_columns, collapse = ", "), call. = FALSE
        )
    }
    if (!group_col %in% colnames(meta)) {
        stop("`group_col = \"", group_col, "\"` was not found in metadata.",
             call. = FALSE)
    }

    if (is.null(sample_id_col)) {
        if (is.null(rownames(meta))) {
            stop("Metadata needs row names or `sample_id_col`.", call. = FALSE)
        }
        meta$.SampleID <- rownames(meta)
    } else {
        if (!sample_id_col %in% colnames(meta)) {
            stop("sample_id_col '", sample_id_col, "' not found in metadata.",
                 call. = FALSE)
        }
        meta$.SampleID <- as.character(meta[[sample_id_col]])
    }
    if (anyNA(meta$.SampleID) || any(!nzchar(meta$.SampleID)) ||
        anyDuplicated(meta$.SampleID)) {
        stop("Metadata sample IDs must be unique, non-empty, and non-missing.",
             call. = FALSE)
    }
    missing_meta_samples <- setdiff(colnames(expr), meta$.SampleID)
    if (length(missing_meta_samples) > 0L) {
        stop(
            "Some expression samples are absent from metadata: ",
            paste(head(missing_meta_samples, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    meta <- meta[match(colnames(expr), meta$.SampleID), , drop = FALSE]
    rownames(meta) <- meta$.SampleID

    observed_groups <- as.character(meta[[group_col]])
    if (anyNA(observed_groups) || any(!nzchar(observed_groups))) {
        stop("The selected group column contains missing or empty values.",
             call. = FALSE)
    }
    observed_groups <- unique(observed_groups)
    if (is.null(group_order)) {
        group_order <- observed_groups
    } else {
        if (!is.character(group_order) || anyNA(group_order) ||
            any(!nzchar(group_order)) || anyDuplicated(group_order)) {
            stop("`group_order` must contain unique non-empty strings.",
                 call. = FALSE)
        }
        group_order <- c(group_order, setdiff(observed_groups, group_order))
        group_order <- group_order[group_order %in% observed_groups]
    }
    meta[[group_col]] <- factor(
        as.character(meta[[group_col]]), levels = group_order
    )

    if (!is.character(side_labels) ||
        !all(c("negative", "positive") %in% names(side_labels))) {
        stop("`side_labels` must contain names 'negative' and 'positive'.",
             call. = FALSE)
    }
    negative_label <- unname(side_labels[["negative"]])
    positive_label <- unname(side_labels[["positive"]])
    if (!scalar_text(negative_label) || !scalar_text(positive_label) ||
        identical(negative_label, positive_label)) {
        stop("Negative and positive side labels must be distinct strings.",
             call. = FALSE)
    }

    de_data <- de_original
    de_data$.GeneKey <- as.character(de_data[[gene_key_col]])
    de_data$.GeneLabel <- if (!is.null(gene_label_col) &&
        gene_label_col %in% colnames(de_data)) {
        as.character(de_data[[gene_label_col]])
    } else {
        de_data$.GeneKey
    }
    bad_labels <- is.na(de_data$.GeneLabel) | !nzchar(de_data$.GeneLabel)
    de_data$.GeneLabel[bad_labels] <- de_data$.GeneKey[bad_labels]
    de_data$.logFC <- suppressWarnings(as.numeric(de_data[[logfc_col]]))
    de_data$.padj <- suppressWarnings(as.numeric(de_data[[padj_col]]))
    usable <- !is.na(de_data$.GeneKey) & nzchar(de_data$.GeneKey) &
        is.finite(de_data$.logFC) & is.finite(de_data$.padj) &
        de_data$.padj >= 0 & de_data$.padj <= 1
    if (verbose && any(!usable)) {
        message(sum(!usable), " unusable DE rows were ignored.")
    }
    de_data <- de_data[usable, , drop = FALSE]
    if (nrow(de_data) == 0L) {
        stop("No usable rows remain in the DE table.", call. = FALSE)
    }
    de_order <- order(de_data$.padj, -abs(de_data$.logFC))
    de_data <- de_data[de_order, , drop = FALSE]
    de_data <- de_data[!duplicated(de_data$.GeneKey), , drop = FALSE]

    available <- de_data$.GeneKey %in% rownames(expr)
    if (verbose && any(!available)) {
        message(sum(!available), " DE genes were not present in expr and were ignored.")
    }
    de_data <- de_data[available, , drop = FALSE]
    if (nrow(de_data) == 0L) {
        stop("No DE genes were present in the expression matrix.", call. = FALSE)
    }
    de_data$.Direction <- NA_character_
    de_data$.Direction[de_data$.logFC < 0] <- negative_label
    de_data$.Direction[de_data$.logFC > 0] <- positive_label

    rank_de <- function(data) {
        if (nrow(data) == 0L) {
            return(data)
        }
        row_order <- if (rank_by == "padj") {
            order(data$.padj, -abs(data$.logFC))
        } else {
            order(-abs(data$.logFC), data$.padj)
        }
        data[row_order, , drop = FALSE]
    }
    bind_unique <- function(...) {
        inputs <- list(...)
        inputs <- inputs[vapply(inputs, nrow, integer(1)) > 0L]
        if (length(inputs) == 0L) {
            return(de_data[FALSE, , drop = FALSE])
        }
        output <- do.call(rbind, inputs)
        output[!duplicated(output$.GeneKey), , drop = FALSE]
    }

    force_genes <- if (is.null(force_genes)) character(0) else
        unique(as.character(force_genes))
    force_genes <- force_genes[!is.na(force_genes) & nzchar(force_genes)]
    force_match <- de_data$.GeneKey %in% force_genes |
        de_data$.GeneLabel %in% force_genes
    force_info <- de_data[force_match & !is.na(de_data$.Direction), , drop = FALSE]
    if (length(force_genes) > 0L && verbose) {
        matched_force <- unique(c(
            force_info$.GeneKey[force_info$.GeneKey %in% force_genes],
            force_info$.GeneLabel[force_info$.GeneLabel %in% force_genes]
        ))
        missing_force <- setdiff(force_genes, matched_force)
        if (length(missing_force) > 0L) {
            warning(
                "Forced genes not found or lacking a non-zero direction: ",
                paste(missing_force, collapse = ", "), call. = FALSE
            )
        }
    }
    force_negative <- force_info[
        force_info$.Direction == negative_label, , drop = FALSE
    ]
    force_positive <- force_info[
        force_info$.Direction == positive_label, , drop = FALSE
    ]

    significant_negative <- de_data[
        de_data$.padj < fdr_cutoff & de_data$.logFC < -logfc_cutoff &
            !de_data$.GeneKey %in% force_info$.GeneKey,
        , drop = FALSE
    ]
    significant_positive <- de_data[
        de_data$.padj < fdr_cutoff & de_data$.logFC > logfc_cutoff &
            !de_data$.GeneKey %in% force_info$.GeneKey,
        , drop = FALSE
    ]
    significant_negative <- rank_de(significant_negative)
    significant_positive <- rank_de(significant_positive)
    negative_slots <- max(0L, as.integer(max_genes_per_side) -
        nrow(force_negative))
    positive_slots <- max(0L, as.integer(max_genes_per_side) -
        nrow(force_positive))
    selected_negative <- bind_unique(
        head(significant_negative, negative_slots), force_negative
    )
    selected_positive <- bind_unique(
        head(significant_positive, positive_slots), force_positive
    )
    selected_de <- bind_unique(selected_negative, selected_positive)
    if (nrow(selected_de) == 0L) {
        stop("No DE genes passed the requested thresholds.", call. = FALSE)
    }

    expr_visual <- expr
    if (correction_requested) {
        correction_columns <- c(batch_col, batch2_col, covariate_cols)
        correction_columns <- correction_columns[!is.na(correction_columns)]
        missing_correction <- setdiff(correction_columns, colnames(meta))
        if (length(missing_correction) > 0L) {
            stop(
                "Missing visualization-correction columns: ",
                paste(missing_correction, collapse = ", "), call. = FALSE
            )
        }
        for (batch_name in c(batch_col, batch2_col)) {
            if (!is.null(batch_name)) {
                batch_values <- meta[[batch_name]]
                if (anyNA(batch_values) || length(unique(batch_values)) < 2L) {
                    stop("Batch columns require complete values and two levels.",
                         call. = FALSE)
                }
            }
        }
        covariate_matrix <- NULL
        if (!is.null(covariate_cols)) {
            if (!is.character(covariate_cols) || anyDuplicated(covariate_cols)) {
                stop("`covariate_cols` must contain unique column names.",
                     call. = FALSE)
            }
            numeric_covariates <- vapply(
                meta[covariate_cols], is.numeric, logical(1)
            )
            if (!all(numeric_covariates) ||
                anyNA(meta[covariate_cols]) ||
                any(!is.finite(as.matrix(meta[covariate_cols])))) {
                stop("Visualization covariates must be complete finite numerics.",
                     call. = FALSE)
            }
            covariate_matrix <- as.matrix(meta[covariate_cols])
        }

        preserve_design <- NULL
        if (!is.null(preserve_cols)) {
            if (!is.character(preserve_cols) || anyNA(preserve_cols) ||
                any(!nzchar(preserve_cols)) || anyDuplicated(preserve_cols)) {
                stop("`preserve_cols` must contain unique non-empty names.",
                     call. = FALSE)
            }
            missing_preserve <- setdiff(preserve_cols, colnames(meta))
            if (length(missing_preserve) > 0L) {
                stop(
                    "Missing preserve_cols: ",
                    paste(missing_preserve, collapse = ", "), call. = FALSE
                )
            }
            preserve_design <- stats::model.matrix(
                stats::reformulate(preserve_cols), data = meta
            )
        }
        expr_visual <- limma::removeBatchEffect(
            expr,
            batch = if (is.null(batch_col)) NULL else meta[[batch_col]],
            batch2 = if (is.null(batch2_col)) NULL else meta[[batch2_col]],
            covariates = covariate_matrix,
            design = preserve_design
        )
        if (anyNA(expr_visual) || any(!is.finite(expr_visual))) {
            stop("Visualization correction produced non-finite values.",
                 call. = FALSE)
        }
    }

    heatmap_matrix <- expr_visual[selected_de$.GeneKey, , drop = FALSE]
    if (is.null(order_by)) {
        order_by <- group_col
    }
    if (!is.character(order_by) || length(order_by) == 0L ||
        anyNA(order_by) || any(!nzchar(order_by)) || anyDuplicated(order_by)) {
        stop("`order_by` must contain unique non-empty metadata columns.",
             call. = FALSE)
    }
    missing_order_columns <- setdiff(order_by, colnames(meta))
    if (length(missing_order_columns) > 0L) {
        stop(
            "Missing order_by columns: ",
            paste(missing_order_columns, collapse = ", "), call. = FALSE
        )
    }
    if (anyNA(meta[order_by])) {
        stop("`order_by` columns cannot contain missing values.", call. = FALSE)
    }
    column_order_index <- do.call(order, c(meta[order_by], list(na.last = TRUE)))
    column_order <- meta$.SampleID[column_order_index]
    heatmap_matrix <- heatmap_matrix[, column_order, drop = FALSE]
    column_metadata <- meta[column_order, , drop = FALSE]

    if (row_zscore) {
        row_sd <- apply(heatmap_matrix, 1L, stats::sd)
        zero_variance <- !is.finite(row_sd) | row_sd == 0
        if (any(zero_variance)) {
            if (verbose) {
                warning(
                    sum(zero_variance),
                    " selected genes had zero variance and were removed ",
                    "before z-scoring.", call. = FALSE
                )
            }
            heatmap_matrix <- heatmap_matrix[!zero_variance, , drop = FALSE]
            selected_de <- selected_de[
                selected_de$.GeneKey %in% rownames(heatmap_matrix), , drop = FALSE
            ]
        }
        if (nrow(heatmap_matrix) == 0L) {
            stop("No variable selected genes remain for row z-scoring.",
                 call. = FALSE)
        }
        heatmap_scaled <- t(scale(t(heatmap_matrix)))
    } else {
        heatmap_scaled <- heatmap_matrix
    }
    if (!is.null(z_cap)) {
        heatmap_scaled <- pmax(pmin(heatmap_scaled, z_cap), -z_cap)
    }

    selected_de <- selected_de[
        match(rownames(heatmap_scaled), selected_de$.GeneKey), , drop = FALSE
    ]
    display_labels <- make.unique(selected_de$.GeneLabel)
    rownames(heatmap_scaled) <- display_labels
    selected_de$.DisplayLabel <- display_labels
    row_split <- factor(
        selected_de$.Direction, levels = c(negative_label, positive_label)
    )

    final_label_genes <- if (is.null(label_genes)) character(0) else
        unique(as.character(label_genes))
    final_label_genes <- final_label_genes[
        !is.na(final_label_genes) & nzchar(final_label_genes)
    ]
    if (label_top_per_side > 0L) {
        top_negative <- head(
            rank_de(selected_de[selected_de$.Direction == negative_label,
                                , drop = FALSE]),
            as.integer(label_top_per_side)
        )
        top_positive <- head(
            rank_de(selected_de[selected_de$.Direction == positive_label,
                                , drop = FALSE]),
            as.integer(label_top_per_side)
        )
        final_label_genes <- unique(c(
            final_label_genes,
            top_negative$.GeneKey,
            top_negative$.GeneLabel,
            top_positive$.GeneKey,
            top_positive$.GeneLabel
        ))
    }
    if (label_force_genes && length(force_genes) > 0L) {
        final_label_genes <- unique(c(final_label_genes, force_genes))
    }
    label_keep <- selected_de$.GeneKey %in% final_label_genes |
        selected_de$.GeneLabel %in% final_label_genes
    label_indices <- which(label_keep)
    label_text <- selected_de$.GeneLabel[label_indices]
    gene_mark <- NULL
    if (length(label_indices) > 0L) {
        gene_mark <- ComplexHeatmap::rowAnnotation(
            Gene = ComplexHeatmap::anno_mark(
                at = label_indices,
                labels = label_text,
                which = "row",
                side = label_side,
                labels_gp = grid::gpar(
                    fontsize = label_fontsize,
                    fontface = if (label_italic) "italic" else "plain"
                ),
                link_gp = grid::gpar(col = "black", lwd = 0.6),
                link_width = grid::unit(label_link_width_mm, "mm"),
                padding = grid::unit(label_padding_mm, "mm"),
                extend = grid::unit(c(2, 2), "mm")
            )
        )
    }

    if (is.null(annotation_cols)) {
        annotation_cols <- group_col
    }
    if (!is.character(annotation_cols) || length(annotation_cols) == 0L ||
        anyNA(annotation_cols) || any(!nzchar(annotation_cols)) ||
        anyDuplicated(annotation_cols)) {
        stop("`annotation_cols` must contain unique non-empty metadata columns.",
             call. = FALSE)
    }
    missing_annotation_columns <- setdiff(annotation_cols, colnames(column_metadata))
    if (length(missing_annotation_columns) > 0L) {
        stop(
            "Missing annotation columns: ",
            paste(missing_annotation_columns, collapse = ", "), call. = FALSE
        )
    }
    if (!is.null(annotation_colors) &&
        (!is.list(annotation_colors) || is.null(names(annotation_colors)) ||
         anyNA(names(annotation_colors)) || any(!nzchar(names(annotation_colors))))) {
        stop("`annotation_colors` must be NULL or a named list.", call. = FALSE)
    }
    annotation_data <- column_metadata[, annotation_cols, drop = FALSE]
    annotation_arguments <- list(
        df = annotation_data,
        show_annotation_name = show_annotation_name
    )
    if (!is.null(annotation_colors)) {
        annotation_arguments$col <- annotation_colors
    }
    top_annotation <- do.call(
        ComplexHeatmap::HeatmapAnnotation, annotation_arguments
    )

    if (!is.null(z_cap)) {
        color_limit <- z_cap
    } else {
        color_limit <- max(abs(heatmap_scaled))
        if (!is.finite(color_limit) || color_limit == 0) {
            color_limit <- 1
        }
    }
    color_function <- circlize::colorRamp2(
        c(-color_limit, 0, color_limit), heatmap_colors
    )
    heatmap <- ComplexHeatmap::Heatmap(
        heatmap_scaled,
        name = heatmap_name,
        col = color_function,
        top_annotation = top_annotation,
        split = row_split,
        row_title = NULL,
        row_gap = grid::unit(row_gap_mm, "mm"),
        cluster_rows = cluster_rows,
        cluster_row_slices = cluster_row_slices,
        clustering_distance_rows = clustering_distance_rows,
        clustering_method_rows = clustering_method_rows,
        show_row_dend = show_row_dend,
        cluster_columns = cluster_columns,
        show_column_dend = show_column_dend,
        show_column_names = show_column_names,
        show_row_names = show_row_names,
        border = border,
        rect_gp = grid::gpar(col = cell_border, lwd = cell_border_lwd),
        use_raster = use_raster,
        heatmap_legend_param = list(
            title = heatmap_legend_title,
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 9)
        )
    )
    heatmap_final <- if (is.null(gene_mark)) {
        heatmap
    } else if (label_side == "left") {
        gene_mark + heatmap
    } else {
        heatmap + gene_mark
    }

    draw_heatmap <- function() {
        ComplexHeatmap::draw(
            heatmap_final,
            heatmap_legend_side = "right",
            annotation_legend_side = "right"
        )
    }
    if (!is.null(save)) {
        save_heatmap <- function() {
            if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
                grDevices::cairo_pdf(filename = if (is.null(save)) NA else save, width = width, height = height)
            } else {
                grDevices::png(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300
                )
            }
            on.exit(grDevices::dev.off(), add = TRUE)
            draw_heatmap()
        }
        save_heatmap()
    }
    if (draw_plot) {
        draw_heatmap()
    }

    if (verbose) {
        message("")
        message("==========================================")
        message(" DE heatmap")
        message("==========================================")
        message(
            negative_label, ": ",
            sum(selected_de$.Direction == negative_label), " genes"
        )
        message(
            positive_label, ": ",
            sum(selected_de$.Direction == positive_label), " genes"
        )
        message("Total genes: ", nrow(heatmap_scaled))
        corrected_for <- c(batch_col, batch2_col, covariate_cols)
        corrected_for <- corrected_for[!is.na(corrected_for)]
        if (length(corrected_for) > 0L) {
            message(
                "Visualization corrected for: ",
                paste(corrected_for, collapse = ", ")
            )
        }
        if (!is.null(preserve_cols)) {
            message("Preserved effects: ", paste(preserve_cols, collapse = ", "))
        }
        message("==========================================")
    }

    invisible(list(
        heatmap = heatmap_final,
        matrix = heatmap_scaled,
        expression_corrected = expr_visual,
        selected_de = selected_de,
        metadata = column_metadata,
        row_split = row_split,
        labels = label_text,
        col_fun = color_function
    ))
}


#' Plot signed Manhattan-style summaries of pairwise differential expression
#'
#' Combine one or more pairwise differential-expression results and display
#' every tested gene by an x-axis ordering and signed adjusted-p-value score.
#' Positive values denote higher expression in the comparison group; negative
#' values denote higher expression in the reference group. Comparisons are
#' shown in columns and can additionally be faceted by sorting or cell group.
#'
#' @param de_list A non-empty list of [edgeR_pairwise()] objects or result data
#'   frames. List names identify analyses. Alternatively, a data frame returned
#'   by [collect_edgeR_pairwise()] can be supplied directly.
#' @param name_sep Non-empty character scalar separating fields in analysis
#'   names. Text before the first separator is used as `Sorting` when collecting
#'   a list.
#' @param gene_col Character scalar naming the preferred gene-label column.
#'   Missing or empty labels fall back to the gene key.
#' @param gene_id_col Character scalar naming the preferred stable gene-ID
#'   column. When absent, original result row names are used.
#' @param logfc_col Character scalar naming the numeric log2-fold-change column.
#' @param padj_col Character scalar naming the numeric adjusted-p-value column.
#'   Observed values must be finite and between zero and one.
#' @param sorting_order `NULL` or a unique character vector ordering sorting
#'   groups. Observed groups omitted from the vector are appended in first-seen
#'   order.
#' @param comparison_order `NULL` or a unique character vector ordering raw
#'   comparison identifiers. Omitted observed comparisons are appended.
#' @param comparison_labels `NULL`, an unnamed character vector aligned to
#'   `comparison_order`, or a named character vector mapping raw comparisons to
#'   unique display labels.
#' @param gene_order Character scalar selecting x-axis order: `"mean_logFC"`
#'   uses each gene's mean effect across all panels; `"alphabetical"` sorts by
#'   gene label; `"input"` keeps first appearance; and `"panel_logFC"` orders
#'   effects independently within each sorting-by-comparison panel.
#' @param fdr_cutoff Numeric scalar in `(0, 1]` defining adjusted-p-value
#'   significance and the optional horizontal threshold lines.
#' @param logfc_cutoff Non-negative numeric scalar defining the minimum absolute
#'   log2 fold change for significant-direction labels. A value of zero uses
#'   strictly positive or negative effects.
#' @param colors `NULL` or a character colour vector. Unnamed values are assigned
#'   in display-label order. Named values may use either raw comparison names or
#'   display labels. The default uses `SIGNED_MANHATTAN_MACARON_COLORS` for up
#'   to 12 comparisons and a pastel HCL palette thereafter.
#' @param ns_alpha,sig_alpha Numeric scalars in `[0, 1]` controlling opacity for
#'   non-significant and significant genes.
#' @param point_size Positive numeric scalar controlling point size.
#' @param rasterized Logical scalar. Whether to use `ggrastr` point layers;
#'   when `ggrastr` is unavailable, a warning is issued and vector points are
#'   used.
#' @param raster_dpi Positive integer passed to `ggrastr` as raster resolution.
#' @param label_top_up,label_top_down Non-negative integers giving the maximum
#'   number of positive- and negative-effect genes labelled per panel, ranked
#'   first by smallest adjusted p-value and then by larger absolute effect in
#'   the requested direction.
#' @param label_genes `NULL` or a character vector of additional gene labels or
#'   stable IDs to label in every matching panel.
#' @param label_only_significant Logical scalar. Whether automatic top-gene
#'   labels are restricted to genes meeting both cutoffs. Manually requested
#'   genes are not restricted.
#' @param label_size Positive numeric scalar controlling gene-label text size.
#' @param italic_gene_labels Logical scalar. Whether gene labels use italic text.
#' @param cap_y `NULL` or a positive numeric scalar. Scores outside
#'   `[-cap_y, cap_y]` are clipped for display and marked in the returned
#'   `Capped` column; uncapped values remain in `SignedLogFDR`.
#' @param show_capped_triangles Logical scalar. When `cap_y` is used, whether
#'   upward and downward open triangles should mark genes clipped at the upper
#'   and lower display boundaries.
#' @param facet_by_sorting Logical scalar. When multiple sorting groups are
#'   present, whether to use sorting rows and comparison columns.
#' @param facet_ncol `NULL` or a positive integer controlling comparison columns
#'   when `facet_wrap()` is used. The default places all comparisons in one row.
#' @param facet_scales One of `"fixed"`, `"free"`, `"free_x"`, or `"free_y"`.
#' @param show_zero_line,show_fdr_line Logical scalars controlling the neutral
#'   zero reference and symmetric adjusted-p-value cutoff lines.
#' @param title `NULL` or a character scalar used as the plot title.
#' @param x_title Character scalar used as the x-axis title.
#' @param width Positive numeric scalar giving saved figure width in inches.
#' @param height `NULL` or a positive numeric scalar giving saved figure height
#'   in inches. The default is 4 inches for a single row or at least 2.4
#'   inches per sorting row.
#' @param save `NULL` or a character scalar giving an output figure path. PDF
#'   files use `grDevices::cairo_pdf`; other formats are written at 300 dpi.
#' @param verbose Logical scalar. Whether to print analysis counts, cutoffs, and
#'   a per-panel significance summary.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{plot}{The assembled `ggplot` object.}
#'     \item{data}{The complete plotting data with gene keys, display labels,
#'       x positions, signed scores, significance classes, and cap indicators.}
#'     \item{labels}{Rows selected for gene labelling.}
#'     \item{summary}{Per-sorting/per-comparison counts of tested, significant,
#'       comparison-up, reference-up, and capped genes.}
#'     \item{label_summary}{Per-panel counts of selected labels, or an empty
#'       data frame when no genes are labelled.}
#'     \item{colors}{Named colours used for displayed comparisons.}
#'     \item{comparison_labels}{Named raw-to-display comparison mapping.}
#'     \item{fdr_cutoff}{Adjusted-p-value cutoff used by the plot.}
#'     \item{logfc_cutoff}{Absolute log2-fold-change cutoff used by the plot.}
#'     \item{cap_y}{The requested signed-score display cap, or `NULL`.}
#'   }
#'
#' @details
#' Exact adjusted p-values of zero are replaced only for plotting by one tenth
#' of the smallest positive value in the combined data, bounded by machine
#' precision. Original adjusted p-values and fold changes are retained. The
#' signed display is a compact exploratory overview across analyses; effect
#' sizes, uncertainty, replication, model design, and full result tables remain
#' necessary for biological interpretation. Required packages are `ggplot2`
#' and, when labels are requested, `ggrepel`.
#'
#' @examples
#' \dontrun{
#' signed_plot <- plot_signed_manhattan(
#'     list(
#'         Myeloid__Treated_vs_Control = myeloid_result,
#'         Lymphoid__Treated_vs_Control = lymphoid_result
#'     ),
#'     comparison_order = "Treated_vs_Control",
#'     comparison_labels = c(Treated_vs_Control = "Treated vs Control"),
#'     label_top_up = 5,
#'     label_top_down = 5,
#'     cap_y = 25,
#'     save = "outputs/signed_manhattan.pdf"
#' )
#' signed_plot$plot
#' signed_plot$summary
#' }
#'
#' @export
plot_signed_manhattan <- function(
    de_list,
    name_sep = "__",
    gene_col = "gene_name",
    gene_id_col = "gene_id",
    logfc_col = "logFC",
    padj_col = "FDR",
    sorting_order = NULL,
    comparison_order = NULL,
    comparison_labels = NULL,
    gene_order = c("alphabetical", "mean_logFC", "input", "panel_logFC"),
    fdr_cutoff = 0.05,
    logfc_cutoff = 0,
    colors = NULL,
    ns_alpha = 0.20,
    sig_alpha = 0.85,
    point_size = 0.5,
    rasterized = FALSE,
    raster_dpi = 300,
    label_top_up = 5,
    label_top_down = 5,
    label_genes = NULL,
    label_only_significant = TRUE,
    label_size = 2.4,
    italic_gene_labels = TRUE,
    cap_y = NULL,
    show_capped_triangles = TRUE,
    facet_by_sorting = TRUE,
    facet_ncol = NULL,
    facet_scales = "fixed",
    show_zero_line = TRUE,
    show_fdr_line = TRUE,
    title = NULL,
    x_title = "Genes",
    width = 11,
    height = NULL,
    save = NULL,
    verbose = TRUE
) {
    gene_order <- match.arg(gene_order)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    scalar_character <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }
    scalar_logical <- function(x) {
        is.logical(x) && length(x) == 1L && !is.na(x)
    }
    column_arguments <- list(
        name_sep = name_sep,
        gene_col = gene_col,
        gene_id_col = gene_id_col,
        logfc_col = logfc_col,
        padj_col = padj_col,
        x_title = x_title
    )
    invalid_columns <- names(column_arguments)[
        !vapply(column_arguments, scalar_character, logical(1))
    ]
    if (length(invalid_columns) > 0L) {
        stop(
            paste(invalid_columns, collapse = ", "),
            " must be non-empty character scalars."
        )
    }
    if (!is.null(title) && !scalar_character(title)) {
        stop("title must be NULL or one non-empty character string.")
    }
    if (!is.null(save) && !scalar_character(save)) {
        stop("save must be NULL or one non-empty character path.")
    }
    logical_arguments <- list(
        rasterized = rasterized,
        label_only_significant = label_only_significant,
        italic_gene_labels = italic_gene_labels,
        show_capped_triangles = show_capped_triangles,
        facet_by_sorting = facet_by_sorting,
        show_zero_line = show_zero_line,
        show_fdr_line = show_fdr_line,
        verbose = verbose
    )
    invalid_logical <- names(logical_arguments)[
        !vapply(logical_arguments, scalar_logical, logical(1))
    ]
    if (length(invalid_logical) > 0L) {
        stop(
            paste(invalid_logical, collapse = ", "),
            " must be non-missing logical scalars."
        )
    }
    validate_positive <- function(value, argument_name, allow_null = FALSE) {
        if (allow_null && is.null(value)) {
            return(invisible(NULL))
        }
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) || value <= 0) {
            stop(argument_name, " must be a positive finite numeric scalar.")
        }
        invisible(NULL)
    }
    validate_positive(point_size, "point_size")
    validate_positive(label_size, "label_size")
    validate_positive(width, "width")
    validate_positive(height, "height", allow_null = TRUE)
    validate_positive(cap_y, "cap_y", allow_null = TRUE)
    for (alpha_name in c("ns_alpha", "sig_alpha")) {
        alpha <- get(alpha_name)
        if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
            !is.finite(alpha) || alpha < 0 || alpha > 1) {
            stop(alpha_name, " must be a finite numeric scalar in [0, 1].")
        }
    }
    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        is.na(fdr_cutoff) || !is.finite(fdr_cutoff) ||
        fdr_cutoff <= 0 || fdr_cutoff > 1) {
        stop("fdr_cutoff must be a finite numeric scalar in (0, 1].")
    }
    if (!is.numeric(logfc_cutoff) || length(logfc_cutoff) != 1L ||
        is.na(logfc_cutoff) || !is.finite(logfc_cutoff) || logfc_cutoff < 0) {
        stop("logfc_cutoff must be a finite non-negative numeric scalar.")
    }
    integer_arguments <- list(
        raster_dpi = raster_dpi,
        label_top_up = label_top_up,
        label_top_down = label_top_down
    )
    for (argument_name in names(integer_arguments)) {
        value <- integer_arguments[[argument_name]]
        minimum <- if (argument_name == "raster_dpi") 1 else 0
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) || value < minimum || value != round(value)) {
            stop(argument_name, " has an invalid integer value.")
        }
    }
    if (!is.null(facet_ncol) &&
        (!is.numeric(facet_ncol) || length(facet_ncol) != 1L ||
            is.na(facet_ncol) || !is.finite(facet_ncol) ||
            facet_ncol < 1 || facet_ncol != round(facet_ncol))) {
        stop("facet_ncol must be NULL or a positive integer.")
    }
    facet_scales <- match.arg(
        facet_scales, c("fixed", "free", "free_x", "free_y")
    )
    if (!is.null(label_genes) &&
        (!is.character(label_genes) || anyNA(label_genes) ||
            any(!nzchar(label_genes)))) {
        stop("label_genes must be NULL or non-empty character identifiers.")
    }

    use_labels <- label_top_up > 0 || label_top_down > 0 ||
        !is.null(label_genes)
    if (use_labels && !requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required for gene labels.")
    }
    if (rasterized && !requireNamespace("ggrastr", quietly = TRUE)) {
        warning("ggrastr is not installed; using normal ggplot2 points instead.")
        rasterized <- FALSE
    }

    if (is.data.frame(de_list)) {
        df <- de_list
        required_collection <- c(
            "Analysis", "Sorting", "Reference", "ComparisonGroup", "Comparison"
        )
        missing_collection <- setdiff(required_collection, colnames(df))
        if (length(missing_collection) > 0L) {
            stop(
                "A collected data frame is missing columns: ",
                paste(missing_collection, collapse = ", ")
            )
        }
        analysis_count <- length(unique(df$Analysis))
    } else {
        df <- .collect_pairwise_results(
            de_list = de_list,
            name_sep = name_sep,
            logfc_col = logfc_col,
            padj_col = padj_col,
            allow_data_frames = TRUE,
            require_names = FALSE
        )
        analysis_count <- length(de_list)
    }
    missing_statistics <- setdiff(c(logfc_col, padj_col), colnames(df))
    if (length(missing_statistics) > 0L) {
        stop(
            "Combined results are missing columns: ",
            paste(missing_statistics, collapse = ", ")
        )
    }
    if (!is.numeric(df[[logfc_col]]) || !is.numeric(df[[padj_col]])) {
        stop("logfc_col and padj_col must identify numeric columns.")
    }
    observed_padj <- df[[padj_col]][!is.na(df[[padj_col]])]
    if (any(!is.finite(observed_padj)) ||
        any(observed_padj < 0 | observed_padj > 1)) {
        stop("Observed adjusted p-values must be finite and in [0, 1].")
    }
    df <- df[
        is.finite(df[[logfc_col]]) & !is.na(df[[padj_col]]),
        ,
        drop = FALSE
    ]
    if (nrow(df) == 0L) {
        stop("No rows with finite log fold changes and observed FDR remain.")
    }

    gene_key <- if (gene_id_col %in% colnames(df)) {
        as.character(df[[gene_id_col]])
    } else if ("GeneRowName" %in% colnames(df)) {
        as.character(df$GeneRowName)
    } else {
        rownames(df)
    }
    bad_key <- is.na(gene_key) | !nzchar(gene_key)
    gene_key[bad_key] <- paste0(
        as.character(df$Analysis[bad_key]), "__row_", which(bad_key)
    )
    gene_label <- if (gene_col %in% colnames(df)) {
        as.character(df[[gene_col]])
    } else {
        gene_key
    }
    bad_label <- is.na(gene_label) | !nzchar(gene_label)
    gene_label[bad_label] <- gene_key[bad_label]
    df$.GeneKey <- gene_key
    df$.GeneLabel <- gene_label

    positive_fdr <- df[[padj_col]][df[[padj_col]] > 0]
    replacement_fdr <- if (length(positive_fdr) > 0L) {
        max(min(positive_fdr) / 10, .Machine$double.xmin)
    } else {
        .Machine$double.xmin
    }
    df$.PlotFDR <- df[[padj_col]]
    df$.PlotFDR[df$.PlotFDR <= 0] <- replacement_fdr
    df$SignedLogFDR <- -log10(df$.PlotFDR) * sign(df[[logfc_col]])
    df$Significance <- "NS"
    significant <- df[[padj_col]] < fdr_cutoff
    if (logfc_cutoff == 0) {
        up <- significant & df[[logfc_col]] > 0
        down <- significant & df[[logfc_col]] < 0
    } else {
        up <- significant & df[[logfc_col]] >= logfc_cutoff
        down <- significant & df[[logfc_col]] <= -logfc_cutoff
    }
    df$Significance[up] <- "Up"
    df$Significance[down] <- "Down"

    validate_order <- function(requested, observed, argument_name) {
        if (is.null(requested)) {
            return(observed)
        }
        if (!is.character(requested) || anyNA(requested) ||
            any(!nzchar(requested)) || anyDuplicated(requested)) {
            stop(argument_name, " must contain unique, non-empty labels.")
        }
        c(requested, setdiff(observed, requested))
    }
    observed_sorting <- unique(as.character(df$Sorting))
    observed_comparison <- unique(as.character(df$Comparison))
    if (anyNA(observed_sorting) || any(!nzchar(observed_sorting)) ||
        anyNA(observed_comparison) || any(!nzchar(observed_comparison))) {
        stop("Sorting and Comparison identifiers must be non-missing and non-empty.")
    }
    sorting_order <- validate_order(
        sorting_order, observed_sorting, "sorting_order"
    )
    comparison_order <- validate_order(
        comparison_order, observed_comparison, "comparison_order"
    )
    df$Sorting <- factor(df$Sorting, levels = sorting_order)
    df$Comparison <- factor(df$Comparison, levels = comparison_order)

    if (is.null(comparison_labels)) {
        label_map <- stats::setNames(comparison_order, comparison_order)
    } else if (is.null(names(comparison_labels))) {
        if (!is.character(comparison_labels) || anyNA(comparison_labels) ||
            length(comparison_labels) != length(comparison_order) ||
            any(!nzchar(comparison_labels))) {
            stop(
                "Unnamed comparison_labels must provide one non-empty label ",
                "per comparison_order value."
            )
        }
        label_map <- stats::setNames(comparison_labels, comparison_order)
    } else {
        if (!is.character(comparison_labels) || anyNA(comparison_labels) ||
            any(!nzchar(comparison_labels)) ||
            anyNA(names(comparison_labels)) ||
            any(!nzchar(names(comparison_labels))) ||
            anyDuplicated(names(comparison_labels))) {
            stop("Named comparison_labels must contain valid unique names.")
        }
        label_map <- comparison_labels
        missing_labels <- setdiff(comparison_order, names(label_map))
        label_map[missing_labels] <- missing_labels
        label_map <- label_map[comparison_order]
    }
    comparison_display_order <- unname(label_map[comparison_order])
    if (anyDuplicated(comparison_display_order)) {
        stop("comparison_labels must produce unique display labels.")
    }
    df$ComparisonLabel <- factor(
        unname(label_map[as.character(df$Comparison)]),
        levels = comparison_display_order
    )
    display_levels <- levels(df$ComparisonLabel)

    if (is.null(colors)) {
        if (length(display_levels) <= length(SIGNED_MANHATTAN_MACARON_COLORS)) {
            colors <- SIGNED_MANHATTAN_MACARON_COLORS[seq_along(display_levels)]
        } else {
            colors <- grDevices::hcl.colors(
                length(display_levels), palette = "Pastel 1"
            )
        }
        names(colors) <- display_levels
    } else {
        if (!is.character(colors) || length(colors) == 0L || anyNA(colors)) {
            stop("colors must be a non-empty character vector without NA.")
        }
        if (is.null(names(colors))) {
            if (length(colors) < length(display_levels)) {
                stop("Not enough colors supplied. Need ", length(display_levels), ".")
            }
            colors <- colors[seq_along(display_levels)]
            names(colors) <- display_levels
        } else {
            if (anyNA(names(colors)) || any(!nzchar(names(colors))) ||
                anyDuplicated(names(colors))) {
                stop("Named colors must have unique, non-empty names.")
            }
            if (all(comparison_order %in% names(colors)) &&
                !all(display_levels %in% names(colors))) {
                colors <- colors[comparison_order]
                names(colors) <- comparison_display_order
            }
            missing_colors <- setdiff(display_levels, names(colors))
            if (length(missing_colors) > 0L) {
                stop("Missing colors for: ", paste(missing_colors, collapse = ", "))
            }
            colors <- colors[display_levels]
        }
    }

    if (gene_order == "alphabetical") {
        gene_labels <- tapply(df$.GeneLabel, df$.GeneKey, function(x) x[1])
        gene_key_order <- names(sort(tolower(gene_labels)))
        df$GeneIndex <- match(df$.GeneKey, gene_key_order)
    } else if (gene_order == "mean_logFC") {
        gene_means <- tapply(df[[logfc_col]], df$.GeneKey, mean, na.rm = TRUE)
        gene_key_order <- names(sort(gene_means, decreasing = FALSE))
        df$GeneIndex <- match(df$.GeneKey, gene_key_order)
    } else if (gene_order == "input") {
        gene_key_order <- unique(df$.GeneKey)
        df$GeneIndex <- match(df$.GeneKey, gene_key_order)
    } else {
        df$GeneIndex <- NA_integer_
        panel_indices <- split(
            seq_len(nrow(df)),
            interaction(df$Sorting, df$ComparisonLabel, drop = TRUE)
        )
        for (indices in panel_indices) {
            ordered_indices <- indices[
                order(df[[logfc_col]][indices], decreasing = FALSE)
            ]
            df$GeneIndex[ordered_indices] <- seq_along(ordered_indices)
        }
    }

    df$SignedLogFDR_plot <- df$SignedLogFDR
    df$Capped <- FALSE
    if (!is.null(cap_y)) {
        df$Capped <- abs(df$SignedLogFDR) > cap_y
        df$SignedLogFDR_plot <- pmax(pmin(df$SignedLogFDR, cap_y), -cap_y)
    }

    df$.row_id <- seq_len(nrow(df))
    panel_indices <- split(
        seq_len(nrow(df)),
        interaction(df$Sorting, df$ComparisonLabel, drop = TRUE)
    )
    label_rows <- integer(0)
    for (indices in panel_indices) {
        candidates <- if (label_only_significant) {
            indices[df$Significance[indices] != "NS"]
        } else {
            indices
        }
        positive <- candidates[df[[logfc_col]][candidates] > 0]
        negative <- candidates[df[[logfc_col]][candidates] < 0]
        if (label_top_up > 0 && length(positive) > 0L) {
            positive <- positive[
                order(df[[padj_col]][positive], -df[[logfc_col]][positive])
            ]
            label_rows <- c(label_rows, head(positive, label_top_up))
        }
        if (label_top_down > 0 && length(negative) > 0L) {
            negative <- negative[
                order(df[[padj_col]][negative], df[[logfc_col]][negative])
            ]
            label_rows <- c(label_rows, head(negative, label_top_down))
        }
    }
    if (!is.null(label_genes)) {
        label_rows <- c(
            label_rows,
            df$.row_id[
                df$.GeneLabel %in% label_genes | df$.GeneKey %in% label_genes
            ]
        )
    }
    label_rows <- unique(label_rows[!is.na(label_rows)])
    label_df <- df[df$.row_id %in% label_rows, , drop = FALSE]

    plot <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = GeneIndex, y = SignedLogFDR_plot)
    )
    if (show_zero_line) {
        plot <- plot + ggplot2::geom_hline(
            yintercept = 0, linewidth = 0.4, colour = "#777777"
        )
    }
    if (show_fdr_line) {
        fdr_line <- -log10(fdr_cutoff)
        plot <- plot + ggplot2::geom_hline(
            yintercept = c(-fdr_line, fdr_line),
            linewidth = 0.3,
            linetype = "dashed",
            colour = "#A0A0A0"
        )
    }
    non_significant <- df[df$Significance == "NS", , drop = FALSE]
    significant_data <- df[df$Significance != "NS", , drop = FALSE]
    if (rasterized) {
        plot <- plot +
            ggrastr::geom_point_rast(
                data = non_significant,
                ggplot2::aes(colour = ComparisonLabel),
                size = point_size,
                alpha = ns_alpha,
                raster.dpi = raster_dpi
            ) +
            ggrastr::geom_point_rast(
                data = significant_data,
                ggplot2::aes(colour = ComparisonLabel),
                size = point_size * 1.2,
                alpha = sig_alpha,
                raster.dpi = raster_dpi
            )
    } else {
        plot <- plot +
            ggplot2::geom_point(
                data = non_significant,
                ggplot2::aes(colour = ComparisonLabel),
                size = point_size,
                alpha = ns_alpha,
                stroke = 0
            ) +
            ggplot2::geom_point(
                data = significant_data,
                ggplot2::aes(colour = ComparisonLabel),
                size = point_size * 1.2,
                alpha = sig_alpha,
                stroke = 0
            )
    }
    if (!is.null(cap_y) && show_capped_triangles) {
        capped_up <- df[df$Capped & df$SignedLogFDR > 0, , drop = FALSE]
        capped_down <- df[df$Capped & df$SignedLogFDR < 0, , drop = FALSE]
        if (nrow(capped_up) > 0L) {
            plot <- plot + ggplot2::geom_point(
                data = capped_up,
                mapping = ggplot2::aes(
                    x = GeneIndex,
                    y = SignedLogFDR_plot,
                    colour = ComparisonLabel
                ),
                inherit.aes = FALSE,
                shape = 24,
                size = point_size * 2.2,
                stroke = 0.35,
                fill = "white"
            )
        }
        if (nrow(capped_down) > 0L) {
            plot <- plot + ggplot2::geom_point(
                data = capped_down,
                mapping = ggplot2::aes(
                    x = GeneIndex,
                    y = SignedLogFDR_plot,
                    colour = ComparisonLabel
                ),
                inherit.aes = FALSE,
                shape = 25,
                size = point_size * 2.2,
                stroke = 0.35,
                fill = "white"
            )
        }
    }

    n_sorting <- length(unique(as.character(df$Sorting)))
    if (facet_by_sorting && n_sorting > 1L) {
        plot <- plot + ggplot2::facet_grid(
            Sorting ~ ComparisonLabel, scales = facet_scales
        )
    } else {
        if (is.null(facet_ncol)) {
            facet_ncol <- length(display_levels)
        }
        plot <- plot + ggplot2::facet_wrap(
            ~ComparisonLabel,
            ncol = as.integer(facet_ncol),
            scales = facet_scales
        )
    }
    plot <- plot + ggplot2::scale_colour_manual(
        values = colors, drop = FALSE
    )
    if (nrow(label_df) > 0L) {
        plot <- plot + ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(
                x = GeneIndex,
                y = SignedLogFDR_plot,
                label = .GeneLabel
            ),
            inherit.aes = FALSE,
            colour = "black",
            size = label_size,
            fontface = if (italic_gene_labels) "italic" else "plain",
            box.padding = 0.35,
            point.padding = 0.20,
            force = 1.5,
            force_pull = 0.3,
            min.segment.length = 0,
            segment.linewidth = 0.25,
            segment.colour = "#777777",
            max.overlaps = Inf,
            seed = 123,
            show.legend = FALSE
        )
    }
    plot <- plot +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.12, 0.12))
        ) +
        ggplot2::coord_cartesian(clip = "off")

    y_axis_title <- if (identical(padj_col, "FDR") &&
        identical(logfc_col, "logFC")) {
        expression(-log[10](FDR) %*% sign(logFC))
    } else {
        paste0("-log10(", padj_col, ") x sign(", logfc_col, ")")
    }
    plot <- plot +
        ggplot2::labs(
            x = x_title,
            y = y_axis_title,
            colour = NULL,
            title = title
        ) +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
            axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
            axis.ticks.length = grid::unit(1.5, "mm"),
            axis.text.y = ggplot2::element_text(colour = "black", size = 7.5),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(colour = "black", size = 9),
            strip.background = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(
                colour = "black", face = "bold", size = 8.5
            ),
            strip.text.y = ggplot2::element_text(
                colour = "black", face = "bold", size = 8.5
            ),
            plot.title = ggplot2::element_text(
                colour = "black", face = "bold", size = 10, hjust = 0
            ),
            legend.position = "none",
            plot.margin = ggplot2::margin(8, 14, 8, 8)
        )

    if (is.null(height)) {
        height <- if (facet_by_sorting && n_sorting > 1L) {
            max(4, n_sorting * 2.4)
        } else {
            4
        }
    }
    if (!is.null(save)) {
        if (grepl("\\.pdf$", save, ignore.case = TRUE)) {
            ggplot2::ggsave(
                filename = save,
                plot = plot,
                width = width,
                height = height,
                units = "in",
                device = grDevices::cairo_pdf
            )
        } else {
            ggplot2::ggsave(
                filename = save,
                plot = plot,
                width = width,
                height = height,
                units = "in",
                dpi = 300
            )
        }
    }

    summary_indices <- split(
        seq_len(nrow(df)),
        interaction(df$Sorting, df$ComparisonLabel, drop = TRUE)
    )
    summary_df <- do.call(rbind, lapply(summary_indices, function(indices) {
        data.frame(
            Sorting = as.character(df$Sorting[indices[1]]),
            Comparison = as.character(df$Comparison[indices[1]]),
            Genes_tested = length(indices),
            Significant = sum(df$Significance[indices] != "NS"),
            Up_comparison = sum(df$Significance[indices] == "Up"),
            Up_reference = sum(df$Significance[indices] == "Down"),
            Capped = sum(df$Capped[indices]),
            stringsAsFactors = FALSE
        )
    }))
    rownames(summary_df) <- NULL

    if (nrow(label_df) > 0L) {
        label_summary_indices <- split(
            seq_len(nrow(label_df)),
            interaction(
                label_df$Sorting,
                label_df$ComparisonLabel,
                drop = TRUE
            )
        )
        label_summary <- do.call(
            rbind,
            lapply(label_summary_indices, function(indices) {
                data.frame(
                    Sorting = as.character(label_df$Sorting[indices[1]]),
                    ComparisonLabel = as.character(
                        label_df$ComparisonLabel[indices[1]]
                    ),
                    Labels = length(indices),
                    stringsAsFactors = FALSE
                )
            })
        )
        rownames(label_summary) <- NULL
    } else {
        label_summary <- data.frame()
    }

    if (verbose) {
        message("")
        message("============================================")
        message(" Signed Manhattan plot")
        message("============================================")
        message("Analyses:            ", analysis_count)
        message("Sorting populations: ", n_sorting)
        message("Comparisons:         ", length(unique(df$Comparison)))
        message("FDR cutoff:          ", fdr_cutoff)
        message("logFC cutoff:        ", logfc_cutoff)
        message("Gene order:          ", gene_order)
        if (!is.null(cap_y)) {
            message("Y-axis cap:          +/-", cap_y)
        }
        message("Genes labeled:       ", nrow(label_df))
        message("")
        print(summary_df)
        if (nrow(label_summary) > 0L) {
            message("")
            message("Labels per panel:")
            print(label_summary)
        }
        message("============================================")
    }

    list(
        plot = plot,
        data = df,
        labels = label_df,
        summary = summary_df,
        label_summary = label_summary,
        colors = colors,
        comparison_labels = label_map,
        fdr_cutoff = fdr_cutoff,
        logfc_cutoff = logfc_cutoff,
        cap_y = cap_y
    )
}


# =============================================================================
# Sample-level TF activity and differential activity for bulk RNA-seq
# =============================================================================
# Requires: edgeR, decoupleR, limma; statmod for robust limma moderation.
# No package installation, network downloads or file writes occur on source().
# References:
# https://saezlab.github.io/decoupleR/reference/run_ulm.html
# https://saezlab.github.io/decoupleR/articles/tf_bk.html
#
# Intended input: prepare_bulk_expression() output (logCPM, metadata, genes, dge).
# Scores use one shared gene universe and network across the selected samples.
# Gene filtering subsets the existing normalized logCPM; it does not renormalize.
# ULM enrichment p-values are not differential-activity p-values.


#' Infer sample-level TF activity from prepared bulk expression
#'
#' @param bulk List with logCPM, metadata and genes; dge is needed for filtering.
#' @param network Signed TF-target data frame, e.g. mouse/human CollecTRI.
#' @param group_cols Metadata columns defining groups for filterByExpr.
#'   For factorial experiments use all biological grouping columns.
#' @param gene_col Gene-symbol column of bulk$genes; NULL uses gene row names.
#' @param source_col,target_col,weight_col Network column names.
#' @param filter_genes Use edgeR::filterByExpr on the supplied dge.
#' @param min_count,min_total_count Passed to edgeR::filterByExpr.
#' @param duplicate_genes Keep highest mean-logCPM row per symbol, or stop.
#' @param center Center each gene across ALL selected samples before scoring.
#'   FALSE follows the direct normalized-expression workflow. TRUE produces
#'   relative activities; use one fixed centering reference for all comparisons.
#' @param min_targets Minimum distinct measured targets per regulator.
#' @param verbose Print processing counts.
#' @return List containing activity_matrix (TF x sample), activity_long,
#'   metadata, network, coverage, logCPM, gene_map, QC and provenance.
#' @export
infer_bulk_tf_activity <- function(
    bulk,
    network,
    group_cols = NULL,
    gene_col = "gene_name",
    source_col = "source",
    target_col = "target",
    weight_col = "mor",
    filter_genes = TRUE,
    min_count = 10,
    min_total_count = 15,
    duplicate_genes = c("highest_mean", "error"),
    center = FALSE,
    min_targets = 10L,
    verbose = TRUE
) {
    # ------------------------------------------------------------
    # Dependencies and input checks
    # ------------------------------------------------------------
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        stop("decoupleR is required.", call. = FALSE)
    }
    duplicate_genes <- match.arg(duplicate_genes)
    for (value in list(filter_genes, center, verbose)) {
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop("filter_genes, center and verbose must be logical scalars.")
        }
    }
    if (!is.numeric(min_targets) || length(min_targets) != 1L ||
        !is.finite(min_targets) || min_targets < 2 || min_targets != floor(min_targets)) {
        stop("min_targets must be an integer >= 2.")
    }
    for (value in list(min_count, min_total_count)) {
        if (!is.numeric(value) || length(value) != 1L ||
            !is.finite(value) || value < 0) stop("Count thresholds must be >= 0.")
    }
    if (!is.list(bulk) || !all(c("logCPM", "metadata", "genes") %in% names(bulk))) {
        stop("bulk must contain logCPM, metadata and genes.")
    }
    mat <- as.matrix(bulk$logCPM)
    meta <- as.data.frame(bulk$metadata)
    anno <- as.data.frame(bulk$genes)
    if (!is.numeric(mat) || nrow(mat) < 3L || ncol(mat) < 3L ||
        any(!is.finite(mat))) stop("logCPM must be finite and numeric, with >=3 genes/samples.")
    if (is.null(rownames(mat)) || is.null(colnames(mat)) ||
        anyNA(rownames(mat)) || anyNA(colnames(mat)) ||
        any(rownames(mat) == "") || any(colnames(mat) == "") ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat))) {
        stop("logCPM must have unique, nonmissing gene and sample names.")
    }
    if (anyDuplicated(rownames(meta)) || anyDuplicated(rownames(anno)) ||
        !setequal(colnames(mat), rownames(meta)) ||
        !setequal(rownames(mat), rownames(anno))) {
        stop("Expression, metadata and annotation identifiers do not match uniquely.")
    }
    meta <- meta[colnames(mat), , drop = FALSE]
    anno <- anno[rownames(mat), , drop = FALSE]
    n_input <- nrow(mat)

    # ------------------------------------------------------------
    # One shared expression filter across all samples
    # ------------------------------------------------------------
    keep <- rep(TRUE, nrow(mat))
    names(keep) <- rownames(mat)
    if (filter_genes) {
        if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is required.")
        if (is.null(bulk$dge) || !inherits(bulk$dge, "DGEList")) {
            stop("Filtering requires bulk$dge as an edgeR DGEList.")
        }
        if (is.null(group_cols) || !is.character(group_cols) ||
            !length(group_cols) || anyNA(group_cols) ||
            anyDuplicated(group_cols) || !all(group_cols %in% colnames(meta))) {
            stop("Provide valid group_cols for expression filtering.")
        }
        group_df <- meta[, group_cols, drop = FALSE]
        if (anyNA(group_df)) stop("Filtering group columns contain missing values.")
        if (any(vapply(group_df, function(x) any(trimws(as.character(x)) == ""), logical(1)))) {
            stop("Filtering group columns contain blank values.")
        }
        # Integer codes avoid collisions when group labels contain separators.
        filter_group <- do.call(interaction, c(lapply(group_df, function(x) as.integer(factor(x))),
                                               list(drop = TRUE, lex.order = TRUE)))
        dge <- bulk$dge
        if (!identical(rownames(dge), rownames(mat)) ||
            !identical(colnames(dge), colnames(mat))) {
            stop("bulk$dge and logCPM must have identical gene/sample order.")
        }
        keep <- edgeR::filterByExpr(dge, group = filter_group,
                                    min.count = min_count,
                                    min.total.count = min_total_count)
        names(keep) <- rownames(mat)
        if (!any(keep)) stop("No genes passed expression filtering.")
    }
    mat <- mat[keep, , drop = FALSE]
    anno <- anno[rownames(mat), , drop = FALSE]
    n_filtered <- nrow(mat)

    # ------------------------------------------------------------
    # Map to unique gene symbols without selecting by significance
    # ------------------------------------------------------------
    if (is.null(gene_col)) {
        symbols <- rownames(mat)
    } else {
        if (!is.character(gene_col) || length(gene_col) != 1L ||
            is.na(gene_col) || !gene_col %in% colnames(anno)) {
            stop("gene_col is absent from bulk$genes.")
        }
        symbols <- trimws(as.character(anno[[gene_col]]))
    }
    valid <- !is.na(symbols) & nzchar(symbols)
    gene_map <- data.frame(GeneID = rownames(mat), Symbol = symbols,
                           MeanLogCPM = rowMeans(mat), stringsAsFactors = FALSE)
    gene_map <- gene_map[valid, , drop = FALSE]
    n_annotated <- nrow(gene_map)
    if (!n_annotated) stop("No valid gene symbols remain.")
    if (duplicate_genes == "error" && anyDuplicated(gene_map$Symbol)) {
        stop("Duplicated gene symbols; use duplicate_genes='highest_mean' or resolve manually.")
    }
    gene_map <- gene_map[order(-gene_map$MeanLogCPM, gene_map$GeneID), , drop = FALSE]
    gene_map <- gene_map[!duplicated(gene_map$Symbol), , drop = FALSE]
    mat <- mat[gene_map$GeneID, , drop = FALSE]
    rownames(mat) <- gene_map$Symbol
    if (nrow(mat) < 3L) stop("Too few uniquely annotated genes for ULM.")

    # ------------------------------------------------------------
    # Validate network and measure regulon coverage
    # ------------------------------------------------------------
    net_cols <- c(source_col, target_col, weight_col)
    if (!is.data.frame(network) || length(net_cols) != 3L ||
        anyNA(net_cols) || anyDuplicated(net_cols) ||
        !all(net_cols %in% colnames(network))) stop("Invalid network columns.")
    raw_weight <- network[[weight_col]]
    if (!is.numeric(raw_weight)) stop("Network weights must be numeric.")
    net <- data.frame(source = trimws(as.character(network[[source_col]])),
                      target = trimws(as.character(network[[target_col]])),
                      mor = raw_weight, stringsAsFactors = FALSE)
    if (anyNA(net) || any(!is.finite(net$mor)) ||
        any(net$source == "") || any(net$target == "")) {
        stop("Network contains missing/blank identifiers or nonfinite weights.")
    }
    net <- unique(net[net$mor != 0, , drop = FALSE])
    if (!nrow(net)) stop("No nonzero network edges.")
    if (anyDuplicated(net[c("source", "target")])) {
        stop("Conflicting weights for a TF-target pair; resolve before scoring.")
    }
    tf_names <- unique(net$source)
    coverage <- data.frame(
        TF = tf_names,
        Targets_network = as.integer(table(factor(net$source, levels = tf_names))),
        Targets_measured = as.integer(table(factor(
            net$source[net$target %in% rownames(mat)], levels = tf_names))),
        stringsAsFactors = FALSE
    )
    coverage$Coverage_fraction <- coverage$Targets_measured / coverage$Targets_network
    coverage$Eligible <- coverage$Targets_measured >= min_targets
    if (!any(coverage$Eligible)) {
        stop("No TF has sufficient targets. Check species, gene symbols and filtering.")
    }
    net_used <- net[net$source %in% coverage$TF[coverage$Eligible] &
                    net$target %in% rownames(mat), , drop = FALSE]

    # ------------------------------------------------------------
    # Infer all samples using the SAME matrix background and network
    # ------------------------------------------------------------
    scoring_mat <- if (center) sweep(mat, 1L, rowMeans(mat), "-") else mat
    if (any(apply(scoring_mat, 2L, stats::sd) == 0)) {
        stop("A sample has no expression variation after preprocessing.")
    }
    acts <- as.data.frame(decoupleR::run_ulm(
        mat = scoring_mat, network = net_used,
        .source = "source", .target = "target", .mor = "mor",
        center = FALSE, minsize = as.integer(min_targets)
    ))
    if (!all(c("source", "condition", "score", "p_value") %in% names(acts)) ||
        any(!is.finite(acts$score))) stop("ULM returned invalid scores.")
    acts$source <- as.character(acts$source)
    acts$condition <- as.character(acts$condition)
    eligible_tfs <- coverage$TF[coverage$Eligible]
    if (!setequal(acts$source, eligible_tfs) ||
        !setequal(acts$condition, colnames(mat)) ||
        anyDuplicated(acts[c("source", "condition")]) ||
        nrow(acts) != length(eligible_tfs) * ncol(mat)) {
        stop("ULM did not return one score for every eligible TF/sample pair.")
    }
    activity_matrix <- matrix(NA_real_, length(eligible_tfs), ncol(mat),
                              dimnames = list(eligible_tfs, colnames(mat)))
    index <- cbind(match(acts$source, eligible_tfs), match(acts$condition, colnames(mat)))
    activity_matrix[index] <- acts$score
    activity_long <- data.frame(TF = acts$source, Sample = acts$condition,
                                Activity = acts$score, ULM_PValue = acts$p_value)
    qc <- data.frame(Input_genes = n_input, Filtered_genes = n_filtered,
                      Annotated_rows = n_annotated, Unique_symbols = nrow(mat),
                      Duplicate_symbols_removed = n_annotated - nrow(mat),
                      Genes_matching_network = sum(rownames(mat) %in% net$target),
                      TFs_scored = nrow(activity_matrix), Samples = ncol(mat))
    if (verbose) {
        message("\n====================================")
        message("Sample-level TF activity complete")
        message("====================================")
        print(qc, row.names = FALSE)
    }
    result <- list(
        method = "ULM", activity_matrix = activity_matrix,
        activity_long = activity_long, metadata = meta,
        coverage = coverage, network = net, network_used = net_used,
        logCPM = mat, gene_map = gene_map, keep_genes = keep, qc = qc,
        settings = list(group_cols = group_cols, gene_col = gene_col,
                         filter_genes = filter_genes, min_count = min_count,
                         min_total_count = min_total_count,
                         duplicate_genes = duplicate_genes,
                         center = center, min_targets = min_targets,
                         normalization = bulk$normalization),
        session = utils::sessionInfo()
    )
    class(result) <- c("bulk_tf_activity", "list")
    result
}


#' Test differential activity between two groups using limma
#'
#' Positive DeltaActivity means higher activity in group2 than group1.
#' Scores are continuous activities, so coefficients are NOT log fold changes.
#' @param tf Output from infer_bulk_tf_activity.
#' @param group_col,group1,group2 Metadata group column, reference and comparison.
#' @param subset Optional expression evaluated in metadata and caller scope.
#' @param covariates Character vector of batch/other adjustment columns.
#' @param block_col Verified subject-ID column for duplicateCorrelation.
#' @param robust,trend Passed to limma::eBayes; trend=FALSE is the default.
#' @param fdr_cutoff BH FDR threshold for Direction annotation.
#' @param delta_cutoff Additional descriptive effect threshold; this does not
#'   change the null hypothesis of the limma test (zero activity difference).
#' @param verbose Logical scalar. Whether to report sample counts, tested TFs,
#'   excluded constant TFs, and the number passing the FDR cutoff.
#' @return results, fitted model, design, contrast, selected metadata, dropped
#'   constant TFs, groups, sample counts, covariates and block correlation.
#' @export
limma_tf_pairwise <- function(
    tf,
    group_col,
    group1,
    group2,
    subset = NULL,
    covariates = NULL,
    block_col = NULL,
    robust = TRUE,
    trend = FALSE,
    fdr_cutoff = 0.05,
    delta_cutoff = 0,
    verbose = TRUE
) {
    if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required.")
    if (!is.list(tf) || !all(c("activity_matrix", "metadata") %in% names(tf))) {
        stop("tf must contain activity_matrix and metadata.")
    }
    for (value in list(robust, trend, verbose)) {
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop("robust, trend and verbose must be logical scalars.")
        }
    }
    if (robust && !requireNamespace("statmod", quietly = TRUE)) {
        stop("statmod is required when robust=TRUE.")
    }
    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1) stop("Invalid fdr_cutoff.")
    if (!is.numeric(delta_cutoff) || length(delta_cutoff) != 1L ||
        !is.finite(delta_cutoff) || delta_cutoff < 0) stop("Invalid delta_cutoff.")
    for (value in list(group_col, group1, group2)) {
        if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
            stop("Group column and labels must be nonmissing scalars.")
        }
    }
    group_col <- as.character(group_col)
    group1 <- as.character(group1)
    group2 <- as.character(group2)
    if (group1 == group2) stop("group1 and group2 must differ.")
    meta <- as.data.frame(tf$metadata)
    mat <- as.matrix(tf$activity_matrix)
    if (!is.numeric(mat) || any(!is.finite(mat)) ||
        is.null(rownames(mat)) || is.null(colnames(mat)) ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) ||
        anyDuplicated(rownames(meta)) || !setequal(colnames(mat), rownames(meta))) {
        stop("Activity matrix and metadata must be finite, unique and aligned.")
    }
    meta <- meta[colnames(mat), , drop = FALSE]
    if (!group_col %in% names(meta)) stop("group_col is absent from metadata.")
    subset_expr <- substitute(subset)
    selected <- rep(TRUE, nrow(meta))
    if (!identical(subset_expr, quote(NULL))) {
        selected <- eval(subset_expr, envir = meta, enclos = parent.frame())
        if (!is.logical(selected) || length(selected) != nrow(meta) || anyNA(selected)) {
            stop("subset must return one nonmissing logical value per sample.")
        }
    }
    group_values <- as.character(meta[[group_col]])
    selected <- selected & !is.na(group_values) & group_values %in% c(group1, group2)
    meta <- meta[selected, , drop = FALSE]
    mat <- mat[, rownames(meta), drop = FALSE]
    group_values <- as.character(meta[[group_col]])
    n1 <- sum(group_values == group1)
    n2 <- sum(group_values == group2)
    if (min(n1, n2) < 2L) stop("At least two biological samples per group are required.")

    # ------------------------------------------------------------
    # Build model; internal covariate names also support spaces in metadata names
    # ------------------------------------------------------------
    if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates) ||
        anyDuplicated(covariates) || !all(covariates %in% names(meta)))) {
        stop("Invalid covariate columns.")
    }
    if (group_col %in% covariates) stop("Do not include group_col as a covariate.")
    model_data <- data.frame(.group = factor(
        ifelse(group_values == group1, "reference", "comparison"),
        levels = c("reference", "comparison")), row.names = rownames(meta))
    internal_covariates <- character()
    for (i in seq_along(covariates)) {
        value <- meta[[covariates[i]]]
        if (anyNA(value)) stop("Missing covariate values: ", covariates[i])
        if (is.character(value) || is.factor(value) || is.logical(value)) {
            if (any(trimws(as.character(value)) == "")) stop("Blank covariate values.")
            value <- droplevels(factor(value))
        } else if (!is.numeric(value) || any(!is.finite(value))) {
            stop("Covariates must be finite numeric or categorical columns.")
        }
        if (length(unique(value)) < 2L) {
            stop("Covariate is constant in selected samples: ", covariates[i],
                 ". Omit it from this comparison.")
        }
        nm <- paste0(".covariate", i)
        model_data[[nm]] <- value
        internal_covariates <- c(internal_covariates, nm)
    }
    design <- stats::model.matrix(
        stats::reformulate(c(".group", internal_covariates)), data = model_data
    )
    if (qr(design)$rank < ncol(design)) {
        stop("Design is rank deficient: check group/batch confounding.")
    }
    if (nrow(design) <= ncol(design)) stop("No residual degrees of freedom.")
    comparison_coef <- match(".groupcomparison", colnames(design))
    if (is.na(comparison_coef)) stop("Internal error: group coefficient not found.")
    contrast <- matrix(0, ncol(design), 1L,
                        dimnames = list(colnames(design), paste0(group2, "_vs_", group1)))
    contrast[comparison_coef, 1L] <- 1
    variable <- apply(mat, 1L, stats::sd) > 0
    dropped_tfs <- rownames(mat)[!variable]
    mat <- mat[variable, , drop = FALSE]
    if (nrow(mat) < 2L) stop("Fewer than two variable TFs remain for limma moderation.")

    # ------------------------------------------------------------
    # Optional repeated-subject correlation
    # ------------------------------------------------------------
    block <- NULL
    correlation <- NULL
    if (!is.null(block_col)) {
        if (!is.character(block_col) || length(block_col) != 1L ||
            is.na(block_col) || !block_col %in% names(meta)) stop("Invalid block_col.")
        if (block_col %in% covariates) stop("Do not model the same column as covariate and block.")
        value <- as.character(meta[[block_col]])
        if (anyNA(value) || any(trimws(value) == "")) stop("Missing/blank subject IDs.")
        block <- factor(value)
        if (sum(table(block) >= 2L) < 2L) stop("At least two repeated subjects are required.")
        if (nrow(design) < ncol(design) + 2L) {
            stop("duplicateCorrelation requires at least two residual degrees of freedom.")
        }
        if (!requireNamespace("statmod", quietly = TRUE)) stop("statmod is required for blocking.")
        correlation <- limma::duplicateCorrelation(mat, design, block = block)$consensus.correlation
        if (!is.finite(correlation)) stop("Could not estimate a finite block correlation.")
        fit <- limma::lmFit(mat, design, block = block, correlation = correlation)
    } else {
        fit <- limma::lmFit(mat, design)
    }
    fit <- limma::contrasts.fit(fit, contrasts = contrast)
    fit <- limma::eBayes(fit, robust = robust, trend = trend)
    x <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "P")
    results <- data.frame(TF = rownames(x), DeltaActivity = x$logFC,
                          MeanActivity = x$AveExpr, t = x$t,
                          PValue = x$P.Value, FDR = x$adj.P.Val, B = x$B,
                          stringsAsFactors = FALSE)
    rownames(results) <- results$TF
    results$Direction <- "NS"
    sig <- !is.na(results$FDR) & results$FDR < fdr_cutoff
    results$Direction[sig & results$DeltaActivity > delta_cutoff] <- paste0("Up_", group2)
    results$Direction[sig & results$DeltaActivity < -delta_cutoff] <- paste0("Up_", group1)
    if (!is.null(tf$coverage)) {
        coverage_cols <- c("Targets_network", "Targets_measured", "Coverage_fraction")
        if (all(c("TF", coverage_cols) %in% names(tf$coverage))) {
            results <- cbind(results, tf$coverage[
                match(results$TF, tf$coverage$TF), coverage_cols, drop = FALSE])
        }
    }
    if (verbose) {
        message("\nComparison: ", group2, " vs ", group1)
        message(group1, ": n = ", n1, "; ", group2, ": n = ", n2)
        message("TFs tested: ", nrow(results))
        message("Constant TFs excluded: ", length(dropped_tfs))
        message("FDR < ", fdr_cutoff, ": ", sum(sig))
    }
    list(method = "limma_TF_activity", comparison = paste0(group2, "_vs_", group1),
         results = results, fit = fit, design = design, contrast = contrast,
         metadata = meta, activity_matrix = mat, dropped_tfs = dropped_tfs,
         group1 = group1, group2 = group2, n_group1 = n1, n_group2 = n2,
         covariates = covariates,
         covariate_map = if (length(covariates) > 0L) {
             stats::setNames(covariates, internal_covariates)
         } else {
             NULL
         },
         block_col = block_col, block_correlation = correlation,
         thresholds = list(fdr_cutoff = fdr_cutoff, delta_cutoff = delta_cutoff),
         session = utils::sessionInfo())
}


#' Plot all significant TFs with subtype/condition and significance annotations
#'
#' TFs are selected from de_results using the requested FDR column, regardless
#' of direction. A subtype's row annotation is Yes if ANY of its comparisons
#' passes that threshold. Scores, not DE coefficients, form the heatmap.
#' Means are calculated on raw sample scores, then each TF is z-scored across
#' the displayed columns. Columns remain in the requested biological order.
#' Rows use either correlation distance 1-r or Euclidean distance.
#' Distance is calculated from unclipped row z-scores; linkage is configurable.
#'
#' @param tf Output of infer_bulk_tf_activity.
#' @param de_results Combined TF DE data frame, e.g. DE_TF_ATM_all.
#' @param subtype_col,condition_col Metadata grouping columns.
#' @param de_subtype_col Subtype column in de_results.
#' @param subtype_order,condition_order Explicit orders; no missing groups allowed.
#' @param subtype_colors,condition_colors Named annotation colour vectors.
#' @param fdr_col,fdr_cutoff Significance criterion, e.g. FDR or FDR_global.
#' @param column_mode "mean" for condition means or "samples" for all replicates.
#' @param feature_label Character label used in the heatmap title.
#' @param clustering_distance_rows Correlation (1-r) or Euclidean distance.
#' @param clustering_method Linkage passed to hclust; average by default.
#' @param correlation_method Pearson or Spearman correlation.
#' @param color Vector of heatmap colours; any length >=2 is supported.
#' @param z_limit Colour saturation threshold only; clustering uses unclipped data.
#' @param fontsize_row,fontsize_col Label sizes.
#' @param width,height Output dimensions in inches; NULL height adapts to TF count.
#' @param save Optional PDF/PNG/TIFF filename, written by pheatmap.
#' @return Plot, raw and z-score matrices, annotations, selected TFs, row tree,
#'   per-group sample counts and dropped TFs with constant group means.
#' @export
plot_bulk_tf_heatmap <- function(
    tf,
    de_results,
    subtype_col = "ATM_subtype",
    condition_col = "Condition",
    de_subtype_col = subtype_col,
    subtype_order = c("pATM", "cATM", "sATM"),
    condition_order = c("ND", "2W", "1M", "2M"),
    subtype_colors,
    condition_colors,
    fdr_col = "FDR",
    fdr_cutoff = 0.05,
    column_mode = c("mean", "samples"),
    feature_label = "TF activities",
    clustering_distance_rows = c("correlation", "euclidean"),
    clustering_method = "average",
    correlation_method = c("pearson", "spearman"),
    color = BULK_ACTIVITY_DIVERGING_COLORS,
    z_limit = 2.5,
    fontsize_row = 9,
    fontsize_col = 10,
    width = 9,
    height = NULL,
    save = NULL
) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) stop("pheatmap is required.")
    column_mode <- match.arg(column_mode)
    if (!is.character(feature_label) || length(feature_label) != 1L ||
        is.na(feature_label) || feature_label == "") stop("Invalid feature_label.")
    clustering_distance_rows <- match.arg(clustering_distance_rows)
    correlation_method <- match.arg(correlation_method)
    allowed_linkages <- c("ward.D", "ward.D2", "single", "complete", "average",
                          "mcquitty", "median", "centroid")
    if (!is.character(clustering_method) || length(clustering_method) != 1L ||
        is.na(clustering_method) || !clustering_method %in% allowed_linkages) {
        stop("Invalid clustering_method for hclust.")
    }
    if (!is.character(color) || length(color) < 2L || anyNA(color)) {
        stop("color must be a vector of at least two valid R colours.")
    }
    tryCatch(grDevices::col2rgb(color), error = function(e) {
        stop("Invalid heatmap color: ", conditionMessage(e), call. = FALSE)
    })
    if (!is.list(tf) || !all(c("activity_matrix", "metadata") %in% names(tf))) {
        stop("tf must contain activity_matrix and metadata.")
    }
    if (!is.data.frame(de_results) ||
        !all(c("TF", de_subtype_col, fdr_col) %in% names(de_results))) {
        stop("de_results must contain TF, subtype and requested FDR columns.")
    }
    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1) stop("Invalid FDR cutoff.")
    if (!is.numeric(z_limit) || length(z_limit) != 1L ||
        !is.finite(z_limit) || z_limit <= 0) stop("z_limit must be positive.")
    for (x in list(subtype_order, condition_order)) {
        if (!is.character(x) || !length(x) || anyNA(x) ||
            any(x == "") || anyDuplicated(x)) stop("Orders must contain unique nonblank labels.")
    }
    mat <- as.matrix(tf$activity_matrix)
    meta <- as.data.frame(tf$metadata)
    if (!is.numeric(mat) || any(!is.finite(mat)) ||
        is.null(rownames(mat)) || is.null(colnames(mat)) ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) ||
        anyDuplicated(rownames(meta)) || !setequal(colnames(mat), rownames(meta))) {
        stop("Activity matrix and metadata must have matching unique identifiers.")
    }
    meta <- meta[colnames(mat), , drop = FALSE]
    if (!all(c(subtype_col, condition_col) %in% names(meta))) stop("Missing metadata columns.")
    subtype <- as.character(meta[[subtype_col]])
    condition <- as.character(meta[[condition_col]])
    if (anyNA(subtype) || anyNA(condition)) stop("Missing subtype/condition labels.")
    selected_samples <- subtype %in% subtype_order & condition %in% condition_order
    meta <- meta[selected_samples, , drop = FALSE]
    mat <- mat[, rownames(meta), drop = FALSE]
    subtype <- as.character(meta[[subtype_col]])
    condition <- as.character(meta[[condition_col]])
    if (!ncol(mat)) stop("No samples match the requested groups.")
    for (i in seq_along(subtype_order)) {
        if (!subtype_order[i] %in% names(subtype_colors)) stop("Missing subtype colour.")
    }
    if (!all(condition_order %in% names(condition_colors))) stop("Missing condition colours.")

    # ------------------------------------------------------------
    # Union of all significant up/down TFs across comparisons
    # ------------------------------------------------------------
    fdr <- de_results[[fdr_col]]
    if (!is.numeric(fdr) || any(!is.na(fdr) & (!is.finite(fdr) | fdr < 0 | fdr > 1))) {
        stop("The FDR column must be numeric, with values between zero and one.")
    }
    sig <- de_results[!is.na(fdr) & fdr < fdr_cutoff &
                       de_results[[de_subtype_col]] %in% subtype_order, , drop = FALSE]
    selected_tfs <- unique(as.character(sig$TF))
    if (!length(selected_tfs)) stop("No significant TFs at this cutoff.")
    if (anyNA(selected_tfs) || any(selected_tfs == "")) stop("Invalid significant TF names.")
    if (!all(selected_tfs %in% rownames(mat))) stop("Significant TFs are absent from the score matrix.")
    mat <- mat[selected_tfs, , drop = FALSE]

    # ------------------------------------------------------------
    # Build columns in subtype-major, condition-minor order
    # ------------------------------------------------------------
    pieces <- list()
    column_annotations <- list()
    counts <- list()
    for (ct in subtype_order) {
        for (cond in condition_order) {
            ids <- rownames(meta)[subtype == ct & condition == cond]
            if (!length(ids)) stop("Missing group: ", ct, " / ", cond)
            tag <- paste(ct, cond, sep = "__")
            counts[[tag]] <- data.frame(Subtype = ct, Condition = cond, N = length(ids))
            if (column_mode == "mean") {
                block <- matrix(rowMeans(mat[, ids, drop = FALSE]), ncol = 1L,
                                dimnames = list(rownames(mat), tag))
            } else {
                block <- mat[, ids, drop = FALSE]
            }
            pieces[[tag]] <- block
            column_annotations[[tag]] <- data.frame(
                Subtype = rep(ct, ncol(block)), Condition = rep(cond, ncol(block)),
                row.names = colnames(block), stringsAsFactors = FALSE
            )
        }
    }
    raw <- do.call(cbind, pieces)
    annotation_col <- do.call(rbind, unname(column_annotations))
    annotation_col <- annotation_col[colnames(raw), , drop = FALSE]
    names(annotation_col)[1] <- subtype_col
    annotation_col[[subtype_col]] <- factor(annotation_col[[subtype_col]], levels = subtype_order)
    annotation_col$Condition <- factor(annotation_col$Condition, levels = condition_order)

    row_sd <- apply(raw, 1L, stats::sd)
    dropped_tfs <- rownames(raw)[!is.finite(row_sd) | row_sd == 0]
    if (length(dropped_tfs)) {
        warning("Constant displayed profiles omitted: ", paste(dropped_tfs, collapse = ", "))
        raw <- raw[!rownames(raw) %in% dropped_tfs, , drop = FALSE]
    }
    if (!nrow(raw)) stop("No variable significant TF profiles remain.")
    z <- sweep(sweep(raw, 1L, rowMeans(raw), "-"), 1L, apply(raw, 1L, stats::sd), "/")

    # ------------------------------------------------------------
    # Independent row bars: a TF may be significant in several subtypes
    # ------------------------------------------------------------
    annotation_row <- data.frame(row.names = rownames(raw))
    annotation_colors <- list()
    annotation_colors[[subtype_col]] <- subtype_colors[subtype_order]
    annotation_colors$Condition <- condition_colors[condition_order]
    for (ct in subtype_order) {
        tf_ct <- as.character(sig$TF[sig[[de_subtype_col]] == ct])
        nm <- paste0("Sig_", ct)
        annotation_row[[nm]] <- factor(ifelse(rownames(raw) %in% tf_ct, "Yes", "No"),
                                      levels = c("No", "Yes"))
        annotation_colors[[nm]] <- c(No = "#F0F0F0", Yes = unname(subtype_colors[ct]))
    }
    tree <- FALSE
    if (nrow(z) > 1L) {
        if (clustering_distance_rows == "correlation") {
            r <- stats::cor(t(z), method = correlation_method)
            if (any(!is.finite(r))) stop("Nonfinite row correlations.")
            distance <- 1 - r
            distance[distance < 0] <- 0
            diag(distance) <- 0
            row_distance <- stats::as.dist(distance)
        } else {
            row_distance <- stats::dist(z, method = "euclidean")
        }
        tree <- stats::hclust(row_distance, method = clustering_method)
    }
    displayed <- z
    displayed[displayed > z_limit] <- z_limit
    displayed[displayed < -z_limit] <- -z_limit
    group_sizes <- as.integer(table(factor(annotation_col[[subtype_col]], levels = subtype_order)))
    gaps_col <- head(cumsum(group_sizes), -1L)
    if (is.null(height)) height <- max(5, 2.3 + 0.16 * nrow(raw))
    labels_col <- as.character(annotation_col$Condition)
    if (column_mode == "samples") {
        labels_col <- paste(labels_col, ave(seq_len(ncol(raw)),
                           interaction(annotation_col[[subtype_col]], annotation_col$Condition),
                           FUN = seq_along), sep = "_")
    }
    args <- list(
        mat = displayed, scale = "none", cluster_rows = tree, cluster_cols = FALSE,
        annotation_row = annotation_row, annotation_col = annotation_col,
        annotation_colors = annotation_colors, drop_levels = FALSE,
        color = color,
        breaks = seq(-z_limit, z_limit, length.out = length(color) + 1L), border_color = NA,
        labels_col = labels_col, gaps_col = gaps_col,
        fontsize_row = fontsize_row, fontsize_col = fontsize_col,
        angle_col = "45", treeheight_row = 45,
        main = paste0("Significant ", feature_label, ": ", fdr_col, " < ", fdr_cutoff,
                       "\n", if (column_mode == "mean") "Condition means" else "Individual samples",
                       "; row z-score"),
        legend_breaks = c(-z_limit, 0, z_limit),
        legend_labels = c(paste0("<= -", z_limit), "0", paste0(">= ", z_limit)),
        silent = TRUE
    )
    if (!is.null(save)) {
        if (!is.character(save) || length(save) != 1L || is.na(save) || save == "") {
            stop("save must be a filename or NULL.")
        }
        if (!dir.exists(dirname(save))) dir.create(dirname(save), recursive = TRUE)
        args$filename <- save
        args$width <- width
        args$height <- height
    }
    plot <- do.call(pheatmap::pheatmap, args)
    grid::grid.newpage()
    grid::grid.draw(plot$gtable)
    list(plot = plot, raw_matrix = raw, z_matrix = z,
         annotation_row = annotation_row, annotation_col = annotation_col,
         significant_tfs = selected_tfs, plotted_tfs = rownames(raw),
         dropped_tfs = dropped_tfs, row_tree = tree,
         sample_counts = do.call(rbind, unname(counts)),
         settings = list(fdr_col = fdr_col, fdr_cutoff = fdr_cutoff,
                          column_mode = column_mode, z_limit = z_limit,
                          clustering_distance_rows = clustering_distance_rows,
                          clustering_method = clustering_method,
                          correlation_method = correlation_method,
                          color = color))
}


# =============================================================================
# Gene profiles, enrichment and condition-response analysis for bulk RNA-seq
# =============================================================================

#' Collect annotated edgeR comparisons with an explicit comparison manifest
#' @param de_list Named list of single-coefficient edgeR pairwise results.
#' @param stratum_col Metadata column identifying the analysed cell type;
#'   NULL uses a single stratum called All.
#' @param gene_col Annotation column containing symbols.
#' @param analysis_col Name of the analysis-identifier column retained from
#'   each result table when available.
#' @param condition_col,reference_col Names used for the comparison and
#'   reference-group fields in the collected output.
#' @param gene_id_col Name used for the original stable gene-row identifier.
#' @return table, comparison_info and QC. Duplicate symbols are resolved using
#'   highest logCPM, independently of significance, separately in each contrast.
#' @export
prepare_bulk_gene_de <- function(
    de_list, stratum_col = "ATM_subtype", gene_col = "gene_name",
    analysis_col = "Analysis", condition_col = "ComparisonGroup",
    reference_col = "Reference", gene_id_col = "GeneRowName"
) {
    # collect_edgeR_pairwise() retains contrast labels and gene identifiers,
    # so saved unfiltered long tables do not require the original fitted objects.
    # Sorting in that table is the first component of each analysis name.
    if (is.data.frame(de_list)) {
        long <- as.data.frame(de_list)
        required <- c(analysis_col, condition_col, reference_col, gene_id_col,
                      gene_col, "logFC", "logCPM", "F", "PValue", "FDR", stratum_col)
        if (!all(required %in% names(long))) {
            stop("Long DE table is missing: ",
                 paste(setdiff(required, names(long)), collapse = ", "),
                 ". collect_edgeR_pairwise() stores population labels in Sorting.")
        }
        if (!nrow(long)) stop("Long DE table is empty.")
        label_cols <- c(analysis_col, condition_col, reference_col, gene_id_col, stratum_col)
        for (column in label_cols) {
            labels <- as.character(long[[column]])
            if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
                stop("Missing or empty long-table identifiers in ", column)
            }
        }
        analyses <- unique(as.character(long[[analysis_col]]))
        de_list <- stats::setNames(vector("list", length(analyses)), analyses)
        for (nm in analyses) {
            x <- long[as.character(long[[analysis_col]]) == nm, , drop = FALSE]
            single <- function(column) {
                value <- unique(as.character(x[[column]]))
                if (length(value) != 1L) stop("Inconsistent ", column, " within ", nm)
                value
            }
            group1 <- single(reference_col)
            group2 <- single(condition_col)
            if (group1 == group2) stop("Reference equals comparison group in ", nm)
            gene_ids <- as.character(x[[gene_id_col]])
            if (anyDuplicated(gene_ids)) stop("Repeated gene IDs within ", nm)
            rownames(x) <- gene_ids
            metadata <- if (is.null(stratum_col)) NULL else
                stats::setNames(data.frame(single(stratum_col), stringsAsFactors = FALSE), stratum_col)
            # This temporary metadata contains contrast labels only, not samples.
            de_list[[nm]] <- list(results = x, metadata = metadata,
                                  group1 = group1, group2 = group2)
        }
    }
    if (!is.list(de_list) || !length(de_list) || is.null(names(de_list)) ||
        anyNA(names(de_list)) || any(names(de_list) == "") || anyDuplicated(names(de_list))) {
        stop("de_list must be a nonempty, uniquely named list.")
    }
    tables <- infos <- qcs <- vector("list", length(de_list))
    for (i in seq_along(de_list)) {
        nm <- names(de_list)[i]
        obj <- de_list[[i]]
        x <- as.data.frame(obj$results)
        needed <- c(gene_col, "logFC", "logCPM", "F", "PValue", "FDR")
        if (!all(needed %in% names(x))) stop("Missing edgeR columns in ", nm)
        if (!is.null(obj$test$df.test) && any(obj$test$df.test != 1)) {
            stop("Directional enrichment requires one-coefficient contrasts: ", nm)
        }
        if (!all(vapply(x[c("logFC", "logCPM", "F", "PValue", "FDR")], is.numeric, logical(1)))) {
            stop("edgeR statistics must be numeric: ", nm)
        }
        if (any(!is.finite(x$F)) || any(x$F < 0) ||
            any(!is.finite(x$logFC)) || any(!is.finite(x$logCPM))) stop("Invalid statistics in ", nm)
        for (p in c("PValue", "FDR")) {
            if (any(!is.na(x[[p]]) & (!is.finite(x[[p]]) | x[[p]] < 0 | x[[p]] > 1))) {
                stop("Invalid probabilities in ", nm)
            }
        }
        if (length(obj$group1) != 1L || length(obj$group2) != 1L ||
            anyNA(c(obj$group1, obj$group2))) stop("Missing group labels in ", nm)
        if (is.null(stratum_col)) {
            stratum <- "All"
        } else {
            if (is.null(obj$metadata) || !stratum_col %in% names(obj$metadata)) {
                stop("Missing stratum metadata in ", nm)
            }
            stratum <- unique(as.character(obj$metadata[[stratum_col]]))
            if (length(stratum) != 1L || anyNA(stratum) || stratum == "") {
                stop("Each comparison must belong to exactly one stratum: ", nm)
            }
        }
        n_tested <- nrow(x)
        x$GeneID <- rownames(x)
        x$Gene <- trimws(as.character(x[[gene_col]]))
        x <- x[!is.na(x$Gene) & x$Gene != "", , drop = FALSE]
        n_annotated <- nrow(x)
        x <- x[order(-x$logCPM, x$GeneID), , drop = FALSE]
        x <- x[!duplicated(x$Gene), , drop = FALSE]
        if (!nrow(x)) stop("No annotated genes in ", nm)
        x$ResultName <- nm
        x$Subtype <- stratum
        x$Condition <- as.character(obj$group2)
        x$Reference <- as.character(obj$group1)
        x$Comparison <- paste0(obj$group2, "_vs_", obj$group1)
        rownames(x) <- NULL
        tables[[i]] <- x
        infos[[i]] <- data.frame(ResultName = nm, Subtype = stratum,
                                Condition = as.character(obj$group2),
                                Reference = as.character(obj$group1),
                                Comparison = paste0(obj$group2, "_vs_", obj$group1))
        qcs[[i]] <- data.frame(ResultName = nm, Tested_rows = n_tested,
                               Annotated_rows = n_annotated, Unique_symbols = nrow(x),
                               Duplicate_symbols_removed = n_annotated - nrow(x))
    }
    result <- list(table = do.call(rbind, tables), comparison_info = do.call(rbind, infos),
                   qc = do.call(rbind, qcs), settings = list(gene_col = gene_col, stratum_col = stratum_col))
    result$table$FDR_global <- stats::p.adjust(result$table$PValue, "BH")
    class(result) <- c("bulk_gene_de", "list")
    result
}


#' Gene-expression heatmap for the union of significant edgeR genes
#' Uses all prepared logCPM rows, avoiding selection of expression rows by DE
#' significance. Duplicate symbols use highest mean expression across samples.
#' Requires plot_bulk_tf_heatmap from the shared helper, used internally as a
#' general feature-matrix renderer. The displayed values are gene logCPM z-scores.
#'
#' @param bulk A prepared expression object returned by
#'   [prepare_bulk_expression()].
#' @param gene_de A collected differential-expression object returned by
#'   [prepare_bulk_gene_de()].
#' @param gene_col Gene-symbol column in `bulk$genes`.
#' @param subtype_col,condition_col Metadata columns defining the row strata
#'   and experimental conditions used in the display.
#' @param subtype_order,condition_order Explicit biological display orders.
#' @param subtype_colors,condition_colors Named colour vectors covering the
#'   requested subtype and condition levels.
#' @param fdr_col,fdr_cutoff Differential-expression significance column and
#'   threshold used to form the union of displayed genes.
#' @param column_mode Either `"mean"` for grouped means or `"samples"` for
#'   individual biological samples.
#' @param top_n `NULL` to retain the complete significant union or a positive
#'   integer limiting genes by their best FDR across comparisons.
#' @param ... Additional arguments forwarded to [plot_bulk_tf_heatmap()].
#'
#' @return A heatmap result list containing the plot, raw and z-score matrices,
#'   annotations, selected/plotted/dropped genes, clustering information,
#'   sample counts, settings, and the resolved `gene_map`.
#' @export
plot_bulk_gene_expression_heatmap <- function(
    bulk, gene_de, gene_col = "gene_name",
    subtype_col = "ATM_subtype", condition_col = "Condition",
    subtype_order = c("pATM", "cATM", "sATM"),
    condition_order = c("ND", "2W", "1M", "2M"),
    subtype_colors, condition_colors, fdr_col = "FDR", fdr_cutoff = 0.05,
    column_mode = "mean", top_n = NULL, ...
) {
    if (!all(c("logCPM", "genes", "metadata") %in% names(bulk)) ||
        !is.list(gene_de) || !"table" %in% names(gene_de)) stop("Invalid prepared inputs.")
    mat <- as.matrix(bulk$logCPM)
    anno <- as.data.frame(bulk$genes)
    if (!is.numeric(mat) || any(!is.finite(mat)) || is.null(rownames(mat)) ||
        anyDuplicated(rownames(mat)) || !setequal(rownames(mat), rownames(anno)) ||
        anyDuplicated(rownames(anno)) || !gene_col %in% names(anno)) stop("Invalid gene expression/annotation.")
    anno <- anno[rownames(mat), , drop = FALSE]
    symbols <- trimws(as.character(anno[[gene_col]]))
    valid <- !is.na(symbols) & symbols != ""
    map <- data.frame(GeneID = rownames(mat)[valid], Gene = symbols[valid],
                       MeanLogCPM = rowMeans(mat)[valid])
    map <- map[order(-map$MeanLogCPM, map$GeneID), , drop = FALSE]
    map <- map[!duplicated(map$Gene), , drop = FALSE]
    mat <- mat[map$GeneID, , drop = FALSE]
    rownames(mat) <- map$Gene
    de <- gene_de$table
    if (!fdr_col %in% names(de)) stop("Missing requested FDR column.")
    if (!is.null(top_n)) {
        if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
            top_n < 1 || top_n != floor(top_n)) stop("Invalid top_n.")
        eligible <- de[!is.na(de[[fdr_col]]) & de[[fdr_col]] < fdr_cutoff &
                         de$Subtype %in% subtype_order, , drop = FALSE]
        ids <- unique(eligible$Gene)
        best <- vapply(ids, function(id) min(eligible[[fdr_col]][eligible$Gene == id]), numeric(1))
        ids <- head(ids[order(best, ids)], top_n)
        de <- de[de$Gene %in% ids, , drop = FALSE]
    }
    de$TF <- de$Gene
    result <- plot_bulk_tf_heatmap(
        tf = list(activity_matrix = mat, metadata = bulk$metadata),
        de_results = de, subtype_col = subtype_col, condition_col = condition_col,
        de_subtype_col = "Subtype", subtype_order = subtype_order,
        condition_order = condition_order, subtype_colors = subtype_colors,
        condition_colors = condition_colors, fdr_col = fdr_col,
        fdr_cutoff = fdr_cutoff, column_mode = column_mode,
        feature_label = "gene expression", ...
    )
    names(result)[names(result) == "significant_tfs"] <- "significant_genes"
    names(result)[names(result) == "plotted_tfs"] <- "plotted_genes"
    names(result)[names(result) == "dropped_tfs"] <- "dropped_genes"
    result$gene_map <- map
    result
}


#' GSEA and optional up/down ORA for every edgeR contrast
#' @param gene_de Output of prepare_bulk_gene_de.
#' @param OrgDb AnnotationDbi organism database; required for GO.
#' @param TERM2GENE Optional two-column custom term-to-SYMBOL mapping. When
#'   supplied, use generic GSEA/enricher; otherwise use gseGO/enrichGO.
#' @param TERM2NAME Optional two-column custom term names.
#' @param ont Gene Ontology namespace passed to `clusterProfiler` when `OrgDb`
#'   is used.
#' @param ranking Ranking statistic: signed square-root F statistic or log-fold
#'   change multiplied by negative log10 p-value.
#' @param minGSSize,maxGSSize Minimum and maximum gene-set sizes tested.
#' @param exponent Weighting exponent used by the GSEA algorithm.
#' @param run_ora Also run separate up/down ORA with each contrast's eligible
#'   measured genes as its own universe. Empty DEG sets are logged, not errors.
#' @param gene_fdr_cutoff,gene_logfc_cutoff FDR and absolute log-fold-change
#'   thresholds used to define up- and downregulated ORA gene lists.
#' @param seed Integer random seed used for reproducible GSEA execution.
#' @param verbose Logical scalar. Whether to report progress for each contrast.
#' @return All tested enrichment results, objects, rankings, universes and QC.
#'   pvalueCutoff=1 retains nonsignificant tests for cross-comparison heatmaps.
#' @export
run_bulk_gene_enrichment <- function(
    gene_de, OrgDb = NULL, ont = "BP", TERM2GENE = NULL, TERM2NAME = NULL,
    ranking = c("signed_sqrt_F", "logFC_logP"),
    minGSSize = 15L, maxGSSize = 500L, exponent = 1,
    run_ora = TRUE, gene_fdr_cutoff = 0.05, gene_logfc_cutoff = 0,
    seed = 123L, verbose = TRUE
) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("clusterProfiler is required.")
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("AnnotationDbi is required.")
    ranking <- match.arg(ranking)
    if (!is.list(gene_de) || !all(c("table", "comparison_info") %in% names(gene_de))) stop("Invalid gene_de.")
    if (!ont %in% c("BP", "MF", "CC", "ALL")) stop("Invalid GO ontology.")
    if (!is.numeric(minGSSize) || length(minGSSize) != 1L || !is.finite(minGSSize) ||
        minGSSize < 2 || minGSSize != floor(minGSSize) ||
        !is.numeric(maxGSSize) || length(maxGSSize) != 1L || !is.finite(maxGSSize) ||
        maxGSSize < minGSSize || maxGSSize != floor(maxGSSize)) stop("Invalid gene-set sizes.")
    if (!is.numeric(exponent) || length(exponent) != 1L || !is.finite(exponent) || exponent < 0) stop("Invalid exponent.")
    for (x in list(run_ora, verbose)) {
        if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("run_ora/verbose must be logical scalars.")
    }
    if (!is.numeric(gene_fdr_cutoff) || length(gene_fdr_cutoff) != 1L ||
        !is.finite(gene_fdr_cutoff) || gene_fdr_cutoff < 0 || gene_fdr_cutoff > 1 ||
        !is.numeric(gene_logfc_cutoff) || length(gene_logfc_cutoff) != 1L ||
        !is.finite(gene_logfc_cutoff) || gene_logfc_cutoff < 0) stop("Invalid DEG thresholds.")
    custom <- !is.null(TERM2GENE)
    if (custom) {
        if (!is.data.frame(TERM2GENE) || ncol(TERM2GENE) != 2L || !nrow(TERM2GENE) || anyNA(TERM2GENE)) {
            stop("TERM2GENE must be a nonempty two-column term-to-gene-SYMBOL data frame.")
        }
        TERM2GENE[] <- lapply(TERM2GENE, as.character)
        if (any(TERM2GENE[[1]] == "") || any(TERM2GENE[[2]] == "")) stop("Blank term/gene IDs.")
        TERM2GENE <- unique(TERM2GENE)
        if (!is.null(TERM2NAME) && (!is.data.frame(TERM2NAME) || ncol(TERM2NAME) != 2L || anyNA(TERM2NAME))) {
            stop("TERM2NAME must have two columns with no missing values.")
        }
        eligible_symbols <- unique(TERM2GENE[[2]])
    } else {
        if (is.null(OrgDb) || !inherits(OrgDb, "OrgDb")) stop("Supply an organism OrgDb for GO enrichment.")
        eligible_symbols <- AnnotationDbi::keys(OrgDb, keytype = "SYMBOL")
    }
    objects <- rankings <- universes <- ora_objects <- list()
    tables <- ora_tables <- qcs <- list()
    info <- gene_de$comparison_info
    add_info <- function(x, i) {
        if (!nrow(x)) return(NULL)
        for (nm in names(info)) x[[nm]] <- info[[nm]][i]
        x
    }
    for (i in seq_len(nrow(info))) {
        nm <- info$ResultName[i]
        if (verbose) message("\nEnrichment: ", nm)
        x <- gene_de$table[gene_de$table$ResultName == nm, , drop = FALSE]
        n_input <- nrow(x)
        x <- x[x$Gene %in% eligible_symbols, , drop = FALSE]
        if (!nrow(x)) stop("No genes match annotation in ", nm, "; check species/namespace.")
        if (anyDuplicated(x$Gene)) stop("Duplicated ranked symbols in ", nm)
        if (ranking == "signed_sqrt_F") {
            stat <- sign(x$logFC) * sqrt(x$F)
        } else {
            if (any(!is.finite(x$PValue))) stop("Nonfinite P values for ranking.")
            stat <- x$logFC * -log10(pmax(x$PValue, .Machine$double.xmin))
        }
        if (any(!is.finite(stat)) || stats::sd(stat) == 0) stop("Invalid ranking in ", nm)
        geneList <- stats::setNames(stat, x$Gene)
        # Deterministic tie order; tied statistics are retained without jitter.
        geneList <- geneList[order(-geneList, names(geneList))]
        if (!any(geneList > 0) || !any(geneList < 0)) {
            warning("Ranking is one-sided in ", nm, "; inspect this comparison.")
        }
        rankings[[nm]] <- geneList
        universes[[nm]] <- names(geneList)
        base_args <- list(geneList = geneList, exponent = exponent,
                          minGSSize = minGSSize, maxGSSize = maxGSSize,
                          pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                          verbose = verbose, seed = FALSE, by = "fgsea")
        # Preserve the caller's RNG state while making this analysis reproducible.
        had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
        restore_rng <- function() {
            if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }
        set.seed(seed)
        obj <- tryCatch({
            if (custom) {
                do.call(clusterProfiler::GSEA, c(base_args, list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)))
            } else {
                do.call(clusterProfiler::gseGO, c(base_args, list(OrgDb = OrgDb, keyType = "SYMBOL", ont = ont)))
            }
        }, finally = restore_rng())
        if (is.null(obj)) stop("No eligible gene sets for GSEA in ", nm, "; adjust set-size bounds.")
        objects[[nm]] <- obj
        g <- as.data.frame(obj)
        tables[[nm]] <- add_info(g, i)
        up <- x$Gene[!is.na(x$FDR) & x$FDR < gene_fdr_cutoff & x$logFC > gene_logfc_cutoff]
        down <- x$Gene[!is.na(x$FDR) & x$FDR < gene_fdr_cutoff & x$logFC < -gene_logfc_cutoff]
        qcs[[nm]] <- data.frame(ResultName = nm, Tested_symbols = n_input,
                                Annotated_symbols = nrow(x), Unmatched_symbols = n_input - nrow(x),
                                Ranking_ties = sum(duplicated(stat)), Up_genes = length(up),
                                Down_genes = length(down), GSEA_terms_returned = nrow(g))
        if (run_ora) {
            ora_objects[[nm]] <- list()
            for (direction in c("Up", "Down")) {
                genes <- if (direction == "Up") up else down
                if (!length(genes)) next
                args <- list(gene = genes, universe = universes[[nm]],
                             minGSSize = minGSSize, maxGSSize = maxGSSize,
                             pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
                if (custom) {
                    o <- do.call(clusterProfiler::enricher, c(args, list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)))
                } else {
                    o <- do.call(clusterProfiler::enrichGO, c(args, list(OrgDb = OrgDb, keyType = "SYMBOL", ont = ont)))
                }
                if (is.null(o)) next
                ora_objects[[nm]][[direction]] <- o
                t <- add_info(as.data.frame(o), i)
                if (!is.null(t)) {
                    t$Direction <- direction
                    ora_tables[[paste(nm, direction, sep = "__")]] <- t
                }
            }
        }
    }
    gsea <- do.call(rbind, tables)
    ora <- do.call(rbind, ora_tables)
    if (is.null(gsea)) gsea <- data.frame()
    if (is.null(ora)) ora <- data.frame()
    if (nrow(gsea)) gsea$FDR_global <- stats::p.adjust(gsea$pvalue, "BH")
    if (nrow(ora)) ora$FDR_global <- stats::p.adjust(ora$pvalue, "BH")
    list(gsea = gsea, ora = ora, gsea_objects = objects, ora_objects = ora_objects,
         rankings = rankings, universes = universes, qc = do.call(rbind, qcs),
         comparison_info = info, annotation = if (custom) list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
             else list(OrgDb_metadata = AnnotationDbi::metadata(OrgDb), ont = ont),
         settings = list(ranking = ranking, minGSSize = minGSSize, maxGSSize = maxGSSize,
                          exponent = exponent, run_ora = run_ora, seed = seed,
                          gene_fdr_cutoff = gene_fdr_cutoff, gene_logfc_cutoff = gene_logfc_cutoff),
         session = utils::sessionInfo())
}


#' Plot gene logFC or pathway NES over an explicit comparison manifest
#' Missing/untested cells remain NA (grey), never zero. Selection uses the union
#' of significant features; top_n=NULL keeps all. Row clustering uses 1-Pearson-r
#' with average linkage. If some distances are undefined, all rows are preserved
#' and clustering is disabled with an explicit warning.
#'
#' @param results Data frame containing one row per feature and comparison,
#'   including effect values and multiple-testing results.
#' @param comparison_info Data frame defining the expected subtype-condition
#'   comparison grid and its display order.
#' @param feature_col,value_col,fdr_col Columns containing feature identifiers,
#'   displayed effect values, and adjusted p-values.
#' @param fdr_cutoff Numeric significance threshold used to select the union of
#'   displayed features.
#' @param label_col `NULL` or a column containing alternative row labels.
#' @param subtype_order,condition_order Explicit column-annotation orders.
#' @param subtype_colors,condition_colors Named annotation colour vectors.
#' @param top_n `NULL` for all significant features or a positive integer
#'   limiting rows by their best FDR.
#' @param row_zscore Logical scalar. Whether each feature is standardized across
#'   comparison columns before display.
#' @param color Heatmap colour vector with at least two valid colours.
#' @param value_limit `NULL` or a positive symmetric display limit.
#' @param label_width Positive integer used to wrap long row labels.
#' @param fontsize_row Row-label font size.
#' @param width,height Saved figure dimensions in inches; `NULL` height adapts
#'   to the number of selected features.
#' @param save `NULL` or a PDF/PNG/TIFF output filename.
#' @param subtype_col,condition_col Column names in `comparison_info` containing
#'   subtype and condition labels.
#'
#' @return A list containing the `pheatmap` object, raw and displayed matrices,
#'   FDR matrix, row and column annotations, optional row tree, selected
#'   features, and aligned comparison manifest.
#' @export
plot_bulk_comparison_heatmap <- function(
    results,
    comparison_info,
    feature_col = "Gene",
    value_col = "logFC",
    fdr_col = "FDR",
    fdr_cutoff = 0.05,
    label_col = NULL,

    subtype_order = c("pATM", "cATM", "sATM"),
    condition_order = c("2W", "1M", "2M"),
    subtype_colors,
    condition_colors,

    top_n = NULL,
    row_zscore = FALSE,

    color = BULK_ACTIVITY_DIVERGING_COLORS,

    value_limit = NULL,
    label_width = 55L,
    fontsize_row = 9,
    width = 11,
    height = NULL,
    save = NULL,

    subtype_col = "Subtype",
    condition_col = "Condition"
) {

    # -------------------------------------------------------------------------
    # Validate inputs
    # -------------------------------------------------------------------------

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package `pheatmap` is required.")
    }

    for (column in list(subtype_col, condition_col)) {
        if (!is.character(column) ||
            length(column) != 1L ||
            is.na(column) ||
            !nzchar(column)) {
            stop("subtype_col and condition_col must be column names.")
        }
    }

    if (subtype_col == condition_col) {
        stop("subtype_col and condition_col must be different.")
    }

    needed_results <- c(
        feature_col, value_col, fdr_col, "ResultName"
    )

    needed_info <- c(
        "ResultName", subtype_col, condition_col, "Reference"
    )

    if (!is.data.frame(results) ||
        !all(needed_results %in% names(results))) {
        stop(
            "Missing result columns: ",
            paste(setdiff(needed_results, names(results)), collapse = ", ")
        )
    }

    if (!is.data.frame(comparison_info) ||
        !all(needed_info %in% names(comparison_info))) {
        stop(
            "Missing comparison manifest columns: ",
            paste(setdiff(needed_info, names(comparison_info)), collapse = ", ")
        )
    }

    if (!is.numeric(results[[value_col]]) ||
        !is.numeric(results[[fdr_col]])) {
        stop("Values and FDR must be numeric.")
    }

    observed_fdr <- results[[fdr_col]][!is.na(results[[fdr_col]])]

    if (any(!is.finite(observed_fdr)) ||
        any(observed_fdr < 0 | observed_fdr > 1)) {
        stop("Observed FDR values must be finite and between 0 and 1.")
    }

    if (!is.numeric(fdr_cutoff) ||
        length(fdr_cutoff) != 1L ||
        !is.finite(fdr_cutoff) ||
        fdr_cutoff < 0 ||
        fdr_cutoff > 1) {
        stop("Invalid FDR cutoff.")
    }

    if (!is.logical(row_zscore) ||
        length(row_zscore) != 1L ||
        is.na(row_zscore)) {
        stop("Invalid row_zscore.")
    }

    subtype_order <- as.character(subtype_order)
    condition_order <- as.character(condition_order)

    if (!length(subtype_order) ||
        !length(condition_order) ||
        anyDuplicated(subtype_order) ||
        anyDuplicated(condition_order) ||
        anyNA(c(subtype_order, condition_order)) ||
        any(!nzchar(c(subtype_order, condition_order)))) {
        stop("Invalid subtype or condition orders.")
    }

    if (!is.character(color) ||
        length(color) < 2L ||
        anyNA(color)) {
        stop("Invalid heatmap color vector.")
    }

    invisible(grDevices::col2rgb(color))

    # Allow named vectors or named lists of annotation colors.
    subtype_colors <- unlist(subtype_colors)
    condition_colors <- unlist(condition_colors)

    if (anyDuplicated(names(subtype_colors)) ||
        anyDuplicated(names(condition_colors)) ||
        !all(subtype_order %in% names(subtype_colors)) ||
        !all(condition_order %in% names(condition_colors))) {
        stop("Missing or duplicated annotation color names.")
    }

    invisible(
        grDevices::col2rgb(
            c(
                subtype_colors[subtype_order],
                condition_colors[condition_order]
            )
        )
    )

    # -------------------------------------------------------------------------
    # Select and order the comparison manifest
    # -------------------------------------------------------------------------

    keep_info <- (
        as.character(comparison_info[[subtype_col]]) %in% subtype_order &
        as.character(comparison_info[[condition_col]]) %in% condition_order
    )

    info <- comparison_info[keep_info, , drop = FALSE]

    if (!nrow(info) ||
        anyNA(info[, needed_info, drop = FALSE]) ||
        anyDuplicated(info$ResultName) ||
        anyDuplicated(info[, c(subtype_col, condition_col), drop = FALSE])) {
        stop("Invalid or ambiguous comparison manifest.")
    }

    info <- info[
        order(
            match(as.character(info[[subtype_col]]), subtype_order),
            match(as.character(info[[condition_col]]), condition_order)
        ),
        ,
        drop = FALSE
    ]

    if (nrow(info) != length(subtype_order) * length(condition_order) ||
        length(unique(info$Reference)) != 1L) {
        stop("Require a full comparison grid with one reference.")
    }

    # -------------------------------------------------------------------------
    # Select the union of significant features
    # -------------------------------------------------------------------------

    x <- results[
        results$ResultName %in% info$ResultName,
        ,
        drop = FALSE
    ]

    if (anyDuplicated(x[, c(feature_col, "ResultName"), drop = FALSE]) ||
        anyNA(x[[feature_col]]) ||
        any(!nzchar(as.character(x[[feature_col]])))) {
        stop("Duplicated or missing feature-comparison identifiers.")
    }

    sig <- x[
        !is.na(x[[fdr_col]]) & x[[fdr_col]] < fdr_cutoff,
        ,
        drop = FALSE
    ]

    ids <- unique(as.character(sig[[feature_col]]))

    if (!length(ids)) {
        stop("No significant features for the selected comparisons.")
    }

    best_fdr <- vapply(
        ids,
        function(id) {
            min(
                sig[[fdr_col]][
                    as.character(sig[[feature_col]]) == id
                ]
            )
        },
        numeric(1)
    )

    ids <- ids[order(best_fdr, ids)]

    if (!is.null(top_n)) {

        if (!is.numeric(top_n) ||
            length(top_n) != 1L ||
            !is.finite(top_n) ||
            top_n < 1 ||
            top_n != floor(top_n)) {
            stop("Invalid top_n.")
        }

        ids <- head(ids, top_n)
    }

    x <- x[
        as.character(x[[feature_col]]) %in% ids,
        ,
        drop = FALSE
    ]

    # -------------------------------------------------------------------------
    # Assemble value and FDR matrices
    # Missing/untested combinations remain NA.
    # -------------------------------------------------------------------------

    raw <- matrix(
        NA_real_,
        nrow = length(ids),
        ncol = nrow(info),
        dimnames = list(ids, as.character(info$ResultName))
    )

    fdr_mat <- raw

    indices <- cbind(
        match(as.character(x[[feature_col]]), ids),
        match(as.character(x$ResultName), as.character(info$ResultName))
    )

    raw[indices] <- x[[value_col]]
    fdr_mat[indices] <- x[[fdr_col]]

    raw[!is.finite(raw)] <- NA_real_

    display <- raw

    if (row_zscore) {

        means <- rowMeans(raw, na.rm = TRUE)
        sds <- apply(raw, 1L, stats::sd, na.rm = TRUE)

        if (any(!is.finite(sds) | sds == 0)) {
            stop("Cannot z-score a constant or sparsely observed feature.")
        }

        display <- sweep(
            sweep(raw, 1L, means, "-"),
            1L,
            sds,
            "/"
        )
    }

    # -------------------------------------------------------------------------
    # Row clustering: 1 - Pearson correlation, average linkage
    # Require at least three shared observations for every row pair.
    # -------------------------------------------------------------------------

    row_tree <- FALSE

    if (nrow(display) > 1L) {

        correlation_matrix <- suppressWarnings(
            stats::cor(
                t(display),
                use = "pairwise.complete.obs"
            )
        )

        overlap <- tcrossprod(
            1L * is.finite(display)
        )

        valid <- is.finite(correlation_matrix) & overlap >= 3L
        diag(valid) <- TRUE

        if (all(valid)) {

            distance_matrix <- 1 - correlation_matrix
            distance_matrix[distance_matrix < 0] <- 0
            diag(distance_matrix) <- 0

            row_tree <- stats::hclust(
                stats::as.dist(distance_matrix),
                method = "average"
            )

        } else {

            warning(
                "Undefined or insufficient row correlations; ",
                "showing all selected features without a dendrogram. ",
                "Missing values remain NA.",
                call. = FALSE
            )
        }
    }

    # -------------------------------------------------------------------------
    # Dynamic column annotations
    # -------------------------------------------------------------------------

    ann_col <- data.frame(
        row.names = as.character(info$ResultName)
    )

    ann_col[[subtype_col]] <- factor(
        as.character(info[[subtype_col]]),
        levels = subtype_order
    )

    ann_col[[condition_col]] <- factor(
        as.character(info[[condition_col]]),
        levels = condition_order
    )

    ann_row <- data.frame(row.names = ids)

    ann_colors <- list()

    ann_colors[[subtype_col]] <- subtype_colors[subtype_order]
    ann_colors[[condition_col]] <- condition_colors[condition_order]

    # Binary row annotation for significance in each subtype.
    for (subtype in subtype_order) {

        subtype_columns <- as.character(
            info$ResultName[
                as.character(info[[subtype_col]]) == subtype
            ]
        )

        significant_in_subtype <- rowSums(
            fdr_mat[, subtype_columns, drop = FALSE] < fdr_cutoff,
            na.rm = TRUE
        ) > 0

        annotation_name <- paste0("Sig_", subtype)

        ann_row[[annotation_name]] <- factor(
            ifelse(significant_in_subtype, "Yes", "No"),
            levels = c("No", "Yes")
        )

        ann_colors[[annotation_name]] <- c(
            No = "#F0F0F0",
            Yes = unname(subtype_colors[subtype])
        )
    }

    # -------------------------------------------------------------------------
    # Row labels
    # -------------------------------------------------------------------------

    labels <- ids

    if (!is.null(label_col)) {

        if (!label_col %in% names(x)) {
            stop("Missing label_col.")
        }

        labels <- as.character(
            x[[label_col]][
                match(ids, as.character(x[[feature_col]]))
            ]
        )

        # Fall back to the stable feature ID for missing descriptions.
        missing_labels <- is.na(labels) | !nzchar(labels)
        labels[missing_labels] <- ids[missing_labels]

        labels <- vapply(
            labels,
            function(label) {
                paste(
                    strwrap(label, width = label_width),
                    collapse = "\n"
                )
            },
            character(1)
        )
    }

    # -------------------------------------------------------------------------
    # Color limits and significance stars
    # -------------------------------------------------------------------------

    if (is.null(value_limit)) {

        if (!any(is.finite(display))) {
            stop("No finite values are available for the selected features.")
        }

        value_limit <- max(abs(display), na.rm = TRUE)
    }

    if (!is.numeric(value_limit) ||
        length(value_limit) != 1L ||
        !is.finite(value_limit) ||
        value_limit < 0) {
        stop("Invalid value_limit.")
    }

    if (value_limit == 0) {
        value_limit <- 1
    }

    clipped <- display

    clipped[
        !is.na(clipped) & clipped > value_limit
    ] <- value_limit

    clipped[
        !is.na(clipped) & clipped < -value_limit
    ] <- -value_limit

    stars <- ifelse(
        !is.na(fdr_mat) & fdr_mat < fdr_cutoff,
        "*",
        ""
    )

    if (is.null(height)) {
        height <- max(
            5,
            2.5 + 0.18 * sum(
                vapply(
                    strsplit(labels, "\n"),
                    length,
                    integer(1)
                )
            )
        )
    }

    subtype_counts <- as.integer(
        table(
            factor(
                as.character(info[[subtype_col]]),
                levels = subtype_order
            )
        )
    )

    gaps_col <- if (length(subtype_counts) > 1L) {
        head(cumsum(subtype_counts), -1L)
    } else {
        NULL
    }

    # -------------------------------------------------------------------------
    # Plot and optionally save
    # -------------------------------------------------------------------------

    args <- list(
        mat = clipped,
        scale = "none",

        cluster_rows = row_tree,
        cluster_cols = FALSE,

        annotation_col = ann_col,
        annotation_row = ann_row,
        annotation_colors = ann_colors,
        drop_levels = FALSE,

        color = color,
        na_col = "#D9D9D9",
        border_color = NA,

        breaks = seq(
            -value_limit,
            value_limit,
            length.out = length(color) + 1L
        ),

        labels_row = labels,
        labels_col = paste(
            as.character(info[[condition_col]]),
            "vs",
            as.character(info$Reference)
        ),

        gaps_col = gaps_col,

        display_numbers = stars,
        number_color = "black",
        fontsize_number = 8,
        fontsize_row = fontsize_row,
        fontsize_col = 10,
        angle_col = "45",

        main = paste0(
            "Significant features: ",
            value_col,
            if (row_zscore) " (row z-score)" else "",
            "\n* ",
            fdr_col,
            " < ",
            fdr_cutoff
        ),

        silent = TRUE
    )

    if (!is.null(save)) {

        if (!is.character(save) ||
            length(save) != 1L ||
            is.na(save) ||
            !nzchar(save)) {
            stop("Invalid save filename.")
        }

        if (!dir.exists(dirname(save))) {
            dir.create(dirname(save), recursive = TRUE)
        }

        args$filename <- save
        args$width <- width
        args$height <- height
    }

    plot <- do.call(
        pheatmap::pheatmap,
        args
    )

    grid::grid.newpage()
    grid::grid.draw(plot$gtable)

    list(
        plot = plot,
        raw_matrix = raw,
        display_matrix = display,
        fdr_matrix = fdr_mat,
        annotation_row = ann_row,
        annotation_col = ann_col,
        row_tree = row_tree,
        selected_features = ids,
        comparison_info = info
    )
}


#' Test an overall condition response or a numerical-duration trend with edgeR
#' All selected samples are fitted together. Omnibus F statistics are unsigned;
#' never pass them to a signed enrichment ranking without an actual contrast.
#' @param test Omnibus tests all non-reference condition coefficients jointly;
#'   linear tests one numerical duration coefficient.
#' @param duration_col Required for linear mode, with verified numerical units.
#' @param covariates Optional categorical/numeric adjustment columns. A verified
#'   subject column can be supplied as a fixed-effect covariate when estimable.
#' @param dat A featureCounts project containing aligned raw counts, metadata,
#'   and gene annotations.
#' @param group_col Metadata column defining the experimental conditions.
#' @param group_order Explicit factor order with the reference condition first.
#' @param subset Optional unquoted expression evaluated within sample metadata.
#' @param robust Logical scalar passed to edgeR dispersion and quasi-likelihood
#'   fitting functions.
#' @param verbose Logical scalar. Whether to print group counts and result totals.
#'
#' @return A list containing the method label, complete edgeR result table,
#'   filtered normalized `DGEList`, fit, quasi-likelihood test, design matrix,
#'   tested coefficients, selected metadata, condition order, duration and
#'   covariate settings, gene filter, and session information.
#' @export
edgeR_condition_response <- function(
    dat, group_col = "Condition", group_order = c("ND", "2W", "1M", "2M"),
    subset = NULL, test = c("omnibus", "linear"), duration_col = NULL,
    covariates = NULL, robust = TRUE, verbose = TRUE
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is required.")
    test <- match.arg(test)
    if (!is.list(dat) || !all(c("counts", "metadata", "genes") %in% names(dat))) stop("Invalid dat.")
    meta <- as.data.frame(dat$metadata)
    counts <- as.matrix(dat$counts)
    anno <- as.data.frame(dat$genes)
    if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) ||
        anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
        !identical(colnames(counts), rownames(meta)) || !identical(rownames(counts), rownames(anno))) {
        stop("Counts, metadata and genes must be finite, nonnegative and aligned.")
    }
    if (!group_col %in% names(meta) || anyDuplicated(group_order) || anyNA(group_order) || length(group_order) < 2L) stop("Invalid condition specification.")
    for (value in list(robust, verbose)) {
        if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("robust/verbose must be logical scalars.")
    }
    expr <- substitute(subset)
    selected <- rep(TRUE, nrow(meta))
    if (!identical(expr, quote(NULL))) {
        selected <- eval(expr, envir = meta, enclos = parent.frame())
        if (!is.logical(selected) || length(selected) != nrow(meta) || anyNA(selected)) stop("Invalid subset expression.")
    }
    selected <- selected & !is.na(meta[[group_col]]) & as.character(meta[[group_col]]) %in% group_order
    meta <- meta[selected, , drop = FALSE]
    counts <- counts[, rownames(meta), drop = FALSE]
    group <- factor(as.character(meta[[group_col]]), levels = group_order)
    if (any(table(group) < 2L)) stop("At least two samples are required in every requested condition.")
    if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates) ||
        anyDuplicated(covariates) || !all(covariates %in% names(meta)))) stop("Invalid covariates.")
    if (group_col %in% covariates || (!is.null(duration_col) && duration_col %in% covariates)) stop("Do not duplicate the main predictor as a covariate.")
    model_data <- data.frame(row.names = rownames(meta))
    if (test == "omnibus") {
        model_data$.group <- group
        predictor <- ".group"
    } else {
        if (is.null(duration_col) || length(duration_col) != 1L || !duration_col %in% names(meta)) {
            stop("Linear trends require a verified numerical duration_col.")
        }
        duration <- meta[[duration_col]]
        if (!is.numeric(duration) || any(!is.finite(duration)) || length(unique(duration)) < 2L) stop("Invalid numerical durations.")
        model_data$.duration <- duration
        predictor <- ".duration"
    }
    adjustment <- character()
    for (i in seq_along(covariates)) {
        value <- meta[[covariates[i]]]
        if (anyNA(value) || length(unique(value)) < 2L) stop("Missing or constant covariate: ", covariates[i])
        if (is.character(value) || is.factor(value) || is.logical(value)) {
            value <- droplevels(factor(value))
        } else if (!is.numeric(value) || any(!is.finite(value))) {
            stop("Invalid covariate type.")
        }
        nm <- paste0(".cov", i)
        model_data[[nm]] <- value
        adjustment <- c(adjustment, nm)
    }
    design <- stats::model.matrix(stats::reformulate(c(predictor, adjustment)), model_data)
    if (qr(design)$rank < ncol(design)) stop("Rank-deficient design: check condition/batch/subject confounding.")
    if (nrow(design) <= ncol(design)) stop("No residual degrees of freedom.")
    dge <- edgeR::DGEList(counts = counts, genes = anno)
    keep <- edgeR::filterByExpr(dge, design = design)
    if (!any(keep)) stop("No genes passed expression filtering.")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge <- edgeR::estimateDisp(dge, design = design, robust = robust)
    fit <- edgeR::glmQLFit(dge, design = design, robust = robust)
    coef <- if (test == "omnibus") which(startsWith(colnames(design), ".group")) else match(".duration", colnames(design))
    qlf <- edgeR::glmQLFTest(fit, coef = coef)
    results <- edgeR::topTags(qlf, n = Inf)$table
    if (verbose) {
        message("\nCondition analysis: ", test)
        print(table(group))
        message("Tested genes: ", nrow(results), "; FDR < 0.05: ", sum(results$FDR < 0.05, na.rm = TRUE))
    }
    list(method = paste0("edgeR_condition_", test), results = results,
         dge = dge, fit = fit, qlf = qlf, design = design, tested_coefficients = coef,
         metadata = meta, group_order = group_order, duration_col = duration_col,
         covariates = covariates, keep_genes = keep, session = utils::sessionInfo())
}


#' Plot bulk expression values with optional within-subject pairing
#'
#' Display observations across an ordered categorical axis, optionally connect
#' repeated measurements belonging to the same verified subject, and facet by
#' an experimental variable. This function expects tidy, already-normalized
#' expression data and does not perform normalization or statistical testing.
#'
#' @param df Data frame containing the plotting, grouping, and optional pairing
#'   columns.
#' @param x,y Column names used for the categorical x-axis and numeric
#'   expression y-axis.
#' @param pair_by `NULL` for unconnected observations or a column containing
#'   subject/pair identifiers to connect within facets.
#' @param facet_by `NULL` or a column used to create facet panels.
#' @param x_order `NULL` or an explicit character order for x-axis levels.
#' @param colors `NULL`, an unnamed colour vector assigned in level order, or a
#'   named vector covering every displayed x-axis level.
#' @param point_size,point_alpha,point_stroke Point size, opacity, and outline
#'   width.
#' @param line_width,line_alpha,line_color Appearance of optional paired lines.
#' @param facet_nrow Positive integer giving the number of facet rows.
#' @param facet_scales Facet scale mode passed to [ggplot2::facet_wrap()].
#' @param rotate_x Numeric x-axis label angle in degrees.
#' @param y_title,x_title,title ggplot2-compatible axis and plot titles.
#' @param show_mean Logical scalar. Whether to overlay group mean markers.
#' @param mean_size Positive numeric size for mean markers.
#' @param width,height Saved figure dimensions in inches.
#' @param save `NULL` or a figure output path. PDF output uses Cairo and other
#'   formats use 300 dpi.
#'
#' @return A [ggplot2::ggplot()] object containing the paired or unpaired
#'   expression display.
#' @export
plot_bulk_paired_expression <- function(
    df,
    x = "Sorting",
    y = "Expression",
    pair_by = NULL,
    facet_by = "Condition",
    x_order = NULL,
    colors = NULL,
    point_size = 3.2,
    point_alpha = 0.95,
    point_stroke = 0.45,
    line_width = 0.55,
    line_alpha = 0.55,
    line_color = "grey55",
    facet_nrow = 1,
    facet_scales = "fixed",
    rotate_x = 35,
    y_title = expression(TMM ~ normalized ~ log[2] ~ CPM),
    x_title = NULL,
    title = NULL,
    show_mean = FALSE,
    mean_size = 4,
    width = 8,
    height = 4,
    save = NULL
) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package `ggplot2` is required.", call. = FALSE)
    }

    # ============================================================
    # Check columns
    # ============================================================

    required_cols <- c(
        x,
        y,
        pair_by,
        facet_by
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(df)
    )

    if (length(missing_cols) > 0) {
        stop(
            "Missing columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    df_plot <- df

    # ============================================================
    # X-axis order
    # ============================================================

    if (!is.null(x_order)) {

        df_plot[[x]] <- factor(
            as.character(df_plot[[x]]),
            levels = x_order
        )

    } else if (!is.factor(df_plot[[x]])) {

        df_plot[[x]] <- factor(
            df_plot[[x]],
            levels = unique(df_plot[[x]])
        )
    }

    x_levels <- levels(
        droplevels(df_plot[[x]])
    )

    # ============================================================
    # Macaron colors
    # ============================================================

    if (is.null(colors)) {

        macaron_colors <- BULK_PAIRED_EXPRESSION_COLORS

        if (length(x_levels) > length(macaron_colors)) {
            stop(
                "More groups than available default colors."
            )
        }

        colors <- macaron_colors[
            seq_along(x_levels)
        ]

        names(colors) <- x_levels

    } else {

        if (!is.character(colors) || anyNA(colors) || any(!nzchar(colors))) {
            stop("colors must contain valid non-missing colour strings.")
        }

        if (is.null(names(colors))) {

            if (length(colors) < length(x_levels)) {
                stop("Not enough colors supplied for the displayed groups.")
            }

            colors <- colors[
                seq_along(x_levels)
            ]

            names(colors) <- x_levels

        } else {

            missing_colors <- setdiff(x_levels, names(colors))

            if (length(missing_colors) > 0L) {
                stop(
                    "Missing colors for group(s): ",
                    paste(missing_colors, collapse = ", ")
                )
            }
        }

        colors <- colors[x_levels]
    }

    tryCatch(
        grDevices::col2rgb(colors),
        error = function(error) {
            stop("colors contains an invalid R colour value.", call. = FALSE)
        }
    )

    # ============================================================
    # Plot
    # ============================================================

    p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(
            x = .data[[x]],
            y = .data[[y]]
        )
    )

    # ------------------------------------------------------------
    # Pairing lines
    # ------------------------------------------------------------

    if (!is.null(pair_by)) {

        if (!is.character(pair_by) || length(pair_by) != 1L ||
            is.na(pair_by) || pair_by == "") {
            stop("pair_by must be one subject-ID column name or NULL.")
        }
        if (anyNA(df_plot[[pair_by]]) ||
            any(trimws(as.character(df_plot[[pair_by]])) == "")) {
            stop("Missing or blank pairing IDs.")
        }
        if (anyDuplicated(df_plot[unique(c(facet_by, pair_by, x))])) {
            stop("Multiple observations for one subject, facet and x group; resolve replicates before pairing.")
        }

        line_mapping <- if (is.null(facet_by)) {
            ggplot2::aes(group = .data[[pair_by]])
        } else {
            ggplot2::aes(
                group = interaction(
                    .data[[facet_by]],
                    .data[[pair_by]]
                )
            )
        }

        p <- p + ggplot2::geom_line(
            mapping = line_mapping,
            colour = line_color,
            linewidth = line_width,
            alpha = line_alpha
        )
    }

    # ------------------------------------------------------------
    # Individual mice
    # ------------------------------------------------------------

    p <- p +
        ggplot2::geom_point(
            ggplot2::aes(
                fill = .data[[x]]
            ),
            shape = 21,
            size = point_size,
            stroke = point_stroke,
            colour = "black",
            alpha = point_alpha
        )

    # ============================================================
    # Optional group mean
    # ============================================================

    if (show_mean) {

        p <- p +
            ggplot2::stat_summary(
                ggplot2::aes(
                    group = .data[[x]]
                ),
                fun = mean,
                geom = "point",
                shape = 23,
                size = mean_size,
                stroke = 0.7,
                fill = "white",
                colour = "black"
            )
    }

    # ============================================================
    # Facets
    # ============================================================

    if (!is.null(facet_by)) {
        facet_formula <- stats::reformulate(facet_by)
        p <- p + ggplot2::facet_wrap(
            facet_formula,
            nrow = facet_nrow,
            scales = facet_scales
        )
    }

    # ============================================================
    # Colors
    # ============================================================

    p <- p +
        ggplot2::scale_fill_manual(
            values = colors,
            drop = FALSE
        )

    # ============================================================
    # Labels
    # ============================================================

    p <- p +
        ggplot2::labs(
            x = x_title,
            y = y_title,
            title = title
        )

    # ============================================================
    # Cell-style theme
    # ============================================================

    p <- p +
        ggplot2::theme_classic(
            base_size = 10
        ) +
        ggplot2::theme(

            panel.grid =
                ggplot2::element_blank(),

            axis.line =
                ggplot2::element_line(
                    colour = "black",
                    linewidth = 0.5
                ),

            axis.ticks =
                ggplot2::element_line(
                    colour = "black",
                    linewidth = 0.45
                ),

            axis.ticks.length =
                grid::unit(
                    1.6,
                    "mm"
                ),

            axis.text =
                ggplot2::element_text(
                    colour = "black",
                    size = 9
                ),

            axis.text.x =
                ggplot2::element_text(
                    angle = rotate_x,
                    hjust = ifelse(
                        rotate_x == 0,
                        0.5,
                        1
                    ),
                    vjust = 1
                ),

            axis.title.x =
                ggplot2::element_text(
                    colour = "black",
                    size = 10
                ),

            axis.title.y =
                ggplot2::element_text(
                    colour = "black",
                    size = 10
                ),

            strip.background =
                ggplot2::element_blank(),

            strip.text =
                ggplot2::element_text(
                    colour = "black",
                    face = "bold",
                    size = 10
                ),

            panel.spacing.x =
                grid::unit(
                    5,
                    "mm"
                ),

            legend.position =
                "none",

            plot.title =
                ggplot2::element_text(
                    colour = "black",
                    face = "bold",
                    size = 11,
                    hjust = 0.5
                ),

            plot.margin =
                ggplot2::margin(
                    5.5,
                    8,
                    5.5,
                    5.5
                )
        )

    # ============================================================
    # Save
    # ============================================================

    if (!is.null(save)) {

        save_args <- list(
            filename = save,
            plot = p,
            width = width,
            height = height,
            units = "in"
        )

        if (
            grepl(
                "\\.pdf$",
                save,
                ignore.case = TRUE
            )
        ) {
            save_args$device <-
                grDevices::cairo_pdf
        } else {
            save_args$dpi <- 300
        }

        do.call(
            ggplot2::ggsave,
            save_args
        )
    }

    return(p)
}


#' Prepare bulk expression for ligand/receptor availability checks
#' CPM has no log prior count. Duplicate symbols use the same representative
#' row in CPM and logCPM, selected by mean logCPM before consulting DE results.
#'
#' @param bulk A prepared expression object returned by
#'   [prepare_bulk_expression()] containing `dge`, `logCPM`, `genes`, and
#'   `metadata` components in matching orders.
#' @param gene_col Gene-symbol column in `bulk$genes`.
#'
#' @return A list containing symbol-indexed `logCPM` and unlogged `cpm`
#'   matrices, aligned sample metadata, and a stable-ID-to-symbol mapping.
#' @export
prepare_bulk_signaling_expression <- function(bulk, gene_col = "gene_name") {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package edgeR is required.")
    if (!all(c("dge", "logCPM", "genes", "metadata") %in% names(bulk))) {
        stop("Use the output of prepare_bulk_expression().")
    }
    x <- bulk$logCPM
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
        is.null(rownames(x)) || is.null(colnames(x)) ||
        anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop("Invalid logCPM matrix.")
    if (!gene_col %in% names(bulk$genes)) stop("Missing gene annotation: ", gene_col)
    if (is.null(rownames(bulk$metadata)) || anyDuplicated(rownames(bulk$metadata)) ||
        !all(colnames(x) %in% rownames(bulk$metadata))) stop("Sample metadata do not match expression.")
    if (!all(rownames(x) %in% rownames(bulk$genes))) stop("Gene annotation does not match expression.")
    cpm <- edgeR::cpm(bulk$dge, log = FALSE)
    if (!all(rownames(x) %in% rownames(cpm)) || !all(colnames(x) %in% colnames(cpm))) {
        stop("DGEList does not match logCPM.")
    }
    cpm <- cpm[rownames(x), colnames(x), drop = FALSE]
    if (any(!is.finite(cpm)) || any(cpm < 0)) stop("Invalid CPM values.")
    symbols <- trimws(as.character(bulk$genes[rownames(x), gene_col]))
    keep <- which(!is.na(symbols) & nzchar(symbols))
    keep <- keep[order(-rowMeans(x)[keep], rownames(x)[keep])]
    keep <- keep[!duplicated(symbols[keep])]
    if (!length(keep)) stop("No valid gene symbols.")
    mapping <- data.frame(GeneID = rownames(x)[keep], Gene = symbols[keep])
    x <- x[keep, , drop = FALSE]
    cpm <- cpm[keep, , drop = FALSE]
    rownames(x) <- rownames(cpm) <- mapping$Gene
    list(logCPM = x, cpm = cpm,
         metadata = bulk$metadata[colnames(x), , drop = FALSE], mapping = mapping)
}


#' Bulk NicheNet: receiver responses followed by separate sender evidence
#' receiver_de/sender_de are outputs of prepare_bulk_gene_de(). Contrast
#' orientation must be Condition versus Reference. All identifiers must use
#' symbols matching the supplied species-specific NicheNet model.
#' Activity scores have no significance/FDR interpretation. The prior is unsigned;
#' a down-response analysis does not establish ligand-mediated repression.
#'
#' @param sender_bulk,receiver_bulk Prepared sender and receiver expression
#'   objects returned by [prepare_bulk_expression()].
#' @param sender_de,receiver_de Collected sender and receiver differential-
#'   expression objects returned by [prepare_bulk_gene_de()].
#' @param ligand_target_matrix Species-matched NicheNet ligand-target matrix.
#' @param lr_network Species-matched ligand-receptor network containing `from`
#'   and `to` identifiers or equivalent recognized columns.
#' @param sender_col,receiver_col Metadata columns identifying sender and
#'   receiver cell populations.
#' @param condition_col Metadata and DE-table column identifying conditions.
#' @param gene_col Gene-symbol annotation column in both prepared objects.
#' @param reference Reference condition used to validate contrast orientation.
#' @param directions Character subset of `"up"` and `"down"` receiver
#'   responses to analyse.
#' @param cpm_cutoff,min_fraction,min_replicates Expression-availability
#'   threshold, required sample fraction, and minimum replicates per context.
#' @param fdr_cutoff,logfc_cutoff Receiver and sender DE thresholds used to
#'   define responsive genes and annotate evidence.
#' @param min_targets,min_background Minimum responsive target-set and measured
#'   background sizes required for a context.
#' @param top_ligands Maximum ligand rankings retained per context.
#' @param prior_top_n Number of prior-weighted targets requested from NicheNet
#'   before final truncation.
#' @param targets_per_ligand Maximum weighted target links retained per ligand.
#' @param verbose Logical scalar. Whether to report progress and context QC.
#'
#' @return A list containing ligand activities, sender-expression evidence,
#'   receptor evidence, weighted targets, QC, receiver gene sets and
#'   backgrounds, sender/receiver gene mappings, and complete settings.
#' @export
run_bulk_nichenet <- function(
    sender_bulk, receiver_bulk, sender_de, receiver_de,
    ligand_target_matrix, lr_network,
    sender_col = "ATM_subtype", receiver_col = "Sorting",
    condition_col = "Condition", gene_col = "gene_name", reference = "ND",
    directions = c("up", "down"), cpm_cutoff = 1, min_fraction = 0.5,
    min_replicates = 2L, fdr_cutoff = 0.05, logfc_cutoff = 0,
    min_targets = 10L, min_background = 100L, top_ligands = 30L,
    prior_top_n = 250L, targets_per_ligand = 20L, verbose = TRUE
) {
    if (!requireNamespace("nichenetr", quietly = TRUE)) stop("Package nichenetr is required.")
    scalar <- function(x, lower, upper = Inf, integer = FALSE) {
        is.numeric(x) && length(x) == 1L && is.finite(x) &&
            x >= lower && x <= upper && (!integer || x == floor(x))
    }
    if (!scalar(cpm_cutoff, 0) || cpm_cutoff == 0 || !scalar(min_fraction, 0, 1) ||
        min_fraction == 0 || !scalar(fdr_cutoff, 0, 1) || !scalar(logfc_cutoff, 0) ||
        !all(vapply(list(min_replicates, min_targets, min_background, top_ligands,
                        prior_top_n, targets_per_ligand), scalar, logical(1), lower = 1, integer = TRUE))) {
        stop("Invalid expression, DE, or size thresholds.")
    }
    if (length(reference) != 1L || is.na(reference) || !nzchar(reference) ||
        !length(directions) || anyNA(directions) || !all(directions %in% c("up", "down"))) {
        stop("Supply a reference and directions chosen from up/down.")
    }
    directions <- unique(directions)
    if (!is.matrix(ligand_target_matrix) || !is.numeric(ligand_target_matrix) ||
        is.null(rownames(ligand_target_matrix)) || is.null(colnames(ligand_target_matrix)) ||
        anyDuplicated(rownames(ligand_target_matrix)) || anyDuplicated(colnames(ligand_target_matrix))) {
        stop("ligand_target_matrix must be numeric, with unique targets in rows and ligands in columns.")
    }
    if (!all(c("from", "to") %in% names(lr_network))) stop("lr_network needs from/to columns.")
    lr <- unique(data.frame(Ligand = as.character(lr_network$from),
                            Receptor = as.character(lr_network$to)))
    lr <- lr[!is.na(lr$Ligand) & !is.na(lr$Receptor) &
             nzchar(lr$Ligand) & nzchar(lr$Receptor), , drop = FALSE]
    s <- prepare_bulk_signaling_expression(sender_bulk, gene_col)
    r <- prepare_bulk_signaling_expression(receiver_bulk, gene_col)
    validate_metadata <- function(x, population_col) {
        cols <- c(population_col, condition_col)
        if (!all(cols %in% names(x$metadata))) stop("Missing metadata: ", paste(cols, collapse = ", "))
        if (anyNA(x$metadata[, cols, drop = FALSE]) ||
            any(vapply(x$metadata[, cols, drop = FALSE], function(z) any(!nzchar(as.character(z))), logical(1)))) {
            stop("Population/condition metadata contain missing labels.")
        }
    }
    validate_metadata(s, sender_col)
    validate_metadata(r, receiver_col)
    validate_de <- function(z) {
        needed <- c("Gene", "Subtype", "Condition", "Reference", "logFC", "FDR", "ResultName")
        if (!is.list(z) || !all(c("table", "comparison_info") %in% names(z)) ||
            !all(needed %in% names(z$table)) ||
            !all(c("Subtype", "Condition", "Reference", "ResultName") %in% names(z$comparison_info))) {
            stop("DE inputs must come from prepare_bulk_gene_de().")
        }
        key <- z$comparison_info[, c("Subtype", "Condition", "Reference"), drop = FALSE]
        if (anyDuplicated(key)) stop("Duplicate population/condition/reference contrasts.")
        if (anyDuplicated(z$table[, c("ResultName", "Gene"), drop = FALSE])) stop("Duplicate DE gene symbols.")
    }
    validate_de(sender_de)
    validate_de(receiver_de)
    info <- receiver_de$comparison_info
    info <- info[info$Reference == reference & info$Condition != reference, , drop = FALSE]
    if (!nrow(info)) stop("No receiver contrasts against ", reference)
    senders <- unique(as.character(s$metadata[[sender_col]]))
    # Validate every requested group before making any predictions.
    group_summary <- function(x, population_col, population, condition) {
        take <- as.character(x$metadata[[population_col]]) == population &
            as.character(x$metadata[[condition_col]]) == condition
        n <- sum(take)
        if (n < min_replicates) stop("Insufficient samples: ", population, " / ", condition)
        fraction <- rowMeans(x$cpm[, take, drop = FALSE] >= cpm_cutoff)
        data.frame(Gene = rownames(x$cpm), MeanLogCPM = rowMeans(x$logCPM[, take, drop = FALSE]),
                   Fraction = fraction, Expressed = fraction >= min_fraction, N = n,
                   row.names = rownames(x$cpm))
    }
    sender_summaries <- list()
    for (condition in unique(c(reference, as.character(info$Condition)))) {
        sender_summaries[[condition]] <- setNames(lapply(senders, function(population) {
            group_summary(s, sender_col, population, condition)
        }), senders)
    }
    for (i in seq_len(nrow(info))) {
        group_summary(r, receiver_col, as.character(info$Subtype[i]), reference)
        group_summary(r, receiver_col, as.character(info$Subtype[i]), as.character(info$Condition[i]))
        for (sender in senders) {
            si <- sender_de$comparison_info
            if (sum(si$Subtype == sender & si$Condition == info$Condition[i] & si$Reference == reference) != 1L) {
                stop("Missing sender DE contrast: ", sender, " / ", info$Condition[i], " vs ", reference)
            }
        }
    }
    activities <- sender_evidence <- receptor_evidence <- targets <- qcs <- list()
    gene_sets <- backgrounds <- list()
    bind <- function(z) if (length(z)) { out <- do.call(rbind, z); rownames(out) <- NULL; out } else data.frame()
    add_de <- function(x, de, prefix) {
        j <- match(x$Gene, de$Gene)
        x[[paste0(prefix, "_logFC")]] <- de$logFC[j]
        x[[paste0(prefix, "_FDR")]] <- de$FDR[j]
        x
    }
    for (i in seq_len(nrow(info))) {
        receiver <- as.character(info$Subtype[i])
        condition <- as.character(info$Condition[i])
        result_name <- as.character(info$ResultName[i])
        rn <- group_summary(r, receiver_col, receiver, reference)
        rc <- group_summary(r, receiver_col, receiver, condition)
        rd <- receiver_de$table[receiver_de$table$ResultName == result_name, , drop = FALSE]
        background <- intersect(rd$Gene, rownames(rn)[rn$Expressed | rc$Expressed])
        background <- intersect(background, rownames(ligand_target_matrix))
        receptors <- rownames(rn)[rn$Expressed | rc$Expressed]
        expressed_ligands <- unique(unlist(lapply(senders, function(sender) {
            sn <- sender_summaries[[reference]][[sender]]
            sc <- sender_summaries[[condition]][[sender]]
            rownames(sn)[sn$Expressed | sc$Expressed]
        }), use.names = FALSE))
        potential <- intersect(intersect(lr$Ligand[lr$Receptor %in% receptors], expressed_ligands),
                               colnames(ligand_target_matrix))
        for (direction in directions) {
            key <- paste(result_name, direction, sep = "::")
            significant <- !is.na(rd$FDR) & rd$FDR < fdr_cutoff &
                (if (direction == "up") rd$logFC > logfc_cutoff else rd$logFC < -logfc_cutoff)
            selected <- rd$Gene[significant]
            geneset <- intersect(selected, background)
            gene_sets[[key]] <- geneset
            backgrounds[[key]] <- background
            reason <- if (length(background) < min_background) "insufficient_background" else
                if (length(geneset) < min_targets) "insufficient_targets" else
                if (!length(potential)) "no_potential_ligands" else "ok"
            qcs[[length(qcs) + 1L]] <- data.frame(Receiver = receiver, Condition = condition,
                Reference = reference, Response = direction, ResultName = result_name,
                Selected_DE = length(selected), Targets = length(geneset),
                Excluded_DE = length(setdiff(selected, geneset)), Background = length(background),
                Potential_ligands = length(potential), Status = reason)
            if (verbose) message(receiver, " / ", condition, " / ", direction, ": ", reason,
                                 " (", length(geneset), " targets; ", length(potential), " ligands)")
            if (reason != "ok") next
            if (any(!is.finite(ligand_target_matrix[background, potential, drop = FALSE]))) {
                stop("Non-finite ligand-target weights in requested model subset.")
            }
            a <- as.data.frame(nichenetr::predict_ligand_activities(
                geneset = geneset, background_expressed_genes = background,
                ligand_target_matrix = ligand_target_matrix, potential_ligands = potential, single = TRUE))
            if (!all(c("test_ligand", "aupr_corrected") %in% names(a))) stop("Unexpected NicheNet activity columns.")
            a <- a[order(-a$aupr_corrected, as.character(a$test_ligand), na.last = TRUE), , drop = FALSE]
            a$Ligand <- as.character(a$test_ligand)
            a$Rank <- seq_len(nrow(a))
            a$Receiver <- receiver
            a$Condition <- condition
            a$Reference <- reference
            a$Response <- direction
            a$ResultName <- result_name
            activities[[length(activities) + 1L]] <- a
            for (sender in senders) {
                sn <- sender_summaries[[reference]][[sender]]
                sc <- sender_summaries[[condition]][[sender]]
                present <- a$Ligand %in% rownames(sn)[sn$Expressed | sc$Expressed]
                se <- a[present, , drop = FALSE]
                if (!nrow(se)) next
                se$Sender <- sender
                se$Gene <- se$Ligand
                j <- match(se$Ligand, sn$Gene)
                se$Ligand_mean_logCPM_reference <- sn$MeanLogCPM[j]
                se$Ligand_mean_logCPM_condition <- sc$MeanLogCPM[j]
                se$Ligand_fraction_reference <- sn$Fraction[j]
                se$Ligand_fraction_condition <- sc$Fraction[j]
                se$Ligand_expressed_reference <- sn$Expressed[j]
                se$Ligand_expressed_condition <- sc$Expressed[j]
                sd <- sender_de$table[sender_de$table$Subtype == sender &
                    sender_de$table$Condition == condition & sender_de$table$Reference == reference, , drop = FALSE]
                se <- add_de(se, sd, "Ligand")
                sender_evidence[[length(sender_evidence) + 1L]] <- se
            }
            links <- lr[lr$Ligand %in% a$Ligand & lr$Receptor %in% receptors, , drop = FALSE]
            links$Gene <- links$Receptor
            j <- match(links$Receptor, rn$Gene)
            links$Receptor_mean_logCPM_reference <- rn$MeanLogCPM[j]
            links$Receptor_mean_logCPM_condition <- rc$MeanLogCPM[j]
            links$Receptor_fraction_reference <- rn$Fraction[j]
            links$Receptor_fraction_condition <- rc$Fraction[j]
            links$Receptor_expressed_reference <- rn$Expressed[j]
            links$Receptor_expressed_condition <- rc$Expressed[j]
            links <- add_de(links, rd, "Receptor")
            links$Receiver <- receiver
            links$Condition <- condition
            links$Reference <- reference
            links$Response <- direction
            links$ResultName <- result_name
            receptor_evidence[[length(receptor_evidence) + 1L]] <- links
            for (ligand in head(a$Ligand, top_ligands)) {
                lt <- as.data.frame(nichenetr::get_weighted_ligand_target_links(
                    ligand = ligand, geneset = geneset, ligand_target_matrix = ligand_target_matrix,
                    n = prior_top_n))
                if (!all(c("target", "weight") %in% names(lt))) stop("Unexpected NicheNet target columns.")
                lt <- lt[!is.na(lt$target) & is.finite(lt$weight), , drop = FALSE]
                lt <- head(lt[order(-lt$weight), , drop = FALSE], targets_per_ligand)
                if (!nrow(lt)) next
                lt$Ligand <- ligand
                lt$Gene <- as.character(lt$target)
                lt <- add_de(lt, rd, "Target")
                lt$Receiver <- receiver
                lt$Condition <- condition
                lt$Reference <- reference
                lt$Response <- direction
                lt$ResultName <- result_name
                targets[[length(targets) + 1L]] <- lt
            }
        }
    }
    list(activities = bind(activities), sender_evidence = bind(sender_evidence),
         receptor_evidence = bind(receptor_evidence), targets = bind(targets), qc = bind(qcs),
         gene_sets = gene_sets, backgrounds = backgrounds,
         gene_mapping = list(sender = s$mapping, receiver = r$mapping),
         settings = list(sender_col = sender_col, receiver_col = receiver_col,
             condition_col = condition_col, reference = reference, directions = directions,
             cpm_cutoff = cpm_cutoff, min_fraction = min_fraction, min_replicates = min_replicates,
             fdr_cutoff = fdr_cutoff, logfc_cutoff = logfc_cutoff, min_targets = min_targets,
             min_background = min_background, top_ligands = top_ligands,
             prior_top_n = prior_top_n, targets_per_ligand = targets_per_ligand,
             nichenetr_version = as.character(utils::packageVersion("nichenetr"))))
}


#' Heatmap of receiver-level ligand rankings; missing contexts remain NA
#' No FDR stars: corrected AUPR is a prediction score, not a DE test.
#'
#' @param result Output returned by [run_bulk_nichenet()].
#' @param response Receiver response direction, either `"up"` or `"down"`.
#' @param receiver_order `NULL` or an explicit receiver display order.
#' @param condition_order Explicit condition display order.
#' @param top_n Positive integer limiting ligands by their best corrected AUPR.
#' @param color Heatmap colour vector with at least two colours.
#' @param receiver_colors,condition_colors `NULL` or named annotation colour
#'   vectors covering the displayed receiver and condition levels.
#' @param save `NULL` or a filename written by `pheatmap`.
#' @param width,height Saved heatmap dimensions in inches.
#'
#' @return Invisibly returns a list containing the activity matrix, column
#'   annotation data frame, and `pheatmap` plot object.
#' @export
plot_bulk_nichenet_activity <- function(
    result, response = c("up", "down"), receiver_order = NULL,
    condition_order = c("2W", "1M", "2M"), top_n = 25L,
    color = BULK_NICHENET_ACTIVITY_COLORS,
    receiver_colors = NULL, condition_colors = NULL, save = NULL,
    width = 11, height = 8
) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Package pheatmap is required.")
    response <- match.arg(response)
    a <- result$activities
    if (!all(c("Receiver", "Condition", "Response", "Ligand", "aupr_corrected") %in% names(a))) {
        stop("No activity results; inspect result$qc.")
    }
    if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != floor(top_n)) stop("Invalid top_n.")
    if (!is.character(color) || length(color) < 2L) stop("color must be a vector of colors.")
    a <- a[a$Response == response, , drop = FALSE]
    if (is.null(receiver_order)) receiver_order <- unique(as.character(result$qc$Receiver))
    if (anyDuplicated(receiver_order) || anyDuplicated(condition_order)) stop("Duplicate display labels.")
    info <- expand.grid(Condition = condition_order, Receiver = receiver_order,
                        stringsAsFactors = FALSE)
    info$Key <- paste(info$Receiver, info$Condition, sep = " / ")
    a$Key <- paste(a$Receiver, a$Condition, sep = " / ")
    a <- a[a$Key %in% info$Key & is.finite(a$aupr_corrected), , drop = FALSE]
    if (!nrow(a)) stop("No finite activity scores for requested contexts.")
    best <- tapply(a$aupr_corrected, a$Ligand, max)
    ligands <- head(names(sort(best, decreasing = TRUE)), top_n)
    mat <- matrix(NA_real_, length(ligands), nrow(info), dimnames = list(ligands, info$Key))
    take <- a$Ligand %in% ligands
    if (anyDuplicated(a[take, c("Ligand", "Key"), drop = FALSE])) stop("Duplicate activity contexts.")
    mat[cbind(match(a$Ligand[take], ligands), match(a$Key[take], info$Key))] <- a$aupr_corrected[take]
    annotation <- info[, c("Receiver", "Condition"), drop = FALSE]
    rownames(annotation) <- info$Key
    anno_colors <- list()
    if (!is.null(receiver_colors)) anno_colors$Receiver <- receiver_colors
    if (!is.null(condition_colors)) anno_colors$Condition <- condition_colors
    if (!length(anno_colors)) anno_colors <- NA
    p <- pheatmap::pheatmap(mat, color = color, cluster_rows = FALSE, cluster_cols = FALSE,
        annotation_col = annotation, annotation_colors = anno_colors, na_col = "#DDDDDD",
        main = paste0("NicheNet corrected AUPR: ASC ", response, " response"),
        border_color = NA, filename = if (is.null(save)) NA else save, width = width, height = height)
    invisible(list(matrix = mat, annotation = annotation, plot = p))
}


#' General gene-name heatmap from prepare_bulk_expression() output
#' Uses existing normalized logCPM; never recalculates normalization or DE.
#' Means are means of logCPM, calculated before optional row standardization.
#' Duplicate symbols use highest mean logCPM across all input samples.
#' color accepts a color vector or a numeric-to-color function (colorRamp2).
#'
#' @param bulk A prepared expression object returned by
#'   [prepare_bulk_expression()].
#' @param genes Character vector of requested gene symbols in display order.
#' @param gene_col Gene-symbol column in `bulk$genes`.
#' @param samples `NULL` for all samples or a character vector of sample IDs.
#' @param column_mode `"samples"` for individual columns or `"mean"` for
#'   metadata-group means.
#' @param group_by Metadata columns defining mean-expression groups.
#' @param column_order Named list defining level order for grouping variables.
#' @param annotation_cols Metadata columns displayed above sample columns.
#' @param annotation_colors `NULL` or a named list of annotation colour maps.
#' @param annotation_row `NULL` or a data frame aligned to resolved gene rows.
#' @param scale Either `"row"` for gene-wise z-scores or `"none"` for logCPM.
#' @param cluster_rows,cluster_cols Logical scalars controlling hierarchical
#'   clustering of rows and columns.
#' @param distance Correlation or Euclidean clustering distance.
#' @param correlation_method Pearson or Spearman method used for correlation
#'   distance.
#' @param clustering_method Linkage method passed to [stats::hclust()].
#' @param color A colour vector or numeric-to-colour mapping function.
#' @param value_limit `NULL` or a positive symmetric colour-scale limit.
#' @param show_colnames Logical scalar. Whether to display column labels.
#' @param fontsize_row,fontsize_col Row and column label sizes.
#' @param title `NULL` or a heatmap title.
#' @param width,height Saved figure dimensions in inches; `NULL` height adapts
#'   to the number of genes.
#' @param save `NULL` or a PDF, PNG, TIFF, JPEG, or BMP output path.
#' @param legend_extra_width_mm Non-negative extra width reserved for legends.
#' @param genes_of_interest `NULL` or requested genes to emphasize in row labels.
#' @param show_rownames Logical scalar. Whether to display gene labels.
#' @param show_row_dendrogram,show_col_dendrogram Logical scalars controlling
#'   dendrogram visibility independently of clustering.
#'
#' @return Invisibly returns a list containing the `pheatmap` object and, when
#'   requested and computable, row and column `hclust` trees.
#' @export
plot_bulk_expression_heatmap <- function(
    bulk,
    genes,
    gene_col = "gene_name",
    samples = NULL,
    column_mode = c("samples", "mean"),
    group_by = NULL,
    column_order = list(),
    annotation_cols = group_by,
    annotation_colors = NULL,
    annotation_row = NULL,
    scale = c("row", "none"),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    distance = c("correlation", "euclidean"),
    correlation_method = c("pearson", "spearman"),
    clustering_method = "ward.D2",
    color = BULK_ACTIVITY_DIVERGING_COLORS,
    value_limit = NULL,
    show_colnames = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    title = NULL,
    width = 10,
    height = NULL,
    save = NULL,
    legend_extra_width_mm = 8,
    genes_of_interest = NULL,
    show_rownames = TRUE,
    show_row_dendrogram = TRUE,
    show_col_dendrogram = TRUE
) {

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package pheatmap is required.")
    }

    column_mode <- match.arg(column_mode)
    scale <- match.arg(scale)
    distance <- match.arg(distance)
    correlation_method <- match.arg(correlation_method)

    if (!clustering_method %in% c(
        "ward.D", "ward.D2", "single", "complete",
        "average", "mcquitty", "median", "centroid"
    )) {
        stop("Invalid clustering_method.")
    }

    for (flag in list(
        cluster_rows, cluster_cols, show_colnames, show_rownames,
        show_row_dendrogram, show_col_dendrogram
    )) {
        if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
            stop("Clustering/display flags must be TRUE/FALSE.")
        }
    }

    if (
        !is.character(genes) ||
        !length(genes) ||
        anyNA(genes) ||
        any(!nzchar(trimws(genes)))
    ) {
        stop("genes must be a nonempty character vector of gene names.")
    }

    genes <- unique(trimws(genes))

    # Optional highlighting; these genes must also be included in `genes`.
    # A character vector or a list of character vectors is accepted.
    if (!is.null(genes_of_interest)) {
        genes_of_interest <- unlist(genes_of_interest, use.names = FALSE)
        if (!is.character(genes_of_interest) ||
            anyNA(genes_of_interest) ||
            any(!nzchar(trimws(genes_of_interest)))) {
            stop("genes_of_interest must contain nonempty gene names.")
        }
        genes_of_interest <- unique(trimws(genes_of_interest))
    }

    if (
        !is.numeric(width) ||
        length(width) != 1L ||
        !is.finite(width) ||
        width <= 0
    ) {
        stop("width must be a positive number in inches.")
    }

    if (!is.null(height)) {
        if (
            !is.numeric(height) ||
            length(height) != 1L ||
            !is.finite(height) ||
            height <= 0
        ) {
            stop("height must be a positive number in inches.")
        }
    }

    if (
        !is.numeric(legend_extra_width_mm) ||
        length(legend_extra_width_mm) != 1L ||
        !is.finite(legend_extra_width_mm) ||
        legend_extra_width_mm < 0
    ) {
        stop("legend_extra_width_mm must be a nonnegative number.")
    }

    if (!is.null(save)) {
        if (
            !is.character(save) ||
            length(save) != 1L ||
            is.na(save) ||
            !nzchar(save)
        ) {
            stop("save must be a single file path.")
        }

        extension <- tolower(tools::file_ext(save))

        if (!extension %in% c(
            "pdf", "png", "tiff", "tif", "jpeg", "jpg", "bmp"
        )) {
            stop(
                "save must end in .pdf, .png, .tiff, ",
                ".tif, .jpeg, .jpg or .bmp."
            )
        }
    }

    # Validate expression data and annotation
    if (!all(c("logCPM", "metadata", "genes") %in% names(bulk))) {
        stop("Use prepare_bulk_expression() output.")
    }

    x <- bulk$logCPM

    if (
        !is.matrix(x) ||
        !is.numeric(x) ||
        any(!is.finite(x)) ||
        is.null(rownames(x)) ||
        is.null(colnames(x)) ||
        anyDuplicated(rownames(x)) ||
        anyDuplicated(colnames(x))
    ) {
        stop("Invalid logCPM matrix.")
    }

    if (
        !gene_col %in% names(bulk$genes) ||
        !all(rownames(x) %in% rownames(bulk$genes)) ||
        anyDuplicated(rownames(bulk$genes))
    ) {
        stop("Invalid gene annotation alignment.")
    }

    meta <- as.data.frame(bulk$metadata)

    if (
        anyDuplicated(rownames(meta)) ||
        !all(colnames(x) %in% rownames(meta))
    ) {
        stop("Metadata sample IDs do not match.")
    }

    # Select and order samples
    if (is.null(samples)) {
        samples <- colnames(x)
    }

    if (
        !is.character(samples) ||
        !length(samples) ||
        anyNA(samples) ||
        anyDuplicated(samples) ||
        !all(samples %in% colnames(x))
    ) {
        stop("samples must contain unique, available sample IDs.")
    }

    meta <- meta[samples, , drop = FALSE]

    if (column_mode == "mean" && !length(group_by)) {
        stop("Mean columns require group_by.")
    }

    if (
        !is.list(column_order) ||
        (
            length(column_order) &&
            (
                is.null(names(column_order)) ||
                anyDuplicated(names(column_order)) ||
                any(!nzchar(names(column_order)))
            )
        )
    ) {
        stop("column_order must be a named list of metadata level orders.")
    }

    needed <- unique(c(
        group_by,
        annotation_cols,
        names(column_order)
    ))

    if (!all(needed %in% names(meta))) {
        stop(
            "Missing metadata columns: ",
            paste(setdiff(needed, names(meta)), collapse = ", ")
        )
    }

    for (column in unique(c(group_by, names(column_order)))) {
        if (anyNA(meta[[column]])) {
            stop("Missing grouping/ordering labels in ", column)
        }
    }

    for (column in names(column_order)) {

        levels <- as.character(column_order[[column]])

        if (
            !length(levels) ||
            anyNA(levels) ||
            anyDuplicated(levels) ||
            !all(as.character(meta[[column]]) %in% levels)
        ) {
            stop("Incomplete or invalid column_order for ", column)
        }

        meta[[column]] <- factor(
            as.character(meta[[column]]),
            levels = levels
        )
    }

    order_cols <- unique(c(names(column_order), group_by))

    if (length(order_cols)) {
        ord <- do.call(
            order,
            c(
                unname(as.list(meta[, order_cols, drop = FALSE])),
                list(na.last = TRUE)
            )
        )

        meta <- meta[ord, , drop = FALSE]
    }

    # Match gene symbols
    # For duplicated symbols, retain the gene ID with highest mean logCPM
    symbols <- trimws(
        as.character(bulk$genes[rownames(x), gene_col])
    )

    idx <- which(!is.na(symbols) & symbols %in% genes)

    idx <- idx[
        order(-rowMeans(x)[idx], rownames(x)[idx])
    ]

    idx <- idx[!duplicated(symbols[idx])]

    missing <- setdiff(genes, symbols[idx])

    if (length(missing)) {
        warning(
            "Genes not found: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    if (!length(idx)) {
        stop("None of the requested gene names were found.")
    }

    idx <- idx[
        match(
            genes[genes %in% symbols[idx]],
            symbols[idx]
        )
    ]

    mapping <- data.frame(
        Gene = symbols[idx],
        GeneID = rownames(x)[idx]
    )

    mat <- x[idx, rownames(meta), drop = FALSE]
    rownames(mat) <- symbols[idx]

    sample_groups <- data.frame(
        Sample = rownames(meta),
        Column = rownames(meta)
    )

    # Optionally average samples within metadata groups
    if (column_mode == "mean") {

        codes <- lapply(
            meta[, group_by, drop = FALSE],
            function(z) {
                match(as.character(z), unique(as.character(z)))
            }
        )

        key <- do.call(
            interaction,
            c(codes, list(drop = TRUE, lex.order = TRUE))
        )

        groups <- unique(as.character(key))

        positions <- lapply(
            groups,
            function(z) which(as.character(key) == z)
        )

        mat <- do.call(
            cbind,
            lapply(
                positions,
                function(j) {
                    rowMeans(mat[, j, drop = FALSE])
                }
            )
        )

        rownames(mat) <- mapping$Gene

        labels <- vapply(
            positions,
            function(j) {
                paste(
                    as.character(
                        unlist(meta[j[1], group_by, drop = FALSE])
                    ),
                    collapse = " / "
                )
            },
            character(1)
        )

        labels <- make.unique(labels)
        colnames(mat) <- labels

        for (column in annotation_cols) {

            varies <- vapply(
                positions,
                function(j) {
                    length(unique(meta[[column]][j])) != 1L
                },
                logical(1)
            )

            if (any(varies)) {
                stop(
                    "Annotation ", column,
                    " varies within groups; include it in ",
                    "group_by or remove it."
                )
            }
        }

        sample_groups$Column <- labels[
            match(as.character(key), groups)
        ]

        meta <- meta[
            vapply(positions, `[`, integer(1), 1L),
            ,
            drop = FALSE
        ]

        rownames(meta) <- labels
    }

    expression_matrix <- mat

    # Remove constant genes when scaling/correlation requires variation
    constant <- if (ncol(mat) > 1L) {
        rownames(mat)[apply(mat, 1, stats::sd) == 0]
    } else {
        character()
    }

    if (
        scale == "row" ||
        (cluster_rows && distance == "correlation")
    ) {
        if (length(constant)) {
            warning(
                "Constant genes excluded: ",
                paste(constant, collapse = ", "),
                call. = FALSE
            )

            mat <- mat[
                !rownames(mat) %in% constant,
                ,
                drop = FALSE
            ]
        }
    }

    if (!nrow(mat)) {
        stop("No variable genes remain for this heatmap.")
    }

    # Row z-scores
    if (scale == "row") {

        if (ncol(mat) < 2L) {
            stop("Row scaling requires at least two columns.")
        }

        mat <- t(base::scale(t(mat)))
    }

    # Construct clustering trees
    make_tree <- function(z, enabled) {

        if (!enabled || nrow(z) < 2L) {
            return(FALSE)
        }

        if (distance == "euclidean") {

            d <- stats::dist(z)

        } else {

            if (
                ncol(z) < 2L ||
                any(apply(z, 1, stats::sd) == 0)
            ) {
                warning(
                    "Correlation clustering disabled for profiles ",
                    "with insufficient variation.",
                    call. = FALSE
                )

                return(FALSE)
            }

            rho <- stats::cor(
                t(z),
                method = correlation_method
            )

            dmat <- 1 - rho
            dmat[dmat < 0] <- 0
            diag(dmat) <- 0

            d <- stats::as.dist(dmat)
        }

        if (any(!is.finite(d))) {
            stop("Non-finite clustering distances.")
        }

        stats::hclust(d, method = clustering_method)
    }

    row_tree <- make_tree(mat, cluster_rows)
    col_tree <- make_tree(t(mat), cluster_cols)

    # Color limits
    if (!is.null(value_limit)) {

        if (
            !is.numeric(value_limit) ||
            length(value_limit) != 1L ||
            !is.finite(value_limit) ||
            value_limit <= 0
        ) {
            stop("Invalid value_limit.")
        }

        limits <- c(-value_limit, value_limit)

    } else {

        limits <- range(mat)
    }

    if (diff(limits) == 0) {
        limits <- limits + c(-0.5, 0.5)
    }

    if (is.function(color)) {
        color <- color(
            seq(limits[1], limits[2], length.out = 100)
        )
    }

    if (
        !is.character(color) ||
        length(color) < 2L ||
        anyNA(color)
    ) {
        stop("Supply a color vector or color function.")
    }

    grDevices::col2rgb(color)

    # Clip displayed values only; clustering uses the unclipped matrix
    displayed <- mat
    displayed[] <- pmax(
        limits[1],
        pmin(limits[2], mat)
    )

    # Align row annotations
    if (!is.null(annotation_row)) {

        if (
            !is.data.frame(annotation_row) ||
            anyDuplicated(rownames(annotation_row)) ||
            !all(rownames(mat) %in% rownames(annotation_row))
        ) {
            stop(
                "annotation_row must have gene-name rownames ",
                "covering plotted genes."
            )
        }

        annotation_row <- annotation_row[
            rownames(mat),
            ,
            drop = FALSE
        ]
    }

    anno <- if (length(annotation_cols)) {
        meta[, annotation_cols, drop = FALSE]
    } else {
        NA
    }

    if (is.null(height)) {
        height <- max(4, 2 + 0.18 * nrow(mat))
    }

    # An interest list overrides show_rownames: only matching genes are labelled.
    interest_provided <- !is.null(genes_of_interest)
    row_labels_to_show <- rownames(mat)
    if (interest_provided) {
        row_labels_to_show[!rownames(mat) %in% genes_of_interest] <- ""
    }

    # Build without drawing/exporting, so the legend can be adjusted first
    p <- pheatmap::pheatmap(
        displayed,
        scale = "none",
        color = color,
        breaks = seq(
            limits[1],
            limits[2],
            length.out = length(color) + 1L
        ),
        cluster_rows = row_tree,
        cluster_cols = col_tree,
        treeheight_row = if (show_row_dendrogram) 50 else 0,
        treeheight_col = if (show_col_dendrogram) 50 else 0,
        show_rownames = if (interest_provided) TRUE else show_rownames,
        labels_row = row_labels_to_show,
        annotation_col = anno,
        annotation_row = if (is.null(annotation_row)) {
            NA
        } else {
            annotation_row
        },
        annotation_colors = if (is.null(annotation_colors)) {
            NA
        } else {
            annotation_colors
        },
        border_color = NA,
        show_colnames = show_colnames,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        main = if (is.null(title)) NA else title,
        filename = NA,
        silent = TRUE
    )

    # Colour the labels in their FINAL displayed order, after clustering.
    # Editing the gtable ensures the colours are preserved in saved files.
    gt <- p$gtable
    if (interest_provided) {
        row_label_idx <- which(gt$layout$name == "row_names")
        if (length(row_label_idx) != 1L) {
            stop("Could not locate the heatmap row labels.")
        }
        row_labels <- gt$grobs[[row_label_idx]]
        row_labels$gp$col <- ifelse(
            as.character(row_labels$label) %in% genes_of_interest,
            "red",
            "black"
        )
        gt$grobs[[row_label_idx]] <- row_labels
    }

    # Reserve additional width for the annotation legend

    legend_col <- unique(
        gt$layout$l[gt$layout$name == "annotation_legend"]
    )

    matrix_col <- unique(
        gt$layout$l[gt$layout$name == "matrix"]
    )

    if (
        legend_extra_width_mm > 0 &&
        length(legend_col) == 1L &&
        length(matrix_col) == 1L &&
        legend_col != matrix_col
    ) {

        extra_space <- grid::unit(
            legend_extra_width_mm,
            "mm"
        )

        gt$widths[legend_col] <-
            gt$widths[legend_col] + extra_space

        gt$widths[matrix_col] <-
            gt$widths[matrix_col] - extra_space
    }

    p$gtable <- gt

    # Export the adjusted layout
    if (!is.null(save)) {

        export_plot <- function() {

            switch(
                extension,
                pdf = grDevices::pdf(
                    file = save,
                    width = width,
                    height = height,
                    useDingbats = FALSE
                ),
                png = grDevices::png(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300
                ),
                tiff =,
                tif = grDevices::tiff(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300,
                    compression = "lzw"
                ),
                jpeg =,
                jpg = grDevices::jpeg(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300,
                    quality = 100
                ),
                bmp = grDevices::bmp(
                    filename = save,
                    width = width,
                    height = height,
                    units = "in",
                    res = 300
                )
            )

            on.exit(grDevices::dev.off(), add = TRUE)

            grid::grid.newpage()
            grid::grid.draw(gt)
        }

        export_plot()

    } else {

        grid::grid.newpage()
        grid::grid.draw(gt)
    }

    # Return only the plot and the requested clustering results.
    # Trees are hclust objects: use cutree() to choose the number of groups.
    # NULL means clustering was requested but could not be computed.
    out <- list(plot = p)
    if (cluster_rows) {
        out["row_tree"] <- list(
            if (inherits(row_tree, "hclust")) row_tree else NULL
        )
    }
    if (cluster_cols) {
        out["column_tree"] <- list(
            if (inherits(col_tree, "hclust")) col_tree else NULL
        )
    }
    invisible(out)
}

# =============================================================================
# Prepare ligand–receptor evidence for plotting
# =============================================================================

#' Prepare NicheNet ligand-receptor evidence for plotting
#'
#' Join ranked ligand activities, sender availability, and receiver receptor
#' evidence for one receiver and response direction across selected senders and
#' conditions. Only links involving the top-ranked ligands are retained.
#'
#' @param result Output returned by [run_bulk_nichenet()].
#' @param receiver Single receiver population to retain.
#' @param response Receiver response direction, either `"up"` or `"down"`.
#' @param senders `NULL` for all available senders or a character vector of
#'   sender populations.
#' @param conditions Character vector of experimental conditions.
#' @param top_n Positive integer number of ranked ligands retained per context.
#'
#' @return A long data frame containing ligand-receptor pairs, context keys,
#'   activity rank and corrected AUPR, sender/receiver expression availability,
#'   log-fold changes, adjusted p-values, and a printable `Pair` label.
#' @export
prepare_bulk_nichenet_lr_plot_data <- function(
    result,
    receiver,
    response = c("up", "down"),
    senders = NULL,
    conditions = c("2W", "1M", "2M"),
    top_n = 10
) {

    response <- match.arg(response)

    if (!is.character(receiver) || length(receiver) != 1L ||
        is.na(receiver) || !nzchar(receiver)) {
        stop("Specify one ASC receiver label.")
    }

    if (!is.numeric(top_n) || length(top_n) != 1L ||
        !is.finite(top_n) || top_n < 1 || top_n != floor(top_n)) {
        stop("top_n must be a positive integer.")
    }

    if (!length(conditions) || anyNA(conditions) ||
        anyDuplicated(conditions)) {
        stop("Provide unique condition labels.")
    }

    keys <- c(
        "Ligand", "Receiver", "Condition",
        "Reference", "Response", "ResultName"
    )

    sender_columns <- c(
        keys, "Sender", "Rank", "aupr_corrected",
        "Ligand_mean_logCPM_condition",
        "Ligand_expressed_reference",
        "Ligand_expressed_condition",
        "Ligand_logFC", "Ligand_FDR"
    )

    receptor_columns <- c(
        keys, "Receptor",
        "Receptor_mean_logCPM_condition",
        "Receptor_expressed_reference",
        "Receptor_expressed_condition",
        "Receptor_logFC", "Receptor_FDR"
    )

    se <- result$sender_evidence
    re <- result$receptor_evidence

    if (!is.data.frame(se) ||
        !all(sender_columns %in% names(se))) {
        stop("Missing fields in result$sender_evidence.")
    }

    if (!is.data.frame(re) ||
        !all(receptor_columns %in% names(re))) {
        stop("Missing fields in result$receptor_evidence.")
    }

    if (is.null(senders)) {
        senders <- unique(as.character(se$Sender))
    }

    if (!length(senders) || anyNA(senders) || anyDuplicated(senders)) {
        stop("Provide unique sender labels.")
    }

    se <- se[
        se$Receiver == receiver &
            se$Response == response &
            se$Condition %in% conditions &
            se$Sender %in% senders &
            !is.na(se$Rank) &
            se$Rank <= top_n,
        sender_columns,
        drop = FALSE
    ]

    re <- re[
        re$Receiver == receiver &
            re$Response == response &
            re$Condition %in% conditions,
        receptor_columns,
        drop = FALSE
    ]

    if (!nrow(se) || !nrow(re)) {
        stop(
            "No matching evidence. Check receiver, response, ",
            "sender labels, conditions, and top_n."
        )
    }

    if (anyDuplicated(se[, c(keys, "Sender"), drop = FALSE])) {
        stop("Duplicate sender evidence identifiers.")
    }

    if (anyDuplicated(re[, c(keys, "Receptor"), drop = FALSE])) {
        stop("Duplicate receptor evidence identifiers.")
    }

    # Each ligand can connect to several receptors and several ATM senders.
    out <- merge(
        se,
        re,
        by = keys,
        all = FALSE,
        sort = FALSE
    )

    if (!nrow(out)) {
        stop("No ligand–receptor links matched the selected sender evidence.")
    }

    numeric_columns <- c(
        "Rank", "aupr_corrected",
        "Ligand_mean_logCPM_condition",
        "Receptor_mean_logCPM_condition",
        "Ligand_logFC", "Receptor_logFC",
        "Ligand_FDR", "Receptor_FDR"
    )

    if (!all(vapply(out[numeric_columns], is.numeric, logical(1)))) {
        stop("Scores, expression, logFC, and FDR fields must be numeric.")
    }

    out$Pair <- paste(out$Ligand, out$Receptor, sep = " → ")

    out
}


# =============================================================================
# Dot plot: separate ligand and receptor expression evidence over HFD
# =============================================================================

#' Plot NicheNet ligand-receptor expression evidence as paired dot panels
#'
#' Display sender ligand and receiver receptor expression evidence for the
#' top-ranked ligand-receptor pairs, with dot size representing expression and
#' colour representing differential-expression effect.
#'
#' @param result Output returned by [run_bulk_nichenet()].
#' @param receiver Single receiver population to plot.
#' @param response Receiver response direction, either `"up"` or `"down"`.
#' @param sender_order,condition_order Explicit sender and condition orders.
#' @param top_n Positive integer number of ranked ligands retained per context.
#' @param fdr_cutoff Numeric FDR threshold used to distinguish supported DE
#'   evidence in the plotting data.
#' @param logfc_limit `NULL` or a positive symmetric colour-scale limit.
#' @param dot_size_range Numeric length-two vector giving minimum and maximum
#'   dot sizes.
#' @param color Three-colour vector for negative, neutral, and positive logFC.
#' @param fontsize Base ggplot2 font size.
#' @param width,height Saved figure dimensions in inches; `NULL` height adapts
#'   to the number of pairs and senders.
#' @param save `NULL` or a figure output path.
#' @param draw Logical scalar. Whether to print the plot on the active device.
#'
#' @return Invisibly returns a list containing the ggplot object, long plotting
#'   data, joined ligand-receptor evidence, and displayed pair order.
#' @export
plot_bulk_nichenet_lr_dotplot <- function(
    result,
    receiver,
    response = c("up", "down"),
    sender_order = c("pATM", "cATM", "sATM"),
    condition_order = c("2W", "1M", "2M"),
    top_n = 10,

    fdr_cutoff = 0.05,
    logfc_limit = NULL,
    dot_size_range = c(1.5, 6),

    color = BULK_SIGNALING_DIVERGING_COLORS,

    fontsize = 10,
    width = 11,
    height = NULL,
    save = NULL,
    draw = TRUE
) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package ggplot2 is required.")
    }

    response <- match.arg(response)

    if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1) {
        stop("Invalid fdr_cutoff.")
    }

    if (length(dot_size_range) != 2L ||
        any(!is.finite(dot_size_range)) ||
        any(dot_size_range <= 0) ||
        dot_size_range[1] > dot_size_range[2]) {
        stop("Invalid dot_size_range.")
    }

    if (!is.character(color) || length(color) < 2L) {
        stop("Provide a vector of heatmap colors.")
    }

    invisible(grDevices::col2rgb(color))

    lr <- prepare_bulk_nichenet_lr_plot_data(
        result = result,
        receiver = receiver,
        response = response,
        senders = sender_order,
        conditions = condition_order,
        top_n = top_n
    )

    # Order pairs by their highest ligand activity across selected conditions.
    pair_score <- tapply(
        lr$aupr_corrected,
        lr$Pair,
        function(z) {
            z <- z[is.finite(z)]
            if (length(z)) max(z) else -Inf
        }
    )

    pair_order <- names(pair_score)[
        order(-pair_score, names(pair_score))
    ]

    # Two evidence panels, retaining the same ligand–receptor row identifiers.
    ligand_data <- lr
    ligand_data$Evidence <- "ATM ligand"
    ligand_data$LogFC <- ligand_data$Ligand_logFC
    ligand_data$FDR <- ligand_data$Ligand_FDR
    ligand_data$MeanLogCPM <- ligand_data$Ligand_mean_logCPM_condition

    receptor_data <- lr
    receptor_data$Evidence <- "ASC receptor"
    receptor_data$LogFC <- receptor_data$Receptor_logFC
    receptor_data$FDR <- receptor_data$Receptor_FDR
    receptor_data$MeanLogCPM <- receptor_data$Receptor_mean_logCPM_condition

    plot_data <- rbind(ligand_data, receptor_data)

    plot_data$Pair <- factor(
        plot_data$Pair,
        levels = rev(pair_order)
    )

    plot_data$Condition <- factor(
        plot_data$Condition,
        levels = condition_order
    )

    plot_data$Sender <- factor(
        plot_data$Sender,
        levels = sender_order
    )

    plot_data$Evidence <- factor(
        plot_data$Evidence,
        levels = c("ATM ligand", "ASC receptor")
    )

    plot_data$Significant <- factor(
        ifelse(
            !is.na(plot_data$FDR) &
                plot_data$FDR < fdr_cutoff,
            "Yes",
            "No / not tested"
        ),
        levels = c("No / not tested", "Yes")
    )

    plot_data$LogFC[!is.finite(plot_data$LogFC)] <- NA_real_

    if (is.null(logfc_limit)) {
        finite_fc <- plot_data$LogFC[is.finite(plot_data$LogFC)]
        logfc_limit <- if (length(finite_fc)) max(abs(finite_fc)) else 1
        if (logfc_limit == 0) logfc_limit <- 1
    }

    if (!is.numeric(logfc_limit) || length(logfc_limit) != 1L ||
        !is.finite(logfc_limit) || logfc_limit <= 0) {
        stop("logfc_limit must be positive.")
    }

    # Clip color values only; original logFC remains in returned data.
    plot_data$PlotLogFC <- pmax(
        -logfc_limit,
        pmin(logfc_limit, plot_data$LogFC)
    )

    p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = Condition,
            y = Pair
        )
    ) +
        ggplot2::geom_point(
            ggplot2::aes(
                fill = PlotLogFC,
                size = MeanLogCPM,
                colour = Significant
            ),
            shape = 21,
            stroke = 0.65
        ) +
        ggplot2::facet_grid(
            Sender ~ Evidence,
            drop = FALSE
        ) +
        ggplot2::scale_fill_gradientn(
            colours = color,
            limits = c(-logfc_limit, logfc_limit),
            na.value = "#D9D9D9",
            name = "Gene logFC"
        ) +
        ggplot2::scale_size_continuous(
            range = dot_size_range,
            name = "Mean log2 CPM\n(HFD condition)"
        ) +
        ggplot2::scale_colour_manual(
            values = c(
                "No / not tested" = "#B8B8B8",
                "Yes" = "#000000"
            ),
            name = paste0("Gene FDR < ", fdr_cutoff)
        ) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::labs(
            x = NULL,
            y = "Ligand → receptor",
            title = paste0(
                receiver, ": ", response, " response"
            ),
            subtitle = paste0(
                "Top ", top_n,
                " ligands per condition by NicheNet rank; ",
                "color shows gene expression change"
            )
        ) +
        ggplot2::theme_classic(base_size = fontsize) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA,
                linewidth = 0.4
            ),
            strip.background = ggplot2::element_rect(
                fill = "#F3F3F3",
                colour = "black"
            ),
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            axis.text.y = ggplot2::element_text(
                colour = "black"
            ),
            legend.position = "right"
        )

    if (is.null(height)) {
        height <- max(
            6,
            2 + length(sender_order) * (0.18 * length(pair_order) + 0.7)
        )
    }

    if (!is.null(save)) {

        dir.create(
            dirname(save),
            recursive = TRUE,
            showWarnings = FALSE
        )

        ggplot2::ggsave(
            filename = save,
            plot = p,
            width = width,
            height = height,
            units = "in",
            dpi = 300,
            limitsize = FALSE
        )
    }

    if (draw) print(p)

    invisible(
        list(
            plot = p,
            data = plot_data,
            lr_evidence = lr,
            pair_order = pair_order
        )
    )
}


# =============================================================================
# Connectivity heatmap for one ASC subtype and one HFD condition
# =============================================================================

#' Plot NicheNet ligand-receptor prior connectivity
#'
#' Build a binary ligand-by-receptor matrix for one receiver, condition, and
#' response direction, with row annotations showing which sender populations
#' meet the requested ligand-availability rule.
#'
#' @param result Output returned by [run_bulk_nichenet()].
#' @param receiver Single receiver population to display.
#' @param condition Single experimental condition to display.
#' @param response Receiver response direction, either `"up"` or `"down"`.
#' @param sender_order Explicit sender-population order.
#' @param sender_colors Named colour vector covering `sender_order`.
#' @param top_n Positive integer number of ranked ligands retained.
#' @param availability `"condition"` to require condition-specific sender
#'   expression or `"either"` to allow reference or condition expression.
#' @param color Two-colour vector for absent and present prior links.
#' @param fontsize_row,fontsize_col Heatmap row and column label sizes.
#' @param width,height Saved heatmap dimensions in inches; `NULL` height adapts
#'   to the number of ligand rows.
#' @param save `NULL` or a filename written by `pheatmap`.
#' @param draw Logical scalar. Whether to draw the heatmap on the active device.
#'
#' @return Invisibly returns a list containing the `pheatmap` object, binary
#'   connectivity matrix, sender-availability row annotation, joined evidence,
#'   and the selected availability rule.
#' @export
plot_bulk_nichenet_lr_connectivity <- function(
    result,
    receiver,
    condition,
    response = c("up", "down"),

    sender_order = c("pATM", "cATM", "sATM"),
    sender_colors = BULK_NICHENET_SENDER_COLORS,

    top_n = 10,
    availability = c("condition", "either"),

    color = BULK_NICHENET_CONNECTIVITY_COLORS,
    fontsize_row = 10,
    fontsize_col = 10,
    width = 10,
    height = NULL,
    save = NULL,
    draw = TRUE
) {

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package pheatmap is required.")
    }

    response <- match.arg(response)
    availability <- match.arg(availability)

    if (!is.character(condition) || length(condition) != 1L ||
        is.na(condition) || !nzchar(condition)) {
        stop("Specify one condition.")
    }

    sender_colors <- unlist(sender_colors)

    if (anyDuplicated(names(sender_colors)) ||
        !all(sender_order %in% names(sender_colors))) {
        stop("sender_colors must have unique names covering sender_order.")
    }

    if (!is.character(color) || length(color) != 2L || anyNA(color)) {
        stop("Provide two colors: no link, prior link.")
    }

    invisible(grDevices::col2rgb(c(color, sender_colors)))

    lr <- prepare_bulk_nichenet_lr_plot_data(
        result = result,
        receiver = receiver,
        response = response,
        senders = sender_order,
        conditions = condition,
        top_n = top_n
    )

    as_flag <- function(z) {
        z <- as.logical(z)
        !is.na(z) & z
    }

    if (availability == "condition") {

        ligand_available <- as_flag(
            lr$Ligand_expressed_condition
        )

        receptor_available <- as_flag(
            lr$Receptor_expressed_condition
        )

    } else {

        ligand_available <- (
            as_flag(lr$Ligand_expressed_condition) |
                as_flag(lr$Ligand_expressed_reference)
        )

        receptor_available <- (
            as_flag(lr$Receptor_expressed_condition) |
                as_flag(lr$Receptor_expressed_reference)
        )
    }

    lr <- lr[
        ligand_available & receptor_available,
        ,
        drop = FALSE
    ]

    if (!nrow(lr)) {
        stop(
            "No links satisfy the selected availability rule. ",
            "Try availability = 'either' or increase top_n."
        )
    }

    ligand_rank <- tapply(
        lr$Rank,
        lr$Ligand,
        min
    )

    ligands <- names(ligand_rank)[
        order(ligand_rank, names(ligand_rank))
    ]

    receptors <- sort(unique(as.character(lr$Receptor)))

    mat <- matrix(
        0,
        nrow = length(ligands),
        ncol = length(receptors),
        dimnames = list(ligands, receptors)
    )

    pairs <- unique(lr[, c("Ligand", "Receptor")])

    mat[
        cbind(
            match(pairs$Ligand, ligands),
            match(pairs$Receptor, receptors)
        )
    ] <- 1

    # Ligand activity is one value per ligand/receiver response,
    # shared across the possible ATM sources.
    activity <- vapply(
        ligands,
        function(ligand) {
            values <- unique(
                lr$aupr_corrected[lr$Ligand == ligand]
            )

            if (length(values) != 1L) {
                stop("Inconsistent ligand activity values for ", ligand)
            }

            values
        },
        numeric(1)
    )

    ann_row <- data.frame(
        Activity = activity,
        row.names = ligands
    )

    ann_colors <- list()

    for (sender in sender_order) {

        available_ligands <- unique(
            lr$Ligand[lr$Sender == sender]
        )

        annotation_name <- paste0("Source_", sender)

        ann_row[[annotation_name]] <- factor(
            ifelse(
                ligands %in% available_ligands,
                "Yes",
                "No"
            ),
            levels = c("No", "Yes")
        )

        ann_colors[[annotation_name]] <- c(
            No = "#F0F0F0",
            Yes = unname(sender_colors[sender])
        )
    }

    if (is.null(height)) {
        height <- max(4, 2 + 0.25 * length(ligands))
    }

    args <- list(
        mat = mat,
        color = color,
        breaks = c(-0.5, 0.5, 1.5),

        cluster_rows = FALSE,
        cluster_cols = FALSE,

        annotation_row = ann_row,
        annotation_colors = ann_colors,
        drop_levels = FALSE,

        legend_breaks = c(0, 1),
        legend_labels = c(
            "No prior link",
            "Prior link"
        ),

        border_color = "#FFFFFF",
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        angle_col = "45",

        main = paste0(
            receiver, " / ", condition,
            "\nLigand–receptor connectivity: ", response, " response"
        ),

        silent = TRUE
    )

    if (!is.null(save)) {

        dir.create(
            dirname(save),
            recursive = TRUE,
            showWarnings = FALSE
        )

        args$filename <- save
        args$width <- width
        args$height <- height
    }

    p <- do.call(
        pheatmap::pheatmap,
        args
    )

    if (draw) {
        grid::grid.newpage()
        grid::grid.draw(p$gtable)
    }

    invisible(
        list(
            plot = p,
            matrix = mat,
            annotation_row = ann_row,
            lr_evidence = lr,
            availability = availability
        )
    )
}
