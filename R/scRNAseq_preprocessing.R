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
#' Calculate TMM-normalized log2 counts per million for selected genes and
#' display their sample-level distributions as faceted violin plots. Optional
#' boxplots and jittered sample points retain distribution summaries and
#' individual biological-replicate context.
#'
#' @param dat A list-like featureCounts project containing `counts`, `metadata`,
#'   and `genes`. Count-matrix columns must be identical to metadata row names,
#'   and count-matrix rows must be identical to gene-annotation row names.
#' @param genes Non-empty character vector of gene symbols or stable gene IDs.
#'   Requested order is preserved. Missing genes produce warnings; when one
#'   request matches multiple annotation rows, the row with the highest mean
#'   log2 CPM is used.
#' @param group_by Character scalar naming the metadata column displayed on the
#'   x-axis and used for violin colours.
#' @param split_by `NULL` or a character scalar naming a metadata column used
#'   for facet columns. Genes form facet rows when this is supplied.
#' @param subset Optional unquoted logical expression evaluated within
#'   `dat$metadata`, with access to the calling environment. Missing results
#'   are treated as `FALSE`.
#' @param gene_col Character scalar naming the preferred gene-symbol column in
#'   `dat$genes`. Stable IDs are used when this column is absent or blank.
#' @param gene_id_col Character scalar naming the preferred stable-ID column in
#'   `dat$genes`. Annotation row names are used when this column is absent.
#' @param group_order `NULL` or a unique character vector defining x-axis group
#'   order. Observed groups omitted from the vector are appended in first-seen
#'   order; unobserved supplied values are dropped from the plotted levels.
#' @param split_order `NULL` or a unique character vector defining split-panel
#'   order. Observed values omitted from the vector are appended.
#' @param colors `NULL` or a character colour vector. Unnamed colours are
#'   assigned in plotted group order. Named colours must cover every observed
#'   group. The default uses `BULK_VIOLIN_MACARON_COLORS` for up to 12 groups
#'   and a pastel HCL palette for larger sets.
#' @param prior_count Non-negative numeric scalar passed to [edgeR::cpm()] when
#'   calculating log2 CPM values.
#' @param violin_width Positive numeric scalar controlling violin width.
#' @param violin_alpha Numeric scalar in `[0, 1]` controlling violin opacity.
#' @param violin_linewidth Non-negative numeric scalar controlling violin
#'   outline width.
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
#'     \item{logCPM}{The complete selected-sample TMM-normalized log2-CPM
#'       matrix before gene selection.}
#'     \item{dge}{The normalized [edgeR::DGEList()] object.}
#'     \item{colors}{The named group-colour mapping used by the plot.}
#'   }
#'
#' @details
#' This function expects raw sample-level bulk RNA-seq counts or
#' replicate-aware pseudobulk counts. TMM normalization is recalculated after
#' optional sample subsetting. The plot is descriptive and does not replace a
#' design-aware differential-expression model, effect estimates, or biological
#' replication checks. Required packages are `edgeR` and `ggplot2`.
#'
#' @examples
#' \dontrun{
#' violin_result <- plot_bulk_violin(
#'     project,
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
    dat,
    genes,
    group_by = "Sorting",
    split_by = NULL,
    subset = NULL,
    gene_col = "gene_name",
    gene_id_col = "gene_id",
    group_order = NULL,
    split_order = NULL,
    colors = NULL,
    prior_count = 2,
    violin_width = 0.85,
    violin_alpha = 0.65,
    violin_linewidth = 0.45,
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
    required_packages <- c("edgeR", "ggplot2")
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

    required_objects <- c("counts", "metadata", "genes")
    if (!is.list(dat)) {
        stop("`dat` must be a list-like featureCounts project.", call. = FALSE)
    }
    missing_objects <- setdiff(required_objects, names(dat))
    if (length(missing_objects) > 0L) {
        stop(
            "`dat` is missing: ", paste(missing_objects, collapse = ", "),
            call. = FALSE
        )
    }

    counts <- dat$counts
    meta <- dat$metadata
    anno <- dat$genes
    if ((!is.matrix(counts) && !inherits(counts, "Matrix")) ||
        !is.numeric(counts)) {
        stop("`dat$counts` must be a numeric matrix-like object.", call. = FALSE)
    }
    if (!is.data.frame(meta) || !is.data.frame(anno)) {
        stop("`dat$metadata` and `dat$genes` must be data frames.", call. = FALSE)
    }
    if (is.null(colnames(counts)) || is.null(rownames(counts)) ||
        is.null(rownames(meta)) || is.null(rownames(anno))) {
        stop("Counts, metadata, and gene annotations require row/column names.",
             call. = FALSE)
    }
    if (anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts)) ||
        anyDuplicated(rownames(meta)) || anyDuplicated(rownames(anno))) {
        stop("Sample and gene row identifiers must be unique.", call. = FALSE)
    }
    if (!identical(colnames(counts), rownames(meta))) {
        stop(
            "colnames(dat$counts) and rownames(dat$metadata) must be identical ",
            "and in the same order.", call. = FALSE
        )
    }
    if (!identical(rownames(counts), rownames(anno))) {
        stop(
            "rownames(dat$counts) and rownames(dat$genes) must be identical ",
            "and in the same order.", call. = FALSE
        )
    }
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
        stop("`dat$counts` must contain finite, non-negative values.",
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
        stop("`group_by = \"", group_by, "\"` was not found in dat$metadata.",
             call. = FALSE)
    }
    if (!is.null(split_by) && !split_by %in% colnames(meta)) {
        stop("`split_by = \"", split_by, "\"` was not found in dat$metadata.",
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
        prior_count = prior_count,
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
    if (prior_count < 0 || violin_width <= 0 || violin_linewidth < 0 ||
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

    subset_expr <- substitute(subset)
    if (!identical(subset_expr, quote(NULL))) {
        keep_samples <- eval(subset_expr, envir = meta, enclos = parent.frame())
        if (!is.logical(keep_samples) || length(keep_samples) != nrow(meta)) {
            stop(
                "`subset` must evaluate to one TRUE/FALSE value for every ",
                "metadata row.", call. = FALSE
            )
        }
        keep_samples[is.na(keep_samples)] <- FALSE
        if (!any(keep_samples)) {
            stop("`subset` retained zero samples.", call. = FALSE)
        }
        meta <- meta[keep_samples, , drop = FALSE]
        counts <- counts[, rownames(meta), drop = FALSE]
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

    library_sizes <- colSums(counts)
    if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
        stop("At least one selected sample has library size <= 0.", call. = FALSE)
    }

    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    logcpm <- edgeR::cpm(
        dge,
        log = TRUE,
        prior.count = prior_count,
        normalized.lib.sizes = TRUE
    )

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
            linewidth = violin_linewidth
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
        message("TMM normalization + log2 CPM")
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
        colors = colors
    )
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
