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
