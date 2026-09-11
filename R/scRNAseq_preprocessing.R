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
#'     \item{original_sample_names}{Sample column names as read from the
#'       featureCounts count file.}
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
    colnames(counts) <- clean_sample_names(
        original_sample_names,
        "featureCounts"
    )

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
        original_sample_names = original_sample_names,
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


# =============================================================================
# Normalization and feature selection
# =============================================================================


# =============================================================================
# Dimensionality reduction and clustering
# =============================================================================


# =============================================================================
# Visualization and export
# =============================================================================
