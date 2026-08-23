# https://github.com/tidyverse/dtplyr/issues/426
.datatable.aware <- TRUE

#' CellAlleleCounts class
#'
#' S4 class storing cell-by-locus allele counts as a single sparse matrix.
#' Behaves like a named list of per-locus matrices (alleles x cells).
#'
#' @slot counts RsparseMatrix with features (locus-allele pairs) as rows and
#'     cells as columns. CSR format makes per-locus row selection efficient.
#' @slot locus_map DataFrame with columns: locus, allele, feature_name.
#' @slot cell_data DataFrame with cell metadata. Rownames are cell identifiers.
#' @slot locus_index Named list mapping each locus name to its row indices in
#'     `counts`. Precomputed at construction for O(1) locus lookup.
#' @param x A CellAlleleCounts object.
#' @param i Index for loci: numeric, character (locus name), or logical.
#' @param j Index for cells: numeric, character (cell name), or logical.
#' @param ... Not used.
#' @param drop Not used.
#' @param object A CellAlleleCounts object.
#' @param X A CellAlleleCounts object (for lapply).
#' @param FUN Function to apply to each locus matrix.
#' @param value A DataFrame of cell metadata to assign.
#' @param do.NULL Not used.
#' @param prefix Not used.
#'
#' @include getFeatureAnno.R
#' @name CellAlleleCounts-class
#' @rdname CellAlleleCounts-class
#' @exportClass CellAlleleCounts
#' @importFrom methods setClass new validObject is as
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom Matrix RsparseMatrix
#' @importFrom Matrix t
setClass("CellAlleleCounts",
    slots = list(
        counts      = "RsparseMatrix",
        locus_map   = "DataFrame",
        cell_data   = "DataFrame",
        locus_index = "list",
        locus_1hot = "CsparseMatrix"
    ),
    validity = function(object) {
        errors <- character()
        if (nrow(object@counts) != nrow(object@locus_map)) {
            errors <- c(errors,
                "nrow(counts) must equal nrow(locus_map)")
        }
        if (ncol(object@counts) != nrow(object@cell_data)) {
            errors <- c(errors,
                "ncol(counts) must equal nrow(cell_data)")
        }
        required_cols <- c("locus", "allele", "feature_name")
        missing <- setdiff(required_cols, colnames(object@locus_map))
        if (length(missing) > 0) {
            errors <- c(errors, sprintf(
                "locus_map missing columns: %s",
                paste(missing, collapse = ", ")))
        }
        if (length(errors) == 0) TRUE else errors
    }
)

.build_locus_index <- function(locus_map) {
    locus_chr <- as.character(locus_map$locus)
    idx <- split(seq_along(locus_chr), locus_chr)
    # preserve first-occurrence order rather than alphabetical
    idx[unique(locus_chr)]
}

#' @importFrom Matrix fac2sparse
.build_locus_1hot <- function(locus_map) {
    mat <- Matrix::fac2sparse(locus_map$locus)
    colnames(mat) <- locus_map$feature_name
    mat
}

.new_cac <- function(counts, locus_map, cell_data) {
    new("CellAlleleCounts",
        counts      = counts,
        locus_map   = locus_map,
        cell_data   = cell_data,
        locus_index = .build_locus_index(locus_map),
        locus_1hot = .build_locus_1hot(locus_map))
}

#' Construct a CellAlleleCounts object
#'
#' @param counts A matrix (dense or sparse) of allele counts with features
#'     (locus-allele pairs) as rows and cells as columns. Feature names should
#'     follow the format `chr{chr}_{position}_{allele}`.
#' @param locus_map A DataFrame with columns `locus`, `allele`, and
#'     `feature_name`, one row per feature (row of `counts`). If NULL, it is
#'     derived from `rownames(counts)` using the standard feature name format.
#' @param cell_data A DataFrame of cell metadata whose rownames are cell
#'     identifiers matching `colnames(counts)`. If NULL, an empty DataFrame
#'     is created from `colnames(counts)`.
#'
#' @returns A `CellAlleleCounts` object.
#'
#' @export
#' @importFrom S4Vectors DataFrame
CellAlleleCounts <- function(counts, locus_map = NULL, cell_data = NULL) {
    if (!is(counts, "RsparseMatrix")) {
        counts <- as(as(counts, "generalMatrix"), "RsparseMatrix")
    }

    if (is.null(locus_map)) {
        feature_names <- rownames(counts)
        if (is.null(feature_names)) {
            stop("counts must have rownames when locus_map is not provided")
        }
        parsed <- .getFeatureAnno(feature_names)
        locus_map <- DataFrame(
            locus = parsed$chr_pos,
            allele = parsed$allele,
            feature_name = feature_names
        )
    }

    if (is.null(cell_data)) {
        cell_names <- colnames(counts)
        if (is.null(cell_names)) {
            stop("counts must have colnames when cell_data is not provided")
        }
        cell_data <- DataFrame(row.names = cell_names)
    }

    .new_cac(counts, locus_map, cell_data)
}

#' Convert a named list of count matrices to a CellAlleleCounts object
#'
#' @param mat_list A list of matrices (alleles x cells), one per
#'     locus.  The list names should be the locus names. All matrices
#'     must have the same column names (cells).
#' @param cell_data Optional DataFrame of cell metadata. Rownames must
#'     match column names of the matrices. If NULL, an empty DataFrame
#'     is created.
#'
#' @returns A `CellAlleleCounts` object.
#'
#' @export
#' @importFrom Matrix t
#' @importFrom S4Vectors DataFrame
CellAlleleCounts_from_list <- function(mat_list, cell_data = NULL) {
    cell_names <- colnames(mat_list[[1]])

    if (is.null(cell_data)) {
        cell_data <- DataFrame(row.names = cell_names)
    }

    locus_names <- names(mat_list)
    feature_names <- unlist(lapply(mat_list, rownames), use.names = FALSE)
    locus_vec <- rep(locus_names,
                     times = vapply(mat_list, nrow, integer(1)))

    parsed <- .getFeatureAnno(feature_names)

    locus_map <- DataFrame(
        locus = locus_vec,
        allele = parsed$allele,
        feature_name = feature_names
    )

    stacked <- do.call(rbind, mat_list)
    counts <- as(stacked, "RsparseMatrix")

    .new_cac(counts, locus_map, cell_data)
}

#' Convert a dataframe produced by sclantern-nf into CellAlleleCounts
#'
#' @param df_counts A long dataframe or data.table, as produced by
#'     sclantern-nf. The table should have columns "cell", "umi",
#'     "allele", "locus", "count", corresponding respectively to the
#'     cell barcode, UMI barcode, allele, locus, reads.
#' @param cell_data Optional DataFrame of cell metadata. Rownames must
#'     match column names of the matrices. If NULL, an empty DataFrame
#'     is created.
#'
#' @returns A `CellAlleleCounts` object.
#'
#' @export
#' @importFrom data.table as.data.table
#' @importFrom dtplyr lazy_dt
#' @importFrom dplyr group_by mutate ungroup filter summarize
#' @importFrom Matrix sparseMatrix
CellAlleleCounts_from_dataframe <- function(df_counts,
                                            min_cells_per_feature = 2) {
    df_counts <- as.data.table(df_counts)

    lazy_dt(df_counts) |>
        group_by(cell, umi, locus) |>
        mutate(frac=count / sum(count)) |>
        ungroup() |>
        filter(frac > .5) |>
        group_by(cell, locus, allele) |>
        summarize(umis=n(), .groups='drop') |>
        mutate(feature=paste0(locus, "_", allele)) |>
        mutate(cell=factor(cell, levels=sort(unique(cell))),
               feature=factor(feature, levels=sort(unique(feature)))) |>
        as.data.frame() ->
        df_counts_summed


    cnt_mat <- with(
        df_counts_summed,
        sparseMatrix(
            i = as.integer(feature),
            j = as.integer(cell),
            x = umis,
            dimnames = list(levels(feature), levels(cell)),
            index1 = TRUE
        )
    )

    ret <- CellAlleleCounts(cnt_mat)

    if (min_cells_per_feature > 0) {
        ret <- filter_cell_allele_counts(
            ret,
            min_cells_per_feature = min_cells_per_feature
        )
    }

    ret
}
