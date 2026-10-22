#' Filter cells and allele features by count support
#'
#' @param cac A `CellAlleleCounts` object.
#' @param min_counts_per_cell Minimum total count required to retain a cell.
#' @param min_cells_per_feature Minimum number of retained cells with nonzero
#'     counts required to retain a feature.
#' @param keep_ref Logical; retain features labelled `REF` in
#'     `locus_map$allele_type` regardless of `min_cells_per_feature`.
#'     Default TRUE
#'
#' @details Cells are filtered first, then features. Loci with only one
#'     retained allele are removed.
#'
#' @returns A `CellAlleleCounts` object with retained counts and metadata.
#'
#' @export
filter_cell_allele_counts <- function(cac,
                                      min_counts_per_cell=0,
                                      min_cells_per_feature=2,
                                      keep_ref=TRUE) {
    if (!is.logical(keep_ref) || length(keep_ref) != 1L || is.na(keep_ref)) {
        stop("keep_ref must be TRUE or FALSE")
    }
    keep_cells <- colSums(cac) >= min_counts_per_cell
    cac <- cac[, keep_cells]

    cnt_mat <- as(cac, 'CsparseMatrix')
    keep_features <- Matrix::rowSums(cnt_mat > 0) >= min_cells_per_feature
    if (keep_ref && "allele_type" %in% colnames(cac@locus_map)) {
        keep_features <- keep_features | cac@locus_map$allele_type %in% "REF"
    }
    cac <- .new_cac(as(cnt_mat[keep_features, , drop = FALSE], "RsparseMatrix"),
                    cac@locus_map[keep_features, , drop = FALSE],
                    cac@cell_data)
    
    features_per_loc <- lengths(cac@locus_index)
    if (length(features_per_loc) > 0L) {
        cac <- cac[features_per_loc > 1]
    }

    cac
}
