
#' @export
filter_cell_allele_counts <- function(cac,
                                      min_counts_per_cell=0,
                                      min_cells_per_feature=2) {
    keep_cells <- colSums(cac) >= min_counts_per_cell
    cac <- cac[, keep_cells]

    cnt_mat <- as(cac, 'CsparseMatrix')
    keep_features <- Matrix::rowSums(cnt_mat > 0) >= min_cells_per_feature
    cac <- CellAlleleCounts(cnt_mat[keep_features,])
    
    features_per_loc <- Matrix::rowSums(cac@locus_1hot)
    cac <- cac[features_per_loc > 1]

    cac
}
