#' Compute cell inner products from feature frequencies
#'
#' Compute the total or within-locus cell-by-cell inner product used by
#' [fitFiniteMixtureModel()] for hierarchical clustering initialization,
#' optionally retaining only its leading eigencomponents.
#'
#' @param cell_allele_counts A \code{\linkS4class{CellAlleleCounts}} object.
#' @param covariance One of \code{"total"} (default) or
#'     \code{"within_locus"}. Both normalize each cell's counts by its total
#'     count across all features. Empty cells have zero frequencies.
#' @param rank Optional positive integer giving the maximum retained rank.
#'     \code{NULL} (default) leaves the inner product unchanged. A rank at
#'     least as large as the number of cells also leaves it unchanged.
#'
#' @details
#' Let \eqn{X} be the cell-by-feature matrix of globally normalized frequencies.
#' The total inner product is \eqn{XX^T}. For within-locus inner products,
#' locus \eqn{l} has \eqn{m_l} features and cell-wise mean vector \eqn{\mu_l},
#' and the result is \eqn{XX^T - \sum_l m_l \mu_l \mu_l^T}.
#' #'
#' @returns A dense symmetric numeric matrix with cell names as both row and
#'     column names, in the input cell order.
#'
#' @examples
#' \dontrun{
#' K <- computeCellInnerProduct(counts, covariance = "within_locus", rank = 10)
#' f2 <- pmax(computeF2dist(K), 0)
#' euclidean <- stats::as.dist(sqrt(f2))
#' }
#'
#' @seealso [computeF2dist()], [fitFiniteMixtureModel()]
#' @export
computeCellInnerProduct <- function(cell_allele_counts,
                                    covariance = c("total", "within_locus"),
                                    rank = NULL) {
    if (!is(cell_allele_counts, "CellAlleleCounts")) {
        stop("cell_allele_counts must be a CellAlleleCounts object")
    }
    covariance <- match.arg(covariance)
    .validate_inner_product_rank(rank)

    Y_t <- MatrixExtra::t_shallow(as(cell_allele_counts, "CsparseMatrix"))
    locus_1hot <- if (covariance == "within_locus") {
        cell_allele_counts@locus_1hot
    } else {
        NULL
    }
    .cell_inner_product(Y_t, locus_1hot, rank)
}

.validate_inner_product_rank <- function(rank, arg = "rank") {
    if (!is.null(rank) &&
        (!is.numeric(rank) || is.complex(rank) || length(rank) != 1L || !is.finite(rank) ||
         rank < 1 || rank != floor(rank))) {
        stop(arg, " must be NULL or a positive integer")
    }
    invisible(NULL)
}

# Shared implementation for the public function and subsampled FMM initialization.
.cell_inner_product <- function(Y_t, locus_1hot = NULL, rank = NULL) {
    rs <- Matrix::rowSums(Y_t)
    rs[rs == 0] <- 1
    freq_mat <- Matrix::Diagonal(x = 1 / rs) %*% Y_t
    rownames(freq_mat) <- rownames(Y_t)

    .frequency_inner_product(freq_mat, locus_1hot, rank)
}

# Assumes frequencies are already normalized. May contain empty
# cells/clusters, avoid renormalizing those.
.frequency_inner_product <- function(freq_mat, locus_1hot = NULL, rank = NULL) {
    crossp <- as.matrix(Matrix::tcrossprod(freq_mat))
    if (!is.null(locus_1hot)) {
        n_features <- Matrix::rowSums(locus_1hot)
        # Z[, l] = sqrt(m_l) * mean frequency in locus l. Then
        # XX' - ZZ' equals the Gram matrix of within-locus residuals.
        locus_freq <- freq_mat %*% Matrix::t(locus_1hot) %*%
            Matrix::Diagonal(x = 1 / sqrt(n_features))
        crossp <- crossp - as.matrix(Matrix::tcrossprod(locus_freq))
    }

    if (!is.null(rank) && rank < nrow(crossp)) {
        eig <- eigen(crossp, symmetric = TRUE)
        keep <- seq_len(rank)
        coords <- sweep(eig$vectors[, keep, drop = FALSE], 2,
                        sqrt(pmax(eig$values[keep], 0)), "*")
        crossp <- tcrossprod(coords)
    }
    dimnames(crossp) <- list(rownames(freq_mat), rownames(freq_mat))
    crossp
}
