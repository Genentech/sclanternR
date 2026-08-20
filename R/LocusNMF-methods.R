
#' @include LocusNMF-class.R
#' @include LocusModel-methods.R

#' Predict from a LocusNMF model
#'
#' @param object A LocusNMF model.
#' @param newx Optional CellAlleleCounts for new cells. Required for
#'     `type="loadings"` and for `pseudocount` shrinkage.
#' @param type Either `"frequencies"` (default) to return predicted allele
#'     frequencies per cell, or `"loadings"` to infer topic proportions for
#'     new cells.
#' @param pseudocount Optional shrinkage parameter. When provided with newx,
#'     shrinks empirical frequencies toward model-predicted frequencies.
#' @param maxiter Maximum iterations for loadings inference (default 50).
#' @param tol Convergence tolerance for loadings inference (default 1e-4).
#' @param ... Not used.
#'
#' @returns For `type="frequencies"`: a named list of A_l x N matrices of
#'     predicted allele frequencies. For `type="loadings"`: an N x K matrix.
#'
#' @rdname LocusNMF-class
#' @exportMethod predict
setMethod("predict", "LocusNMF", function(object, newx = NULL,
    type = c("frequencies", "loadings"), pseudocount = NULL,
    maxiter = 50, tol = 1e-4, ...) {
    type <- match.arg(type)
    if (type == "frequencies") {
        .predict_frequencies(object, newx, pseudocount)
    } else {
        if (is.null(newx)) stop("predict(type='loadings') requires newx")
        .predict_loadings_nmf(object, newx, maxiter = maxiter, tol = tol)
    }
})

.predict_loadings.LocusNMF <- function(object, newx) {
    .predict_loadings_nmf(object, newx)
}

.predict_loadings_nmf <- function(object, newx, maxiter = 50, tol = 1e-4,
                                  init_L = NULL) {
    Y_t <- t(as(newx, 'CsparseMatrix'))
    Fmat <- object@Fmat
    K <- ncol(Fmat)
    N <- nrow(Y_t)
    cell_pseudocount <- object@cell_pseudocount

    if (!is.null(init_L)) {
        L <- init_L
    } else {
        L <- matrix(1/K, nrow = N, ncol = K)
    }

    for (iter in seq_len(maxiter)) {
        expected_sddmm <- sddmm_csc(Y_t, L, t(Fmat))
        expected_sddmm@x <- pmax(expected_sddmm@x, 1e-12)

        denom_sddmm <- expected_sddmm
        denom_sddmm@x <- 1 / denom_sddmm@x
        R_t <- Y_t * denom_sddmm

        L_stat <- L * as.matrix(R_t %*% Fmat)
        L_stat <- L_stat + cell_pseudocount
        L_new <- sweep(L_stat, 1, rowSums(L_stat), '/')

        if (iter > 1 && max(abs(L_new - L)) < tol) break
        L <- L_new
    }

    rownames(L) <- colnames(newx)
    colnames(L) <- colnames(object@Fmat)
    L
}
