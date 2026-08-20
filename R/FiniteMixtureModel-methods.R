
#' @include FiniteMixtureModel-class.R
#' @include LocusModel-methods.R

#' Predict from a FiniteMixtureModel
#'
#' @param object A FiniteMixtureModel.
#' @param newx Optional CellAlleleCounts for new cells. Required for
#'     `type="loadings"` and for `pseudocount` shrinkage.
#' @param type Either `"frequencies"` (default) to return predicted allele
#'     frequencies per cell, or `"loadings"` to infer clone responsibilities
#'     for new cells.
#' @param pseudocount Optional shrinkage parameter. When provided with newx,
#'     shrinks empirical frequencies toward model-predicted frequencies.
#' @param loadings For `type="frequencies"`, overrides the loadings in the
#'     model object (when `newx` is NULL) or the loadings inferred from
#'     non-NULL `newx`
#' @param ... Not used.
#'
#' @returns For `type="frequencies"`: a named list of A_l x N matrices of
#'     predicted allele frequencies. For `type="loadings"`: an N x K matrix.
#'
#' @rdname FiniteMixtureModel-class
#' @exportMethod predict
setMethod("predict", "FiniteMixtureModel", function(object, newx = NULL,
    type = c("frequencies", "loadings"), pseudocount = NULL, loadings = NULL, ...) {
    type <- match.arg(type)
    if (type == "frequencies") {
        .predict_frequencies(object, newx, pseudocount, loadings)
    } else {
        if (is.null(newx)) stop("predict(type='loadings') requires newx")
        .predict_loadings_fmm(object, newx)
    }
})

.predict_loadings.FiniteMixtureModel <- function(object, newx) {
    .predict_loadings_fmm(object, newx)
}

.fmm_e_step <- function(Y_t, Fmat, pi) {
    Fmat <- pmax(Fmat, 1e-12)
    N <- nrow(Y_t)
    K <- length(pi)

    log_L <- matrix(log(pi), nrow = N, ncol = K, byrow = TRUE) +
             as.matrix(Y_t %*% log(Fmat))

    m <- apply(log_L, 1, max)
    exp_shifted <- exp(log_L - m)
    exp_shifted / rowSums(exp_shifted)
}

.predict_loadings_fmm <- function(object, newx) {
    Y_t <- t(as(newx, 'CsparseMatrix'))
    L <- .fmm_e_step(Y_t, object@Fmat, object@pi)
    rownames(L) <- colnames(newx)
    colnames(L) <- colnames(object@Fmat)
    L
}

#' Get mixing proportions from a FiniteMixtureModel
#'
#' @param object A FiniteMixtureModel.
#' @returns Numeric vector of global clone proportions (sums to 1).
#' @rdname FiniteMixtureModel-class
#' @export
mixing_proportions <- function(object) {
    if (!is(object, "FiniteMixtureModel"))
        stop("mixing_proportions requires a FiniteMixtureModel")
    object@pi
}
