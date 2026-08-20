
#' Estimate Dirichlet pseudocount from data
#'
#' Estimates the Dirichlet concentration parameter (alpha_0) by fitting a
#' Dirichlet-Multinomial model to observed allele counts, with model-predicted
#' allele frequencies as the fixed mean direction. Uses Cox-Reid adjusted
#' profile likelihood (Newton-Raphson) to reduce the upward bias of the
#' standard MLE, initialized by a method-of-moments estimate.
#'
#' The estimated alpha_0 can be used as the \code{pseudocount} argument to
#' \code{\link{fitCellDistancePipeline}} or to the \code{predict} method of
#' a \code{\linkS4class{LocusModel}}, controlling how strongly cell-level
#' frequencies are shrunk toward model predictions.
#'
#' @param model A \code{\linkS4class{LocusModel}} (FiniteMixtureModel or
#'     LocusNMF).
#' @param cell_allele_counts A \code{\linkS4class{CellAlleleCounts}} object,
#'     matched to the same loci as \code{model}.
#' @param init Initial estimate for alpha_0, or \code{NULL} (default) to use
#'     a method-of-moments estimate based on the Pearson chi-squared statistic.
#' @param maxiter Maximum number of Newton-Raphson iterations (default 100).
#' @param tol Convergence tolerance on relative change in alpha_0
#'     (default 1e-6).
#' @param min_alpha Minimum allowed alpha_0 (default 1e-4).
#' @param max_alpha Maximum allowed alpha_0 (default 1e6).
#' @param loadings If non-NULL, use these loadings instead of inferring them
#'     from `cell_allele_counts`
#' @param verbose Verbosity level (default \code{FALSE}). Set to \code{TRUE}
#'     (or a positive integer) for progress messages.
#'
#' @returns An S3 object of class \code{"DirichletPseudocountEstimate"} with:
#'   \describe{
#'     \item{pseudocount}{The estimated alpha_0.}
#'     \item{mom_estimate}{The method-of-moments initialization (or the
#'         user-supplied \code{init} value).}
#'     \item{iterations}{Number of Newton-Raphson iterations run.}
#'     \item{converged}{Logical; whether the algorithm converged within
#'         \code{maxiter} iterations.}
#'     \item{trace}{Numeric vector of alpha_0 values at each iteration.}
#'   }
#'
#' @seealso \code{\link{fitCellDistancePipeline}}
#'
#' @examples
#' \dontrun{
#' fmm <- fitFiniteMixtureModel(counts, n_clusters = 20)
#' est <- estimateDirichletPseudocount(fmm, counts)
#' est$pseudocount
#'
#' # Use in the pipeline:
#' result <- fitCellDistancePipeline(counts, n_loci = 100, pseudocount = "auto")
#' }
#'
#' @include LocusModel-class.R
#' @export
estimateDirichletPseudocount <- function(
    model,
    cell_allele_counts,
    init = NULL,
    maxiter = 100,
    tol = 1e-6,
    min_alpha = 1e-4,
    max_alpha = 1e6,
    loadings = NULL,
    verbose = FALSE
) {
    if (!is(model, "LocusModel"))
        stop("model must be a LocusModel (FiniteMixtureModel or LocusNMF)")
    if (!is(cell_allele_counts, "CellAlleleCounts"))
        stop("cell_allele_counts must be a CellAlleleCounts object")

    .log_msg(1, verbose, "Estimating Dirichlet pseudocount (%d loci, %d cells)",
             length(model), ncol(model))

    if (is.null(loadings)) {
        L <- .predict_loadings(model, cell_allele_counts)
    } else {
        L <- loadings
    }

    Fmat <- model@Fmat
    P <- t(L %*% t(Fmat))

    P_denom <- as.matrix(
        Matrix::t(model@locus_1hot) %*% model@locus_1hot %*% P)
    P <- pmax(P / P_denom, 1e-12)

    X <- as.matrix(as(cell_allele_counts, 'CsparseMatrix'))
    N_cl <- as.matrix(model@locus_1hot %*% X)

    if (is.null(init)) {
        alpha <- .mom_dirichlet_pseudocount(X, P, N_cl, model@locus_1hot,
                                            min_alpha, max_alpha)
        .log_msg(1, verbose, "MOM initialization: %.4g", alpha)
    } else {
        alpha <- max(min_alpha, min(max_alpha, init))
        .log_msg(1, verbose, "User initialization: %.4g", alpha)
    }
    mom_est <- alpha

    P2 <- P * P
    P3 <- P2 * P

    trace <- numeric(maxiter)
    converged <- FALSE

    for (iter in seq_len(maxiter)) {
        trace[iter] <- alpha
        alpha_P <- alpha * P

        score <- sum(P * (digamma(X + alpha_P) - digamma(alpha_P))) -
                 sum(digamma(N_cl + alpha) - digamma(alpha))

        hessian <- sum(P2 * (trigamma(X + alpha_P) - trigamma(alpha_P))) +
                   sum(trigamma(alpha) - trigamma(N_cl + alpha))

        third <- sum(P3 * (psigamma(X + alpha_P, 2) - psigamma(alpha_P, 2))) +
                 sum(psigamma(alpha, 2) - psigamma(N_cl + alpha, 2))

        if (hessian >= 0) {
            .log_msg(1, verbose,
                     "Stopping: non-negative Hessian at iter %d (alpha=%.4g)",
                     iter, alpha)
            break
        }

        adj_score <- score - 0.5 * third / hessian

        alpha_new <- alpha - adj_score / hessian
        alpha_new <- max(min_alpha, min(max_alpha, alpha_new))

        .log_msg(2, verbose, "  iter %d: alpha=%.4g -> %.4g (score=%.4g)",
                 iter, alpha, alpha_new, adj_score)

        if (abs(alpha_new - alpha) / (abs(alpha) + 1e-12) < tol) {
            alpha <- alpha_new
            converged <- TRUE
            break
        }
        alpha <- alpha_new
    }

    trace <- trace[seq_len(iter)]
    .log_msg(1, verbose, "Converged: %s after %d iterations (alpha=%.4g)",
             converged, iter, alpha)

    result <- list(
        pseudocount = alpha,
        mom_estimate = mom_est,
        iterations = iter,
        converged = converged,
        trace = trace
    )
    class(result) <- "DirichletPseudocountEstimate"
    result
}

# Method of moments via Pearson chi-squared
.mom_dirichlet_pseudocount <- function(X, P, N_cl, locus_1hot,
                                       min_alpha, max_alpha) {
    A <- as.integer(Matrix::rowSums(locus_1hot))

    N_feat <- as.matrix(Matrix::t(locus_1hot) %*% N_cl)
    E <- N_feat * P

    chi2 <- (X - E)^2 / pmax(E, 1e-12)
    chi2[N_feat == 0] <- 0

    S <- sum(chi2)

    Am1 <- A - 1L
    D <- sum(Am1 * (N_cl > 0))
    W <- sum(Am1 * N_cl)

    if (S <= D || W <= S) return(max_alpha)

    alpha <- (W - S) / (S - D)
    max(min_alpha, min(max_alpha, alpha))
}

#' @export
print.DirichletPseudocountEstimate <- function(x, ...) {
    cat("DirichletPseudocountEstimate (Cox-Reid adjusted)\n")
    cat(sprintf("  pseudocount: %.4g\n", x$pseudocount))
    cat(sprintf("  MOM init: %.4g\n", x$mom_estimate))
    cat(sprintf("  iterations: %d (converged: %s)\n",
                x$iterations, x$converged))
    invisible(x)
}
