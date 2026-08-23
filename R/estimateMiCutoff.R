
#' Estimate a mutual information cutoff for locus selection
#'
#' Computes per-locus mutual information and applies one of several
#' automatic threshold methods to select informative loci.
#'
#' @param model A \code{\linkS4class{LocusModel}} (LocusNMF or
#'     FiniteMixtureModel).
#' @param x Optional \code{\linkS4class{CellAlleleCounts}}; passed to
#'     \code{\link{mutualinfo}}.
#' @param method Cutoff method. One of:
#'   \describe{
#'     \item{\code{"fdr"}}{(Default) BH-adjusted p-values from a G-test
#'         (chi-squared approximation).}
#'     \item{\code{"mad"}}{Loci with MI exceeding
#'         \code{median + mad_threshold * MAD}. Conservative and robust
#'         to outliers.}
#'     \item{\code{"elbow"}}{Knee-point detection on the sorted MI curve
#'         via maximum distance to the baseline (Kneedle algorithm).}
#'   }
#' @param alpha Significance level for the \code{"fdr"} method (default 0.05).
#' @param mad_threshold Number of MADs above the median for the \code{"mad"}
#'     method (default 3).
#' @param min_loci Minimum number of loci to keep (default 10)
#' @param max_loci Maximum number of loci to keep (default 300)
#'
#' @returns An S3 object of class \code{"MiCutoffEstimate"} (a named list)
#'   with:
#'   \describe{
#'     \item{mutualinfo}{Named numeric vector of per-locus MI scores,
#'         unsorted (preserving the model's locus order).}
#'     \item{cutoff}{Numeric scalar: MI threshold. Loci with
#'         \code{mi >= cutoff} are selected.}
#'     \item{n_loci}{Number of loci selected.}
#'     \item{loci}{Character vector of selected locus names.}
#'     \item{method}{Character string naming the method used.}
#'   }
#'   Additional method-specific elements:
#'   \describe{
#'     \item{pvalues}{(\code{fdr}) BH-adjusted p-values, unsorted.}
#'     \item{df}{(\code{fdr}) Per-locus degrees of freedom.}
#'     \item{mad_stats}{(\code{mad}) List with \code{median} and \code{mad}.}
#'     \item{curvature}{(\code{elbow}) Difference-curve values (sorted).}
#'   }
#'
#' @examples
#' \dontrun{
#' fmm <- fitFiniteMixtureModel(counts, n_clusters = 10)
#'
#' # Default: MAD-based cutoff
#' est <- estimateMiCutoff(fmm, counts)
#' plot(est)
#'
#' # Compare methods:
#' est_fdr   <- estimateMiCutoff(fmm, counts, method = "fdr")
#' est_elbow <- estimateMiCutoff(fmm, counts, method = "elbow")
#' }
#'
#' @seealso \code{\link{mutualinfo}}, \code{\link{fitCellDistancePipeline}}
#'
#' @importFrom stats pchisq p.adjust mad median quantile
#' @include LocusModel-class.R
#' @export
estimateMiCutoff <- function(
    model,
    x = NULL,
    method = c("fdr", "mad", "elbow"),
    alpha = 0.05,
    mad_threshold = 3,
    min_loci = 10,
    max_loci = 300
) {
    if (!is(model, "LocusModel"))
        stop("model must be a LocusModel (FiniteMixtureModel or LocusNMF)")
    method <- match.arg(method)

    mi <- mutualinfo(model, x)

    n_alleles <- as.integer(Matrix::rowSums(model@locus_1hot))
    K <- ncol(model@Fmat)
    method_result <- .apply_mi_cutoff(mi, method, n_alleles, K,
                                      alpha, mad_threshold)
    cutoff <- method_result$cutoff
    passes <- mi >= cutoff

    n_loci <- pmax(min_loci, pmin(max_loci, sum(passes)))

    mi_rank <- rank(-mi, ties.method="first")
    passes <- mi_rank <= n_loci

    loci <- names(mi)[passes]

    result <- c(
        list(
            mutualinfo = mi,
            n_loci = sum(passes),
            loci = loci
        ),
        method_result
    )

    class(result) <- "MiCutoffEstimate"
    result
}

# --- Internal: dispatch to method helpers (non-permutation) ---

.apply_mi_cutoff <- function(mi, method, n_alleles, K,
                             alpha = 0.05, mad_threshold = 3) {
    switch(method,
        fdr   = .estimate_mi_fdr(mi, n_alleles, K, alpha),
        elbow = .estimate_mi_elbow(mi),
        mad   = .estimate_mi_mad(mi, mad_threshold)
    )
}

# --- FDR via G-test chi-squared approximation ---

.estimate_mi_fdr <- function(mi, n_alleles, K, alpha) {
    # current MI values are G/2 = sum(O * log(O/E))
    G <- 2 * mi
    df <- (K - 1L) * (n_alleles - 1L)
    names(df) <- names(mi)

    pvals_raw <- ifelse(df > 0, pchisq(G, df = df, lower.tail = FALSE), 1)
    pvals_adj <- p.adjust(pvals_raw, method = "BH")
    names(pvals_adj) <- names(mi)

    passes <- pvals_adj < alpha
    cutoff <- if (any(passes)) min(mi[passes]) else Inf

    list(cutoff = cutoff, pvalues = pvals_adj, df = df)
}

# --- Elbow detection via Kneedle algorithm ---

.estimate_mi_elbow <- function(mi) {
    mi_sorted <- sort(mi, decreasing = TRUE)
    n <- length(mi_sorted)

    if (n <= 2L) {
        return(list(cutoff = min(mi_sorted),
                    curvature = setNames(rep(0, n), names(mi_sorted))))
    }

    x_norm <- seq(0, 1, length.out = n)
    y_range <- range(mi_sorted)
    if (diff(y_range) < .Machine$double.eps) {
        return(list(cutoff = mi_sorted[1],
                    curvature = setNames(rep(0, n), names(mi_sorted))))
    }
    y_norm <- (mi_sorted - y_range[1]) / diff(y_range)

    y_line <- y_norm[1] + (y_norm[n] - y_norm[1]) * x_norm
    diff_curve <- y_norm - y_line
    names(diff_curve) <- names(mi_sorted)

    elbow_idx <- which.max(diff_curve)
    cutoff <- mi_sorted[elbow_idx]

    list(cutoff = cutoff, curvature = diff_curve)
}

# --- MAD-based outlier detection ---

.estimate_mi_mad <- function(mi, mad_threshold) {
    med <- median(mi)
    mad_val <- mad(mi)

    cutoff <- med + mad_threshold * mad_val

    list(cutoff = cutoff, mad_stats = list(median = med, mad = mad_val))
}

# --- S3 methods ---

#' @export
print.MiCutoffEstimate <- function(x, ...) {
    cat("MiCutoffEstimate \n")
    cat(sprintf("  %d / %d loci selected (cutoff = %.4g)\n",
                x$n_loci, length(x$mutualinfo), x$cutoff))
    if (x$n_loci > 0 && x$n_loci <= 10) {
        cat(sprintf("  loci: %s\n", paste(x$loci, collapse = ", ")))
    } else if (x$n_loci > 10) {
        cat(sprintf("  loci: %s, ... (%d more)\n",
                    paste(x$loci[seq_len(5)], collapse = ", "), x$n_loci - 5L))
    }
    invisible(x)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline scale_y_continuous scale_x_log10 theme_bw theme
#' @export
plot.MiCutoffEstimate <- function(x) {
    df <- data.frame(
        mutualinfo = x$mutualinfo,
        rank = rank(-x$mutualinfo))
    ggplot(df, aes(x = rank, y = mutualinfo)) +
        geom_vline(aes(xintercept = x$n_loci + 0.5),
                   lty = "dashed") +
        geom_point(shape = 1) +
        scale_y_continuous(trans = "log1p", breaks = c(0, 10^(0:6))) +
        scale_x_log10() +
        theme_bw(base_size=16) +
        theme(legend.position = "bottom")
}
