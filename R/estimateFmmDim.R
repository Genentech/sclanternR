
#' Estimate FMM dimensionality via PCA
#'
#' Performs CLR-weighted PCA on cluster allele frequencies and applies a
#' dimensionality selection criterion.
#'
#' @param model A \code{\linkS4class{FiniteMixtureModel}}.
#' @param method Dimensionality selection method. One of:
#'   \describe{
#'     \item{\code{"gavish_donoho"}}{(Default) Optimal hard threshold on
#'         singular values (Gavish & Donoho 2014).}
#'     \item{\code{"marchenko_pastur"}}{Retain eigenvalues exceeding
#'         the Marchenko-Pastur upper edge.}
#'     \item{\code{"kneedle"}}{Kneedle algorithm (Satopaa et al. 2011): finds
#'         the elbow in the screeplot by maximizing the distance from each
#'         eigenvalue to the baseline connecting the first and last points.}
#'   }
#'
#' @returns An S3 object of class \code{"FmmDimEstimate"} (a named list) with:
#'   \describe{
#'     \item{eigenvalues}{Numeric vector of PCA eigenvalues.}
#'     \item{cutoff_index}{Number of components passing the criterion.}
#'     \item{method}{Character string indicating the method used.}
#'     \item{threshold}{(\code{gavish_donoho}, \code{marchenko_pastur})
#'         Eigenvalue threshold used for the cutoff.}
#'     \item{curvature}{(\code{kneedle}) Difference-curve values.}
#'   }
#'
#' @export
#' @include FiniteMixtureModel-class.R
estimateFmmDimensionality <- function(
    model,
    method = c("gavish_donoho", "marchenko_pastur", "kneedle")
) {
    if (!is(model, "FiniteMixtureModel"))
        stop("model must be a FiniteMixtureModel")
    method <- match.arg(method)

    Fmat <- model@Fmat
    pi_vec <- model@pi
    locus_1hot <- model@locus_1hot

    pca <- .clr_weighted_pca(Fmat, pi_vec, locus_1hot)

    result <- switch(method,
        gavish_donoho = .estimate_gavish_donoho(pca),
        marchenko_pastur = .estimate_marchenko_pastur(pca),
        kneedle = .estimate_kneedle(pca$eigenvalues)
    )

    class(result) <- "FmmDimEstimate"
    result
}

# --- CLR + weighted PCA ---

.clr_weighted_pca <- function(Fmat, pi_vec, locus_1hot) {
    keep <- rowSums(Fmat) > 0
    if (!all(keep)) {
        Fmat <- Fmat[keep, , drop = FALSE]
        locus_1hot <- locus_1hot[, keep, drop = FALSE]
        loci_alive <- Matrix::rowSums(locus_1hot) > 0
        locus_1hot <- locus_1hot[loci_alive, , drop = FALSE]
    }

    Fmat_safe <- pmax(Fmat, 1e-12)
    log_Fmat <- log(Fmat_safe)

    n_alleles <- as.vector(Matrix::rowSums(locus_1hot))
    locus_sum_log <- as.matrix(locus_1hot %*% log_Fmat)
    locus_mean_log <- locus_sum_log / n_alleles
    mean_log_broadcast <- as.matrix(Matrix::t(locus_1hot) %*% locus_mean_log)
    clr_Fmat <- log_Fmat - mean_log_broadcast

    K <- length(pi_vec)
    mu <- as.vector(clr_Fmat %*% pi_vec)
    clr_centered <- clr_Fmat - mu
    clr_weighted <- sweep(clr_centered, 2, sqrt(pi_vec), '*')

    sv <- svd(clr_weighted)
    s <- sv$d
    eigenvalues <- s^2
    n_ev <- min(K - 1, length(eigenvalues))
    eigenvalues <- eigenvalues[seq_len(n_ev)]

    m <- nrow(clr_weighted)
    n <- ncol(clr_weighted)

    list(eigenvalues = eigenvalues, singular_values = s[seq_len(n_ev)],
         m = m, n = n, locus_1hot = locus_1hot,
         U = sv$u, V = sv$v, d = s, mu = mu, pi_vec = pi_vec,
         keep = keep)
}

# --- Gavish-Donoho optimal hard threshold ---

.estimate_gavish_donoho <- function(pca) {
    eigenvalues <- pca$eigenvalues
    sv <- pca$singular_values
    m <- pca$m
    n <- pca$n

    beta <- min(m, n) / max(m, n)

    lambda_star <- sqrt(2 * (beta + 1) +
                        8 * beta / ((beta + 1) + sqrt(beta^2 + 14 * beta + 1)))

    mp_med <- .marchenko_pastur_median(beta)
    omega <- lambda_star / sqrt(mp_med)

    tau_sv <- omega * median(sv)
    threshold <- tau_sv^2

    cutoff_index <- sum(eigenvalues > threshold)

    list(
        eigenvalues = eigenvalues,
        cutoff_index = cutoff_index,
        method = "gavish_donoho",
        threshold = threshold
    )
}

.marchenko_pastur_median <- function(beta) {
    a <- (1 - sqrt(beta))^2
    b <- (1 + sqrt(beta))^2

    n_grid <- 2000L
    x_grid <- seq(a, b, length.out = n_grid + 2L)[2:(n_grid + 1L)]
    y <- sqrt((b - x_grid) * (x_grid - a)) / (2 * pi * beta * x_grid)

    dx <- x_grid[2] - x_grid[1]
    cdf <- cumsum(y) * dx
    cdf <- cdf / cdf[n_grid]

    approx(cdf, x_grid, xout = 0.5)$y
}

# --- Marchenko-Pastur upper edge ---

.estimate_marchenko_pastur <- function(pca) {
    eigenvalues <- pca$eigenvalues
    sv <- pca$singular_values
    m <- pca$m
    n <- pca$n

    beta <- min(m, n) / max(m, n)
    mp_med <- .marchenko_pastur_median(beta)

    sigma_hat <- median(sv) / sqrt(mp_med)
    upper_edge <- sigma_hat^2 * (1 + sqrt(beta))^2

    cutoff_index <- sum(eigenvalues > upper_edge)

    list(
        eigenvalues = eigenvalues,
        cutoff_index = cutoff_index,
        method = "marchenko_pastur",
        threshold = upper_edge
    )
}

# --- Kneedle (elbow detection) ---

.estimate_kneedle <- function(eigenvalues) {
    p <- length(eigenvalues)

    if (p <= 2L) {
        return(list(
            eigenvalues = eigenvalues,
            cutoff_index = max(1L, p - 1L),
            method = "kneedle",
            curvature = rep(0, p)
        ))
    }

    x_norm <- seq(0, 1, length.out = p)
    y_range <- range(eigenvalues)
    if (diff(y_range) < .Machine$double.eps) {
        return(list(
            eigenvalues = eigenvalues,
            cutoff_index = 1L,
            method = "kneedle",
            curvature = rep(0, p)
        ))
    }
    y_norm <- (eigenvalues - y_range[1]) / diff(y_range)

    y_line <- y_norm[1] + (y_norm[p] - y_norm[1]) * x_norm
    curvature <- y_norm - y_line

    cutoff_index <- which.max(curvature)

    list(
        eigenvalues = eigenvalues,
        cutoff_index = cutoff_index,
        method = "kneedle",
        curvature = curvature
    )
}

# --- S3 methods ---

#' @export
print.FmmDimEstimate <- function(x, ...) {
    cat(sprintf("FmmDimEstimate (method: %s)\n", x$method))
    cat(sprintf("  cutoff_index: %d\n", x$cutoff_index))
    cat(sprintf("  eigenvalues: %s\n",
                paste(format(x$eigenvalues[seq_len(min(10, length(x$eigenvalues)))],
                             digits = 3), collapse = ", ")))
    if (length(x$eigenvalues) > 10)
        cat("  ...\n")
    invisible(x)
}

#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline theme_bw theme
#' @export
plot.FmmDimEstimate <- function(x, manual_cutoff = NULL, ...) {
    df <- data.frame(rank = seq_along(x$eigenvalues),
                     eigenvalues = x$eigenvalues)
    cutoff_df <- data.frame(
        cutoff = x$cutoff_index,
        method = x$method)
    if (!is.null(manual_cutoff)) {
        cutoff_df <- rbind(cutoff_df,
            data.frame(cutoff = manual_cutoff, method = "manual"))
    }
    p <- ggplot(df, aes(x = rank, y = eigenvalues)) +
        geom_line() + geom_point() +
        geom_vline(aes(xintercept = cutoff + 0.5, color = method),
                   lty = "dashed", data = cutoff_df) +
        theme_bw(base_size=16) +
        theme(legend.position = "bottom")
    p
}
