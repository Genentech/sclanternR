
#' Compute F2-distance from (imputed) allele frequencies
#'
#' @param freqs A list of per-locus allele frequency matrices (A_l x N) as
#'     returned by \code{predict()}, or a precomputed N x N inner
#'     product matrix.
#'
#' @returns A matrix of the F2 distances between every pair of cells.
#'
#' @export
computeF2dist <- function(freqs) {
    if (is.matrix(freqs)) {
        d <- freqs
    } else {
        mat <- t(do.call(rbind, freqs))
        d <- mat %*% t(mat)
    }
    d_diag <- diag(d)
    d <- -2 * d
    d <- sweep(d, 1, d_diag, '+')
    d <- sweep(d, 2, d_diag, '+')
    d
}
