
#' Compute per-locus mutual information for a fitted model
#'
#' Computes the scaled mutual information (G-test statistic) between
#' topics/clones and alleles at each locus. Can be called on a cell subset
#' by first subsetting both the model and data.
#'
#' @param object A LocusModel (LocusNMF or FiniteMixtureModel).
#' @param x A CellAlleleCounts providing the count data, or NULL (the
#'   default) to compute mutual information from the model's frequency
#'   matrix weighted by topic proportions (pi for FMM, colSums(L) for
#'   LocusNMF).
#' @param ... Not used.
#'
#' @returns Named numeric vector of per-locus mutual information scores.
#'
#' @examples
#' \dontrun{
#' model <- fitLocusNMF(counts)
#' mi <- mutualinfo(model, counts)
#'
#' # Without data, uses model frequencies weighted by topic proportions:
#' mi_model <- mutualinfo(model)
#'
#' # On a cell subset:
#' cells <- c("cell1", "cell2")
#' mi_sub <- mutualinfo(model[, cells], counts[, cells])
#' }
#'
#' @include LocusNMF-class.R
#' @include FiniteMixtureModel-class.R
#' @include LocusModel-methods.R
#' @export
setGeneric("mutualinfo", function(object, x = NULL, ...) standardGeneric("mutualinfo"))

#' @rdname mutualinfo
#' @exportMethod mutualinfo
setMethod("mutualinfo", signature(object = "LocusNMF", x = "NULL"),
    function(object, x = NULL, ...) {
        Fmat_stat <- sweep(object@Fmat, 2, colSums(object@L), `*`)
        .compute_vectorized_mutinfo(Fmat_stat, object@locus_1hot)
    }
)

#' @rdname mutualinfo
#' @exportMethod mutualinfo
setMethod("mutualinfo", signature(object = "LocusNMF", x = "CellAlleleCounts"),
    function(object, x, ...) {
        Y_t <- t(as(x, 'CsparseMatrix'))
        L <- object@L
        Fmat <- object@Fmat

        expected_sddmm <- sddmm_csc(Y_t, L, t(Fmat))
        expected_sddmm@x <- pmax(expected_sddmm@x, 1e-12)

        denom_sddmm <- expected_sddmm
        denom_sddmm@x <- 1 / denom_sddmm@x
        R_t <- Y_t * denom_sddmm

        Fmat_stat <- Fmat * as.matrix(t(R_t) %*% L)

        .compute_vectorized_mutinfo(Fmat_stat, object@locus_1hot)
    }
)

#' @rdname mutualinfo
#' @exportMethod mutualinfo
setMethod("mutualinfo", signature(object = "FiniteMixtureModel", x = "NULL"),
    function(object, x = NULL, ...) {
        Fmat_stat <- sweep(object@Fmat, 2, object@pi, `*`)
        .compute_vectorized_mutinfo(Fmat_stat, object@locus_1hot)
    }
)

#' @rdname mutualinfo
#' @exportMethod mutualinfo
setMethod("mutualinfo", signature(object = "FiniteMixtureModel", x = "CellAlleleCounts"),
    function(object, x, ...) {
        Y_t <- t(as(x, 'CsparseMatrix'))
        L <- object@L

        Fmat_stat <- t(as.matrix(t(L) %*% Y_t))

        .compute_vectorized_mutinfo(Fmat_stat, object@locus_1hot)
    }
)

# Returns mutual info scaled by the total weight in x.
# see also: g-test statistic
# https://en.wikipedia.org/wiki/Mutual_information#For_discrete_data
# https://en.wikipedia.org/wiki/G-test
.compute_scaled_mutinfo <- function(x) {
    p1 <- rowSums(x) / sum(x)
    p2 <- colSums(x) / sum(x)
    p <- x / sum(x)
    ret <- p * log(p / outer(p1, p2))
    ret[x == 0] <- 0
    sum(x) * sum(ret)
}

# Vectorized G-test across all loci simultaneously.
# Fmat_stat: (n_features x K) matrix of joint counts/weights.
# locus_1hot: (n_loci x n_features) sparse indicator matrix.
.compute_vectorized_mutinfo <- function(Fmat_stat, locus_1hot) {
    C <- as.matrix(locus_1hot %*% Fmat_stat)
    S <- rowSums(C)
    r <- rowSums(Fmat_stat)

    E <- r * as.matrix(Matrix::crossprod(locus_1hot, C / S))

    G <- Fmat_stat * (log(Fmat_stat) - log(E))
    #G[Fmat_stat == 0] <- 0
    G[Fmat_stat < .Machine$double.eps] <- 0

    mi <- as.numeric(locus_1hot %*% rowSums(G))
    names(mi) <- rownames(locus_1hot)
    mi
}
