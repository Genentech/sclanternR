
#' LocusModel virtual class
#'
#' Virtual parent class for locus-aware models (LocusNMF and
#' FiniteMixtureModel). Stores per-locus allele frequency matrices (Fmat),
#' cell loadings (L), and locus metadata. Supports subsetting by loci and/or
#' cells via `[`.
#'
#' @slot Fmat Dense matrix (features x K) of allele frequencies. Columns sum
#'     to 1 within each locus.
#' @slot L Dense matrix (N x K) of cell loadings/responsibilities. Rows sum
#'     to 1. May be 0-row if constructed without cell data.
#' @slot loglik Numeric vector of log-likelihood values from fitting.
#' @slot locus_map DataFrame with columns: locus, allele, feature_name.
#' @slot locus_index Named list mapping locus names to row indices in Fmat.
#' @slot locus_1hot CsparseMatrix indicator (loci x features).
#' @slot cell_names Character vector of cell identifiers.
#' @slot cell_pseudocount Numeric; Dirichlet prior on loadings used in fitting.
#' @slot loc_pseudocount Numeric; Dirichlet prior on allele frequencies.
#' @slot nobs Numeric; total observation count (sum of count matrix) used for
#'     BIC computation. NA if not set (e.g. manually constructed models).
#'
#' @param x A LocusModel object.
#' @param i Index for loci: numeric, character, or logical.
#' @param j Index for cells: numeric, character, or logical.
#' @param object A LocusModel object.
#' @param do.NULL Not used.
#' @param prefix Not used.
#' @param ... Not used.
#' @param drop Not used.
#'
#' @include getFeatureAnno.R
#' @name LocusModel-class
#' @rdname LocusModel-class
#' @importFrom methods setClass new validObject is as
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom Matrix CsparseMatrix
setClass("LocusModel",
    contains = "VIRTUAL",
    slots = list(
        Fmat            = "matrix",
        L               = "matrix",
        loglik          = "numeric",
        locus_map       = "DataFrame",
        locus_index     = "list",
        locus_1hot      = "CsparseMatrix",
        cell_names      = "character",
        cell_pseudocount = "numeric",
        loc_pseudocount  = "numeric",
        nobs            = "numeric"
    ),
    validity = function(object) {
        errors <- character()
        if (nrow(object@Fmat) != nrow(object@locus_map))
            errors <- c(errors, "nrow(Fmat) must equal nrow(locus_map)")
        if (ncol(object@Fmat) != ncol(object@L) && nrow(object@L) > 0)
            errors <- c(errors, "ncol(Fmat) must equal ncol(L)")
        if (nrow(object@L) != length(object@cell_names))
            errors <- c(errors, "nrow(L) must equal length(cell_names)")
        if (length(errors) == 0) TRUE else errors
    }
)

.build_model_locus_index <- function(locus_map) {
    locus_chr <- as.character(locus_map$locus)
    idx <- split(seq_along(locus_chr), locus_chr)
    idx[unique(locus_chr)]
}

#' @importFrom Matrix fac2sparse
.build_model_locus_1hot <- function(locus_map) {
    Matrix::fac2sparse(locus_map$locus)
}

.make_locus_model_slots <- function(Fmat, L, loglik, locus_map,
                                    cell_names, cell_pseudocount,
                                    loc_pseudocount, nobs = NA_real_) {
    list(
        Fmat             = Fmat,
        L                = L,
        loglik           = loglik,
        locus_map        = locus_map,
        locus_index      = .build_model_locus_index(locus_map),
        locus_1hot       = .build_model_locus_1hot(locus_map),
        cell_names       = cell_names,
        cell_pseudocount = cell_pseudocount,
        loc_pseudocount  = loc_pseudocount,
        nobs             = nobs
    )
}

.flist_to_fmat <- function(F_list) {
    locus_names <- names(F_list)
    Fmat <- do.call(rbind, F_list)
    feature_names <- rownames(Fmat)
    K <- ncol(Fmat)

    parsed <- .getFeatureAnno(feature_names)
    locus_vec <- rep(locus_names, times = vapply(F_list, nrow, integer(1)))
    locus_map <- DataFrame(
        locus = locus_vec,
        allele = parsed$allele,
        feature_name = feature_names
    )

    list(Fmat = Fmat, locus_map = locus_map)
}
