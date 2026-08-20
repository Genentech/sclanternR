
#' FiniteMixtureModel class
#'
#' S4 class representing a fitted Finite Mixture Model (Latent Class Model).
#' Each cell belongs to one of K latent clones; the global clone prevalence
#' vector pi drives assignments, and L represents posterior responsibilities.
#'
#' @slot pi Numeric vector (length K) of global clone proportions (sums to 1).
#'
#' @include LocusModel-class.R
#' @name FiniteMixtureModel-class
#' @rdname FiniteMixtureModel-class
#' @exportClass FiniteMixtureModel
setClass("FiniteMixtureModel",
    contains = "LocusModel",
    slots = list(pi = "numeric"),
    validity = function(object) {
        errors <- character()
        if (ncol(object@L) > 0 && length(object@pi) != ncol(object@L))
            errors <- c(errors, "length(pi) must equal ncol(L) (K)")
        if (abs(sum(object@pi) - 1) > 1e-10)
            errors <- c(errors, "pi must sum to 1")
        if (length(errors) == 0) TRUE else errors
    }
)

#' Construct a FiniteMixtureModel
#'
#' @param F_list Named list of A_l x K matrices (columns sum to 1 within each
#'     locus). Names are locus identifiers; rownames are feature names.
#' @param pi Numeric vector (length K) of global clone proportions (sums to 1).
#' @param L Optional N x K matrix of cell-clone responsibilities (rows sum to 1).
#'     If NULL, the model is constructed without cell data.
#' @param loglik Optional numeric vector of log-likelihood trace.
#' @param cell_pseudocount Pseudocount for clone proportions (default 1).
#' @param loc_pseudocount Pseudocount per locus (default 1).
#'
#' @returns A `FiniteMixtureModel` object.
#' @export
FiniteMixtureModel <- function(F_list, pi, L = NULL, loglik = numeric(0),
                               cell_pseudocount = 1, loc_pseudocount = 1,
                               nobs = NA_real_) {
    parsed <- .flist_to_fmat(F_list)
    Fmat <- parsed$Fmat
    locus_map <- parsed$locus_map

    if (is.null(L)) {
        L <- matrix(numeric(0), nrow = 0, ncol = ncol(Fmat))
        cell_names <- character(0)
    } else {
        cell_names <- rownames(L) %||% paste0("cell", seq_len(nrow(L)))
    }

    slots <- .make_locus_model_slots(
        Fmat = Fmat, L = L, loglik = loglik,
        locus_map = locus_map, cell_names = cell_names,
        cell_pseudocount = cell_pseudocount,
        loc_pseudocount = loc_pseudocount,
        nobs = nobs
    )

    do.call(new, c("FiniteMixtureModel", slots, list(pi = pi)))
}
