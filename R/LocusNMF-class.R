
#' LocusNMF class
#'
#' S4 class representing a fitted locus-aware NMF model. Each cell's allele
#' frequencies are a linear combination of topic-specific frequencies, weighted
#' by the cell's topic proportions (L).
#'
#' @include LocusModel-class.R
#' @name LocusNMF-class
#' @rdname LocusNMF-class
#' @exportClass LocusNMF
setClass("LocusNMF", contains = "LocusModel")

#' Construct a LocusNMF model
#'
#' @param F_list Named list of A_l x K matrices (columns sum to 1 within each
#'     locus). Names are locus identifiers; rownames are feature names in
#'     format `chr{chr}_{position}_{allele}`.
#' @param L Optional N x K matrix of cell-topic proportions (rows sum to 1).
#'     If NULL, the model is constructed without cell data and predict()
#'     requires newx.
#' @param loglik Optional numeric vector of log-likelihood trace.
#' @param cell_pseudocount Pseudocount for topic proportions (default 0.01).
#' @param loc_pseudocount Pseudocount per locus (default 1).
#'
#' @returns A `LocusNMF` object.
#' @export
LocusNMF <- function(F_list, L = NULL, loglik = numeric(0),
                     cell_pseudocount = 0.01, loc_pseudocount = 1,
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

    do.call(new, c("LocusNMF", slots))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
