
#' @include LocusModel-class.R
#' @include LocusNMF-class.R
#' @include FiniteMixtureModel-class.R
#' @importFrom methods setMethod setGeneric is as new
#' @importFrom stats setNames coef predict loadings logLik AIC BIC nobs
#' @importFrom utils tail

#' @rdname LocusModel-class
#' @exportMethod length
setMethod("length", "LocusModel", function(x) {
    length(x@locus_index)
})

#' @rdname LocusModel-class
#' @exportMethod names
setMethod("names", "LocusModel", function(x) {
    names(x@locus_index)
})

#' @rdname LocusModel-class
#' @importFrom BiocGenerics ncol
#' @exportMethod ncol
setMethod("ncol", "LocusModel", function(x) {
    nrow(x@L)
})

#' @rdname LocusModel-class
#' @importFrom BiocGenerics colnames
#' @exportMethod colnames
setMethod("colnames", "LocusModel", function(x, do.NULL = TRUE, prefix = "col") {
    x@cell_names
})

#' @rdname LocusModel-class
#' @exportMethod show
setMethod("show", "LocusModel", function(object) {
    cat(sprintf("%s: %d loci, %d cells, K=%d\n",
                class(object), length(object), ncol(object), ncol(object@Fmat)))
    if (length(object@loglik) > 0) {
        cat(sprintf("  log-lik: %.2f (%d iterations)\n",
                    tail(object@loglik, 1), length(object@loglik)))
    }
})

#' @rdname LocusModel-class
#' @exportMethod coef
setGeneric("coef")
setMethod("coef", "LocusModel", function(object, ...) {
    lapply(object@locus_index,
           function(idx) object@Fmat[idx, , drop = FALSE])
})

#' @rdname LocusModel-class
#' @exportMethod loadings
setGeneric("loadings")
setMethod("loadings", "LocusModel", function(x) {
    x@L
})

# --- Subsetting ---

#' @rdname LocusModel-class
#' @exportMethod [
setMethod("[", signature(x = "LocusModel", i = "ANY", j = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        if (!missing(i)) x <- .subset_model_loci_any(x, i)
        if (!missing(j)) x <- .subset_model_cells_any(x, j)
        x
    }
)

#' @rdname LocusModel-class
#' @exportMethod [
setMethod("[", signature(x = "LocusModel", i = "ANY", j = "missing"),
    function(x, i, j, ..., drop = TRUE) {
        .subset_model_loci_any(x, i)
    }
)

#' @rdname LocusModel-class
#' @exportMethod [
setMethod("[", signature(x = "LocusModel", i = "missing", j = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        .subset_model_cells_any(x, j)
    }
)

.subset_model_loci_any <- function(x, i) {
    all_loci <- names(x@locus_index)
    if (is.logical(i)) {
        i <- all_loci[i]
    } else if (is.numeric(i)) {
        i <- all_loci[i]
    }
    .subset_model_loci(x, i)
}

.subset_model_loci <- function(x, loci) {
    idx <- unlist(x@locus_index[loci], use.names = FALSE)
    sub_Fmat <- x@Fmat[idx, , drop = FALSE]
    sub_locus_map <- x@locus_map[idx, , drop = FALSE]

    slots <- .make_locus_model_slots(
        Fmat = sub_Fmat, L = x@L, loglik = x@loglik,
        locus_map = sub_locus_map, cell_names = x@cell_names,
        cell_pseudocount = x@cell_pseudocount,
        loc_pseudocount = x@loc_pseudocount
    )

    if (is(x, "FiniteMixtureModel")) {
        do.call(new, c("FiniteMixtureModel", slots, list(pi = x@pi)))
    } else {
        do.call(new, c(class(x), slots))
    }
}

.subset_model_cells_any <- function(x, j) {
    if (is.character(j)) {
        cidx <- match(j, x@cell_names)
    } else if (is.logical(j)) {
        cidx <- which(j)
    } else {
        cidx <- j
    }

    sub_L <- x@L[cidx, , drop = FALSE]
    sub_cell_names <- x@cell_names[cidx]

    slots <- .make_locus_model_slots(
        Fmat = x@Fmat, L = sub_L, loglik = x@loglik,
        locus_map = x@locus_map, cell_names = sub_cell_names,
        cell_pseudocount = x@cell_pseudocount,
        loc_pseudocount = x@loc_pseudocount
    )

    if (is(x, "FiniteMixtureModel")) {
        do.call(new, c("FiniteMixtureModel", slots, list(pi = x@pi)))
    } else {
        do.call(new, c(class(x), slots))
    }
}

# --- Shared predict helper for frequencies ---

.predict_frequencies <- function(object, newx, pseudocount, loadings=NULL) {
    if (is.null(newx) && nrow(object@L) == 0) {
        stop("predict(type='frequencies') without newx requires L in model")
    }

    L <- if (!is.null(loadings)) {
        loadings
    } else if (!is.null(newx)) {
        .predict_loadings(object, newx)
    } else {
        object@L
    }

    Fmat <- object@Fmat
    proj <- t(L %*% t(Fmat))

    if (!is.null(pseudocount)) {
        if (is.null(newx)) {
            stop("pseudocount requires newx")
        }
        if (pseudocount < Inf) {
            proj <- pseudocount * proj + as.matrix(as(newx, 'CsparseMatrix'))
        }
        proj <- proj / as.matrix(
            t(object@locus_1hot) %*% object@locus_1hot %*% proj)
    }

    lapply(object@locus_index,
           function(idx) proj[idx, , drop = FALSE])
}

# Dispatches to model-specific loadings inference
.predict_loadings <- function(object, newx) {
    UseMethod(".predict_loadings")
}

# --- logLik, AIC, BIC ---

.nfree <- function(object) {
    K <- ncol(object@Fmat)
    n_features <- nrow(object@Fmat)
    n_loci <- length(object@locus_index)
    N <- nrow(object@L)

    df_Fmat <- (n_features - n_loci) * K
    df_L <- N * (K - 1)
    df <- df_Fmat + df_L

    if (is(object, "FiniteMixtureModel")) {
        df <- df + (K - 1)
    }
    df
}

#' @rdname LocusModel-class
#' @export
logLik.LocusModel <- function(object, ...) {
    ll <- tail(object@loglik, 1)
    if (length(ll) == 0) stop("model has no log-likelihood (not fitted)")
    structure(ll, df = .nfree(object), nobs = object@nobs, class = "logLik")
}

#' @rdname LocusModel-class
#' @export
logLik.LocusNMF <- logLik.LocusModel

#' @rdname LocusModel-class
#' @export
logLik.FiniteMixtureModel <- logLik.LocusModel

#' @rdname LocusModel-class
#' @export
AIC.LocusModel <- function(object, ..., k = 2) {
    stats::AIC(logLik(object), k = k)
}

#' @rdname LocusModel-class
#' @export
AIC.LocusNMF <- AIC.LocusModel

#' @rdname LocusModel-class
#' @export
AIC.FiniteMixtureModel <- AIC.LocusModel

#' @rdname LocusModel-class
#' @export
BIC.LocusModel <- function(object, ...) {
    stats::BIC(logLik(object))
}

#' @rdname LocusModel-class
#' @export
BIC.LocusNMF <- BIC.LocusModel

#' @rdname LocusModel-class
#' @export
BIC.FiniteMixtureModel <- BIC.LocusModel

#' @rdname LocusModel-class
#' @export
nobs.LocusModel <- function(object, ...) {
    object@nobs
}

#' @rdname LocusModel-class
#' @export
nobs.LocusNMF <- nobs.LocusModel

#' @rdname LocusModel-class
#' @export
nobs.FiniteMixtureModel <- nobs.LocusModel
