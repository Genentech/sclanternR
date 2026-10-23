#' Cluster cells and select informative loci
#'
#' From an initial finite mixture model (fitting one if not provided),
#' selects loci with mutual information. Optionally, the model can be
#' refit with the selected loci. The refitting procedure can be
#' iteratively rerun.
#'
#' @inheritParams fitCellDistancePipeline
#' @param n_loci Positive integer or cutoff method (\code{"fdr"},
#'     \code{"mad"}, or \code{"elbow"}). Numeric selection retains the top
#'     loci by mutual information, capped to the available loci.
#' @param batch Optional batch labels, matched to input cells by name or,
#'     if unnamed, by position. Selects loci separately per batch and retains
#'     their union. Stored for downstream batch correction.
#' @param init A previously fitted \code{\linkS4class{FiniteMixtureModel}},
#'     which skips initial fitting, or \code{NULL} to fit an initial FMM.
#' @param reestimate_fmm Non-negative integer giving the number of FMM refits
#'     on selected loci (default 0). Logical \code{FALSE} and \code{TRUE}
#'     are accepted as 0 and 1. With 0, only initial fitting and selection
#'     are performed. With 1, the model is refit on the selected loci. With 2 or
#'     more, all loci are rescored and reselected before each additional refit.
#' @returns A named list of class \code{CellClustering} containing:
#'   \describe{
#'     \item{fmm}{The final FMM, fitted on selected loci when
#'         \code{reestimate_fmm > 0}, otherwise the initial FMM.}
#'     \item{keep_loci}{Names of the selected loci.}
#'     \item{mutualinfo}{Scores used for the final selection, for all input
#'         loci: a named vector, or a loci-by-batches matrix. These are
#'         computed before the final refit.}
#'     \item{mi_cutoff}{The final \code{MiCutoffEstimate}, or \code{NULL}
#'         for numeric or per-batch selection.}
#'     \item{batch}{Resolved batch factor, or \code{NULL}.}
#'     \item{reestimate_fmm}{Number of refits performed.}
#'   }
#' @seealso \code{\link{fitFiniteMixtureModel}}, \code{\link{mutualinfo}}
#' @include FiniteMixtureModel-class.R
#' @export
fitCellClustering <- function(
    cell_allele_counts,
    n_clusters = 20,
    n_loci = "fdr",
    batch = NULL,
    init = NULL,
    fmm_control = list(),
    n_loci_control = list(),
    reestimate_fmm = 0,
    verbose = FALSE
) {
    if (!is(cell_allele_counts, "CellAlleleCounts"))
        stop("cell_allele_counts must be a CellAlleleCounts object")
    reestimate_fmm <- .validate_reestimate_fmm(reestimate_fmm)
    n_loci <- .validate_n_loci(n_loci)
    if (!is.null(init) && !is(init, "FiniteMixtureModel"))
        stop("init must be NULL or a FiniteMixtureModel")

    if (is.null(init)) {
        .log_msg(1, verbose, "Fitting initial FMM with %d clusters", n_clusters)
        fmm_args <- list(cell_allele_cnt_list = cell_allele_counts,
                         n_clusters = n_clusters, verbose = verbose)
        init <- do.call(fitFiniteMixtureModel, modifyList(fmm_args, fmm_control))
    }

    if (!setequal(colnames(init), colnames(cell_allele_counts)))
        stop("init and cell_allele_counts must contain the same cells")
    if (!is.null(batch)) {
        batch <- .resolve_batch(batch, colnames(cell_allele_counts))
        names(batch) <- colnames(cell_allele_counts)
        batch <- .resolve_batch(batch, colnames(init))
    }
    cell_allele_counts <- cell_allele_counts[, colnames(init)]

    fmm <- init
    selection <- .select_fmm_loci(fmm, cell_allele_counts, n_loci,
                                  batch, n_loci_control)

    for (iteration in seq_len(reestimate_fmm)) {
        if (iteration > 1L) {
            selection <- .select_fmm_loci(fmm, cell_allele_counts, n_loci,
                                          batch, n_loci_control)
        }
        .log_msg(1, verbose, "Refitting FMM on %d loci (round %d of %d)",
                 length(selection$keep_loci), iteration, reestimate_fmm)
        fmm_args <- modifyList(list(verbose = verbose), fmm_control)
        # These define the refinement step and must not be overridden by controls.
        fmm_args$cell_allele_cnt_list <- cell_allele_counts[selection$keep_loci]
        fmm_args$init <- fmm@L
        fmm_args$combined_multinomial <- FALSE
        fmm <- do.call(fitFiniteMixtureModel, fmm_args)
    }
    .log_msg(1, verbose, "Selected %d informative loci", length(selection$keep_loci))
    structure(c(list(fmm = fmm), selection,
                list(batch = batch, reestimate_fmm = reestimate_fmm)),
              class = "CellClustering")
}

.select_fmm_loci <- function(fmm, counts, n_loci, batch, control) {
    select_one <- function(model, x, cutoff_control) {
        cutoff <- NULL
        if (is.character(n_loci)) {
            cutoff <- do.call(estimateMiCutoff,
                              c(list(model = model, x = x, method = n_loci),
                                cutoff_control))
            mi <- cutoff$mutualinfo
            selected <- cutoff$loci
        } else {
            mi <- mutualinfo(model, x)
            ranked <- sort(mi, decreasing = TRUE)
            selected <- names(ranked)[seq_len(min(n_loci, length(ranked)))]
        }
        list(keep_loci = selected, mutualinfo = mi, mi_cutoff = cutoff)
    }
    if (is.null(batch)) {
        selection <- select_one(fmm, counts, control)
    } else {
        # Preserve threshold-only selection within each batch by default.
        batch_control <- modifyList(list(min_loci = 0, max_loci = Inf), control)
        by_batch <- lapply(levels(batch), function(b) {
            cells <- colnames(fmm)[batch == b]
            select_one(fmm[, cells], counts[, cells], batch_control)
        })
        mi <- do.call(cbind, lapply(by_batch, `[[`, "mutualinfo"))
        colnames(mi) <- levels(batch)
        selection <- list(
            keep_loci = unique(unlist(lapply(by_batch, `[[`, "keep_loci"),
                                      use.names = FALSE)),
            mutualinfo = mi, mi_cutoff = NULL)
    }
    if (length(selection$keep_loci) == 0L) {
        stop("Automatic MI cutoff selected 0 loci. ",
             "Try a different method or specify n_loci as a number.")
    }
    selection
}

.validate_reestimate_fmm <- function(reestimate_fmm) {
    if (!(is.numeric(reestimate_fmm) || is.logical(reestimate_fmm)) ||
        is.complex(reestimate_fmm) || length(reestimate_fmm) != 1L ||
        is.na(reestimate_fmm) || !is.finite(reestimate_fmm) ||
        reestimate_fmm < 0 || reestimate_fmm > .Machine$integer.max ||
        reestimate_fmm != floor(reestimate_fmm)) {
        stop("reestimate_fmm must be a non-negative integer (or TRUE/FALSE)")
    }
    as.integer(reestimate_fmm)
}

.validate_n_loci <- function(n_loci) {
    if (is.character(n_loci)) {
        return(match.arg(n_loci, c("fdr", "mad", "elbow")))
    }
    if (!is.numeric(n_loci) || is.complex(n_loci) || length(n_loci) != 1L ||
        is.na(n_loci) || !is.finite(n_loci) || n_loci < 1 ||
        n_loci != floor(n_loci)) {
        stop("n_loci must be a positive integer or a method name (character)")
    }
    n_loci
}
