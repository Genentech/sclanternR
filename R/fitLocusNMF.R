
#' Fit a locus-aware NMF model
#'
#' Refines a fastTopics model using an EM algorithm that enforces per-locus
#' simplex constraints: allele frequencies sum to 1 within each locus
#' (rather than across all features as in the fastTopics bag-of-words model).
#' Each cell's allele frequencies at a locus are a linear combination of the
#' topic-specific frequencies, weighted by the cell's topic proportions.
#'
#' Initialized from a fastTopics or LocusModel object by taking L
#' (row-normalized) and F renormalized within each locus.
#'
#' @param cell_allele_cnt_list A CellAlleleCounts object.
#' @param init Initialization for L and F. Can be one of:
#'     \itemize{
#'       \item A \code{\linkS4class{LocusModel}} object (e.g. a previously
#'             fitted \code{LocusNMF} or \code{FiniteMixtureModel}).
#'       \item A fastTopics model (list with \code{$F} and \code{$L} matrices,
#'             as returned by \code{fit_topic_model}).
#'       \item \code{NULL} (the default): an initial fastTopics model is fit
#'             automatically on \code{cell_allele_cnt_list}.
#'     }
#'     All allele names in \code{cell_allele_cnt_list} must appear in the
#'     initialization model's feature names.
#' @param n_topics Number of topics to use when \code{init=NULL} is auto-fit
#'     (default 10). Ignored when \code{init} is supplied.
#' @param loc_pseudocount Pseudocount per locus (default 1). Corresponds to
#'     a Dirichlet prior on alleles with shape `1+loc_pseudocount*global_allele_freq`
#' @param cell_pseudocount Pseudocount for clone proportions (default .01).
#'     Corresponds to a Dirichlet prior on topic-loadings with shape
#'     `1+cell_pseudocount`
#' @param maxiter Maximum number of EM iterations (default 100).
#' @param tol Convergence tolerance: stop when relative change in
#'     log-likelihood falls below this value (default 1e-4).
#' @param update_final_L If TRUE, recompute L (cell-topic loadings) after convergence.
#'     If FALSE, take L from the final EM step.
#' @param verbose Verbosity level: 0/FALSE for silent, 1/TRUE for milestones,
#'     2 for per-iteration log-likelihood, 3 for sub-step timing.
#'
#' @returns A \code{\linkS4class{LocusNMF}} object.
#'
#' @export
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom fastTopics fit_topic_model
#' @importFrom MatrixExtra t_shallow
fitLocusNMF <- function(
    cell_allele_cnt_list,
    init = NULL,
    n_topics = 10,
    loc_pseudocount  = 1,
    cell_pseudocount = .01,
    maxiter = 100,
    tol     = 1e-12,
    update_final_L = FALSE,
    verbose = FALSE
) {
    if (is.null(init)) {
        .log_msg(1, verbose, "Fitting initial fastTopics model")
        mat      <- as(cell_allele_cnt_list, 'CsparseMatrix')
        if (verbose >= 1) {
            ft_verbose <- "progressbar"
        } else {
            ft_verbose <- "none"
        }
        ft_model <- fit_topic_model(
            t(mat), k = n_topics,
            init.method = 'topicscore', numiter.main = 3, numiter.refine = 0,
            verbose = ft_verbose
        )
        .log_msg(1, verbose, "Finished fitting initial fastTopics model")
    } else {
        ft_model <- .as_ft_init(init)
        .log_msg(1, verbose, "Sanitizing features")
        drop_result <- .drop_loci_absent_from_ft(cell_allele_cnt_list, ft_model)
        cell_allele_cnt_list <- drop_result$cac
        .log_msg(1, verbose, "Finished sanitizing features")
    }

    .log_msg(1, verbose, "Initializing model")
    Y_t <- t(as(cell_allele_cnt_list, 'CsparseMatrix'))
    locus_1hot <- cell_allele_cnt_list@locus_1hot

    feature_pseudocount <- Matrix::colSums(Y_t)
    feature_pseudocount_denom <- t(locus_1hot) %*% (locus_1hot %*% feature_pseudocount)
    feature_pseudocount <- loc_pseudocount * feature_pseudocount / feature_pseudocount_denom

    nmf_init <- .init_structure_from_fasttopics(ft_model, cell_allele_cnt_list)
    L    <- nmf_init$L
    Fmat <- nmf_init$Fmat

    N <- nrow(L)
    K <- ncol(L)
    loglik_trace <- numeric(maxiter)

    .log_msg(1, verbose, "Starting EM (%d loci, %d iterations max)",
             length(cell_allele_cnt_list), maxiter)

    for (iter in seq_len(maxiter)) {
        .log_msg(3, verbose, "Iter %d: E-step start", iter)

        expected_sddmm <- sddmm_csc(Y_t, L, t(Fmat))
        expected_sddmm@x <- pmax(expected_sddmm@x, 1e-12)

        denom_sddmm <- expected_sddmm
        denom_sddmm@x <- 1 / denom_sddmm@x
        R_t <- Y_t * denom_sddmm

        log_sddmm <- expected_sddmm
        log_sddmm@x <- log(log_sddmm@x)
        ll <- sum(Y_t * log_sddmm)

        if (cell_pseudocount > 0) {
            ll <- ll + sum(log(L) * cell_pseudocount)
        }

        if (loc_pseudocount > 0) {
            ll <- ll + sum(log(pmax(Fmat, 1e-12)) * feature_pseudocount)
        }

        L_stat <- L * as.matrix(R_t %*% Fmat)

        Fmat_stat <- Fmat * as.matrix(t(R_t) %*% L)
        Fmat_stat <- sweep(Fmat_stat, 1, feature_pseudocount, '+')
        Fmat_stat_denom <- as.matrix(t(locus_1hot) %*% (locus_1hot %*% Fmat_stat))
        Fmat <- Fmat_stat / Fmat_stat_denom

        .log_msg(3, verbose, "Iter %d: E-step done", iter)

        L_stat <- L_stat + cell_pseudocount
        L      <- sweep(L_stat, 1, rowSums(L_stat), '/')
        .log_msg(3, verbose, "Iter %d: M-step done", iter)

        loglik_trace[iter] <- ll
        .log_msg(2, verbose, "Iter %d: log-lik = %.4f", iter, ll)

        if (iter > 1) {
            rel_change <- abs(ll - loglik_trace[iter - 1]) /
                          (abs(loglik_trace[iter - 1]) + 1e-12)
            if (rel_change < tol) {
                loglik_trace <- loglik_trace[seq_len(iter)]
                break
            }
        }
    }

    .log_msg(1, verbose, "Building model object")

    locus_map <- cell_allele_cnt_list@locus_map
    cell_names <- colnames(cell_allele_cnt_list)

    slots <- .make_locus_model_slots(
        Fmat = Fmat, L = L, loglik = loglik_trace,
        locus_map = locus_map, cell_names = cell_names,
        cell_pseudocount = cell_pseudocount,
        loc_pseudocount = loc_pseudocount,
        nobs = sum(Y_t@x)
    )
    model <- do.call(new, c("LocusNMF", slots))

    if (update_final_L) {
        model@L <- .predict_loadings_nmf(model, cell_allele_cnt_list, init_L = L)
    }

    model
}

.as_ft_init <- function(init) {
    if (is(init, "LocusModel")) {
        F_mat <- init@Fmat
        rownames(F_mat) <- as.character(init@locus_map$feature_name)
        L_mat <- init@L
        list(F = F_mat, L = L_mat)
    } else {
        init
    }
}

.init_structure_from_fasttopics <- function(ft_model, cell_allele_cnt_list) {
    L <- ft_model$L
    L <- sweep(L, 1, rowSums(L), '/')

    Fmat <- .init_Fmat_from_fasttopics(ft_model, cell_allele_cnt_list, as_list=FALSE)

    list(L=L, Fmat=Fmat)
}

# Build per-locus F matrix from ft_model$F in one vectorized pass.
# Returns F_features x K matrix (same orientation as fastTopics F).
.init_Fmat_from_fasttopics <- function(ft_model, cell_allele_cnt_list, as_list=TRUE) {
    K           <- ncol(ft_model$F)
    topic_names <- colnames(ft_model$F)
    feature_names <- as.character(cell_allele_cnt_list@locus_map$feature_name)
    n_features <- length(feature_names)

    ft_idx   <- match(feature_names, rownames(ft_model$F))
    has_match <- !is.na(ft_idx)

    Fmat <- matrix(.Machine$double.eps, nrow = n_features, ncol = K,
                   dimnames = list(feature_names, topic_names))
    Fmat[has_match, ] <- ft_model$F[ft_idx[has_match], , drop = FALSE]

    locus_idx <- cell_allele_cnt_list@locus_index
    locus_1hot <- cell_allele_cnt_list@locus_1hot

    locus_colsums <- locus_1hot %*% Fmat  # L x K
    normalizer <- t(locus_1hot) %*% locus_colsums  # F x K

    Fmat <- Fmat / as.matrix(normalizer)

    if (as_list) {
        lapply(locus_idx, function(idx) Fmat[idx, , drop = FALSE])
    } else {
        Fmat
    }
}

# Vectorized check: drop loci where no feature appears in ft_model$F.
.drop_loci_absent_from_ft <- function(cell_allele_cnt_list, ft_model) {
    feature_in_ft <- as.character(cell_allele_cnt_list@locus_map$feature_name) %in%
                     rownames(ft_model$F)
    loci_with_feature <- unique(
        as.character(cell_allele_cnt_list@locus_map$locus[feature_in_ft]))
    all_loci <- names(cell_allele_cnt_list)
    in_model <- all_loci %in% loci_with_feature

    if (!all(in_model)) {
        message(sprintf(
            "Dropping %d loci not present in init (keeping %d). ",
            sum(!in_model), sum(in_model)),
            "Subset cell_allele_cnt_list to the model's loci to suppress this.")
        cell_allele_cnt_list <- cell_allele_cnt_list[all_loci[in_model]]
    }

    list(cac = cell_allele_cnt_list)
}
