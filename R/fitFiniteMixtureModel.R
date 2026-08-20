
#' Fit a Finite Mixture Model (Latent Class Model) for single-cell data
#'
#' Fits a Finite Mixture Model where each cell belongs entirely to one of K
#' latent clones. Unlike the STRUCTURE model (where cells are admixtures of
#' topics), here a single global clone-prevalence vector pi drives assignments.
#'
#' Initialized from an L matrix (cell-topic loadings): pi is the column mean
#' of row-normalized L; F is computed from the data.
#'
#' @param cell_allele_cnt_list A CellAlleleCounts object.
#' @param init Initialization for the L matrix. Can be one of:
#'     \itemize{
#'       \item A matrix (N x K) of cell-topic loadings.
#'       \item A \code{\linkS4class{LocusModel}} object (e.g. a previously
#'             fitted \code{LocusNMF} or \code{FiniteMixtureModel}), from
#'             which the L matrix is extracted.
#'       \item A list with an \code{$L} element (e.g. a fastTopics model),
#'             from which the L matrix is extracted.
#'       \item \code{"random"}: L is initialized from a uniform Dirichlet
#'             distribution on the K-simplex, where K is set by
#'             \code{n_clusters}.
#'       \item \code{"fastTopics"}: an initial fastTopics model is fit
#'             and L is extracted from it.
#'       \item \code{"hclust"} (the default): hierarchical clustering on
#'             Euclidean distances of global allele frequencies is used
#'             to obtain initial cluster assignments, from which L is
#'             derived.
#'     }
#' @param n_clusters Number of clusters (K). Used when \code{init="hclust"},
#'     \code{init="random"}, or \code{init="fastTopics"}
#'     (default 20). Ignored when \code{init} supplies an L matrix.
#' @param loc_pseudocount Pseudocount per locus (default 1). Corresponds to
#'     a Dirichlet prior on alleles with shape `1+loc_pseudocount*global_allele_freq`
#' @param cell_pseudocount Pseudocount for global clone proportions (default 1).
#'     Corresponds to a Dirichlet prior on global proportions with shape
#'     `1+cell_pseudocount`
#' @param combined_multinomial One of \code{"init"} (default), \code{TRUE},
#'     or \code{FALSE}. When \code{"init"}, a combined-multinomial FMM is
#'     first fitted (using \code{init} for its initialization), then used
#'     to initialize the per-locus FMM. When \code{TRUE}, the model treats
#'     all features as one locus (allele frequencies sum to 1 across all
#'     features). When \code{FALSE}, per-locus constraints are used directly
#'     without a combined-multinomial initialization step.
#' @param n_init_cells Maximum number of cells to use for the default
#'     hierarchical clustering initialization (default 2000). If the dataset
#'     has more cells, a random subsample is drawn. Set to \code{NULL} to
#'     use all cells. Only used when \code{init="hclust"}.
#' @param maxiter Maximum number of EM iterations (default 100).
#' @param tol Convergence tolerance: stop when relative change in
#'     log-likelihood falls below this value (default 1e-4).
#' @param verbose Verbosity level: 0/FALSE for silent, 1/TRUE for milestones,
#'     2 for per-iteration log-likelihood, 3 for sub-step timing.
#'
#' @returns A \code{\linkS4class{FiniteMixtureModel}} object.
#'
#' @export
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom fastTopics fit_topic_model
#' @importFrom Matrix colSums Diagonal
#' @importFrom MatrixExtra t_shallow
fitFiniteMixtureModel <- function(
    cell_allele_cnt_list,
    n_clusters = 20,
    init = "hclust",
    loc_pseudocount  = 1,
    cell_pseudocount = 1,
    combined_multinomial = "init",
    n_init_cells = NULL,
    maxiter = 100,
    tol     = 1e-12,
    verbose = FALSE
) {
    if (!identical(combined_multinomial, "init") &&
        !identical(combined_multinomial, TRUE) &&
        !identical(combined_multinomial, FALSE)) {
        stop("combined_multinomial must be \"init\", TRUE, or FALSE")
    }

    if (identical(combined_multinomial, "init")) {
        .log_msg(1, verbose, "Fitting combined multinomial initialization")
        global_fmm <- fitFiniteMixtureModel(
            cell_allele_cnt_list,
            n_clusters = n_clusters, init = init,
            combined_multinomial = TRUE,
            loc_pseudocount = loc_pseudocount,
            cell_pseudocount = cell_pseudocount,
            n_init_cells = n_init_cells,
            maxiter = maxiter, tol = tol,
            verbose = verbose
        )
        init <- global_fmm
        combined_multinomial <- FALSE
    }

    Y_t <- MatrixExtra::t_shallow(as(cell_allele_cnt_list, 'CsparseMatrix'))
    if (combined_multinomial) {
        n_features <- ncol(Y_t)
        locus_1hot <- Matrix::Matrix(1, nrow = 1, ncol = n_features, sparse = TRUE)
        colnames(locus_1hot) <- colnames(Y_t)
        rownames(locus_1hot) <- "combined"
    } else {
        locus_1hot <- cell_allele_cnt_list@locus_1hot
    }

    feature_pseudocount <- Matrix::colSums(Y_t)
    feature_pseudocount_denom <- as.vector(t(locus_1hot) %*%
                                               (locus_1hot %*% feature_pseudocount))
    feature_pseudocount_denom[feature_pseudocount_denom == 0] <- 1
    feature_pseudocount <- loc_pseudocount * feature_pseudocount / feature_pseudocount_denom

    if (is.character(init) && init == "hclust") {
        .log_msg(1, verbose, "Computing hclust initialization")
        N_total <- nrow(Y_t)

        if (!is.null(n_init_cells) && N_total > n_init_cells) {
            sub_idx <- sample(N_total, n_init_cells)
            Y_t_sub <- Y_t[sub_idx, , drop = FALSE]
        } else {
            # distance matrix size grows quadratically, warn if >1GiB expected
            if (nrow(Y_t)^2 * 8 / 2^30 > 1) {
                warning(sprintf("%d x %d matrix to be constructed (consider setting n_init_cells)", nrow(Y_t), nrow(Y_t)))
            }

            Y_t_sub <- Y_t
        }
        N_sub <- nrow(Y_t_sub)

        rs <- Matrix::rowSums(Y_t_sub)
        rs[rs == 0] <- 1
        freq_mat <- Matrix::Diagonal(x = 1 / rs) %*% Y_t_sub

        crossp <- as.matrix(Matrix::tcrossprod(freq_mat))
        f2 <- computeF2dist(crossp)
        euc_dist <- sqrt(pmax(f2, 0))

        hc <- stats::hclust(stats::as.dist(euc_dist), method = "ward.D2")
        clusters <- stats::cutree(hc, k = n_clusters)

        L_sub <- matrix(0, nrow = N_sub, ncol = n_clusters)
        L_sub[cbind(seq_len(N_sub), clusters)] <- 1

        pi_init <- colSums(L_sub) + cell_pseudocount
        pi_init <- pi_init / sum(pi_init)

        Fmat_init <- .update_Fmat(Y_t_sub, L_sub, feature_pseudocount, locus_1hot)

        L_init <- .fmm_e_step(Y_t, Fmat_init, pi_init)

        .log_msg(1, verbose, "Finished hclust initialization")
    } else if (is.character(init) && init == "fastTopics") {
        .log_msg(1, verbose, "Fitting initial fastTopics model")
        mat <- as(cell_allele_cnt_list, 'CsparseMatrix')

        if (verbose >= 1) {
            ft_verbose <- "progressbar"
        } else {
            ft_verbose <- "none"
        }
        ft_model <- fit_topic_model(
            t(mat), k = n_clusters,
            init.method = 'topicscore',
            verbose = ft_verbose
        )
        .log_msg(1, verbose, "Finished fitting initial fastTopics model")
        L_init <- ft_model$L
    } else if (is.character(init) && init == "random") {
        N <- ncol(cell_allele_cnt_list)
        L_init <- .random_simplex_matrix(N, n_clusters)
    } else if (is.matrix(init)) {
        L_init <- init
    } else if (is(init, "LocusModel")) {
        L_init <- init@L
    } else if (is.list(init) && !is.null(init$L)) {
        L_init <- init$L
    } else {
        stop("'init' must be \"hclust\", \"fastTopics\", \"random\", a matrix, a LocusModel, or a list with $L")
    }

    .log_msg(1, verbose, "Initializing FMM")

    L <- sweep(L_init, 1, rowSums(L_init), '/')
    pi <- colMeans(L)
    pi <- pi / sum(pi)
    Fmat <- .update_Fmat(Y_t, L, feature_pseudocount, locus_1hot)

    .log_msg(1, verbose, "Finished initializing FMM")

    N <- nrow(L_init)
    K <- ncol(L_init)
    loglik_trace <- numeric(maxiter)

    locus_names <- names(cell_allele_cnt_list)

    .log_msg(1, verbose, "Starting EM (%d loci, %d iterations max)",
             length(locus_names), maxiter)

    for (iter in seq_len(maxiter)) {
        .log_msg(3, verbose, "Iter %d: E-step start", iter)

        Fmat <- pmax(Fmat, 1e-12)

        log_L <- matrix(log(pi), nrow=N, ncol=K, byrow=TRUE) +
            as.matrix(Y_t %*% log(Fmat))

        .log_msg(3, verbose, "Iter %d: E-step done", iter)

        m           <- apply(log_L, 1, max)
        exp_shifted <- exp(log_L - m)
        row_sums    <- rowSums(exp_shifted)
        ll          <- sum(m + log(row_sums))

        if (cell_pseudocount > 0) {
            ll <- ll + sum(log(pi) * cell_pseudocount)
        }

        if (loc_pseudocount > 0) {
            ll <- ll + sum(log(Fmat) * feature_pseudocount)
        }

        .log_msg(3, verbose, "Iter %d: M-step start", iter)

        L         <- exp_shifted / row_sums

        pi_stat <- colSums(L) + cell_pseudocount
        pi      <- pi_stat / sum(pi_stat)

        Fmat <- .update_Fmat(Y_t, L, feature_pseudocount, locus_1hot)

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
    if (combined_multinomial) {
        slots$locus_1hot <- locus_1hot
    }
    model <- do.call(new, c("FiniteMixtureModel", slots, list(pi = pi)))

    model@L <- predict(model, newx = cell_allele_cnt_list, type = "loadings")

    model
}

.random_simplex_matrix <- function(N, K) {
    L <- matrix(stats::rgamma(N * K, shape = 1, rate = 1), nrow = N, ncol = K)
    sweep(L, 1, rowSums(L), '/')
}

.update_Fmat <- function(Y_t, L, feature_pseudocount, locus_1hot) {
    Fmat_stat <- as.matrix(t(Y_t) %*% L)
    Fmat_stat <- sweep(Fmat_stat, 1, feature_pseudocount, '+')
    Fmat_stat_denom <- as.matrix(t(locus_1hot) %*% (locus_1hot %*% Fmat_stat))
    Fmat_stat_denom[Fmat_stat_denom == 0] <- 1
    Fmat_stat / Fmat_stat_denom
}
