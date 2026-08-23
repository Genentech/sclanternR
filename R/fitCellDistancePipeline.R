
#' Fit a cell distance pipeline via FMM and CLR-PCA
#'
#' Multi-step pipeline: (1) fit a FiniteMixtureModel, (2) perform CLR-PCA on
#' topic frequencies and replace with low-rank projection, (3) select top loci
#' by mutual information, (4) shrink cell-level frequencies toward model
#' predictions, (5) compute PCA on shrunken frequencies, and optionally (6)
#' compute pairwise F2 distances.
#'
#' When batch labels are provided, LDA-based batch correction is applied twice:
#' once to the CLR-PCA topic coordinates (removing batch discriminant directions
#' before reconstructing the denoised frequency matrix), and once to the final
#' cell PCA coordinates. Loci are also selected per-batch (union of top loci).
#'
#' The pipeline supports partial re-entry: pass a previously fitted FMM via
#' \code{init} to skip step 1, manually specify \code{n_pcs} to override
#' Marchenko-Pastur, or change \code{n_loci} on a rerun.
#'
#' @param cell_allele_counts A \code{\linkS4class{CellAlleleCounts}} object.
#' @param n_clusters Number of clusters for the initial FMM (default 20).
#' @param n_loci Integer or character. When numeric, retains the top \code{n_loci}
#'     loci by mutual information (capped to available loci). When a character
#'     string, names a method from \code{\link{estimateMiCutoff}} (\code{"fdr"},
#'     \code{"mad"}, \code{"elbow"}) and retains all loci
#'     exceeding the estimated threshold. Ignored if \code{keep_loci} non-NULL.
#' @param keep_loci Loci to keep for the shrunken frequencies. Overrides `n_loci`
#'     when non-null.
#' @param pseudocount Numeric shrinkage parameter or \code{"auto"} (default 1).
#'     When numeric, cell-level frequencies are shrunk toward model predictions
#'     using this as the Dirichlet prior weight. When \code{"auto"}, the
#'     concentration is estimated from the data via
#'     \code{\link{estimateDirichletPseudocount}}.
#' @param n_pcs Number of principal components, or a character string naming a
#'     method from \code{\link{estimateFmmDimensionality}} (default
#'     \code{NULL} uses \code{"gavish_donoho"}). When automatic, a minimum
#'     of 2 dimensions is enforced (or \code{2 + nlevels(batch) - 1} when
#'     batch correction is active) to prevent LDA from collapsing
#'     dimensionality.
#' @param n_pcs2 Number of principal components for the 2nd round of PCA on
#'     after computing shrunken frequencies. By default, same value is used as n_pcs.
#' @param batch Optional character or factor vector of batch labels, one per
#'     cell. If named, matched to cells by name; if unnamed, matched by
#'     position. When provided, LDA batch correction is applied.
#' @param init A previously fitted \code{\linkS4class{FiniteMixtureModel}},
#'     a previous pipeline result (list with class
#'     \code{"CellDistancePipeline"}), or a fastTopics model (list with
#'     \code{$F} and \code{$L} matrices). When a FiniteMixtureModel or
#'     pipeline result is provided, FMM fitting is skipped entirely. When a
#'     fastTopics model is provided, it is used to initialize
#'     \code{\link{fitFiniteMixtureModel}}.
#' @param compute_dist Logical. If \code{TRUE} (default), compute the full
#'     N x N pairwise distance matrix. Set to \code{FALSE} for large datasets
#'     to return only the cell PCA coordinates.
#' @param fmm_control Named list of additional arguments passed to
#'     \code{\link{fitFiniteMixtureModel}}.
#' @param n_loci_control Named list of additional arguments passed to
#'     \code{\link{estimateMiCutoff}}.
#' @param reestimate_loadings Logical or \code{"auto"} (default \code{TRUE}).
#'     If \code{TRUE}, cell loadings are re-estimated on the informative loci
#'     only before predicting frequencies. If \code{FALSE}, the full-model
#'     loadings are used. If \code{"auto"}, both options are compared via
#'     leave-one-out cross-validation on the loci, with the higher-scoring option selected.
#' @param verbose Verbosity level (default \code{FALSE}).
#'
#' @returns A named list with S3 class \code{"CellDistancePipeline"} containing:
#'   \describe{
#'     \item{fmm}{The fitted \code{\linkS4class{FiniteMixtureModel}}.}
#'     \item{scree}{An \code{FmmDimEstimate} object with all eigenvalues from
#'         the CLR-PCA (for screeplot inspection).}
#'     \item{mutualinfo}{Named numeric vector of per-locus MI scores (unsorted
#'         when \code{n_loci} is character; sorted descending when numeric),
#'         or a matrix (loci x batches) when batch labels are provided.}
#'     \item{mi_cutoff}{An \code{MiCutoffEstimate} object (see
#'         \code{\link{estimateMiCutoff}}), or \code{NULL} when
#'         \code{n_loci} is numeric.}
#'     \item{pseudocount}{The pseudocount value used for shrinkage (numeric
#'         scalar). When \code{pseudocount = "auto"}, this is the estimated
#'         value.}
#'     \item{pseudocount_estimate}{A \code{DirichletPseudocountEstimate} object
#'         (see \code{\link{estimateDirichletPseudocount}}), or \code{NULL} when
#'         a numeric pseudocount was provided.}
#'     \item{cell_pca}{N x n_pcs matrix of final cell PCA coordinates (always
#'         returned). Can be used to compute distances downstream via
#'         \code{\link{computeF2dist}(tcrossprod(cell_pca))}.}
#'     \item{dist}{N x N matrix of squared Euclidean (F2) distances, or
#'         \code{NULL} if \code{compute_dist = FALSE}.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Full pipeline:
#' result <- fitCellDistancePipeline(counts, n_clusters = 20, n_loci = 100)
#'
#' # Inspect screeplot, then rerun with manual n_pcs:
#' result1 <- fitCellDistancePipeline(counts, n_loci = 100)
#' plot(result1)
#' result2 <- fitCellDistancePipeline(counts, n_loci = 100,
#'                                    init = result1, n_pcs = 5)
#'
#' # With batch correction:
#' batch_vec <- setNames(cell_data$batch, cell_data$cell_id)
#' result <- fitCellDistancePipeline(counts, n_loci = 100, batch = batch_vec)
#'
#' # Large dataset: skip distance matrix, compute downstream
#' result <- fitCellDistancePipeline(counts, n_loci = 100, compute_dist = FALSE)
#' dist_mat <- computeF2dist(tcrossprod(result$cell_pca))
#' }
#'
#' @seealso \code{\link{fitFiniteMixtureModel}},
#'     \code{\link{estimateFmmDimensionality}}, \code{\link{mutualinfo}},
#'     \code{\link{computeF2dist}}
#'
#' @include FiniteMixtureModel-class.R
#' @include estimateDirichletPseudocount.R
#' @include estimateFmmDim.R
#' @include estimateMiCutoff.R
#' @export
fitCellDistancePipeline <- function(
    cell_allele_counts,
    n_clusters = 20,
    n_loci = "fdr",
    keep_loci = NULL,
    pseudocount = 1,
    n_pcs = NULL,
    n_pcs2 = NULL,
    batch = NULL,
    init = NULL,
    compute_dist = TRUE,
    fmm_control = list(),
    n_loci_control = list(),
    reestimate_loadings = TRUE,
    verbose = FALSE
) {
    if (!is(cell_allele_counts, "CellAlleleCounts"))
        stop("cell_allele_counts must be a CellAlleleCounts object")
    auto_pseudocount <- identical(pseudocount, "auto")
    if (!auto_pseudocount && (!is.numeric(pseudocount) || length(pseudocount) != 1
                              || pseudocount <= 0)) {
        stop("pseudocount must be a positive number or \"auto\"")
    }
    mi_method <- NULL
    if (is.character(n_loci)) {
        mi_method <- match.arg(n_loci, c("fdr", "mad", "elbow"))
    } else if (!is.null(n_loci) && (!is.numeric(n_loci) || length(n_loci) != 1 || n_loci < 1)) {
        stop("n_loci must be a positive integer or a method name (character)")
    }
    if (is.null(n_loci) && is.null(keep_loci)) {
        stop("n_loci must be specified or keep_loci must be non-NULL")
    }
    auto_loadings <- identical(reestimate_loadings, "auto")
    if (!auto_loadings && !is.logical(reestimate_loadings)) {
        stop("reestimate_loadings must be TRUE, FALSE, or \"auto\"")
    }

    # --- Step 1: Fit FMM ---
    if (is.null(init)) {
        .log_msg(1, verbose, "Step 1: Fitting FMM with %d clusters", n_clusters)
        fmm_args <- modifyList(
            list(cell_allele_cnt_list = cell_allele_counts,
                 n_clusters = n_clusters, verbose = verbose),
            fmm_control
        )
        fmm <- do.call(fitFiniteMixtureModel, fmm_args)
    } else if (inherits(init, "CellDistancePipeline")) {
        fmm <- init$init_fmm
        .log_msg(1, verbose, "Step 1: Using FMM from previous pipeline result")
    } else if (is(init, "FiniteMixtureModel")) {
        fmm <- init
        .log_msg(1, verbose, "Step 1: Using provided FMM")
    } else if (is.list(init) && !is.null(init$F) && !is.null(init$L)) {
        .log_msg(1, verbose,
                 "Step 1: Fitting FMM with %d clusters (fastTopics init)",
                 n_clusters)
        fmm_args <- modifyList(
            list(cell_allele_cnt_list = cell_allele_counts,
                 init = init, n_clusters = n_clusters, verbose = verbose),
            fmm_control
        )
        fmm <- do.call(fitFiniteMixtureModel, fmm_args)
    } else {
        stop("init must be NULL, a FiniteMixtureModel, a CellDistancePipeline result, or a fastTopics model")
    }

    cell_names <- colnames(fmm)

    if (!is.null(batch)) {
        batch <- .resolve_batch(batch, cell_names)
    }

    # --- Step 2: CLR-PCA on topic frequencies ---
    .log_msg(1, verbose, "Step 2: CLR-PCA on topic frequencies")
    pca <- .clr_weighted_pca(fmm@Fmat, fmm@pi, fmm@locus_1hot)

    if (is.null(n_pcs) || is.character(n_pcs)) {
        method <- if (is.character(n_pcs)) n_pcs else "gavish_donoho"
        scree <- estimateFmmDimensionality(fmm, method = method)
        min_pcs <- if (is.null(batch)) 2L else 1L + nlevels(batch)
        n_pcs <- max(min_pcs, scree$cutoff_index)
    } else {
        scree <- NULL
    }
    .log_msg(1, verbose, "Using %d PCs", n_pcs)

    # Topic PCA coordinates: K x n_pcs
    r <- min(n_pcs, length(pca$d))
    M <- pca$V[, seq_len(r), drop = FALSE] *
         rep(pca$d[seq_len(r)], each = nrow(pca$V))
    U_r <- pca$U[, seq_len(r), drop = FALSE]

    # --- Step 2b: LDA correction on CLR-PCA (if batch) ---
    if (!is.null(batch)) {
        .log_msg(1, verbose, "Applying LDA batch correction to CLR-PCA")
        X_cell <- loadings(fmm) %*% M
        Q <- .compute_lda_basis(X_cell, batch)
        if (!is.null(Q)) {
            M <- M - M %*% tcrossprod(Q)
        }
    }

    # Reconstruct Fmat from (possibly corrected) low-rank CLR
    Fmat_proj <- .reconstruct_fmat(U_r, M, pca$mu, pca$pi_vec, pca$locus_1hot,
                                   pca$keep, fmm@locus_1hot)
    projected_fmm <- .replace_fmm_fmat(fmm, Fmat_proj)

    # --- Step 3: Compute MI and select top loci ---

    #mi_fmm <- projected_fmm
    # use original FMM to be consistent with user manually selecting cutoff from it
    mi_fmm <- fmm
    mi_cutoff <- NULL

    if (is.null(batch)) {
        .log_msg(1, verbose, "Step 3: Computing mutual information across %d loci",
                 length(mi_fmm))
        if (is.null(keep_loci) && !is.null(mi_method)) {
            mi_cutoff <- estimateMiCutoff(mi_fmm, cell_allele_counts, method = mi_method)
            mi <- mi_cutoff$mutualinfo
            keep_loci <- mi_cutoff$loci
            .log_msg(1, verbose,
                     "Automatic MI cutoff (%s): %.4g, %d loci selected",
                     mi_method, mi_cutoff$cutoff, length(keep_loci))
        } else {
            mi <- mutualinfo(mi_fmm, cell_allele_counts)
            if (is.null(keep_loci)) {
                mi_sort <- sort(mi, decreasing = TRUE)
                keep_loci <- names(mi_sort)[seq_len(min(n_loci, length(mi_sort)))]
            }
        }
    } else {
        .log_msg(1, verbose,
                 "Step 3: Computing per-batch mutual information (%d batches)",
                 nlevels(batch))

        mi <- sapply(levels(batch), function(b) {
            cells_b <- cell_names[batch == b]
            mutualinfo(mi_fmm[, cells_b], cell_allele_counts[, cells_b])
        })
        if (is.null(keep_loci)) {
            if (!is.null(mi_method)) {
                n_alleles <- as.integer(Matrix::rowSums(mi_fmm@locus_1hot))
                K <- ncol(mi_fmm@Fmat)
                keep_loci <- lapply(colnames(mi), function(b) {
                    mi_col <- mi[, b]
                    res <- .apply_mi_cutoff(mi_col, mi_method,
                                           n_alleles, K)
                    names(mi_col)[mi_col >= res$cutoff]
                })
                keep_loci <- unique(unlist(keep_loci))
            } else {
                keep_loci <- apply(mi, 2, function(x) {
                    names(sort(x, decreasing=TRUE))[seq_len(n_loci)]
                }, simplify=FALSE)
                keep_loci <- unique(unlist(keep_loci))
            }
        }
    }
    if (length(keep_loci) == 0L) {
        stop("Automatic MI cutoff selected 0 loci. ",
             "Try a different method or specify n_loci as a number.")
    }
    .log_msg(1, verbose, "Selected %d informative loci", length(keep_loci))

    # --- Subset model and counts to selected loci ---
    model_sub <- projected_fmm[keep_loci]
    cac_sub <- cell_allele_counts[keep_loci]

    if (auto_loadings) {
        # NOTE: It seems to work better to use the original FMM than projected
        # model when deciding whether the original loadings are informative.
        # (TODO: Check how projected model minus LDA would do?)
        ll_diff <- .loo_loci_loglik(fmm[keep_loci], cac_sub, fmm@L)

        if (ll_diff >= 0) {
            cell_loadings <- .predict_loadings(model_sub, cac_sub)
            loadings_method <- "auto:reestimated"
        } else {
            cell_loadings <- model_sub@L
            loadings_method <- "auto:full"
        }
        .log_msg(1, verbose, "Loadings: %s (mean held-out LL diff = %.2f)",
                 loadings_method, ll_diff)
    } else if (reestimate_loadings) {
        cell_loadings <- .predict_loadings(model_sub, cac_sub)
        loadings_method <- "reestimated"
        ll_diff <- NULL
    } else {
        cell_loadings <- model_sub@L
        loadings_method <- "full"
        ll_diff <- NULL
    }

    # --- Step 3b: Estimate pseudocount (if auto) ---
    pc_estimate <- NULL
    if (auto_pseudocount) {
        pc_estimate <- estimateDirichletPseudocount(model_sub, cac_sub,
                                                    loadings = cell_loadings,
                                                    verbose = verbose)
        pseudocount <- pc_estimate$pseudocount
    }

    # --- Step 4: Predict shrunken frequencies ---
    .log_msg(1, verbose, "Step 4: Predicting frequencies (pseudocount = %g)",
             pseudocount)
    freq_list <- predict(model_sub, newx = cac_sub, pseudocount = pseudocount,
                         loadings = cell_loadings)
    freq_mat <- t(do.call(rbind, freq_list))

    # --- Step 5: PCA on shrunken frequencies ---
    if (is.null(n_pcs2)) {
        n_pcs2 <- n_pcs
    }
    N <- nrow(freq_mat)
    P <- ncol(freq_mat)
    n_pcs_freq <- min(n_pcs2, N - 1L, P)

    col_means <- colMeans(freq_mat)
    freq_centered <- sweep(freq_mat, 2, col_means, '-')

    if (n_pcs_freq < P) {
        .log_msg(1, verbose, "Step 5: PCA: reducing %d features to %d PCs",
                 P, n_pcs_freq)
        sv <- svd(freq_centered, nu = n_pcs_freq, nv = 0)
        coords <- sv$u[, seq_len(n_pcs_freq), drop = FALSE] *
                  rep(sv$d[seq_len(n_pcs_freq)], each = N)
    } else {
        coords <- freq_centered
    }

    # --- Step 5b: LDA batch correction on frequency PCA ---
    if (!is.null(batch)) {
        coords <- .lda_batch_correct(coords, batch, verbose)
    }

    rownames(coords) <- cell_names
    colnames(coords) <- paste0("PC", 1:ncol(coords))

    # --- Step 6: Compute distances ---
    dist_mat <- NULL
    if (compute_dist) {
        .log_msg(1, verbose, "Step 6: Computing pairwise distances")
        if (nrow(coords)^2 * 8 / 2^30 > 1) {
            warning(sprintf("%d x %d matrix to be constructed (consider setting compute_dist=FALSE)", nrow(coords), nrow(coords)))
        }
        inner_prod <- tcrossprod(coords)
        dist_mat <- computeF2dist(inner_prod)
        rownames(dist_mat) <- cell_names
        colnames(dist_mat) <- cell_names
    }

    result <- list(
        init_fmm = fmm,
        final_fmm = model_sub,
        projected_fmm = projected_fmm,
        scree = scree,
        mutualinfo = mi,
        mi_cutoff = mi_cutoff,
        pseudocount = pseudocount,
        pseudocount_estimate = pc_estimate,
        cell_pca = coords,
        dist = dist_mat,
        freq_list = freq_list,
        cluster_loadings = cell_loadings,
        loadings_choice = list(method = loadings_method, ll_diff = ll_diff)
    )
    class(result) <- "CellDistancePipeline"
    result
}

# --- Helper: reconstruct Fmat from CLR-PCA components ---

.reconstruct_fmat <- function(U_r, M, mu, pi_vec, locus_1hot, keep, orig_locus_1hot) {
    clr_weighted <- U_r %*% t(M)
    clr_centered <- sweep(clr_weighted, 2, sqrt(pi_vec), '/')
    clr_approx <- clr_centered + mu
    fmat_sub <- exp(clr_approx)
    locus_sums <- as.matrix(Matrix::t(locus_1hot) %*% locus_1hot %*% fmat_sub)
    fmat_sub <- fmat_sub / locus_sums

    if (all(keep)) {
        fmat <- fmat_sub
    } else {
        P_orig <- length(keep)
        K <- ncol(fmat_sub)
        fmat <- matrix(0, nrow = P_orig, ncol = K)
        fmat[keep, ] <- fmat_sub
    }

    rownames(fmat) <- names(keep)

    fmat
}

# --- Helper: replace Fmat in an existing FMM ---

.replace_fmm_fmat <- function(fmm, new_Fmat) {
    slots <- .make_locus_model_slots(
        Fmat = new_Fmat, L = fmm@L, loglik = fmm@loglik,
        locus_map = fmm@locus_map, cell_names = fmm@cell_names,
        cell_pseudocount = fmm@cell_pseudocount,
        loc_pseudocount = fmm@loc_pseudocount,
        nobs = fmm@nobs
    )
    do.call(new, c("FiniteMixtureModel", slots, list(pi = fmm@pi)))
}

# --- LOO-CV helper for auto loadings selection ---

.loo_loci_loglik <- function(model, x, loadings_full) {
    Y <- as(x, "CsparseMatrix")
    Y_t <- t(Y)
    Fmat <- model@Fmat
    log_Fmat <- log(pmax(Fmat, 1e-12))
    pi_vec <- model@pi
    N <- ncol(Y)
    K <- length(pi_vec)
    log_pi <- matrix(log(pi_vec), nrow = N, ncol = K, byrow = TRUE)

    S_by_locus <- lapply(model@locus_index, function(idx) {
        as.matrix(Y_t[, idx, drop = FALSE] %*% log_Fmat[idx, , drop = FALSE])
    })
    S_full <- Reduce("+", S_by_locus)

    pred_full <- t(loadings_full %*% t(Fmat))

    ll_diff <- 0
    for (locus in names(model@locus_index)) {
        allele_idx <- model@locus_index[[locus]]

        log_L <- log_pi + S_full - S_by_locus[[locus]]
        m <- apply(log_L, 1, max)
        exp_shifted <- exp(log_L - m)
        L_reest <- exp_shifted / rowSums(exp_shifted)

        Fmat_l <- Fmat[allele_idx, , drop = FALSE]
        pred_reest_l <- t(L_reest %*% t(Fmat_l))
        pred_full_l <- pred_full[allele_idx, , drop = FALSE]

        Y_l <- as.matrix(Y[allele_idx,])
        ll_diff <- ll_diff +
            sum(Y_l * log(pmax(pred_reest_l, 1e-300))) -
            sum(Y_l * log(pmax(pred_full_l, 1e-300)))
    }
    ll_diff
}

# --- Batch correction helpers ---

.resolve_batch <- function(batch, cell_names) {
    if (!is.null(names(batch))) {
        idx <- match(cell_names, names(batch))
        if (anyNA(idx))
            stop("batch vector is missing names for some cells")
        batch <- batch[idx]
    } else {
        if (length(batch) != length(cell_names))
            stop("batch must have one entry per cell")
    }
    as.factor(batch)
}

.compute_lda_basis <- function(X, batch) {
    batch <- as.factor(batch)
    groups <- levels(batch)
    G <- length(groups)
    if (G <= 1L) return(NULL)

    p <- ncol(X)
    global_mean <- colMeans(X)

    S_w <- matrix(0, p, p)
    for (g in groups) {
        mask <- batch == g
        X_g <- X[mask, , drop = FALSE]
        X_g_c <- sweep(X_g, 2, colMeans(X_g), '-')
        S_w <- S_w + crossprod(X_g_c)
    }

    S_b <- matrix(0, p, p)
    for (g in groups) {
        mask <- batch == g
        n_g <- sum(mask)
        diff_g <- colMeans(X[mask, , drop = FALSE]) - global_mean
        S_b <- S_b + n_g * tcrossprod(diff_g)
    }

    S_w_reg <- S_w + diag(1e-8, p)
    eig <- eigen(solve(S_w_reg, S_b), symmetric = FALSE)
    eig_vals <- Re(eig$values)
    n_disc <- min(G - 1L, sum(eig_vals > max(eig_vals) * 1e-6))
    if (n_disc == 0L) return(NULL)

    V <- Re(eig$vectors[, seq_len(n_disc), drop = FALSE])
    qr.Q(qr(V))
}

.lda_batch_correct <- function(X, batch, verbose = FALSE) {
    Q <- .compute_lda_basis(X, batch)
    if (is.null(Q)) return(X)

    .log_msg(1, verbose,
             "Batch correction: projecting out %d LDA components (%d batches)",
             ncol(Q), nlevels(as.factor(batch)))
    X - X %*% tcrossprod(Q)
}

# --- S3 methods ---

#' @export
print.CellDistancePipeline <- function(x, ...) {
    cat("CellDistancePipeline result\n")
    cat(sprintf("  FMM: %d clusters, %d cells\n",
                ncol(x$init_fmm@Fmat), ncol(x$init_fmm)))
    cat(sprintf("  cell_pca: %d cells x %d PCs\n",
                nrow(x$cell_pca), ncol(x$cell_pca)))
    cat(sprintf("  Top loci used: %d (of %d scored)\n",
                sum(names(x$mutualinfo) %in% colnames(x$cell_pca) == FALSE),
                length(x$mutualinfo)))
    if (!is.null(x$dist)) {
        cat(sprintf("  dist: %d x %d\n", nrow(x$dist), ncol(x$dist)))
    } else {
        cat("  dist: not computed\n")
    }
    invisible(x)
}
