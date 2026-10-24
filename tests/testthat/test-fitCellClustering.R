test_that("clustering alone preserves zero and one refit behavior", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    mi <- mutualinfo(fmm, cnt)
    selected <- names(sort(mi, decreasing = TRUE))[1:2]
    zero <- fitCellClustering(cnt, init = fmm, n_loci = 2, reestimate_fmm = 0)
    expect_s3_class(zero, "CellClustering")
    expect_identical(zero$fmm, fmm)
    expect_equal(zero$keep_loci, selected)
    expect_equal(zero$mutualinfo, mi)
    expect_null(zero$cell_pca)
    expect_identical(zero, fitCellClustering(cnt, init = fmm, n_loci = 2,
                                            reestimate_fmm = FALSE))

    one <- fitCellClustering(cnt, init = fmm, n_loci = 2, reestimate_fmm = 1,
                             fmm_control = list(maxiter = 10))
    expected <- fitFiniteMixtureModel(cnt[selected], init = fmm@L,
                                      combined_multinomial = FALSE, maxiter = 10)
    expect_identical(one$fmm, expected)
    expect_equal(one$mutualinfo, mi)
    expect_identical(one, fitCellClustering(cnt, init = fmm, n_loci = 2,
                                           reestimate_fmm = TRUE,
                                           fmm_control = list(maxiter = 10)))
})

test_that("multiple refits match manual all-locus rescoring and fitting", {
    cnt <- .make_cell_cnt_large()
    initial <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    model <- initial
    for (iteration in 1:3) {
        # Independent per-locus contingency tables, including excluded loci.
        mi <- vapply(as.list(cnt), function(x) {
            sclanternR:::.compute_scaled_mutinfo(as.matrix(x %*% model@L))
        }, numeric(1))
        selected <- names(sort(mi, decreasing = TRUE))[1:2]
        model <- fitFiniteMixtureModel(cnt[selected], init = model@L,
                                       combined_multinomial = FALSE, maxiter = 10)
        actual <- fitCellClustering(cnt, init = initial, n_loci = 2,
                                     reestimate_fmm = iteration,
                                     fmm_control = list(maxiter = 10))
        expect_equal(actual$mutualinfo, mi)
        expect_equal(actual$keep_loci, selected)
        expect_identical(actual$fmm, model)
        expect_equal(names(actual$mutualinfo), names(cnt))
    }
})

test_that("refits preserve explicit multinomial modes except initialization-only mode", {
    cnt <- .make_cell_cnt_large()
    initial <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    selected <- names(sort(mutualinfo(initial, cnt), decreasing = TRUE))[1:2]
    for (mode in list(FALSE, TRUE, "init")) {
        control <- list(combined_multinomial = mode, maxiter = 10)
        expected <- fitFiniteMixtureModel(
            cnt[selected], init = initial@L, maxiter = 10,
            combined_multinomial = identical(mode, TRUE))
        clustering <- fitCellClustering(
            cnt, init = initial, n_loci = 2, reestimate_fmm = 1,
            fmm_control = control)
        expect_identical(clustering$fmm, expected)
        pipeline <- fitCellDistancePipeline(
            cnt, init = initial, keep_loci = selected, reestimate_fmm = 1,
            n_pcs = 2, fmm_control = control)
        expect_equal(pipeline$projected_fmm@L, expected@L)
        expect_equal(pipeline$projected_fmm@pi, expected@pi)
    }
})

test_that("an excluded locus can re-enter after responsibilities change", {
    cells <- paste0("cell", 1:4)
    loci <- c("chr1_100", "chr2_200")
    groups <- list(c(1, 1, 2, 2), c(1, 2, 1, 2))
    counts <- lapply(seq_along(loci), function(i) {
        rbind(REF = ifelse(groups[[i]] == 1, 20, 0),
              ALT = ifelse(groups[[i]] == 2, 20, 0))
    })
    for (i in seq_along(loci)) {
        dimnames(counts[[i]]) <- list(paste0(loci[i], c("_REF", "_ALT")), cells)
    }
    cnt <- CellAlleleCounts_from_list(setNames(counts, loci))
    L <- lapply(groups, function(g) {
        matrix(as.numeric(outer(g, 1:2, `==`)), nrow = 4,
               dimnames = list(cells, NULL))
    })
    initial <- FiniteMixtureModel(lapply(as.list(cnt), function(x) {
        matrix(0.5, nrow = 2, ncol = 2, dimnames = list(rownames(x), NULL))
    }), pi = c(0.5, 0.5), L = L[[1]])
    fit_loci <- list()
    local_mocked_bindings(fitFiniteMixtureModel = function(cell_allele_cnt_list,
                                                          init, ...) {
        fit_loci[[length(fit_loci) + 1L]] <<- names(cell_allele_cnt_list)
        model <- initial[names(cell_allele_cnt_list)]
        model@L <- L[[2]]
        model
    }, .package = "sclanternR")
    result <- fitCellClustering(cnt, init = initial, n_loci = 1, reestimate_fmm = 2)
    expect_identical(fit_loci, as.list(loci))
    expect_identical(result$keep_loci, loci[2])
    expect_equal(names(result$mutualinfo), loci)
    expect_gt(result$mutualinfo[2], result$mutualinfo[1])
})

test_that("cutoff controls apply across rounds", {
    cnt <- .make_cell_cnt_large()
    initial <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    automatic <- fitCellClustering(cnt, init = initial, n_loci = "fdr", reestimate_fmm = 2,
                                   n_loci_control = list(min_loci = 2, max_loci = 2),
                                   fmm_control = list(maxiter = 10))
    expect_length(automatic$keep_loci, 2)
    expect_s3_class(automatic$mi_cutoff, "MiCutoffEstimate")
    expect_length(automatic$mi_cutoff$df, length(cnt))
    expect_identical(automatic$mi_cutoff$loci, automatic$keep_loci)
})

test_that("batch selection rescans all loci and caps numeric selections", {
    cnt <- .make_cell_cnt_large()
    batch <- setNames(rep(c("A", "B"), each = 10), colnames(cnt))
    initial <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    first <- fitCellClustering(cnt, init = initial, batch = batch, n_loci = 1,
                               reestimate_fmm = 1, fmm_control = list(maxiter = 10))
    result <- fitCellClustering(cnt, init = initial, batch = batch, n_loci = 1,
                                reestimate_fmm = 2, fmm_control = list(maxiter = 10))
    expected <- sapply(c("A", "B"), function(b) {
        cells <- names(batch)[batch == b]
        mutualinfo(first$fmm[, cells], cnt[, cells])
    })
    expect_equal(result$mutualinfo, expected)
    expect_equal(result$keep_loci,
                 unique(rownames(expected)[apply(expected, 2, which.max)]))
    capped <- fitCellClustering(cnt, init = initial, batch = batch, n_loci = 100)
    expect_setequal(capped$keep_loci, names(cnt))
    expect_false(anyNA(capped$keep_loci))
    automatic <- fitCellClustering(cnt, init = initial, batch = batch, n_loci = "fdr",
                                   reestimate_fmm = 2,
                                   n_loci_control = list(min_loci = 1, max_loci = 1),
                                   fmm_control = list(maxiter = 10))
    expect_equal(dim(automatic$mutualinfo), c(length(cnt), 2))
    expect_true(length(automatic$keep_loci) %in% 1:2)
})

test_that("refit count rejects invalid values before fitting", {
    cnt <- .make_cell_cnt_large()
    local_mocked_bindings(fitFiniteMixtureModel = function(...) stop("unexpected fit"),
                          .package = "sclanternR")
    for (value in list(-1, 0.5, NA, NA_integer_, Inf, "2", c(1, 2), NULL, 1i)) {
        expect_error(fitCellClustering(cnt, reestimate_fmm = value),
                     "reestimate_fmm must be a non-negative integer")
        expect_error(fitCellDistancePipeline(cnt, reestimate_fmm = value),
                     "reestimate_fmm must be a non-negative integer")
        expect_error(fitCellDistancePipeline(cnt, keep_loci = names(cnt)[1],
                                             reestimate_fmm = value),
                     "reestimate_fmm must be a non-negative integer")
    }
})

test_that("resuming clustering gives the same pipeline without refitting", {
    cnt <- .make_cell_cnt_large()
    batch <- rep(c("A", "B"), each = 10)
    full <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 2,
                                    batch = batch, reestimate_fmm = 2,
                                    fmm_control = list(maxiter = 10), n_pcs = 2)
    clustering <- fitCellClustering(cnt, n_clusters = 3, n_loci = 2,
                                    batch = batch, reestimate_fmm = 2,
                                    fmm_control = list(maxiter = 10))
    local_mocked_bindings(fitCellClustering = function(...) stop("unexpected fit"),
                          fitFiniteMixtureModel = function(...) stop("unexpected fit"),
                          .package = "sclanternR")
    resumed <- fitCellDistancePipeline(cnt, init = clustering$fmm,
                                       keep_loci = clustering$keep_loci,
                                       batch = clustering$batch, n_pcs = 2)
    downstream <- c("keep_loci", "final_fmm", "projected_fmm", "cell_pca", "dist",
                    "freq_list", "cluster_loadings", "pseudocount", "loadings_choice")
    expect_identical(resumed[downstream], full[downstream])
    reordered <- fitCellDistancePipeline(cnt[, rev(colnames(cnt))],
                                         init = clustering$fmm,
                                         keep_loci = clustering$keep_loci,
                                         batch = clustering$batch, n_pcs = 2)
    expect_identical(reordered, resumed)
})

test_that("selection arguments and initial models are checked before fitting", {
    cnt <- .make_cell_cnt_large()
    local_mocked_bindings(fitFiniteMixtureModel = function(...) stop("unexpected fit"),
                          .package = "sclanternR")
    for (value in list(0, -1, 0.5, NA, Inf, numeric(), c(1, 2), 1i)) {
        expect_error(fitCellClustering(cnt, n_loci = value),
                     "n_loci must be a positive integer")
        expect_error(fitCellDistancePipeline(cnt, n_loci = value),
                     "n_loci must be a positive integer")
    }
    expect_error(fitCellClustering(cnt, init = list()),
                 "init must be NULL or a FiniteMixtureModel")
})

test_that("FMM scores and cutoff degrees of freedom use all supplied loci", {
    cnt <- .make_cell_cnt_unequal_loci()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 5)
    subset <- fmm[1]
    expect_equal(mutualinfo(subset, cnt), mutualinfo(fmm, cnt))
    result <- estimateMiCutoff(subset, cnt)
    expect_equal(unname(result$df), as.numeric(Matrix::rowSums(cnt@locus_1hot) - 1))
    expect_setequal(names(result$mutualinfo), names(cnt))
})
