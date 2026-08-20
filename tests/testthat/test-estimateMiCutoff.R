
test_that("estimateMiCutoff: mad returns valid result", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 5, maxiter = 15)
    result <- estimateMiCutoff(fmm, cnt, method = "mad")

    expect_s3_class(result, "MiCutoffEstimate")
    expect_equal(result$method, "mad")
    expect_true(is.numeric(result$mutualinfo))
    expect_equal(length(result$mutualinfo), length(fmm))
    expect_true(is.numeric(result$cutoff))
    expect_equal(length(result$cutoff), 1L)
    expect_true(result$n_loci >= 0)
    expect_true(result$n_loci <= length(result$mutualinfo))
    expect_equal(length(result$loci), result$n_loci)
    expect_true(!is.null(result$mad_stats))
    expect_true(!is.null(result$mad_stats$median))
    expect_true(!is.null(result$mad_stats$mad))
})

test_that("estimateMiCutoff: elbow returns valid result", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 5, maxiter = 15)
    result <- estimateMiCutoff(fmm, cnt, method = "elbow")

    expect_s3_class(result, "MiCutoffEstimate")
    expect_equal(result$method, "elbow")
    expect_true(result$n_loci >= 1)
    expect_true(!is.null(result$curvature))
})

test_that("estimateMiCutoff: fdr returns valid result", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 5, maxiter = 15)
    result <- estimateMiCutoff(fmm, cnt, method = "fdr")

    expect_s3_class(result, "MiCutoffEstimate")
    expect_equal(result$method, "fdr")
    expect_true(!is.null(result$pvalues))
    expect_equal(length(result$pvalues), length(result$mutualinfo))
    expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
    expect_true(!is.null(result$df))
    expect_true(all(result$df >= 0))
})

test_that("estimateMiCutoff: mutualinfo preserves model locus order", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateMiCutoff(fmm, cnt)

    expect_equal(names(result$mutualinfo), names(fmm))
})

test_that("estimateMiCutoff: works with x = NULL", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateMiCutoff(fmm)

    expect_s3_class(result, "MiCutoffEstimate")
    expect_equal(length(result$mutualinfo), length(fmm))
})

test_that("estimateMiCutoff: works with LocusNMF", {
    cnt <- .make_cell_cnt_large()
    nmf <- fitLocusNMF(cnt, init = NULL, n_topics = 3, maxiter = 5)
    result <- estimateMiCutoff(nmf, cnt)

    expect_s3_class(result, "MiCutoffEstimate")
    expect_equal(length(result$mutualinfo), length(nmf))
})

test_that("estimateMiCutoff: rejects non-LocusModel input", {
    expect_error(estimateMiCutoff(list()), "LocusModel")
})

test_that("estimateMiCutoff: fdr df is (K-1)*(A-1)", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 5, maxiter = 15)
    result <- estimateMiCutoff(fmm, cnt, method = "fdr")

    K <- ncol(fmm@Fmat)
    n_alleles <- as.integer(Matrix::rowSums(fmm@locus_1hot))
    expected_df <- (K - 1L) * (n_alleles - 1L)

    expect_equal(as.integer(result$df), expected_df)
})

test_that("estimateMiCutoff: lower mad_threshold selects more loci", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 5, maxiter = 15)
    result_strict <- estimateMiCutoff(fmm, cnt, method = "mad", mad_threshold = 5)
    result_loose  <- estimateMiCutoff(fmm, cnt, method = "mad", mad_threshold = 1)

    expect_true(result_loose$n_loci >= result_strict$n_loci)
    expect_true(result_loose$cutoff <= result_strict$cutoff)
})

test_that("print.MiCutoffEstimate works", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateMiCutoff(fmm, cnt)

    expect_output(print(result), "MiCutoffEstimate")
    expect_output(print(result), "loci selected")
})

test_that("plot.MiCutoffEstimate works for all methods", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)

    for (m in c("mad", "elbow", "fdr")) {
        result <- estimateMiCutoff(fmm, cnt, method = m)
        expect_no_error(plot(result))
    }
})

# --- Pipeline integration ---

test_that("pipeline: character n_loci = 'elbow' works", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = "elbow",
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_s3_class(result$mi_cutoff, "MiCutoffEstimate")
    expect_equal(result$mi_cutoff$method, "elbow")
    expect_true(is.matrix(result$cell_pca))
})

test_that("pipeline: character n_loci = 'fdr' works", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = "fdr",
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_equal(result$mi_cutoff$method, "fdr")
})

test_that("pipeline: numeric n_loci returns mi_cutoff = NULL", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      fmm_control = list(maxiter = 10))

    expect_null(result$mi_cutoff)
})

test_that("pipeline: character n_loci = 'mad' clamps to min_loci with random data", {
    cnt <- .make_cell_cnt_large()
    expect_warning(
        result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = "mad",
                                          fmm_control = list(maxiter = 10)),
        "clamping to min_loci"
    )
    expect_s3_class(result, "CellDistancePipeline")
})

test_that("pipeline: character n_loci with batch correction", {
    cnt <- .make_cell_cnt_large()
    batch <- rep(c("A", "B"), each = 10)
    names(batch) <- paste0("cell", 1:20)

    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = "elbow",
                                      batch = batch,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_true(is.matrix(result$dist))
})
