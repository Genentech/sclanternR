
.make_fmm <- function(n_clusters = 5) {
    cnt <- .make_cell_cnt_large()
    fitFiniteMixtureModel(cnt, n_clusters = n_clusters, maxiter = 15)
}

test_that("estimateFmmDimensionality: marchenko_pastur (default) returns valid result", {
    fmm <- .make_fmm()
    result <- estimateFmmDimensionality(fmm)

    expect_s3_class(result, "FmmDimEstimate")
    expect_equal(result$method, "gavish_donoho")
    expect_true(all(result$eigenvalues >= 0))
    expect_true(is.numeric(result$threshold))
    expect_true(result$threshold > 0)
    expect_null(result$n_topics)
})

test_that("estimateFmmDimensionality: gavish_donoho returns valid result", {
    fmm <- .make_fmm()
    result <- estimateFmmDimensionality(fmm, method = "gavish_donoho")

    expect_s3_class(result, "FmmDimEstimate")
    expect_equal(result$method, "gavish_donoho")
    expect_true(all(result$eigenvalues >= 0))
    expect_true(is.numeric(result$threshold))
    expect_true(result$threshold > 0)
})

test_that("estimateFmmDimensionality: marchenko_pastur is less conservative than gavish_donoho", {
    fmm <- .make_fmm()
    gd <- estimateFmmDimensionality(fmm, method = "gavish_donoho")
    mp <- estimateFmmDimensionality(fmm, method = "marchenko_pastur")

    expect_true(mp$threshold <= gd$threshold)
    expect_true(mp$cutoff_index >= gd$cutoff_index)
})

test_that("estimateFmmDimensionality: kneedle returns valid result", {
    fmm <- .make_fmm()
    result <- estimateFmmDimensionality(fmm, method = "kneedle")

    expect_s3_class(result, "FmmDimEstimate")
    expect_equal(result$method, "kneedle")
    expect_true(all(result$eigenvalues >= 0))
    expect_equal(length(result$curvature), length(result$eigenvalues))
    expect_true(result$cutoff_index >= 1)
    expect_true(result$cutoff_index <= length(result$eigenvalues))
})

test_that("estimateFmmDimensionality: rejects non-FMM input", {
    cnt <- .make_cell_cnt_large()
    nmf <- fitLocusNMF(cnt, init = NULL, n_topics = 2, maxiter = 5)
    expect_error(estimateFmmDimensionality(nmf), "FiniteMixtureModel")
})

test_that("plot.FmmDimEstimate runs without error for all methods", {
    fmm <- .make_fmm()
    for (method in c("gavish_donoho", "marchenko_pastur", "kneedle")) {
        result <- estimateFmmDimensionality(fmm, method = method)
        expect_no_error(plot(result))
    }
})

test_that("print.FmmDimEstimate runs without error", {
    fmm <- .make_fmm()
    result <- estimateFmmDimensionality(fmm)
    expect_output(print(result), "FmmDimEstimate")
    expect_output(print(result), "cutoff_index")
})

