
test_that("fitCellDistancePipeline runs full pipeline", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_s4_class(result$init_fmm, "FiniteMixtureModel")
    expect_s3_class(result$scree, "FmmDimEstimate")
    expect_true(is.numeric(result$mutualinfo))
    expect_equal(length(result$mutualinfo), 5L)
    expect_true(is.matrix(result$cell_pca))
    expect_equal(nrow(result$cell_pca), 20L)
    expect_true(is.matrix(result$dist))
    expect_equal(dim(result$dist), c(20L, 20L))
    expect_true(all(diag(result$dist) < 1e-10))
    expect_true(all(result$dist >= -1e-10))
})

test_that("compute_dist = FALSE returns NULL dist but valid cell_pca", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      compute_dist = FALSE,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_null(result$dist)
    expect_true(is.matrix(result$cell_pca))
    expect_equal(nrow(result$cell_pca), 20L)
})

test_that("reusing FMM from previous pipeline result", {
    cnt <- .make_cell_cnt_large()
    result1 <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                       fmm_control = list(maxiter = 10))
    result2 <- fitCellDistancePipeline(cnt, n_loci = 3, init = result1)

    expect_identical(result2$init_fmm, result1$init_fmm)
    expect_true(is.matrix(result2$dist))
})

test_that("providing bare FiniteMixtureModel", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- fitCellDistancePipeline(cnt, n_loci = 3, init = fmm)

    expect_identical(result$init_fmm, fmm)
    expect_true(is.matrix(result$cell_pca))
})

test_that("manual n_pcs override", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 5, n_loci = 3,
                                      n_pcs = 2,
                                      fmm_control = list(maxiter = 10))

    expect_equal(ncol(result$cell_pca), 2L)
})

test_that("scree has all eigenvalues", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 5, n_loci = 3,
                                      fmm_control = list(maxiter = 10))

    expect_equal(length(result$scree$eigenvalues), 4L)
    expect_true(all(result$scree$eigenvalues >= 0))
    expect_equal(result$scree$method, "gavish_donoho")
})

test_that("n_loci capped to available loci", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 100,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_true(is.matrix(result$dist))
})

test_that("with batch correction", {
    cnt <- .make_cell_cnt_large()
    batch <- rep(c("A", "B"), each = 10)
    names(batch) <- paste0("cell", 1:20)

    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      batch = batch,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_true(is.matrix(result$dist))
    expect_equal(dim(result$dist), c(20L, 20L))
})

test_that("batch correction with positional batch vector", {
    cnt <- .make_cell_cnt_large()
    batch <- rep(c("A", "B"), each = 10)

    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      batch = batch,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_true(is.matrix(result$dist))
})

test_that("print method works", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      fmm_control = list(maxiter = 10))

    expect_output(print(result), "CellDistancePipeline result")
})

test_that("character n_pcs selects method", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      n_pcs = "gavish_donoho",
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_equal(result$scree$method, "gavish_donoho")
    expect_true(ncol(result$cell_pca) >= 2L)
})

test_that("automatic n_pcs enforces minimum with batch", {
    cnt <- .make_cell_cnt_large()
    batch <- rep(c("A", "B", "C"), length.out = 20)
    names(batch) <- paste0("cell", 1:20)

    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      batch = batch,
                                      fmm_control = list(maxiter = 10))

    expect_true(ncol(result$cell_pca) >= 1L + nlevels(as.factor(batch)))
})

test_that("pseudocount = 'auto' estimates from data", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      pseudocount = "auto",
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_true(is.numeric(result$pseudocount))
    expect_true(result$pseudocount > 0)
    expect_s3_class(result$pseudocount_estimate,
                    "DirichletPseudocountEstimate")
    expect_true(is.matrix(result$cell_pca))
})

test_that("numeric pseudocount returns NULL estimate", {
    cnt <- .make_cell_cnt_large()
    result <- fitCellDistancePipeline(cnt, n_clusters = 3, n_loci = 3,
                                      pseudocount = 5,
                                      fmm_control = list(maxiter = 10))

    expect_s3_class(result, "CellDistancePipeline")
    expect_equal(result$pseudocount, 5)
    expect_null(result$pseudocount_estimate)
})

test_that("explicit loci are refit once for any positive reestimate_fmm", {
    cnt <- .make_cell_cnt_large()
    initial <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    fixed <- names(cnt)[c(3, 1)]
    refitted <- fitFiniteMixtureModel(cnt[fixed], init = initial@L,
                                      combined_multinomial = FALSE, maxiter = 10)
    for (refits in list(FALSE, 0, TRUE, 1, 2, 3)) {
        expected <- if (refits > 0) refitted else initial
        result <- fitCellDistancePipeline(cnt, init = initial, keep_loci = fixed,
                                          n_loci = NULL, reestimate_fmm = refits,
                                          n_pcs = 2, fmm_control = list(maxiter = 10))
        expect_identical(result$keep_loci, fixed)
        expect_identical(names(result$projected_fmm), names(expected))
        expect_equal(result$projected_fmm@L, expected@L)
        expect_equal(result$projected_fmm@pi, expected@pi)
    }
})

test_that("batch labels stay aligned when a supplied FMM reorders cells", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    permutation <- as.vector(rbind(1:10, 11:20))
    shuffled <- cnt[, permutation]
    batch <- rep(c("A", "B"), each = 10)[permutation]
    named_batch <- setNames(batch, colnames(shuffled))
    for (fixed in list(NULL, names(cnt)[1:2])) {
        positional <- fitCellDistancePipeline(shuffled, init = fmm, batch = batch,
                                              n_loci = 2, keep_loci = fixed, n_pcs = 2)
        named <- fitCellDistancePipeline(shuffled, init = fmm, batch = named_batch,
                                         n_loci = 2, keep_loci = fixed, n_pcs = 2)
        expect_identical(positional, named)
    }
})

test_that("explicit locus indices refer to the input counts", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    reordered <- cnt[rev(names(cnt))]
    indexed <- fitCellDistancePipeline(reordered, init = fmm, keep_loci = 1:2, n_pcs = 2)
    named <- fitCellDistancePipeline(reordered, init = fmm,
                                     keep_loci = names(reordered)[1:2], n_pcs = 2)
    expect_identical(indexed, named)
})
