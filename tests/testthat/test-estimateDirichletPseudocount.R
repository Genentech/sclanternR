
test_that("estimateDirichletPseudocount returns valid result", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateDirichletPseudocount(fmm, cnt)

    expect_s3_class(result, "DirichletPseudocountEstimate")
    expect_true(is.numeric(result$pseudocount))
    expect_true(result$pseudocount > 0)
    expect_true(is.numeric(result$mom_estimate))
    expect_true(result$mom_estimate > 0)
    expect_true(is.logical(result$converged))
    expect_true(result$iterations >= 1)
    expect_equal(length(result$trace), result$iterations)
})

test_that("pseudocount grows large for near-multinomial data", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateDirichletPseudocount(fmm, cnt)

    # Synthetic Poisson data has no Dirichlet overdispersion,
    # so the estimated concentration should be large
    expect_true(result$pseudocount > 10)
})

test_that("MOM provides initialization", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateDirichletPseudocount(fmm, cnt)

    # MOM estimate should be the starting point of the trace
    expect_equal(result$mom_estimate, result$trace[1])
})

test_that("manual init overrides MOM", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateDirichletPseudocount(fmm, cnt, init = 5)

    expect_equal(result$mom_estimate, 5)
    expect_equal(result$trace[1], 5)
})

test_that("estimateDirichletPseudocount respects min/max bounds", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)

    result <- estimateDirichletPseudocount(fmm, cnt,
                                           min_alpha = 5, max_alpha = 5)
    expect_equal(result$pseudocount, 5)
})

test_that("estimateDirichletPseudocount works with subsetted model", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    loci <- names(fmm)[1:3]
    result <- estimateDirichletPseudocount(fmm[loci], cnt[loci])

    expect_s3_class(result, "DirichletPseudocountEstimate")
    expect_true(result$pseudocount > 0)
})

test_that("estimateDirichletPseudocount rejects invalid inputs", {
    cnt <- .make_cell_cnt_large()
    expect_error(estimateDirichletPseudocount(list(), cnt), "LocusModel")

    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    expect_error(estimateDirichletPseudocount(fmm, "bad"),
                 "CellAlleleCounts")
})

test_that("different init values both produce large pseudocount", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)

    r1 <- estimateDirichletPseudocount(fmm, cnt, init = 0.1)
    r2 <- estimateDirichletPseudocount(fmm, cnt, init = 10)

    # Both should grow large for near-multinomial synthetic data
    expect_true(r1$pseudocount > 10)
    expect_true(r2$pseudocount > 10)
})

test_that("print.DirichletPseudocountEstimate works", {
    cnt <- .make_cell_cnt_large()
    fmm <- fitFiniteMixtureModel(cnt, n_clusters = 3, maxiter = 10)
    result <- estimateDirichletPseudocount(fmm, cnt)

    expect_output(print(result), "Cox-Reid")
    expect_output(print(result), "pseudocount")
    expect_output(print(result), "MOM init")
})
