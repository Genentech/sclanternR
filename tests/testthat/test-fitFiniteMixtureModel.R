
test_that("fitFiniteMixtureModel returns a FiniteMixtureModel object", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    expect_s4_class(res, "FiniteMixtureModel")
    expect_length(mixing_proportions(res), 2L)
    expect_equal(nrow(loadings(res)), 4L)
    expect_equal(ncol(loadings(res)), 2L)
    expect_equal(length(coef(res)), 2L)
    expect_equal(dim(coef(res)[[1]]), c(2L, 2L))
})

test_that("fitFiniteMixtureModel: pi sums to 1, L rows sum to 1, F_list columns sum to 1", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    expect_equal(sum(mixing_proportions(res)), 1, tolerance = 1e-10)
    expect_equal(rowSums(loadings(res)), rep(1, 4), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitFiniteMixtureModel: log-likelihood is monotonically non-decreasing", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 30, tol = 1e-8)

    expect_true(all(diff(res@loglik) >= -1e-6))
})

test_that("fitFiniteMixtureModel: missing allele within a locus keeps the locus", {
    cnt <- .make_cell_cnt()
    ref_feat <- c("chr1_100_REF", "chr2_200_REF")
    ft <- .make_ft_model(ref_feat, colnames(cnt[[1]]), K = 2)

    expect_no_message(res <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20))

    expect_equal(length(coef(res)), 2L)
    expect_true(all(c("chr1_100", "chr2_200") %in% names(coef(res))))
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitFiniteMixtureModel: init='hclust' (default) returns correct structure", {
    cnt <- .make_cell_cnt_large()
    res <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 5)

    expect_s4_class(res, "FiniteMixtureModel")
    expect_length(mixing_proportions(res), 2L)
    expect_equal(nrow(loadings(res)), 20L)
    expect_equal(ncol(loadings(res)), 2L)
    expect_equal(length(coef(res)), 5L)
    expect_equal(sum(mixing_proportions(res)), 1, tolerance = 1e-10)
    expect_equal(rowSums(loadings(res)), rep(1, 20), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitFiniteMixtureModel: init='fastTopics' returns correct structure", {
    cnt <- .make_cell_cnt_large()
    res <- fitFiniteMixtureModel(cnt, n_clusters = 2, init = "fastTopics",
                                 maxiter = 5)

    expect_s4_class(res, "FiniteMixtureModel")
    expect_length(mixing_proportions(res), 2L)
    expect_equal(nrow(loadings(res)), 20L)
    expect_equal(ncol(loadings(res)), 2L)
    expect_equal(length(coef(res)), 5L)
    expect_equal(sum(mixing_proportions(res)), 1, tolerance = 1e-10)
    expect_equal(rowSums(loadings(res)), rep(1, 20), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitFiniteMixtureModel: init from L matrix only", {
    cnt <- .make_cell_cnt()
    L_mat <- matrix(abs(rnorm(4 * 2)) + 0.1, nrow = 4, ncol = 2)
    res <- fitFiniteMixtureModel(cnt, init = L_mat, maxiter = 5)

    expect_s4_class(res, "FiniteMixtureModel")
    expect_length(mixing_proportions(res), 2L)
    expect_equal(nrow(loadings(res)), 4L)
    expect_equal(length(coef(res)), 2L)
    expect_equal(sum(mixing_proportions(res)), 1, tolerance = 1e-10)
})

test_that("fitFiniteMixtureModel: init='random' works", {
    cnt <- .make_cell_cnt()
    res <- fitFiniteMixtureModel(cnt, n_clusters = 3, init = "random", maxiter = 5)

    expect_s4_class(res, "FiniteMixtureModel")
    expect_length(mixing_proportions(res), 3L)
    expect_equal(nrow(loadings(res)), 4L)
    expect_equal(ncol(loadings(res)), 3L)
    expect_equal(sum(mixing_proportions(res)), 1, tolerance = 1e-10)
})

test_that("FiniteMixtureModel: subsetting by loci", {
    cnt <- .make_cell_cnt_large()
    res <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 5)

    sub <- res[1:2]
    expect_s4_class(sub, "FiniteMixtureModel")
    expect_equal(length(sub), 2L)
    expect_equal(ncol(sub), 20L)
    expect_equal(mixing_proportions(sub), mixing_proportions(res))
})

test_that("FiniteMixtureModel: subsetting by cells", {
    cnt <- .make_cell_cnt_large()
    res <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 5)

    sub <- res[, 1:5]
    expect_s4_class(sub, "FiniteMixtureModel")
    expect_equal(length(sub), 5L)
    expect_equal(ncol(sub), 5L)
    expect_equal(nrow(loadings(sub)), 5L)
})

test_that("FiniteMixtureModel: predict returns frequencies", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    freqs <- predict(res)
    expect_true(is.list(freqs))
    expect_equal(length(freqs), 2L)
    expect_equal(ncol(freqs[[1]]), 4L)
})

test_that("FiniteMixtureModel: predict with newx returns loadings", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    L_pred <- predict(res, newx = cnt, type = "loadings")
    expect_true(is.matrix(L_pred))
    expect_equal(nrow(L_pred), 4L)
    expect_equal(ncol(L_pred), 2L)
    expect_equal(rowSums(L_pred), rep(1, 4), tolerance = 1e-6,
                 ignore_attr = TRUE)
})

test_that("FiniteMixtureModel: mutualinfo works", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    mi <- mutualinfo(res, cnt)
    expect_true(is.numeric(mi))
    expect_equal(length(mi), 2L)
    expect_true(all(mi >= 0))
})

test_that("FiniteMixtureModel: constructor from F_list and pi", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)

    model2 <- FiniteMixtureModel(coef(res), pi = mixing_proportions(res))
    expect_s4_class(model2, "FiniteMixtureModel")
    expect_equal(length(model2), 2L)
    expect_equal(ncol(model2), 0L)

    L_pred <- predict(model2, newx = cnt, type = "loadings")
    expect_equal(nrow(L_pred), 4L)
})

test_that("fitFiniteMixtureModel: init from LocusModel object", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res1 <- fitLocusNMF(cnt, ft, maxiter = 5)

    res2 <- fitFiniteMixtureModel(cnt, init = res1, maxiter = 5)

    expect_s4_class(res2, "FiniteMixtureModel")
    expect_length(mixing_proportions(res2), 2L)
    expect_equal(nrow(loadings(res2)), 4L)
    expect_equal(ncol(loadings(res2)), 2L)
    expect_equal(length(coef(res2)), 2L)
    expect_equal(sum(mixing_proportions(res2)), 1, tolerance = 1e-10)
    expect_equal(rowSums(loadings(res2)), rep(1, 4), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res2)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})
