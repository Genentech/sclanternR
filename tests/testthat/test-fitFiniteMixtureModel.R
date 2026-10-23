
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

# Independent dense reference: explicitly center features within each locus.
.fmm_hclust_reference_dist <- function(Y_t, locus_index = NULL) {
    freq <- as.matrix(Y_t) / pmax(Matrix::rowSums(Y_t), 1)
    for (idx in locus_index) {
        freq[, idx] <- sweep(freq[, idx, drop = FALSE], 1,
                            rowMeans(freq[, idx, drop = FALSE]), "-")
    }
    stats::dist(freq)
}

test_that("hclust distances match explicit within-locus centering", {
    cnt <- .make_cell_cnt_unequal_loci()
    Y_t <- t(cnt@counts)
    total <- as.matrix(sclanternR:::.fmm_hclust_dist(Y_t))
    within <- as.matrix(sclanternR:::.fmm_hclust_dist(Y_t, cnt@locus_1hot))

    expect_equal(total^2, as.matrix(.fmm_hclust_reference_dist(Y_t))^2,
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(within^2,
                 as.matrix(.fmm_hclust_reference_dist(Y_t, cnt@locus_index))^2,
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_true(all(is.finite(within)))
    expect_true(all(within >= 0))
    expect_gt(max(abs(total^2 - within^2)), 0.01)
})

test_that("within-locus distances remove locus means and singleton features", {
    cnt <- .make_cell_cnt_unequal_loci()
    # The second cell adds constants (1, -1, 1) within the three loci;
    # the total count stays fixed, so their within-locus residuals agree.
    Y_t <- Matrix::Matrix(rbind(c(8, 6, 3, 2, 1, 4),
                               c(9, 7, 2, 1, 0, 5),
                               c(6, 8, 3, 2, 1, 4)), sparse = TRUE)
    within <- as.matrix(sclanternR:::.fmm_hclust_dist(Y_t, cnt@locus_1hot))
    expect_equal(within[1, 2]^2, 0, tolerance = 1e-12)
    expect_equal(within[1, 3], within[2, 3], tolerance = 1e-12)
    expect_gt(within[1, 3], 0)

    # A single locus only removes the common global mean for nonempty cells.
    one_locus <- Matrix::Matrix(1, nrow = 1, ncol = ncol(Y_t), sparse = TRUE)
    expect_equal(as.matrix(sclanternR:::.fmm_hclust_dist(Y_t, one_locus))^2,
                 as.matrix(stats::dist(as.matrix(Y_t) / 24))^2,
                 tolerance = 1e-12)

    singletons <- Matrix::Diagonal(ncol(Y_t))
    expect_equal(as.numeric(sclanternR:::.fmm_hclust_dist(Y_t, singletons)),
                 rep(0, 3), tolerance = 1e-12)
})

test_that("within-locus initialization uses original loci in every model mode", {
    cnt <- .make_cell_cnt_unequal_loci()
    original_dist <- sclanternR:::.fmm_hclust_dist
    sampled_sizes <- integer()
    local_mocked_bindings(.fmm_hclust_dist = function(Y_t, locus_1hot = NULL, rank = NULL) {
        actual <- original_dist(Y_t, locus_1hot, rank)
        expected <- .fmm_hclust_reference_dist(Y_t, cnt@locus_index)
        expect_equal(as.matrix(actual)^2, as.matrix(expected)^2,
                     tolerance = 1e-12, ignore_attr = TRUE)
        sampled_sizes <<- c(sampled_sizes, nrow(Y_t))
        actual
    }, .package = "sclanternR")

    for (mode in list(FALSE, TRUE, "init")) {
        for (n_init in list(NULL, 4L)) {
            set.seed(42)
            res <- fitFiniteMixtureModel(
                cnt, n_clusters = 2, combined_multinomial = mode,
                hclust_covariance = "within_locus", n_init_cells = n_init,
                maxiter = 2)
            expect_s4_class(res, "FiniteMixtureModel")
            expect_equal(dim(loadings(res)), c(6L, 2L))
            expect_true(all(is.finite(loadings(res))))
            expect_equal(unname(rowSums(loadings(res))), rep(1, 6))
        }
    }
    expect_equal(sampled_sizes, rep(c(6L, 4L), 3))
})

test_that("hclust covariance defaults to total and validates its value", {
    cnt <- .make_cell_cnt_unequal_loci()
    default <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 2)
    total <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 2,
                                   hclust_covariance = "total")
    expect_equal(default, total)
    expect_error(fitFiniteMixtureModel(cnt, hclust_covariance = "invalid"),
                 "'arg' should be one of")
})

test_that("hclust rank reduction matches the public inner product on sampled cells", {
    cnt <- .make_cell_cnt_unequal_loci()
    original_dist <- sclanternR:::.fmm_hclust_dist
    sampled_sizes <- integer()
    local_mocked_bindings(.fmm_hclust_dist = function(Y_t, locus_1hot = NULL, rank = NULL) {
        expect_identical(rank, 1L)
        actual <- original_dist(Y_t, locus_1hot, rank)
        sampled <- cnt[, rownames(Y_t)]
        expected <- computeF2dist(computeCellInnerProduct(
            sampled, covariance = covariance, rank = 1L))
        expect_equal(as.matrix(actual)^2, pmax(expected, 0), tolerance = 1e-12)
        sampled_sizes <<- c(sampled_sizes, nrow(Y_t))
        actual
    }, .package = "sclanternR")

    for (covariance in c("total", "within_locus")) {
        for (mode in list(FALSE, TRUE, "init")) {
            for (n_init in list(NULL, 4L)) {
                set.seed(42)
                res <- fitFiniteMixtureModel(
                    cnt, n_clusters = 2, combined_multinomial = mode,
                    hclust_covariance = covariance, hclust_rank = 1L,
                    n_init_cells = n_init, maxiter = 2)
                expect_true(all(is.finite(loadings(res))))
                expect_equal(dim(loadings(res)), c(6L, 2L))
            }
        }
    }
    expect_equal(sampled_sizes, rep(c(6L, 4L), 6))
})

test_that("hclust rank defaults to no truncation and validates its value", {
    cnt <- .make_cell_cnt_unequal_loci()
    default <- fitFiniteMixtureModel(cnt, n_clusters = 2, maxiter = 2)
    expect_equal(default, fitFiniteMixtureModel(
        cnt, n_clusters = 2, maxiter = 2, hclust_rank = NULL))
    expect_equal(default, fitFiniteMixtureModel(
        cnt, n_clusters = 2, maxiter = 2, hclust_rank = 100))
    expect_error(fitFiniteMixtureModel(cnt, hclust_rank = 1.5),
                 "hclust_rank must be NULL or a positive integer")
})
