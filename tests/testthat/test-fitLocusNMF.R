
test_that("fitLocusNMF returns a LocusNMF object", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    expect_s4_class(res, "LocusNMF")
    expect_equal(nrow(loadings(res)), 4L)
    expect_equal(ncol(loadings(res)), 2L)
    expect_equal(length(coef(res)), 2L)
    expect_equal(dim(coef(res)[[1]]), c(2L, 2L))
    expect_true(length(res@loglik) > 0)
})

test_that("fitLocusNMF: L rows and F_list columns sum to 1", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    expect_equal(rowSums(loadings(res)), rep(1, 4), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitLocusNMF: log-likelihood is monotonically non-decreasing", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 30, tol = 1e-8)

    expect_true(all(diff(res@loglik) >= -1e-6))
})

test_that("fitLocusNMF: missing allele within a locus keeps the locus", {
    cnt <- .make_cell_cnt()
    ref_feat <- c("chr1_100_REF", "chr2_200_REF")
    ft <- .make_ft_model(ref_feat, colnames(cnt[[1]]), K = 2)

    expect_no_message(res <- fitLocusNMF(cnt, ft, maxiter = 20))

    expect_equal(length(coef(res)), 2L)
    expect_true(all(c("chr1_100", "chr2_200") %in% names(coef(res))))

    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitLocusNMF: init=NULL auto-fits and returns correct structure", {
    cnt <- .make_cell_cnt_large()
    res <- fitLocusNMF(cnt, init = NULL, n_topics = 2, maxiter = 5)

    expect_s4_class(res, "LocusNMF")
    expect_equal(nrow(loadings(res)), 20L)
    expect_equal(ncol(loadings(res)), 2L)
    expect_equal(length(coef(res)), 5L)
    expect_equal(rowSums(loadings(res)), rep(1, 20), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})

test_that("fitLocusNMF: locus entirely absent from init is dropped", {
    cnt <- .make_cell_cnt()
    feat <- c("chr1_100_REF", "chr1_100_ALT")
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)

    expect_message(res <- fitLocusNMF(cnt, ft, maxiter = 5),
                   regexp = "Dropping 1 loci not present in init")

    expect_equal(length(coef(res)), 1L)
    expect_true("chr1_100" %in% names(coef(res)))
    expect_false("chr2_200" %in% names(coef(res)))
})

test_that("LocusNMF: subsetting by loci", {
    cnt <- .make_cell_cnt_large()
    res <- fitLocusNMF(cnt, init = NULL, n_topics = 2, maxiter = 5)

    sub <- res[1:2]
    expect_s4_class(sub, "LocusNMF")
    expect_equal(length(sub), 2L)
    expect_equal(ncol(sub), 20L)

    sub2 <- res["chr1_100"]
    expect_equal(length(sub2), 1L)
    expect_equal(names(sub2), "chr1_100")
})

test_that("LocusNMF: subsetting by cells", {
    cnt <- .make_cell_cnt_large()
    res <- fitLocusNMF(cnt, init = NULL, n_topics = 2, maxiter = 5)

    sub <- res[, 1:5]
    expect_s4_class(sub, "LocusNMF")
    expect_equal(length(sub), 5L)
    expect_equal(ncol(sub), 5L)
    expect_equal(nrow(loadings(sub)), 5L)
})

test_that("LocusNMF: predict returns frequencies", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    freqs <- predict(res)
    expect_true(is.list(freqs))
    expect_equal(length(freqs), 2L)
    expect_equal(ncol(freqs[[1]]), 4L)
})

test_that("LocusNMF: predict with newx returns loadings", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    L_pred <- predict(res, newx = cnt, type = "loadings")
    expect_true(is.matrix(L_pred))
    expect_equal(nrow(L_pred), 4L)
    expect_equal(ncol(L_pred), 2L)
    expect_equal(rowSums(L_pred), rep(1, 4), tolerance = 1e-6,
                 ignore_attr = TRUE)
})

test_that("LocusNMF: mutualinfo works", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    mi <- mutualinfo(res, cnt)
    expect_true(is.numeric(mi))
    expect_equal(length(mi), 2L)
    expect_true(all(mi >= 0))
})

test_that("vectorized mutualinfo matches per-locus sapply version", {
    cnt  <- .make_cell_cnt_large()
    feat <- unlist(lapply(as(cnt, "list"), rownames))
    ft   <- .make_ft_model(feat, colnames(cnt), K = 3)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    # Reference: per-locus sapply using .compute_scaled_mutinfo
    ref_sapply <- function(Fmat_stat, locus_index) {
        sapply(locus_index, function(idx)
            sclanternR:::.compute_scaled_mutinfo(
                t(Fmat_stat[idx, , drop = FALSE])))
    }

    # LocusNMF, x = NULL
    Fmat_stat1 <- sweep(res@Fmat, 2, colSums(res@L), `*`)
    mi_vec1 <- sclanternR:::.compute_vectorized_mutinfo(Fmat_stat1, res@locus_1hot)
    mi_ref1 <- ref_sapply(Fmat_stat1, res@locus_index)
    expect_equal(mi_vec1, mi_ref1, tolerance = 1e-12)

    # LocusNMF, x = data
    mi_data <- mutualinfo(res, cnt)
    expect_equal(length(mi_data), 5L)
    expect_true(all(mi_data >= 0))

    # FiniteMixtureModel, x = NULL
    fmm <- fitFiniteMixtureModel(cnt, init = ft, maxiter = 20)
    Fmat_stat2 <- sweep(fmm@Fmat, 2, fmm@pi, `*`)
    mi_vec2 <- sclanternR:::.compute_vectorized_mutinfo(Fmat_stat2, fmm@locus_1hot)
    mi_ref2 <- ref_sapply(Fmat_stat2, fmm@locus_index)
    expect_equal(mi_vec2, mi_ref2, tolerance = 1e-12)

    # FiniteMixtureModel, x = data
    Y_t <- t(as(cnt, 'CsparseMatrix'))
    Fmat_stat3 <- t(as.matrix(t(fmm@L) %*% Y_t))
    mi_vec3 <- sclanternR:::.compute_vectorized_mutinfo(Fmat_stat3, fmm@locus_1hot)
    mi_ref3 <- ref_sapply(Fmat_stat3, fmm@locus_index)
    expect_equal(mi_vec3, mi_ref3, tolerance = 1e-12)
})

test_that("LocusNMF: constructor from F_list", {
    cnt <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res  <- fitLocusNMF(cnt, ft, maxiter = 20)

    model2 <- LocusNMF(coef(res))
    expect_s4_class(model2, "LocusNMF")
    expect_equal(length(model2), 2L)
    expect_equal(ncol(model2), 0L)

    L_pred <- predict(model2, newx = cnt, type = "loadings")
    expect_equal(nrow(L_pred), 4L)
})

test_that("fitLocusNMF: init from LocusModel object", {
    cnt  <- .make_cell_cnt()
    feat <- unlist(lapply(cnt, rownames))
    ft   <- .make_ft_model(feat, colnames(cnt[[1]]), K = 2)
    res1 <- fitLocusNMF(cnt, ft, maxiter = 5)

    res2 <- fitLocusNMF(cnt, init = res1, maxiter = 5)

    expect_s4_class(res2, "LocusNMF")
    expect_equal(nrow(loadings(res2)), 4L)
    expect_equal(ncol(loadings(res2)), 2L)
    expect_equal(length(coef(res2)), 2L)
    expect_equal(rowSums(loadings(res2)), rep(1, 4), tolerance = 1e-10,
                 ignore_attr = TRUE)
    for (F_l in coef(res2)) {
        expect_equal(colSums(F_l), rep(1, 2), tolerance = 1e-10,
                     ignore_attr = TRUE)
    }
})
