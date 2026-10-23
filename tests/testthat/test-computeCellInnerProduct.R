# Explicit dense reference for testing the sparse covariance calculation.
.inner_product_reference_features <- function(cnt, covariance) {
    X <- t(as.matrix(cnt@counts))
    X <- X / pmax(rowSums(X), 1)
    if (covariance == "within_locus") {
        for (idx in cnt@locus_index) {
            X[, idx] <- sweep(X[, idx, drop = FALSE], 1,
                             rowMeans(X[, idx, drop = FALSE]), "-")
        }
    }
    X
}

test_that("exported cell inner products match explicit feature inner products", {
    cnt <- .make_cell_cnt_unequal_loci()[, c(6, 2, 5, 1, 4, 3)]
    for (covariance in c("total", "within_locus")) {
        X <- .inner_product_reference_features(cnt, covariance)
        K <- sclanternR::computeCellInnerProduct(cnt, covariance = covariance)
        expect_true(is.matrix(K))
        expect_equal(K, tcrossprod(X), tolerance = 1e-12)
        expect_identical(dimnames(K), rep(list(colnames(cnt)), 2))
        expect_equal(computeF2dist(K), as.matrix(stats::dist(X))^2,
                     tolerance = 1e-12)
    }
    expect_equal(computeCellInnerProduct(cnt),
                 computeCellInnerProduct(cnt, covariance = "total"))
})

test_that("truncated inner products match an independent feature SVD", {
    cnt <- .make_cell_cnt_unequal_loci()
    for (covariance in c("total", "within_locus")) {
        X <- .inner_product_reference_features(cnt, covariance)
        sv <- svd(X)
        original <- computeCellInnerProduct(cnt, covariance = covariance)
        for (k in 1:2) {
            coords <- sweep(sv$u[, seq_len(k), drop = FALSE], 2,
                            sv$d[seq_len(k)], "*")
            rownames(coords) <- colnames(cnt)
            K <- computeCellInnerProduct(cnt, covariance = covariance, rank = k)
            expect_equal(K, tcrossprod(coords), tolerance = 1e-12)
            expect_equal(computeF2dist(K), as.matrix(stats::dist(coords))^2,
                         tolerance = 1e-12)
            expect_equal(sum(svd(K, nu = 0, nv = 0)$d > 1e-10), k)
            expect_gte(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values),
                       -1e-12)
            expect_gt(max(abs(K - original)), 1e-4)
        }
        expect_identical(original, computeCellInnerProduct(
            cnt, covariance = covariance, rank = ncol(cnt)))
        expect_identical(original, computeCellInnerProduct(
            cnt, covariance = covariance, rank = ncol(cnt) + 5))
    }
})

test_that("rank reduction handles zero inner products and single cells", {
    cnt <- .make_cell_cnt_unequal_loci()
    empty <- cnt
    empty@counts@x[] <- 0
    for (covariance in c("total", "within_locus")) {
        K <- computeCellInnerProduct(empty, covariance = covariance, rank = 1)
        expect_true(all(is.finite(K)))
        expect_equal(unname(K), matrix(0, ncol(cnt), ncol(cnt)))
        single <- cnt[, 1]
        expect_equal(computeCellInnerProduct(single, covariance = covariance, rank = 1),
                     computeCellInnerProduct(single, covariance = covariance))
    }
    singleton <- cnt["chr3_300"]
    K <- computeCellInnerProduct(singleton, covariance = "within_locus", rank = 1)
    expect_equal(unname(K), matrix(0, ncol(cnt), ncol(cnt)))
})

test_that("cell inner product arguments are validated", {
    cnt <- .make_cell_cnt_unequal_loci()
    expect_error(computeCellInnerProduct(matrix(1, 2, 2)),
                 "cell_allele_counts must be a CellAlleleCounts object")
    expect_error(computeCellInnerProduct(cnt, covariance = "invalid"),
                 "'arg' should be one of")
    for (bad_rank in list(0, -1, 1.5, NA, Inf, NaN, c(1, 2), numeric(), "2", TRUE, 1i)) {
        expect_error(computeCellInnerProduct(cnt, rank = bad_rank),
                     "rank must be NULL or a positive integer")
    }
})
