.pseudobulk_test_clustering <- function(counts, groups = NULL, L = NULL,
                                       keep_loci = names(counts)) {
    if (is.null(L)) {
        L <- outer(groups, seq_len(max(groups)), `==`) * 1
    }
    rownames(L) <- colnames(counts)
    K <- ncol(L)
    F_list <- lapply(as.list(counts), function(x) {
        matrix(1 / nrow(x), nrow(x), K,
               dimnames = list(rownames(x), NULL))
    })
    structure(list(fmm = FiniteMixtureModel(F_list, rep(1 / K, K), L),
                   keep_loci = keep_loci), class = "CellClustering")
}

.pseudobulk_test_coordinates <- function(counts) {
    x <- t(as.matrix(counts@counts))
    totals <- rowSums(x)
    totals[totals == 0] <- 1
    x <- x / totals
    for (idx in counts@locus_index) {
        x[, idx] <- x[, idx, drop = FALSE] - rowMeans(x[, idx, drop = FALSE])
    }
    x
}

test_that("pseudobulk wrapper fits and reuses CellClustering results", {
    cnt <- .make_cell_cnt_large()
    result <- fitPseudobulkClustering(
        cnt, min_cells = 4, n_clusters = 3, n_loci = 2,
        reestimate_fmm = 1, fmm_control = list(maxiter = 5))
    expected <- fitCellClustering(
        cnt, n_clusters = 3, n_loci = 2,
        reestimate_fmm = 1, fmm_control = list(maxiter = 5))
    expect_s3_class(result, "PseudobulkClustering")
    expect_identical(result$clustering, expected)
    expect_identical(result$initial_clusters,
                     setNames(max.col(expected$fmm@L, ties.method = "first"),
                              colnames(cnt)))
    expect_equal(sum(result$cluster_sizes), ncol(cnt))
    expect_true(all(result$cluster_sizes >= 4))

    local_mocked_bindings(fitCellClustering = function(...) stop("unexpected fit"),
                          .package = "sclanternR")
    reused <- fitPseudobulkClustering(cnt, 4, clustering = expected)
    expect_identical(reused, result)
    reordered <- fitPseudobulkClustering(cnt[, rev(colnames(cnt))], 4,
                                        clustering = expected)
    expect_identical(reordered$clusters, result$clusters[rev(colnames(cnt))])
    expect_identical(reordered$initial_clusters,
                     result$initial_clusters[rev(colnames(cnt))])
    expect_equal(reordered$merges, result$merges, tolerance = 1e-12)
    expect_identical(reordered$clustering, expected)
})

test_that("hard assignments resolve exact ties and omit empty clusters", {
    cnt <- .make_cell_cnt()
    L <- rbind(c(0.5, 0.5, 0), c(0.1, 0.9, 0),
               c(0.2, 0.8, 0), c(0.6, 0.4, 0))
    clustering <- .pseudobulk_test_clustering(cnt, L = L)
    set.seed(123)
    rng <- .Random.seed
    result <- fitPseudobulkClustering(cnt, 1, clustering)
    expect_identical(.Random.seed, rng)
    expect_identical(result$clusters, setNames(c(1L, 2L, 2L, 1L), colnames(cnt)))
    expect_identical(result$initial_clusters, result$clusters)
    expect_identical(result$cluster_map, c(`1` = 1L, `2` = 2L))
    expect_identical(result$cluster_sizes, c(`1` = 2L, `2` = 2L))
    expect_equal(nrow(result$merges), 0L)
    expect_identical(result$clustering, clustering)
    expect_identical(fitPseudobulkClustering(cnt, 2, clustering)$clusters,
                     result$clusters)
})

test_that("Ward weighting can select a farther but smaller neighbor", {
    # Group 2 is farther from the singleton than group 3, but its smaller
    # size means absorbing the singleton costs less in total variance.
    ref <- c(50, rep(60, 2), rep(59, 10))
    cnt <- CellAlleleCounts(rbind(chr1_100_REF = ref, chr1_100_ALT = 100 - ref),
                            cell_data = S4Vectors::DataFrame(
                                row.names = paste0("cell", seq_along(ref))))
    clustering <- .pseudobulk_test_clustering(cnt, c(1, rep(2, 2), rep(3, 10)))
    result <- fitPseudobulkClustering(cnt, 2, clustering)
    expect_equal(result$merges$from, 1L)
    expect_equal(result$merges$to, 2L)
    expect_equal(result$merges$ward_cost, (1 * 2 / 3) * 2 * (0.6 - 0.5)^2)
    expect_identical(result$cluster_sizes, c(`2` = 3L, `3` = 10L))
    expect_equal(unname(result$clusters), c(rep(2L, 3), rep(3L, 10)))
})

test_that("within-locus subtraction changes the merge chosen by total variance", {
    # Groups 1 and 2 differ only by locus means (including a singleton locus).
    # Group 3 is closer in total frequency distance, but has different
    # within-locus residuals. Only group 1 is undersized.
    x <- rbind(c(8, 6, 3, 2, 1, 4), c(9, 7, 2, 1, 0, 5),
               c(9, 7, 2, 1, 0, 5), c(9, 5, 3, 2, 1, 4),
               c(9, 5, 3, 2, 1, 4))
    dimnames(x) <- list(paste0("cell", 1:5),
                        rownames(.make_cell_cnt_unequal_loci()@counts))
    cnt <- CellAlleleCounts(t(x))
    clustering <- .pseudobulk_test_clustering(cnt, c(1, 2, 2, 3, 3))
    result <- fitPseudobulkClustering(cnt, 2, clustering)
    total <- as.matrix(stats::dist(x / rowSums(x)))
    expect_lt(total[1, 4], total[1, 2])
    expect_equal(result$merges$from, 1L)
    expect_equal(result$merges$to, 2L)
    expect_equal(result$merges$ward_cost, 0, tolerance = 1e-12)
})

test_that("every merge minimizes the actual increase in within-locus SSE", {
    cnt <- .make_cell_cnt_unequal_loci()
    # Unequal depth and locus sizes, and an empty cell. Centroids must average
    # normalized cells, not pooled counts or renormalized cluster frequencies.
    cnt <- CellAlleleCounts(cnt@counts %*%
                            Matrix::Diagonal(x = c(1, 3, 5, 7, 1, 2)),
                            cell_data = cnt@cell_data)
    clustering <- .pseudobulk_test_clustering(cnt, seq_len(ncol(cnt)))
    result <- fitPseudobulkClustering(cnt, 3, clustering)
    coords <- .pseudobulk_test_coordinates(cnt)
    sse <- function(groups) {
        sum(vapply(split(seq_along(groups), groups), function(idx) {
            residual <- sweep(coords[idx, , drop = FALSE], 2,
                              colMeans(coords[idx, , drop = FALSE]), "-")
            sum(residual^2)
        }, numeric(1)))
    }
    groups <- result$initial_clusters
    expect_gte(nrow(result$merges), 2L)
    for (i in seq_len(nrow(result$merges))) {
        merge <- result$merges[i, ]
        sizes <- table(groups)
        candidates <- combn(as.integer(names(sizes)), 2)
        eligible <- apply(candidates, 2, function(pair) {
            any(sizes[as.character(pair)] < 3)
        })
        increases <- apply(candidates[, eligible, drop = FALSE], 2, function(pair) {
            next_groups <- groups
            next_groups[next_groups == pair[1]] <- pair[2]
            sse(next_groups) - sse(groups)
        })
        expect_equal(merge$ward_cost, min(increases), tolerance = 1e-12)
        expect_equal(merge$n_from, unname(sizes[as.character(merge$from)]),
                     ignore_attr = TRUE)
        expect_equal(merge$n_to, unname(sizes[as.character(merge$to)]),
                     ignore_attr = TRUE)
        next_groups <- groups
        next_groups[next_groups == merge$from] <- merge$to
        expect_equal(merge$ward_cost, sse(next_groups) - sse(groups),
                     tolerance = 1e-12)
        groups <- next_groups
    }
    expect_identical(result$clusters, groups)
    expect_true(all(result$cluster_sizes >= 3))
    expect_identical(unname(result$cluster_map[as.character(result$initial_clusters)]),
                     unname(result$clusters))
})

test_that("merging uses selected loci even when the FMM contains all loci", {
    cnt <- .make_cell_cnt_unequal_loci()
    selected <- names(cnt)[2]
    clustering <- .pseudobulk_test_clustering(cnt, c(1, 2, 2, 3, 3, 4),
                                             keep_loci = selected)
    actual <- fitPseudobulkClustering(cnt, 3, clustering)
    expected <- fitPseudobulkClustering(cnt[selected], 3, clustering)
    expect_identical(actual, expected)

    # Changes outside the selected loci cannot affect the output.
    counts <- as.list(cnt)
    excluded <- setdiff(names(cnt), selected)
    for (locus in excluded) counts[[locus]] <- counts[[locus]] * 1000
    changed <- CellAlleleCounts_from_list(counts)
    expect_identical(fitPseudobulkClustering(changed, 3, clustering), actual)
})

test_that("small groups can combine and equal merge costs are deterministic", {
    cnt <- .make_cell_cnt()
    # All cells have exactly the same globally normalized frequencies.
    cnt@counts <- as(Matrix::Matrix(1, nrow = 4, ncol = 4, sparse = TRUE,
                                   dimnames = dimnames(cnt@counts)), "RsparseMatrix")
    clustering <- .pseudobulk_test_clustering(cnt, 1:4)
    result <- fitPseudobulkClustering(cnt, 3, clustering)
    expect_identical(result$merges$from, c(2L, 3L, 4L))
    expect_identical(result$merges$to, rep(1L, 3))
    expect_equal(result$merges$ward_cost, rep(0, 3), tolerance = 1e-12)
    expect_identical(result$cluster_sizes, c(`1` = 4L))
    expect_identical(fitPseudobulkClustering(cnt, 100, clustering)$clusters,
                     result$clusters)

    single <- cnt[, 1]
    one <- fitPseudobulkClustering(single, 20,
                                  .pseudobulk_test_clustering(single, 1L))
    expect_identical(one$cluster_sizes, c(`1` = 1L))
    expect_equal(nrow(one$merges), 0L)
})

test_that("pseudobulk arguments and reusable results are validated", {
    cnt <- .make_cell_cnt()
    clustering <- .pseudobulk_test_clustering(cnt, c(1, 1, 2, 2))
    local_mocked_bindings(fitCellClustering = function(...) stop("unexpected fit"),
                          .package = "sclanternR")
    expect_error(fitPseudobulkClustering(list()), "must be a CellAlleleCounts")
    for (value in list(0, -1, 0.5, NA, Inf, "2", c(1, 2), NULL, 1i, TRUE)) {
        expect_error(fitPseudobulkClustering(cnt, min_cells = value),
                     "min_cells must be a positive integer")
    }
    expect_error(fitPseudobulkClustering(cnt[, integer()]), "cell identifiers")
    expect_error(fitPseudobulkClustering(cnt, clustering = list()),
                 "must be a CellClustering result")
    expect_error(fitPseudobulkClustering(cnt, clustering = clustering, n_loci = 1),
                 "cannot be supplied")
    expect_error(fitPseudobulkClustering(cnt[, 1:3], clustering = clustering),
                 "must contain the same cells")
    bad <- clustering
    bad$keep_loci <- "missing"
    expect_error(fitPseudobulkClustering(cnt, clustering = bad), "keep_loci")
    bad <- clustering
    bad$fmm@L[1, ] <- NA_real_
    expect_error(fitPseudobulkClustering(cnt, clustering = bad), "responsibilities")
    bad$fmm@L[1, ] <- 0
    expect_error(fitPseudobulkClustering(cnt, clustering = bad), "responsibilities")
})
