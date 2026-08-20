test_that("returns a CellAlleleCounts with correct dimensions", {
    expect_s4_class(list_cnt, "CellAlleleCounts")
    expect_equal(names(list_cnt), c("chr1_100", "chr2_200"))
    expect_equal(length(list_cnt), 2L)
    for (nm in names(list_cnt)) {
        mat <- list_cnt[[nm]]
        expect_true(is(mat, "Matrix") || is.matrix(mat))
        expect_equal(nrow(mat), 2L)  # 2 alleles
        expect_equal(ncol(mat), 2L)  # 2 cells
    }
})

test_that("values are non-negative and column sums equal UMI count per cell", {
    for (nm in names(list_cnt)) {
        mat <- as.matrix(list_cnt[[nm]])
        expect_true(all(mat >= 0))
        # each cell has 2 UMIs; each UMI contributes proportions summing to 1
        expect_equal(unname(colSums(mat)), c(2, 2), tolerance = 1e-10)
    }
})

test_that("[[ by index and by name return same result", {
    expect_equal(list_cnt[[1]], list_cnt[["chr1_100"]])
    expect_equal(list_cnt[[2]], list_cnt[["chr2_200"]])
})

test_that("[ subsetting preserves requested order", {
    sub <- list_cnt[c("chr2_200", "chr1_100")]
    expect_equal(names(sub), c("chr2_200", "chr1_100"))
    expect_equal(sub[["chr2_200"]], list_cnt[["chr2_200"]])
})

test_that("1D [ subsets loci (not cells)", {
    cnt <- .make_cell_cnt()
    sub <- cnt["chr1_100"]
    expect_s4_class(sub, "CellAlleleCounts")
    expect_equal(length(sub), 1L)
    expect_equal(names(sub), "chr1_100")
    expect_equal(ncol(sub), 4L)  # all cells preserved
    sub2 <- cnt[1]
    expect_equal(names(sub2), "chr1_100")
})

test_that("[ subsetting with logical vector works", {
    sub <- list_cnt[c(TRUE, FALSE)]
    expect_equal(length(sub), 1L)
    expect_equal(names(sub), "chr1_100")
})

test_that("coercion to CsparseMatrix returns features x cells sparse matrix", {
    mat <- as(list_cnt, "CsparseMatrix")
    expect_true(is(mat, "CsparseMatrix"))
    expect_equal(nrow(mat), 4L)  # 2 alleles x 2 loci
    expect_equal(ncol(mat), 2L)  # 2 cells
    expect_true(setequal(rownames(mat),
                         c("chr1_100_REF", "chr1_100_ALT",
                           "chr2_200_REF", "chr2_200_ALT")))
})

test_that("as.list roundtrips correctly", {
    lst <- as.list(list_cnt)
    expect_true(is.list(lst))
    expect_equal(names(lst), names(list_cnt))
    for (nm in names(lst)) {
        expect_equal(lst[[nm]], list_cnt[[nm]])
    }
})

test_that("lapply works on CellAlleleCounts", {
    result <- lapply(list_cnt, nrow)
    expect_true(is.list(result))
    expect_equal(names(result), names(list_cnt))
    expect_equal(result[["chr1_100"]], 2L)
    expect_equal(result[["chr2_200"]], 2L)
})

test_that("cell subsetting with [,j] works", {
    sub <- list_cnt[, "cellA"]
    expect_s4_class(sub, "CellAlleleCounts")
    expect_equal(length(sub), 2L)
    expect_equal(colnames(sub), "cellA")
    expect_equal(ncol(sub[[1]]), 1L)
})

test_that("2D subsetting [i,j] subsets both loci and cells", {
    sub <- list_cnt["chr1_100", "cellB"]
    expect_equal(length(sub), 1L)
    expect_equal(names(sub), "chr1_100")
    expect_equal(colnames(sub), "cellB")
})

test_that("cell_data accessor works", {
    cd <- cell_data(list_cnt)
    expect_s4_class(cd, "DataFrame")
    expect_equal(rownames(cd), c("cellA", "cellB"))
})

test_that("cell_data<- setter works", {
    cac <- list_cnt
    cell_data(cac) <- DataFrame(group = c("A", "B"),
                                row.names = c("cellA", "cellB"))
    expect_equal(colnames(cell_data(cac)), "group")
    expect_equal(cell_data(cac)$group, c("A", "B"))
})

test_that("CellAlleleCounts() constructor works with dense matrix", {
    cells <- paste0("cell", 1:3)
    mat <- matrix(c(5, 1, 3, 3, 0, 6), nrow = 2, ncol = 3,
                  dimnames = list(c("chr1_100_REF", "chr1_100_ALT"), cells))
    cac <- CellAlleleCounts(mat)
    expect_s4_class(cac, "CellAlleleCounts")
    expect_equal(length(cac), 1L)
    expect_equal(names(cac), "chr1_100")
    expect_equal(colnames(cac), cells)
    expect_equal(as.matrix(cac[["chr1_100"]]), mat)
})

test_that("CellAlleleCounts() constructor works with sparse CSC matrix", {
    cells <- paste0("cell", 1:3)
    mat <- matrix(c(5, 1, 3, 3, 0, 6), nrow = 2, ncol = 3,
                  dimnames = list(c("chr1_100_REF", "chr1_100_ALT"), cells))
    sparse <- as(mat, "CsparseMatrix")
    cac <- CellAlleleCounts(sparse)
    expect_s4_class(cac, "CellAlleleCounts")
    expect_equal(as.matrix(cac[["chr1_100"]]), mat)
})

test_that("CellAlleleCounts() derives locus_map from rownames", {
    cells <- paste0("cell", 1:2)
    mat <- matrix(c(5, 1, 3, 0, 2, 4), nrow = 3, ncol = 2,
                  dimnames = list(c("chr1_100_REF", "chr1_100_ALT",
                                    "chr2_200_G"), cells))
    cac <- CellAlleleCounts(mat)
    expect_equal(length(cac), 2L)
    expect_equal(names(cac), c("chr1_100", "chr2_200"))
})

test_that("CellAlleleCounts() accepts explicit locus_map and cell_data", {
    cells <- paste0("cell", 1:2)
    mat <- matrix(c(5, 1, 3, 0), nrow = 2, ncol = 2,
                  dimnames = list(c("chr1_100_REF", "chr1_100_ALT"), cells))
    lm <- DataFrame(locus = c("L1", "L1"), allele = c("R", "A"),
                     feature_name = c("chr1_100_REF", "chr1_100_ALT"))
    cd <- DataFrame(batch = c("X", "Y"), row.names = cells)
    cac <- CellAlleleCounts(mat, locus_map = lm, cell_data = cd)
    expect_equal(names(cac), "L1")
    expect_equal(cell_data(cac)$batch, c("X", "Y"))
})

test_that("CellAlleleCounts_from_list accepts cell_data", {
    cells <- paste0("cell", 1:4)
    l1 <- matrix(c(8, 0, 7, 1, 1, 7, 0, 8), nrow = 2, ncol = 4,
                 dimnames = list(c("chr1_100_REF", "chr1_100_ALT"), cells))
    cd <- DataFrame(batch = rep(c("X", "Y"), 2), row.names = cells)
    cac <- CellAlleleCounts_from_list(list(chr1_100 = l1), cell_data = cd)
    expect_equal(cell_data(cac)$batch, c("X", "Y", "X", "Y"))
})
