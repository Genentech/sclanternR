test_that("matrix constructor supports empty axes without weakening name checks", {
    empty <- CellAlleleCounts(matrix(numeric(), 0, 0))
    expect_true(validObject(empty))
    expect_equal(ncol(empty), 0L)
    expect_equal(length(empty), 0L)
    expect_identical(colnames(empty), character())

    no_features <- CellAlleleCounts(matrix(numeric(), 0, 2,
        dimnames = list(NULL, c("cell1", "cell2"))))
    expect_true(validObject(no_features))
    expect_equal(colnames(no_features), c("cell1", "cell2"))
    expect_equal(unname(colSums(no_features)), c(0, 0))

    no_cells <- CellAlleleCounts(matrix(numeric(), 2, 0,
        dimnames = list(c("chr1_10_A", "chr1_10_C"), NULL)))
    expect_true(validObject(no_cells))
    expect_equal(length(no_cells), 1L)
    expect_equal(ncol(no_cells), 0L)
    expect_error(CellAlleleCounts(matrix(0, 1, 1)), "rownames")
    expect_error(CellAlleleCounts(matrix(0, 1, 1,
        dimnames = list("chr1_10_A", NULL))), "colnames")
})

test_that("empty long tables support filtering and RNA-editing collapse", {
    df <- data.frame(cell = character(), umi = character(), allele = character(),
                     locus = character(), count = integer(), allele_type = character())
    for (threshold in c(0, 2)) {
        counts <- CellAlleleCounts_from_dataframe(df, min_cells_per_feature = threshold)
        expect_true(validObject(counts))
        expect_identical(colnames(counts), character())
        expect_equal(dim(counts@counts), c(0L, 0L))
        expect_true("allele_type" %in% colnames(counts@locus_map))
        expect_true(validObject(collapse_rna_editing(counts)))
        expect_true(validObject(filter_cell_allele_counts(counts)))
    }
})

test_that("ambiguous-only cells remain as zero-count columns", {
    df <- data.frame(cell = rep("ambiguous", 2), umi = "AC", allele = c("A", "C"),
                     locus = "chr1_10", count = 1L, allele_type = c("REF", "ALT"))
    counts <- CellAlleleCounts_from_dataframe(df)
    expect_true(validObject(counts))
    expect_equal(colnames(counts), "ambiguous")
    expect_equal(length(counts), 0L)
    expect_equal(unname(colSums(counts)), 0)
    expect_true(validObject(collapse_rna_editing(counts)))

    supported <- data.frame(cell = c("c1", "c2", "c3"), umi = "AC",
        allele = c("A", "C", "C"), locus = "chr1_10", count = 1L,
        allele_type = c("REF", "ALT", "ALT"))
    counts <- CellAlleleCounts_from_dataframe(rbind(df, supported))
    expect_equal(colnames(counts), c("ambiguous", "c1", "c2", "c3"))
    expect_equal(unname(colSums(counts)), c(0, 1, 1, 1))
    expect_equal(length(counts), 1L)
    expect_equal(colnames(counts[, colSums(counts) > 0]), c("c1", "c2", "c3"))
})

test_that("empty counts can be subset repeatedly and all cells can be filtered", {
    empty <- CellAlleleCounts(matrix(numeric(), 0, 0))
    expect_true(validObject(empty[integer(), integer()]))
    counts <- CellAlleleCounts(matrix(c(1, 0, 0, 1), 2, 2,
        dimnames = list(c("chr1_10_A", "chr1_10_C"), c("c1", "c2"))))
    expect_true(validObject(counts[integer(), integer()]))
    for (x in list(empty, counts, counts[integer(), ])) {
        filtered <- filter_cell_allele_counts(x, min_counts_per_cell = 10)
        expect_true(validObject(filtered))
        expect_equal(ncol(filtered), 0L)
        expect_equal(length(filtered), 0L)
        expect_true(validObject(collapse_rna_editing(filtered)))
    }
})
