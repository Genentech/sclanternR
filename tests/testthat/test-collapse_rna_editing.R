.editing_counts <- function(alleles, counts = NULL, allele_type = NULL) {
    locus <- rep(paste0("chr1_", seq_along(alleles)), lengths(alleles))
    allele <- unlist(alleles, use.names = FALSE)
    features <- paste0(locus, "_", allele)
    if (is.null(counts)) {
        counts <- matrix(rep(c(8, 2), length.out = length(allele) * 2),
                         nrow = length(allele))
    }
    dimnames(counts) <- list(features, c("cellB", "cellA"))
    lm <- DataFrame(locus = locus, allele = allele, feature_name = features,
                    extra = seq_along(allele))
    if (!is.null(allele_type)) lm$allele_type <- allele_type
    CellAlleleCounts(counts, lm,
        DataFrame(batch = c("b", "a"), row.names = colnames(counts)))
}

test_that("dataframe import retains aligned REF/ALT labels through filtering", {
    df <- df_counts
    df$allele_type <- df$allele
    df$allele <- ifelse(df$allele_type == "REF", "A", "G")
    # Each allele occurs in only one cell: disable the default two-cell cutoff.
    cac <- CellAlleleCounts_from_dataframe(df[nrow(df):1, ], 1)
    expect_equal(cac@locus_map$allele_type, rep(c("REF", "ALT"), 2))
    expect_equal(cac@locus_map$allele, rep(c("A", "G"), 2))
    expect_equal(length(collapse_rna_editing(cac)), 0L)
    expect_false("allele_type" %in% colnames(list_cnt@locus_map))

    # A missing label does not override an available label for that feature.
    df$allele_type[1] <- NA_character_
    expect_equal(CellAlleleCounts_from_dataframe(df, 0)@locus_map$allele_type,
                 rep(c("REF", "ALT"), 2))
    df$allele_type[1] <- "ALT"
    expect_error(CellAlleleCounts_from_dataframe(df, 0),
                 "Conflicting allele_type")
})

test_that("collapsing sums cell counts into annotated REF and preserves metadata", {
    cac <- .editing_counts(list(c("G", "A", "C"), c("T", "C"), "A"),
        matrix(c(10, 1, 3, 6, 2, 9, 4, 2, 7, 1, 5, 8), nrow = 6),
        c("ALT", "REF", "ALT", "REF", "ALT", "REF"))
    original <- cac
    out <- collapse_rna_editing(cac)
    expect_equal(names(out), "chr1_1")
    expect_equal(as.matrix(out@counts),
        matrix(c(11, 3, 6, 7), nrow = 2,
               dimnames = list(c("chr1_1_A", "chr1_1_C"), c("cellB", "cellA"))))
    expect_equal(out@locus_map, cac@locus_map[c(2, 3), ])
    expect_equal(cell_data(out), cell_data(cac))
    expect_equal(unname(colSums(out)), unname(colSums(cac[["chr1_1"]])))
    expect_equal(out@locus_index, list(chr1_1 = 1:2))
    expect_equal(as.matrix(out@locus_1hot),
        matrix(1L, nrow = 1, ncol = 2,
               dimnames = list("chr1_1", c("chr1_1_A", "chr1_1_C"))))
    expect_s4_class(out@counts, "RsparseMatrix")
    expect_true(validObject(out))
    expect_equal(cac, original)
    expect_equal(collapse_rna_editing(out), out)
})

test_that("default and all-transition modes distinguish every substitution", {
    pairs <- expand.grid(ref = c("A", "C", "G", "T"),
                         alt = c("A", "C", "G", "T"),
                         stringsAsFactors = FALSE)
    pairs <- pairs[pairs$ref != pairs$alt, ]
    alleles <- Map(c, pairs$ref, pairs$alt)
    cac <- .editing_counts(alleles, allele_type = rep(c("REF", "ALT"), 12))
    changes <- paste0(pairs$ref, pairs$alt)
    expect_equal(names(collapse_rna_editing(cac)),
                 names(cac)[!changes %in% c("AG", "TC")])
    expect_equal(names(collapse_rna_editing(cac, all_transitions = TRUE)),
                 names(cac)[!changes %in% c("AG", "TC", "GA", "CT")])
})

test_that("explicit major reference uses total counts with deterministic ties", {
    cac <- .editing_counts(list(c("G", "A", "C"), c("G", "A"), c("T", "C")),
        matrix(c(2, 1, 1, 4, 2, 3, 1, 0, 5, 1, 4, 1, 0, 2), nrow = 7))
    expect_error(collapse_rna_editing(cac),
                 'REF annotation is missing.*use reference = "major"')
    out <- collapse_rna_editing(cac, reference = "major")
    # A wins by total counts at locus 1; G wins at locus 2; T wins the tie at 3.
    expect_equal(names(out), c("chr1_1", "chr1_2"))
    expect_equal(as.matrix(out[["chr1_1"]]),
                 matrix(c(3, 1, 5, 1), nrow = 2,
                     dimnames = list(c("chr1_1_A", "chr1_1_C"), colnames(cac))))
    expect_equal(out[["chr1_2"]], cac[["chr1_2"]])
    expect_equal(names(collapse_rna_editing(cac, TRUE, reference = "major")),
                 "chr1_1")
    expect_false("allele_type" %in% colnames(out@locus_map))
})

test_that("major mode overrides REF annotation without relabelling metadata", {
    cac <- .editing_counts(list(c("A", "G", "C")),
        matrix(c(1, 5, 1, 2, 8, 1), nrow = 3), c("REF", "ALT", "ALT"))
    expect_equal(nrow(collapse_rna_editing(cac)@counts), 2L)
    expect_equal(collapse_rna_editing(cac, reference = "major"), cac)
    out <- collapse_rna_editing(cac, all_transitions = TRUE, reference = "major")
    expect_equal(out@locus_map, cac@locus_map[2:3, ])
    expect_equal(as.matrix(out@counts),
        matrix(c(6, 1, 10, 1), nrow = 2,
               dimnames = list(c("chr1_1_G", "chr1_1_C"), colnames(cac))))
    expect_equal(cell_data(out), cell_data(cac))
})

test_that("sequence alleles collapse only when all differences are allowed", {
    alleles <- list(c("AAT", "GAC"), c("AA", "GG"), c("AT", "GA"),
                    c("A", "AG"), c("AN", "GN"), c("a", "g"),
                    c("GC", "AT"), c("REF", "ALT"))
    cac <- .editing_counts(alleles, allele_type = rep(c("REF", "ALT"), 8))
    expect_equal(names(collapse_rna_editing(cac)), names(cac)[c(3:5, 7:8)])
    expect_equal(names(collapse_rna_editing(cac, TRUE)), names(cac)[c(3:5, 8)])
})

test_that("every locus must have REF unless major mode is explicit", {
    cac <- .editing_counts(list(c("A", "G"), c("A", "G"), c("A", "G")),
        allele_type = c("ALT", "ALT", NA, NA, "REF", NA))
    expect_error(collapse_rna_editing(cac),
                 'REF allele is missing at locus \'chr1_1\'; use reference = "major"')
    expect_error(collapse_rna_editing(cac["chr1_2"]),
                 'REF allele is missing.*use reference = "major"')
    cac@locus_map$allele_type[1] <- "REF"
    expect_error(collapse_rna_editing(cac),
                 'REF allele is missing at locus \'chr1_2\'; use reference = "major"')
    expect_length(collapse_rna_editing(cac, reference = "major"), 0)
    # An unlabelled candidate is not merged when a REF is available.
    expect_equal(collapse_rna_editing(cac["chr1_3"]), cac["chr1_3"])
    cac@locus_map$allele_type[1:2] <- "REF"
    expect_error(collapse_rna_editing(cac), "Multiple REF alleles at locus 'chr1_1'")
})

test_that("empty results remain valid and retain cells and annotations", {
    cac <- .editing_counts(list(c("A", "G")), allele_type = c("REF", "ALT"))
    out <- collapse_rna_editing(cac)
    expect_length(out, 0)
    expect_equal(dim(out@counts), c(0L, 2L))
    expect_equal(dim(out@locus_1hot), c(0L, 0L))
    expect_equal(colnames(out), colnames(cac))
    expect_equal(colnames(out@locus_map), colnames(cac@locus_map))
    expect_equal(cell_data(out), cell_data(cac))
    expect_true(validObject(out))
    expect_equal(collapse_rna_editing(out), out)
    expect_equal(length(collapse_rna_editing(cac[, integer()])), 0L)
})

test_that("feature filtering preserves annotations, cells, and locus order", {
    cac <- .editing_counts(list(c("A", "G"), c("C", "T")),
        matrix(c(1, 2, 3, 0, 4, 5, 6, 0), nrow = 4),
        c("REF", "ALT", "REF", "ALT"))
    cac <- cac[c("chr1_2", "chr1_1")]
    out <- filter_cell_allele_counts(cac, min_cells_per_feature = 1)
    expect_equal(out, cac["chr1_1"])
    empty <- filter_cell_allele_counts(cac, min_cells_per_feature = 3)
    expect_length(empty, 0)
    expect_equal(cell_data(empty), cell_data(cac))
})

test_that("collapse options are validated", {
    expect_error(collapse_rna_editing(list()), "CellAlleleCounts")
    for (arg in list(NA, NULL, c(TRUE, FALSE), 1, "TRUE")) {
        expect_error(collapse_rna_editing(list_cnt, arg), "TRUE or FALSE")
    }
    expect_error(collapse_rna_editing(list_cnt, reference = "unknown"), "arg")
})

test_that("filtering protects rare and zero-count REF alleles by default", {
    cac <- .editing_counts(list(c("A", "G", "C"), c("T", "C")),
        matrix(c(0, 4, 3, 1, 6, 0, 5, 2, 0, 10), nrow = 5),
        c("REF", "ALT", "ALT", "REF", "ALT"))
    expect_equal(filter_cell_allele_counts(cac), cac)
    # Cell filtering can remove all support for REF; its row still survives.
    out <- filter_cell_allele_counts(cac, min_counts_per_cell = 15,
                                     min_cells_per_feature = 1)
    expect_equal(out, cac[, "cellA"])
    legacy <- filter_cell_allele_counts(cac, keep_ref = FALSE)
    expect_equal(names(legacy), "chr1_1")
    expect_equal(legacy@locus_map, cac@locus_map[2:3, ])
    expect_error(collapse_rna_editing(legacy),
                 'REF allele is missing.*use reference = "major"')
    expect_equal(collapse_rna_editing(legacy, reference = "major"), legacy)
    cac@locus_map$allele_type <- NULL
    expect_equal(filter_cell_allele_counts(cac),
                 filter_cell_allele_counts(cac, keep_ref = FALSE))
})

test_that("unknown annotations are filtered and REF-only loci are removed", {
    cac <- .editing_counts(list(c("A", "G", "C"), c("T", "C")),
        matrix(c(0, 4, 0, 2, 0, 0, 5, 0, 3, 0), nrow = 5),
        c("REF", "ALT", NA, "REF", "ALT"))
    out <- filter_cell_allele_counts(cac)
    expect_equal(names(out), "chr1_1")
    expect_equal(out@locus_map, cac@locus_map[1:2, ])
    for (arg in list(NA, NULL, c(TRUE, FALSE), 1, "TRUE")) {
        expect_error(filter_cell_allele_counts(cac, keep_ref = arg),
                     "keep_ref must be TRUE or FALSE")
    }
})

test_that("dataframe import exposes REF retention and legacy filtering", {
    df <- data.frame(cell = c("c1", "c1", "c2", "c1", "c2"),
                     umi = paste0("u", 1:5), locus = "chr1_1",
                     allele = c("A", "G", "G", "C", "C"),
                     allele_type = c("REF", rep("ALT", 4)), count = 1)
    cac <- CellAlleleCounts_from_dataframe(df)
    expect_equal(cac@locus_map$allele_type, c("REF", "ALT", "ALT"))
    expect_equal(nrow(collapse_rna_editing(cac)@counts), 2L)
    legacy <- CellAlleleCounts_from_dataframe(df, keep_ref = FALSE)
    expect_equal(legacy, filter_cell_allele_counts(cac, keep_ref = FALSE))
    expect_error(collapse_rna_editing(legacy),
                 'REF allele is missing.*use reference = "major"')
})
