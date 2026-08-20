test_that("parses standard feature name into chr, pos, allele, chr_pos", {
    anno <- sclanternR:::.getFeatureAnno("chr1_100_REF")
    expect_equal(anno["chr1_100_REF", "chr"],     "chr1")
    expect_equal(anno["chr1_100_REF", "pos"],      100L)
    expect_equal(anno["chr1_100_REF", "allele"],  "REF")
    expect_equal(anno["chr1_100_REF", "chr_pos"], "chr1_100")
})

test_that("handles multi-character alleles and multiple inputs", {
    anno <- sclanternR:::.getFeatureAnno(c("chr2_200_ACGT", "chrX_999_DEL"))
    expect_equal(anno["chr2_200_ACGT", "allele"], "ACGT")
    expect_equal(anno["chrX_999_DEL",  "chr"],    "chrX")
    expect_equal(anno["chrX_999_DEL",  "pos"],    999L)
})
