#' HAP1 cell-line allele counts
#'
#' Allele counts for the top 2,000 loci in a subset of the HAP1 cell-line data,
#' used in the package vignette to demonstrate lineage analysis.
#'
#' @format A sparse matrix of class `dgRMatrix`, with locus-allele pairs in
#'   rows and cells in columns. Row names use the format `CHROM_POS_ALLELE`.
#'   Column names combine the ground-truth clone and cell barcode, separated
#'   by an underscore. Entries are allele counts.
"hap1_subset"
