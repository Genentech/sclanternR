#' @useDynLib sclanternR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats approx
#' @importFrom utils modifyList
#' @importFrom dplyr n
NULL

# Column names evaluated in dplyr/dtplyr data masks and ggplot2 aesthetics.
utils::globalVariables(c(
    ".data", "allele", "cell", "chr", "chr_pos", "count", "cutoff",
    "eigenvalues", "feature", "frac", "locus", "method", "pos", "umi"
))
