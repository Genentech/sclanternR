library(Matrix)

# Minimal synthetic dataset: 2 loci, 2 cells, 2 UMIs per cell per locus,
# 2 alleles per locus. cellA strongly favors REF; cellB strongly favors ALT.
df_counts <- data.frame(
    locus  = c(rep("chr1_100", 8), rep("chr2_200", 8)),
    cell   = c("cellA","cellA","cellA","cellA",
               "cellB","cellB","cellB","cellB",
               "cellA","cellA","cellA","cellA",
               "cellB","cellB","cellB","cellB"),
    umi    = c("u1","u1","u2","u2","u3","u3","u4","u4",
               "u5","u5","u6","u6","u7","u7","u8","u8"),
    allele = c("REF","ALT","REF","ALT","REF","ALT","REF","ALT",
               "REF","ALT","REF","ALT","REF","ALT","REF","ALT"),
    count  = c(10,0,8,2,0,10,1,9,
               10,0,9,1,2,8,0,10)
)

list_cnt <- CellAlleleCounts_from_dataframe(df_counts, min_cells_per_feature = 1)

# Helper: minimal fake ft_model (list with $F and $L) for testing topic models.
# F is features x K (column-normalized); L is cells x K (row-normalized).
.make_ft_model <- function(feature_names, cell_names, K, seed = 42) {
    set.seed(seed)
    nf <- length(feature_names)
    nc <- length(cell_names)
    F_mat <- matrix(abs(rnorm(nf * K)) + 0.1, nrow = nf, ncol = K,
                    dimnames = list(feature_names, paste0("k", seq_len(K))))
    F_mat <- sweep(F_mat, 2, colSums(F_mat), "/")
    L_mat <- matrix(abs(rnorm(nc * K)) + 0.1, nrow = nc, ncol = K,
                    dimnames = list(cell_names, paste0("k", seq_len(K))))
    L_mat <- sweep(L_mat, 1, rowSums(L_mat), "/")
    list(F = F_mat, L = L_mat)
}

# Helper: minimal dense count matrices (2 loci, 4 cells, 2 alleles each).
# cells 1-2 favor REF; cells 3-4 favor ALT.
.make_cell_cnt <- function() {
    cells <- paste0("cell", 1:4)
    l1 <- matrix(c(8, 0, 7, 1, 1, 7, 0, 8), nrow = 2, ncol = 4,
                 dimnames = list(c("chr1_100_REF", "chr1_100_ALT"), cells))
    l2 <- matrix(c(9, 0, 8, 1, 0, 9, 1, 8), nrow = 2, ncol = 4,
                 dimnames = list(c("chr2_200_REF", "chr2_200_ALT"), cells))
    CellAlleleCounts_from_list(list(chr1_100 = l1, chr2_200 = l2))
}

# Helper: larger count matrices suitable for passing directly to fit_topic_model
# (5 loci, 20 cells, 2 alleles each -> 10 features x 20 cells).
.make_cell_cnt_large <- function() {
    set.seed(7)
    loci  <- c("chr1_100", "chr2_200", "chr3_300", "chr4_400", "chr5_500")
    cells <- paste0("cell", 1:20)
    mat_list <- lapply(setNames(loci, loci), function(loc) {
        mat <- matrix(rpois(40, lambda = 5), nrow = 2, ncol = 20,
                      dimnames = list(paste0(loc, c("_REF", "_ALT")), cells))
        mat
    })
    CellAlleleCounts_from_list(mat_list)
}
