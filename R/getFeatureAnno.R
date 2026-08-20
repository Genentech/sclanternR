
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames
.getFeatureAnno <- function(features) {
    df_anno <- str_match(
        features,
        '^(chr[^_]+)_(\\d+)_(.*)$'
    )
    colnames(df_anno) <- c('full', 'chr', 'pos', 'allele')
    df_anno <- as.data.frame(df_anno)
    df_anno$chr_pos <- with(df_anno, paste0(chr, '_', pos))
    df_anno$pos <- as.integer(df_anno$pos)

    df_anno <- column_to_rownames(df_anno, 'full')

    df_anno
}
