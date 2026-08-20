

#' @export
#' @importFrom forcats fct_drop
#' @importFrom dplyr arrange mutate group_by summarize if_else
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom pals polychrome kovesi.rainbow trubetskoy
heatmapLocusMatrixList <- function(
    list_mat,
    cell_meta,
    cell_meta_colors = NULL,
    cluster_rows = TRUE,
    gene_names = NULL,
    colors = c('white','yellow','red','black'),
    color_limits = NULL,
    color_title = NULL,
    max_allele_length = NULL,
    column_names_max_cm=10
) {
    if (!is.null(gene_names) && is.null(names(gene_names))) {
        names(gene_names) <- names(list_mat)
    }

    list_mat <- lapply(list_mat, as.matrix)

    mat_plt <- t(do.call(rbind, list_mat))

    df_colanno <- .getFeatureAnno(colnames(mat_plt))

    df_colanno$chr <- factor(df_colanno$chr,
                             levels=paste0('chr', c(1:22, 'X', 'Y')))
    df_colanno$chr <- fct_drop(df_colanno$chr)

    df_colanno <- arrange(df_colanno, chr, pos)
    df_colanno <- mutate(df_colanno, log10pos=log10(pos), sqrt_pos=sqrt(pos))

    mat_plt <- mat_plt[,rownames(df_colanno)]

    if (is.null(cell_meta_colors)) {
        right_anno <- rowAnnotation(
            df=cell_meta[rownames(mat_plt),,drop=FALSE]
        )
    } else {
        right_anno <- rowAnnotation(
            df=cell_meta[rownames(mat_plt),,drop=FALSE],
            col=cell_meta_colors
        )
    }

    chr_pos_levels <- unique(df_colanno$chr_pos)
    df_colanno$chr_pos <- factor(df_colanno$chr_pos, levels=chr_pos_levels)
    chr_pos_colors <- rep_len(trubetskoy(), length(chr_pos_levels))
    names(chr_pos_colors) <- chr_pos_levels

    if (is.null(color_limits)) {
        color_limits <- c(min(mat_plt, na.rm=T), max(mat_plt, na.rm=T))
    }

    col <- colorRamp2(
        breaks=seq(color_limits[1], color_limits[2], length.out=length(colors)),
        colors=colors
    )

    if (is.null(color_title)) {
        if (color_limits[1] >= 0 & color_limits[2] <= 1) {
            color_title <- "p"
        } else {
            color_title <- "x"
        }
    }

    colsplit <- as.character(df_colanno$chr_pos)
    if (!is.null(gene_names)) {
        colsplit <- paste0(colsplit, "_", gene_names[colsplit])
    }
    colsplit <- as.character(colsplit)
    colsplit <- factor(colsplit, levels=unique(colsplit))

    if (!is.null(max_allele_length)) {
        loc_max_len <- summarize(group_by(df_colanno, chr_pos),
                                 max_len = max(nchar(.data$allele)))
        loc_max_len <- setNames(loc_max_len$max_len, loc_max_len$chr_pos)
        loc_max_len <- loc_max_len[df_colanno$chr_pos]

        df_colanno$max_len <- pmax(
            ceiling(max_allele_length / 2),
            max_allele_length - (loc_max_len - nchar(df_colanno$allele))
        )

        #df_colanno$allele <- str_trunc(df_colanno$allele, df_colanno$max_len)
        df_colanno$allele <- if_else(
            nchar(df_colanno$allele) <= df_colanno$max_len,
            df_colanno$allele,
            paste0(substr(df_colanno$allele, 1, df_colanno$max_len), "...")
        )
    }

    Heatmap(
        mat_plt,
        col=col,
        column_split = colsplit,
        #column_title = FALSE,
        column_title_rot = 45,
        cluster_rows = cluster_rows,
        clustering_method_rows = 'ward.D2',
        cluster_columns=FALSE,
        show_row_names=FALSE,
        right_annotation=right_anno,
        top_annotation=ComplexHeatmap::columnAnnotation(
            df=df_colanno[,c('chr'),drop=FALSE],
            col=list(
                chr=`names<-`(
                    polychrome(length(levels(df_colanno$chr))),
                    levels(df_colanno$chr)
                ),
                log10pos=colorRamp2(
                    breaks=seq(min(df_colanno$log10pos),
                               max(df_colanno$log10pos),length.out=10),
                    colors=kovesi.rainbow(10)
                ),
                chr_pos = chr_pos_colors
            ),
            show_legend=c(TRUE,TRUE,FALSE)
        ),
        column_names_max_height=unit(column_names_max_cm,'cm'),
        column_labels=df_colanno[, 'allele'],
        heatmap_legend_param=list(title=color_title)
    )
}
