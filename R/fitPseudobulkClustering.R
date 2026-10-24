#' Construct pseudobulks clusters from a probabilistic clustering
#'
#' Run [fitCellClustering()], assign each cell to its most probable cluster,
#' and merge small clusters using Ward's minimum variance criterion on
#' the within-locus variance of the allele frequencies.
#'
#' @param cell_allele_counts A \code{\linkS4class{CellAlleleCounts}} object.
#' @param min_cells Positive integer giving the minimum number of cells per
#'     pseudobulk (default 20). With 1, only hard assignment is performed.
#' @param clustering An existing \code{CellClustering} result from
#'     [fitCellClustering()], or \code{NULL} to fit one. Cells must match the
#'     counts, but may be in a different order.
#' @param ... Arguments passed to [fitCellClustering()] when
#'     \code{clustering = NULL}, including \code{n_clusters}, \code{n_loci},
#'     \code{init}, \code{reestimate_fmm}, and \code{fmm_control}. Cannot be
#'     supplied with an existing \code{clustering} result.
#'
#' @details
#' Hard assignments use the largest responsibility in the final FMM's
#' \code{L} matrix; exact ties go to the first model column. Unoccupied
#' clusters are discarded. Merge distances use observed counts at
#' \code{clustering$keep_loci}.
#'
#' As with \code{hclust_covariance = "within_locus"} in
#' [fitFiniteMixtureModel()], the within-locus variance is computed by
#' subtracting the between-locus variance from the total
#' variance. This is equivalent to subtracting out the locus mean from
#' the allele frequencies before computing the variance.
#'
#' At each step, the pair with the smallest increase among pairs containing
#' at least one cluster smaller than \code{min_cells} is merged.
#'
#' The larger cluster's original model column index is retained as the group
#' ID (the smaller ID wins a size tie). Equal merge costs are resolved in
#' increasing order of the larger, then smaller, ID in each pair.
#'
#' @returns A list of class \code{PseudobulkClustering} containing:
#'   \describe{
#'     \item{clusters}{Named integer vector of final group IDs, in input
#'         cell order.}
#'     \item{initial_clusters}{Named integer vector of hard FMM assignments,
#'         in input cell order.}
#'     \item{cluster_sizes}{Named integer vector of final cell counts by ID.}
#'     \item{cluster_map}{Named integer vector mapping occupied initial
#'         cluster IDs to final IDs. Empty model clusters are omitted.}
#'     \item{merges}{Data frame with one row per merge: \code{from},
#'         \code{to}, \code{n_from}, \code{n_to} (sizes before merging),
#'         and \code{ward_cost} (increase in the sum of squares above).}
#'     \item{clustering}{The original \code{CellClustering} result.}
#'     \item{min_cells}{The requested minimum group size.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- fitPseudobulkClustering(counts, min_cells = 30,
#'                                   n_clusters = 20, n_loci = 100)
#' barcodes <- split(names(result$clusters), result$clusters)
#' assignments <- data.frame(cell = names(result$clusters),
#'                           cluster = unname(result$clusters))
#' write.table(assignments, "pseudobulk_clusters.tsv", sep = "\t",
#'             quote = FALSE, row.names = FALSE)
#'
#' # Try a different size threshold without refitting the model.
#' smaller <- fitPseudobulkClustering(counts, min_cells = 10,
#'                                    clustering = result$clustering)
#' }
#'
#' @seealso [fitCellClustering()], [computeCellInnerProduct()]
#' @export
fitPseudobulkClustering <- function(cell_allele_counts, min_cells = 20,
                                   clustering = NULL, ...) {
    if (!is(cell_allele_counts, "CellAlleleCounts")) {
        stop("cell_allele_counts must be a CellAlleleCounts object")
    }
    if (!is.numeric(min_cells) || is.complex(min_cells) ||
        length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 1 ||
        min_cells != floor(min_cells)) {
        stop("min_cells must be a positive integer")
    }
    cells <- colnames(cell_allele_counts)
    if (length(cells) == 0L || anyNA(cells) || any(!nzchar(cells)) ||
        anyDuplicated(cells)) {
        stop("cell_allele_counts must contain nonempty, unique cell identifiers")
    }
    if (is.null(clustering)) {
        clustering <- fitCellClustering(cell_allele_counts, ...)
    } else if (length(list(...)) > 0L) {
        stop("... cannot be supplied with an existing clustering result")
    }
    if (!inherits(clustering, "CellClustering") ||
        !is(clustering$fmm, "FiniteMixtureModel")) {
        stop("clustering must be a CellClustering result")
    }
    model_cells <- colnames(clustering$fmm)
    if (anyDuplicated(model_cells) || !setequal(model_cells, cells)) {
        stop("clustering and cell_allele_counts must contain the same cells")
    }
    loci <- clustering$keep_loci
    if (!is.character(loci) || length(loci) == 0L || anyNA(loci) ||
        anyDuplicated(loci) || !all(loci %in% names(cell_allele_counts))) {
        stop("clustering$keep_loci must name distinct loci present in cell_allele_counts")
    }
    L <- clustering$fmm@L[match(cells, model_cells), , drop = FALSE]
    if (ncol(L) == 0L || any(!is.finite(L)) || any(L < 0) ||
        any(rowSums(L) <= 0)) {
        stop("clustering responsibilities must be finite, nonnegative, and have positive row sums")
    }
    initial <- setNames(max.col(L, ties.method = "first"), cells)
    merged <- .merge_small_cell_clusters(cell_allele_counts[loci], initial,
                                         min_cells)
    structure(c(merged, list(initial_clusters = initial,
                             clustering = clustering, min_cells = min_cells)),
              class = "PseudobulkClustering")
}

.merge_small_cell_clusters <- function(counts, clusters, min_cells) {
    ids <- sort(unique(unname(clusters)))
    group <- match(clusters, ids)
    sizes <- tabulate(group, nbins = length(ids))
    cluster_map <- setNames(ids, ids)
    merges <- data.frame(from = integer(), to = integer(),
                         n_from = integer(), n_to = integer(),
                         ward_cost = numeric())

    if (length(ids) > 1L && any(sizes < min_cells)) {
        Y_t <- MatrixExtra::t_shallow(as(counts, "CsparseMatrix"))
        if (any(!is.finite(Y_t@x)) || any(Y_t@x < 0)) {
            stop("counts at selected loci must be finite and nonnegative")
        }
        totals <- Matrix::rowSums(Y_t)
        totals[totals == 0] <- 1
        freq <- Matrix::Diagonal(x = 1 / totals) %*% Y_t
        average <- Matrix::sparseMatrix(
            i = group, j = seq_along(group), x = 1 / sizes[group],
            dims = c(length(ids), length(group)))
        centroids <- average %*% freq
        gram <- .frequency_inner_product(centroids, counts@locus_1hot)

        while (length(ids) > 1L && any(sizes < min_cells)) {
            costs <- pmax(computeF2dist(gram), 0) *
                (outer(sizes, sizes, "*") / outer(sizes, sizes, "+"))
            eligible <- upper.tri(costs) &
                outer(sizes < min_cells, sizes < min_cells, "|")
            costs[!eligible] <- Inf
            pair <- arrayInd(which.min(costs), dim(costs))
            to <- pair[1L]
            from <- pair[2L]
            if (sizes[from] > sizes[to]) {
                to <- pair[2L]
                from <- pair[1L]
            }
            merges[nrow(merges) + 1L, ] <- list(
                ids[from], ids[to], sizes[from], sizes[to], costs[pair])

            # Update the Gram matrix of centroids by the cell-count-weighted
            # average, using the original diagonal and cross term.
            weight <- sizes[to] / (sizes[to] + sizes[from])
            new_diag <- weight^2 * gram[to, to] +
                (1 - weight)^2 * gram[from, from] +
                2 * weight * (1 - weight) * gram[to, from]
            new_row <- weight * gram[to, ] + (1 - weight) * gram[from, ]
            gram[to, ] <- gram[, to] <- new_row
            gram[to, to] <- new_diag
            cluster_map[cluster_map == ids[from]] <- ids[to]
            sizes[to] <- sizes[to] + sizes[from]
            gram <- gram[-from, -from, drop = FALSE]
            sizes <- sizes[-from]
            ids <- ids[-from]
        }
    }
    final <- setNames(unname(cluster_map[as.character(clusters)]), names(clusters))
    list(clusters = final, cluster_sizes = setNames(sizes, ids),
         cluster_map = cluster_map, merges = merges)
}
