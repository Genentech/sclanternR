
#' @include CellAlleleCounts-class.R
#' @importFrom methods setMethod setAs show is as
#' @importFrom stats setNames
#' @importFrom Matrix t

#' @rdname CellAlleleCounts-class
#' @exportMethod length
setMethod("length", "CellAlleleCounts", function(x) {
    length(x@locus_index)
})

#' @rdname CellAlleleCounts-class
#' @exportMethod names
setMethod("names", "CellAlleleCounts", function(x) {
    names(x@locus_index)
})

#' @rdname CellAlleleCounts-class
#' @importFrom BiocGenerics ncol
#' @exportMethod ncol
setMethod("ncol", "CellAlleleCounts", function(x) {
    nrow(x@cell_data)
})

#' @rdname CellAlleleCounts-class
#' @importFrom BiocGenerics colnames
#' @exportMethod colnames
setMethod("colnames", "CellAlleleCounts", function(x, do.NULL = TRUE, prefix = "col") {
    rownames(x@cell_data)
})

#' @rdname CellAlleleCounts-class
#' @exportMethod colSums
setMethod("colSums", "CellAlleleCounts", function(x) {
    colSums(x@counts)
})

#' @rdname CellAlleleCounts-class
#' @exportMethod [[
setMethod("[[", signature(x = "CellAlleleCounts", i = "numeric"),
    function(x, i) {
        loci <- unique(as.character(x@locus_map$locus))
        if (i < 1 || i > length(loci)) stop("subscript out of bounds")
        .extract_locus(x, loci[i])
    }
)

#' @rdname CellAlleleCounts-class
#' @exportMethod [[
setMethod("[[", signature(x = "CellAlleleCounts", i = "character"),
    function(x, i) {
        .extract_locus(x, i)
    }
)

.extract_locus <- function(x, locus_name) {
    idx <- x@locus_index[[locus_name]]
    if (is.null(idx)) stop(sprintf("locus '%s' not found", locus_name))
    sub <- x@counts[idx, , drop = FALSE]
    rownames(sub) <- as.character(x@locus_map$feature_name[idx])
    colnames(sub) <- rownames(x@cell_data)
    sub
}

#' @rdname CellAlleleCounts-class
#' @exportMethod [
setMethod("[", signature(x = "CellAlleleCounts", i = "ANY", j = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        if (!missing(i)) {
            x <- .subset_loci_any(x, i)
        }
        if (!missing(j)) {
            x <- .subset_cells_any(x, j)
        }
        x
    }
)

#' @rdname CellAlleleCounts-class
#' @exportMethod [
setMethod("[", signature(x = "CellAlleleCounts", i = "ANY", j = "missing"),
    function(x, i, j, ..., drop = TRUE) {
        .subset_loci_any(x, i)
    }
)

#' @rdname CellAlleleCounts-class
#' @exportMethod [
setMethod("[", signature(x = "CellAlleleCounts", i = "missing", j = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        .subset_cells_any(x, j)
    }
)

.subset_loci_any <- function(x, i) {
    if (is.logical(i)) {
        loci <- unique(as.character(x@locus_map$locus))
        i <- loci[i]
    } else if (is.numeric(i)) {
        loci <- unique(as.character(x@locus_map$locus))
        i <- loci[i]
    }
    .subset_loci(x, i)
}

.subset_loci <- function(x, loci) {
    idx <- unlist(x@locus_index[loci], use.names = FALSE)

    sub_counts <- x@counts[idx, , drop = FALSE]
    if (!is(sub_counts, "RsparseMatrix")) {
        sub_counts <- as(sub_counts, "RsparseMatrix")
    }
    .new_cac(sub_counts,
             x@locus_map[idx, , drop = FALSE],
             x@cell_data)
}

.subset_cells_any <- function(x, j) {
    cell_names <- rownames(x@cell_data)
    if (is.character(j)) {
        cidx <- match(j, cell_names)
    } else if (is.logical(j)) {
        cidx <- which(j)
    } else {
        cidx <- j
    }

    sub_counts <- x@counts[, cidx, drop = FALSE]
    if (!is(sub_counts, "RsparseMatrix")) {
        sub_counts <- as(sub_counts, "RsparseMatrix")
    }
    .new_cac(sub_counts,
             x@locus_map,
             x@cell_data[cidx, , drop = FALSE])
}

#' @rdname CellAlleleCounts-class
#' @exportMethod show
setMethod("show", "CellAlleleCounts", function(object) {
    n_loci <- length(object)
    n_cells <- nrow(object@cell_data)
    n_features <- nrow(object@counts)
    cat(sprintf("CellAlleleCounts: %d loci, %d cells, %d features\n",
                n_loci, n_cells, n_features))
    if (ncol(object@cell_data) > 0) {
        cat(sprintf("  cell_data columns: %s\n",
                    paste(colnames(object@cell_data), collapse = ", ")))
    }
})

setAs("CellAlleleCounts", "CsparseMatrix", function(from) {
    out <- as(from@counts, "CsparseMatrix")
    rownames(out) <- as.character(from@locus_map$feature_name)
    colnames(out) <- rownames(from@cell_data)
    out
})

#' @rdname CellAlleleCounts-class
#' @export
as.list.CellAlleleCounts <- function(x, ...) {
    loci <- names(x@locus_index)
    setNames(base::lapply(loci, function(l) x[[l]]), loci)
}

#' @rdname CellAlleleCounts-class
#' @importFrom BiocGenerics lapply
#' @exportMethod lapply
setMethod("lapply", "CellAlleleCounts", function(X, FUN, ...) {
    loci <- names(X@locus_index)
    setNames(base::lapply(loci, function(l) FUN(X[[l]], ...)), loci)
})

#' Get cell metadata from a CellAlleleCounts object
#' @rdname CellAlleleCounts-class
#' @export
cell_data <- function(x) {
    x@cell_data
}

#' Set cell metadata on a CellAlleleCounts object
#' @rdname CellAlleleCounts-class
#' @export
`cell_data<-` <- function(x, value) {
    x@cell_data <- value
    validObject(x)
    x
}
