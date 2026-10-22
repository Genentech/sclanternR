#' Collapse alleles consistent with RNA editing
#'
#' Sum allele counts into a chosen reference allele when the sequences differ
#' only by A>G or T>C substitutions. Loci with only one resulting allele are
#' removed.
#'
#' @param cac A `CellAlleleCounts` object.
#' @param all_transitions Logical; also collapse G>A and C>T substitutions,
#'     so that all transitions (A/G and C/T) are collapsed. Default FALSE.
#' @param reference Which allele to use as the reference: `"REF"` (default)
#'     uses the annotated REF allele; `"major"` uses the allele with the
#'     greatest total count across cells, regardless of REF/ALT annotation.
#'
#' @details
#' With `reference = "REF"`, `locus_map` must contain `allele_type` and every
#' locus must have exactly one allele labelled `REF`. Only alleles labelled
#' `ALT` can be merged into it. Missing annotations or a missing REF at any
#' locus cause an error directing the caller to use `reference = "major"`.
#' Multiple REF alleles at a locus are also an error.
#'
#' With `reference = "major"`, all other alleles are candidates for merging
#' into the major allele. Ties use the first allele in `locus_map`.
#'
#' Sequences must contain only A, C, G, and T (case insensitive). Indels and ambiguous
#' sequences are left uncollapsed.
#'
#' @returns A `CellAlleleCounts` object with collapsed counts and updated
#'     locus indices. The input object is not modified.
#'
#' @examples
#' counts <- matrix(c(5, 2, 1, 3, 4, 2), nrow = 3,
#'     dimnames = list(paste0("chr1_100_", c("A", "G", "C")), c("c1", "c2")))
#' cac <- CellAlleleCounts(counts)
#' collapse_rna_editing(cac, reference = "major")
#' cac@locus_map$allele_type <- c("REF", "ALT", "ALT")
#' collapse_rna_editing(cac)
#' collapse_rna_editing(cac, all_transitions = TRUE)
#'
#' @export
collapse_rna_editing <- function(cac, all_transitions = FALSE,
                                 reference = c("REF", "major")) {
    if (!is(cac, "CellAlleleCounts")) {
        stop("cac must be a CellAlleleCounts object")
    }
    if (!is.logical(all_transitions) || length(all_transitions) != 1L ||
        is.na(all_transitions)) {
        stop("all_transitions must be TRUE or FALSE")
    }
    reference <- match.arg(reference)

    locus_map <- cac@locus_map
    alleles <- toupper(as.character(locus_map$allele))
    if (reference == "REF") {
        if (!"allele_type" %in% colnames(locus_map)) {
            stop('REF annotation is missing (no allele_type column); ',
                 'use reference = "major".')
        }
        allele_type <- as.character(locus_map$allele_type)
    } else {
        totals <- Matrix::rowSums(cac@counts)
    }
    allowed <- c("AG", "TC")
    if (all_transitions) allowed <- c(allowed, "GA", "CT")
    target <- seq_len(nrow(locus_map))

    for (idx in cac@locus_index) {
        if (reference == "REF") {
            ref <- idx[which(allele_type[idx] == "REF")]
            if (length(ref) > 1L) {
                stop(sprintf("Multiple REF alleles at locus '%s'",
                             locus_map$locus[idx[1L]]))
            }
            if (length(ref) == 0L) {
                stop(sprintf('REF allele is missing at locus \'%s\'; use reference = "major".',
                             locus_map$locus[idx[1L]]))
            }
            alt <- idx[which(allele_type[idx] == "ALT")]
        } else {
            ref <- idx[which.max(totals[idx])]
            alt <- setdiff(idx, ref)
        }
        ref_seq <- alleles[ref]
        if (is.na(ref_seq) || !grepl("^[ACGT]+$", ref_seq)) next
        ref_bases <- strsplit(ref_seq, "", fixed = TRUE)[[1L]]
        for (a in alt) {
            alt_seq <- alleles[a]
            if (is.na(alt_seq) || !grepl("^[ACGT]+$", alt_seq) ||
                nchar(alt_seq) != length(ref_bases)) next
            alt_bases <- strsplit(alt_seq, "", fixed = TRUE)[[1L]]
            differs <- ref_bases != alt_bases
            if (any(differs) && all(paste0(ref_bases[differs],
                                           alt_bases[differs]) %in% allowed)) {
                target[a] <- ref
            }
        }
    }

    # Representatives stay in their original feature order, even if a REF
    # follows one of its ALT rows in the input matrix.
    representatives <- which(target == seq_along(target))
    remaining <- .build_locus_index(locus_map[representatives, , drop = FALSE])
    keep_loci <- names(remaining)[lengths(remaining) > 1L]
    representatives <- representatives[
        as.character(locus_map$locus[representatives]) %in% keep_loci]
    output_row <- match(target, representatives)
    keep <- which(!is.na(output_row))
    aggregation <- Matrix::sparseMatrix(
        i = output_row[keep], j = keep, x = rep(1, length(keep)),
        dims = c(length(representatives), nrow(locus_map)))
    counts <- as(aggregation %*% cac@counts, "RsparseMatrix")
    locus_map <- locus_map[representatives, , drop = FALSE]
    dimnames(counts) <- list(as.character(locus_map$feature_name),
                            rownames(cac@cell_data))
    .new_cac(counts, locus_map, cac@cell_data)
}
