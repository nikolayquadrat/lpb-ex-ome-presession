#' Build a per-gene "hubness" table from a signature (pathway) collection
#'
#' Hubness = how many signatures a gene belongs to. This is the covariate you
#' match on when building a null pool for per-gene enrichment-overlap
#' prioritisation: hub genes are members of many pathways and therefore overlap
#' any enriched set by chance, so a candidate is only interesting if it overlaps
#' MORE than genes of equal hubness.
#'
#' Two hubness measures are produced:
#'   n_signatures  -- raw count of signatures containing the gene.
#'   n_clusters    -- (optional) redundancy-collapsed count: signatures are
#'                    grouped into near-duplicate clusters (connected components
#'                    of the Jaccard >= `collapse_jaccard` graph) and this counts
#'                    the number of DISTINCT clusters the gene participates in.
#'                    This is the honest hubness for a redundant union like
#'                    MSigDB, where one biological pathway is mirrored across
#'                    Hallmark / GO / Reactome / etc. Match on THIS, not the raw
#'                    count, when your collection is redundant.
#'
#' @param pathway_list Named list: list("signature name" = c(genes...), ...).
#' @param all_genes   Optional character vector defining the gene universe you
#'                    will sample the null pool from (e.g. exome-testable genes).
#'                    Genes in `all_genes` that are in NO signature are included
#'                    with hubness 0 (a valid, important stratum: they can never
#'                    be members of an enriched pathway). If NULL, only genes
#'                    that appear in >=1 signature are listed. If given, the
#'                    table is restricted to `all_genes`.
#' @param collapse    Compute the redundancy-collapsed `n_clusters` column
#'                    (needs Matrix + igraph). Default FALSE (raw only).
#' @param collapse_jaccard Jaccard cutoff for merging near-duplicate signatures
#'                    into one cluster (default 0.5).
#' @param max_signatures_for_collapse Safety cap: skip collapse (with a warning)
#'                    above this many signatures, because the signature-signature
#'                    similarity can blow up in memory for very large unions.
#'                    Default 8000. Raise at your own risk, or run collapse on a
#'                    reduced collection / the enriched subset.
#' @param add_bins    Add `hub_percentile` and `hub_bin` columns (quantile bins
#'                    of the hubness used for matching). Default TRUE.
#' @param n_bins      Number of quantile bins (default 10). Uses n_clusters when
#'                    available, else n_signatures.
#' @param list_signatures Add a `signatures` column listing each gene's
#'                    signatures ("; "-joined). Off by default (can be huge).
#'
#' @return data.frame, one row per gene, sorted by hubness (desc).
build_hubness_table <- function(pathway_list,
                                all_genes = NULL,
                                collapse = FALSE,
                                collapse_jaccard = 0.5,
                                max_signatures_for_collapse = 8000,
                                add_bins = TRUE,
                                n_bins = 10,
                                list_signatures = FALSE) {

    stopifnot(is.list(pathway_list), length(pathway_list) > 0,
              !is.null(names(pathway_list)))

    # ---- normalise: unique genes per signature, drop empties ----
    sets <- lapply(pathway_list, function(g) unique(as.character(g)))
    sets <- sets[vapply(sets, length, integer(1)) > 0]
    if (length(sets) == 0) stop("All signatures are empty.")
    sig_names <- names(sets)
    set_sizes <- vapply(sets, length, integer(1))

    # ---- long (gene, signature) membership ----
    gene_col <- unlist(sets, use.names = FALSE)
    sig_col  <- rep(sig_names, times = set_sizes)

    # ---- raw hubness ----
    tab <- table(gene_col)
    hub <- data.frame(gene = names(tab),
                      n_signatures = as.integer(tab),
                      stringsAsFactors = FALSE)

    # ---- extend to the requested universe (0-hubness genes matter) ----
    if (!is.null(all_genes)) {
        all_genes <- unique(as.character(all_genes))
        add <- setdiff(all_genes, hub$gene)
        if (length(add) > 0)
            hub <- rbind(hub, data.frame(gene = add, n_signatures = 0L))
        hub <- hub[hub$gene %in% all_genes, , drop = FALSE]
    }

    # ---- optional per-gene signature listing ----
    if (list_signatures) {
        by_gene <- split(sig_col, factor(gene_col, levels = hub$gene))
        hub$signatures <- vapply(by_gene, function(s)
            if (length(s) == 0) "" else paste(sort(unique(s)), collapse = "; "),
            character(1))
    }

    # ---- optional redundancy-collapsed hubness ----
    if (collapse) {
        ok <- requireNamespace("Matrix", quietly = TRUE) &&
              requireNamespace("igraph", quietly = TRUE)
        if (!ok) {
            warning("collapse requested but Matrix/igraph not available; skipping n_clusters.")
        } else if (length(sets) > max_signatures_for_collapse) {
            warning(sprintf(paste0("collapse skipped: %d signatures exceed ",
                    "max_signatures_for_collapse=%d. Run collapse on a reduced ",
                    "collection (e.g. the enriched subset), or raise the cap."),
                    length(sets), max_signatures_for_collapse))
        } else {
            genes_u <- unique(gene_col)
            gi <- match(gene_col, genes_u)
            si <- match(sig_col, sig_names)
            # gene x signature sparse incidence
            I <- Matrix::sparseMatrix(i = gi, j = si, x = 1,
                                      dims = c(length(genes_u), length(sig_names)))
            sizes <- Matrix::colSums(I)
            inter <- Matrix::crossprod(I)                 # sig x sig intersection counts
            inter <- methods::as(inter, "TsparseMatrix")
            ii <- inter@i + 1L; jj <- inter@j + 1L; xx <- inter@x
            up <- ii < jj                                 # upper triangle, no diagonal
            ii <- ii[up]; jj <- jj[up]; xx <- xx[up]
            jac <- xx / (sizes[ii] + sizes[jj] - xx)
            e   <- which(jac >= collapse_jaccard)
            g_graph <- igraph::make_empty_graph(n = length(sig_names), directed = FALSE)
            if (length(e) > 0)
                g_graph <- igraph::add_edges(g_graph, rbind(ii[e], jj[e]))
            sig_cluster <- igraph::components(g_graph)$membership   # length = #signatures
            # collapsed hubness per gene = #distinct clusters among its signatures
            gclust <- sig_cluster[si]
            n_clusters_by_gene <- tapply(gclust, gene_col, function(x) length(unique(x)))
            hub$n_clusters <- as.integer(n_clusters_by_gene[hub$gene])
            hub$n_clusters[is.na(hub$n_clusters)] <- 0L
            attr(hub, "n_signature_clusters") <- max(sig_cluster)
        }
    }

    # ---- sort by hubness (collapsed if present, else raw) ----
    key <- if ("n_clusters" %in% names(hub)) hub$n_clusters else hub$n_signatures
    hub <- hub[order(-key), , drop = FALSE]

    # ---- bins for matched sampling ----
    if (add_bins) {
        bkey <- if ("n_clusters" %in% names(hub)) hub$n_clusters else hub$n_signatures
        hub$hub_percentile <- rank(bkey, ties.method = "average") / length(bkey)
        br <- unique(stats::quantile(bkey, probs = seq(0, 1, length.out = n_bins + 1),
                                     na.rm = TRUE))
        # If hubness is highly tied (many zeros), quantile breaks collapse; guard.
        if (length(br) < 2) {
            hub$hub_bin <- 1L
        } else {
            hub$hub_bin <- as.integer(cut(bkey, breaks = br,
                                          include.lowest = TRUE, labels = FALSE))
        }
    }

    rownames(hub) <- NULL
    hub
}