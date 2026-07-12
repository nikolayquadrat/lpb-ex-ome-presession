#' Enrichment map (emap-style network) directly from an fgsea result
#'
#' Builds the same picture as enrichplot::emapplot -- pathways as nodes,
#' gene-overlap as edges, force-directed layout, communities as redundancy
#' groups -- but straight from an fgsea() / fgseaMultilevel() result, with no
#' re-run through clusterProfiler and no S4 object coercion. Returns a ggplot
#' object (via ggraph), so you can add titles / theming / ggsave() downstream.
#'
#' Edge similarity can be computed two ways:
#'   - "gene_sets"    : Jaccard/overlap over the FULL gene sets (structural
#'                      redundancy). Requires `pathway_list`.
#'   - "leading_edge" : Jaccard/overlap over the fgsea leadingEdge genes
#'                      (data-driven redundancy -- pathways connected only if
#'                      the SAME genes drove their enrichment in THIS ranking).
#'                      Uses the fgsea result's leadingEdge column; no
#'                      pathway_list needed. This is closer to what
#'                      collapsePathways() collapses on.
#'
#' @param fgsea_res   data.frame/data.table from fgsea. Needs column `pathway`;
#'                    uses `padj`, `NES`, `size` and `leadingEdge` when present.
#' @param pathway_list Named list of gene-set character vectors (only needed
#'                    when edge_source = "gene_sets").
#' @param edge_source "gene_sets" or "leading_edge" (see above).
#' @param similarity  "jaccard" (|A∩B|/|A∪B|) or "overlap" (|A∩B|/min(|A|,|B|)).
#' @param padj_cutoff Keep pathways with padj <= this (default 0.05). Set NULL
#'                    to keep all rows passed in (e.g. if you already filtered
#'                    to mainPathways).
#' @param pathways    Optional explicit character vector of pathway names to
#'                    display (e.g. your collapsePathways mainPathways). Applied
#'                    after padj_cutoff; overrides top_n selection order.
#' @param top_n       Keep only the top_n pathways by padj (after padj_cutoff).
#'                    NULL = keep all.
#' @param edge_cutoff Minimum similarity for an edge to be drawn (default 0.3).
#' @param cluster_method "louvain" (default) colours/groups nodes by community;
#'                    "none" skips community detection.
#' @param color_by    "NES" (diverging), "neglog10padj", or "cluster".
#' @param label       Label nodes? (default TRUE)
#' @param max_label   If set, label only this many nodes (the most significant
#'                    by padj); others unlabelled to reduce clutter. NULL =
#'                    label all shown nodes.
#' @param label_size  Text size for node labels (default 2.5).
#' @param layout      ggraph/igraph layout name (default "fr" = Fruchterman-
#'                    Reingold). "kk", "drl", "graphopt" also work.
#' @param node_size_range Point-size range mapped from pathway `size`.
#' @param seed        RNG seed for the (stochastic) layout, for reproducibility.
#'
#' @return A ggplot object.
fgsea_emap <- function(fgsea_res,
                       pathway_list   = NULL,
                       edge_source    = c("gene_sets", "leading_edge"),
                       similarity     = c("jaccard", "overlap"),
                       padj_cutoff    = 0.05,
                       pathways       = NULL,
                       top_n          = NULL,
                       edge_cutoff    = 0.3,
                       cluster_method = c("louvain", "none"),
                       color_by       = c("NES", "neglog10padj", "cluster"),
                       label          = TRUE,
                       max_label      = NULL,
                       label_size     = 2.5,
                       layout         = "fr",
                       node_size_range = c(3, 12),
                       seed           = 1) {
    
    edge_source    <- match.arg(edge_source)
    similarity     <- match.arg(similarity)
    cluster_method <- match.arg(cluster_method)
    color_by       <- match.arg(color_by)
    
    for (pkg in c("igraph", "ggraph", "ggplot2")) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required")
    }
    
    res <- as.data.frame(fgsea_res)
    if (!"pathway" %in% names(res)) stop("fgsea_res must have a 'pathway' column")
    
    # ---- select which pathways to show -------------------------------------
    if (!is.null(padj_cutoff) && "padj" %in% names(res)) {
        res <- res[!is.na(res$padj) & res$padj <= padj_cutoff, , drop = FALSE]
    }
    if (!is.null(pathways)) {
        res <- res[res$pathway %in% pathways, , drop = FALSE]
    }
    if (!is.null(top_n) && nrow(res) > top_n && "padj" %in% names(res)) {
        res <- res[order(res$padj), , drop = FALSE][seq_len(top_n), , drop = FALSE]
    }
    if (nrow(res) == 0) stop("No pathways left to plot after filtering.")
    
    sel <- res$pathway
    
    # ---- assemble the gene sets that define the edges ----------------------
    if (edge_source == "gene_sets") {
        if (is.null(pathway_list))
            stop("edge_source='gene_sets' needs `pathway_list`.")
        missing <- setdiff(sel, names(pathway_list))
        if (length(missing) > 0)
            stop(sprintf("%d selected pathways are not in pathway_list (e.g. %s)",
                         length(missing), paste(head(missing, 3), collapse = ", ")))
        sets <- lapply(pathway_list[sel], function(x) unique(as.character(x)))
    } else {
        if (!"leadingEdge" %in% names(res))
            stop("edge_source='leading_edge' needs a 'leadingEdge' column.")
        le <- res$leadingEdge
        # leadingEdge may be a list-column (usual) or a delimited string
        # (e.g. after read from TSV). Normalise to a list of char vectors.
        if (is.list(le)) {
            sets <- lapply(le, function(x) unique(as.character(x)))
        } else {
            sets <- lapply(as.character(le), function(s)
                unique(strsplit(s, "[/,; ]+")[[1]]))
        }
        names(sets) <- sel
    }
    
    n <- length(sets)
    
    # ---- pairwise similarity matrix ----------------------------------------
    sim_fun <- if (similarity == "jaccard") {
        function(a, b) { u <- length(union(a, b)); if (u == 0) 0 else length(intersect(a, b)) / u }
    } else {
        function(a, b) { m <- min(length(a), length(b)); if (m == 0) 0 else length(intersect(a, b)) / m }
    }
    M <- matrix(0, n, n, dimnames = list(sel, sel))
    if (n >= 2) {
        for (i in seq_len(n - 1)) for (j in (i + 1):n) {
            s <- sim_fun(sets[[i]], sets[[j]])
            M[i, j] <- M[j, i] <- s
        }
    }
    
    # ---- build graph (all selected pathways are nodes; edges above cutoff) --
    adj <- (M >= edge_cutoff) & (M > 0)
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    # carry similarity as edge weight (for edge width)
    if (igraph::ecount(g) > 0) {
        el <- igraph::as_edgelist(g, names = TRUE)
        igraph::E(g)$weight <- vapply(seq_len(nrow(el)),
                                      function(k) M[el[k, 1], el[k, 2]], numeric(1))
    }
    
    # ---- node attributes from the fgsea result ------------------------------
    idx <- match(igraph::V(g)$name, res$pathway)
    if ("NES"  %in% names(res)) igraph::V(g)$NES  <- res$NES[idx]
    if ("padj" %in% names(res)) igraph::V(g)$padj <- res$padj[idx]
    if ("size" %in% names(res)) {
        igraph::V(g)$size <- res$size[idx]
    } else {
        igraph::V(g)$size <- vapply(sets[igraph::V(g)$name], length, integer(1))
    }
    igraph::V(g)$neglog10padj <- if ("padj" %in% names(res))
        -log10(pmax(res$padj[idx], .Machine$double.xmin)) else NA_real_
    
    # ---- community detection (redundancy groups) ---------------------------
    if (cluster_method == "louvain") {
        comm <- igraph::cluster_louvain(g)          # isolated nodes -> own community
        igraph::V(g)$cluster <- as.factor(igraph::membership(comm))
    } else {
        igraph::V(g)$cluster <- factor(rep(1, igraph::vcount(g)))
    }
    
    # ---- which nodes to label ----------------------------------------------
    igraph::V(g)$do_label <- label
    if (label && !is.null(max_label) && "padj" %in% names(res)) {
        ord <- order(res$padj[idx])
        keep_lab <- igraph::V(g)$name[ord[seq_len(min(max_label, igraph::vcount(g)))]]
        igraph::V(g)$do_label <- igraph::V(g)$name %in% keep_lab
    }
    
    # ---- plot ---------------------------------------------------------------
    set.seed(seed)
    p <- ggraph::ggraph(g, layout = layout)
    
    if (igraph::ecount(g) > 0) {
        p <- p + ggraph::geom_edge_link(
            ggplot2::aes(width = weight), alpha = 0.25, colour = "grey55",
            show.legend = FALSE)
    }
    
    node_aes <- switch(color_by,
                       NES          = ggplot2::aes(size = size, colour = NES),
                       neglog10padj = ggplot2::aes(size = size, colour = neglog10padj),
                       cluster      = ggplot2::aes(size = size, colour = cluster))
    p <- p + ggraph::geom_node_point(node_aes)
    
    p <- p + switch(color_by,
                    NES          = ggplot2::scale_colour_gradient2(
                        low = "#3b6fb0", mid = "grey85", high = "#c0392b",
                        midpoint = 0, name = "NES"),
                    neglog10padj = ggplot2::scale_colour_viridis_c(
                        name = expression(-log[10]~padj), option = "C"),
                    cluster      = ggplot2::scale_colour_hue(name = "cluster"))
    
    p <- p + ggraph::scale_edge_width(range = c(0.2, 2.5)) +
        ggplot2::scale_size(range = node_size_range, name = "set size") +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "right")
    
    if (label) {
        p <- p + ggraph::geom_node_text(
            ggplot2::aes(label = ifelse(do_label, name, "")),
            repel = TRUE, size = label_size, max.overlaps = Inf)
    }
    p
}