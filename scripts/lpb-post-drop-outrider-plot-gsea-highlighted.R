#' Plot a GSEA results table with tier-gene highlighting + leading edge
#'
#' Mirrors the layout of fgsea's plotGseaTable.
#'
#' @param pathways      Named list of gene-set vectors (gene symbols matching
#'                      names(stats)). Same as plotGseaTable's first arg.
#' @param stats         Named numeric vector of ranking statistics (e.g.
#'                      z-scores), names = gene symbols.
#' @param fgseaRes      The data.table / data.frame returned by fgsea() or
#'                      fgseaMultilevel(). Must have columns: pathway, NES,
#'                      padj, and (for leading-edge coloring) leadingEdge as
#'                      a list-column of character vectors per pathway. If
#'                      leadingEdge is missing, leading-edge coloring is
#'                      silently disabled.
#' @param highlight_genes Character vector of gene symbols to mark in the
#'                      highlight color and list in the right-hand column.
#'                      Genes not present in `stats` are silently ignored
#'                      (with a count printed).
#' @param highlight_color Color for the highlight-set rank-ticks (default "red").
#' @param highlight_lwd Line width for highlighted ticks (default 0.8).
#' @param leading_edge_color Color for ticks of pathway genes in the leading
#'                      edge (default "#2ca02c", a green). Set to NULL to
#'                      disable leading-edge coloring and use a single color
#'                      for all non-highlight pathway genes.
#' @param leading_edge_lwd Line width for leading-edge ticks (default 0.5,
#'                      slightly thicker than regular pathway-gene ticks).
#' @param highlight_max_listed Maximum number of highlight-gene symbols to list
#'                      in the right-hand column per pathway. Default 8. Genes
#'                      in the leading edge get an asterisk suffix.
#' @param padj_threshold Significance cutoff on the padj column (default 0.05).
#'                      Pathways that do NOT reach significance (padj >=
#'                      padj_threshold, or padj is NA) have ONLY their name in
#'                      the "Pathway" column de-emphasised (drawn in
#'                      `nonsig_color`). Everything else -- rank-ticks, NES,
#'                      padj and GOI text -- keeps its normal colour. Set
#'                      `nonsig_color` to NULL to disable this entirely. To grey
#'                      the OTHER side instead (i.e. the significant pathway
#'                      names), see the one-line note in the row loop where
#'                      `is_sig` is computed -- flip `<` to `>=`.
#' @param nonsig_color  Colour used to grey out non-significant pathway NAMES
#'                      (default "grey70"). NULL disables the feature so every
#'                      pathway name is drawn in black regardless of padj.
#' @param colwidths     Column widths (relative units): c(pathway_name,
#'                      gene_ranks, NES, padj, hits).
#' @param gseaParam     Same as fgsea::plotGseaTable: stats are pre-scaled by
#'                      `sign(x) * |x|^gseaParam` before ranking. The transformed
#'                      stat also becomes the tick height.
#' @param render        If TRUE (default), draw the table on the current graphics
#'                      device. If FALSE, return the gtable without drawing.
#' @param ...           Reserved.
#'
#' @return A gtable object (invisibly when render=TRUE).
plotGseaTableHighlighted <- function(pathways,
                                      stats,
                                      fgseaRes,
                                      highlight_genes,
                                      highlight_color = "#e66c2c",
                                      highlight_lwd = 0.8,
                                      leading_edge_color = "#2ca6e6",
                                      leading_edge_lwd = 0.5,
                                      highlight_max_listed = 8,
                                      padj_threshold = 0.05,
                                      nonsig_color = "grey70",
                                      colwidths = c(5, 4, 0.8, 1.2, 2.5),
                                      gseaParam = 1,
                                      render = TRUE) {

    if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("ggplot2 is required")
    if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra is required")
    if (!requireNamespace("grid",      quietly = TRUE)) stop("grid is required")

    # ---- normalise highlight set to genes actually in the ranking ----
    highlight_genes <- unique(as.character(highlight_genes))
    in_ranking <- highlight_genes %in% names(stats)
    # if (any(!in_ranking)) {
    #     message(sprintf("[plotGseaTableHighlighted] %d/%d highlight genes are not in `stats` and will be ignored",
    #                     sum(!in_ranking), length(highlight_genes)))
    # }
    highlight_genes <- highlight_genes[in_ranking]

    # ---- replicate plotGseaTable's rank transform ----
    # fgsea's plotGseaTable scales ranking stats via |x|^gseaParam keeping sign;
    # then ranks genes by that scaled stat (decreasing) and reports each gene's
    # position in the ranked vector. We also use statsAdj as the TICK HEIGHT,
    # so genes with strong positive z-scores get tall up-ticks, strong negative
    # get tall down-ticks, and middling genes get short ticks.
    statsAdj <- stats
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    # Order from largest to smallest (so the top of the ranked list is at left).
    ord <- order(statsAdj, decreasing = TRUE)
    ranked_genes <- names(statsAdj)[ord]
    rank_of <- setNames(seq_along(ranked_genes), ranked_genes)
    n_genes <- length(ranked_genes)

    # Tick-height y-axis range, shared across all per-pathway panels so the
    # vertical scale is comparable from row to row. Use the symmetric max
    # so positive and negative ticks fan out equally from y=0.
    y_max_abs <- max(abs(statsAdj), na.rm = TRUE)
    # Add a small margin so the tallest ticks aren't flush with the panel edge.
    y_lim <- c(-y_max_abs, y_max_abs) * 1.05

    # ---- leading-edge handling ----
    # fgsea returns leadingEdge as a list-column. Convert to a named list
    # keyed by pathway name for fast lookup. If the column is absent, all
    # leading-edge sets are empty and the leading_edge_color is unused.
    if ("leadingEdge" %in% colnames(fgseaRes)) {
        le_by_pathway <- setNames(
            as.list(fgseaRes$leadingEdge),
            fgseaRes$pathway
        )
    } else {
        message("[plotGseaTableHighlighted] fgseaRes has no 'leadingEdge' column; ",
                "leading-edge coloring disabled.")
        le_by_pathway <- list()
    }

    # Quick lookup: is a gene in the highlight set?
    highlight_set <- setNames(rep(TRUE, length(highlight_genes)), highlight_genes)

    # ---- build one row of grobs per pathway ----
    # Each row has:
    #   1. pathway name (textGrob)
    #   2. rank-bars ggplot with black ticks + red ticks for highlighted genes
    #   3. NES value (textGrob)
    #   4. padj value (textGrob)
    #   5. comma-separated highlight-gene list (textGrob)
    pathway_names <- names(pathways)
    if (is.null(pathway_names) || any(pathway_names == "")) {
        stop("`pathways` must be a NAMED list (names become the row labels)")
    }

    # Match fgsea rows to pathway list by name
    fgsea_idx <- match(pathway_names, fgseaRes$pathway)
    if (any(is.na(fgsea_idx))) {
        stop(sprintf("Pathways missing from fgseaRes: %s",
                     paste(pathway_names[is.na(fgsea_idx)], collapse = ", ")))
    }

    # Helper: build the rank-bars panel for one pathway.
    # `pathway_name` is needed to look up the leading-edge set. Tick colours
    # are ALWAYS the full palette; the padj-based de-emphasis applies only to
    # the pathway-name text in the "Pathway" column (see the row loop below).
    make_ticks_grob <- function(pathway_name, pathway_genes) {
        # Genes from the pathway that are present in the ranking
        present <- pathway_genes[pathway_genes %in% ranked_genes]
        if (length(present) == 0) {
            ggplot_obj <- ggplot2::ggplot() +
                ggplot2::xlim(0, n_genes) +
                ggplot2::ylim(y_lim) +
                ggplot2::theme_void()
            return(ggplot2::ggplotGrob(ggplot_obj))
        }

        # Build a data frame for each gene: rank, height, and category
        # (regular / leading-edge / highlight). Categories are mutually
        # exclusive in the plotting layers: a gene that is BOTH in the
        # leading edge AND in highlight_genes is drawn as a highlight
        # (the user-specified emphasis wins).
        le_set <- le_by_pathway[[pathway_name]]
        if (is.null(le_set)) le_set <- character(0)

        is_hl <- present %in% highlight_genes
        is_le <- (present %in% le_set) & !is_hl
        is_reg <- !is_hl & !is_le

        # Construct per-category data frames. The y-end is statsAdj at this gene
        # (so up-ticks for positively-ranked genes, down-ticks for negatively-
        # ranked ones); y=0 is the baseline.
        mk_df <- function(genes) {
            if (length(genes) == 0) return(NULL)
            data.frame(
                rank = unname(rank_of[genes]),
                yend = unname(statsAdj[genes]),
                stringsAsFactors = FALSE
            )
        }
        df_reg <- mk_df(present[is_reg])
        df_le  <- mk_df(present[is_le])
        df_hl  <- mk_df(present[is_hl])

        # Tick colours: always the full palette (de-emphasis is name-only).
        reg_col <- "black"
        hl_col  <- highlight_color
        # Leading-edge layer: the configured LE colour, or (LE colouring
        # disabled via NULL) fall back to black at the regular line width.
        if (!is.null(leading_edge_color)) {
            le_col <- leading_edge_color;   le_draw_lwd <- leading_edge_lwd
        } else {
            le_col <- "black";              le_draw_lwd <- 0.3
        }

        # Build the plot. theme_void hides axes/grids.
        p <- ggplot2::ggplot() +
            ggplot2::xlim(0, n_genes) +
            ggplot2::ylim(y_lim) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))

        # Layer 1: regular pathway-gene ticks (background)
        if (!is.null(df_reg)) {
            p <- p + ggplot2::geom_segment(
                data = df_reg,
                ggplot2::aes(x = rank, xend = rank, y = 0, yend = yend),
                linewidth = 0.3, color = reg_col
            )
        }
        # Layer 2: leading-edge ticks (middle)
        if (!is.null(df_le)) {
            p <- p + ggplot2::geom_segment(
                data = df_le,
                ggplot2::aes(x = rank, xend = rank, y = 0, yend = yend),
                linewidth = le_draw_lwd, color = le_col
            )
        }
        # Layer 3: highlight-gene ticks (foreground). Drawn LAST so they sit on
        # top of anything else, even when packed densely.
        if (!is.null(df_hl)) {
            p <- p + ggplot2::geom_segment(
                data = df_hl,
                ggplot2::aes(x = rank, xend = rank, y = 0, yend = yend),
                linewidth = highlight_lwd, color = hl_col
            )
        }
        ggplot2::ggplotGrob(p)
    }

    # Helper: format a numeric padj value the way plotGseaTable does
    fmt_padj <- function(p) {
        if (is.na(p)) return("NA")
        if (p == 0)  return("<1e-300")
        sprintf("%.1e", p)
    }
    fmt_nes <- function(n) {
        if (is.na(n)) "NA" else sprintf("%.2f", n)
    }

    # Helper: list the highlight genes present in a pathway, ordered by rank.
    # Genes that are also in the pathway's leading edge get a "*" suffix.
    highlight_listing <- function(pathway_name, pathway_genes) {
        present <- pathway_genes[pathway_genes %in% highlight_genes]
        if (length(present) == 0) return("")
        # Order by rank in the ranked list (top of the list first)
        present <- present[order(rank_of[present])]

        # Mark leading-edge tier genes
        le_set <- le_by_pathway[[pathway_name]]
        if (is.null(le_set)) le_set <- character(0)
        labels <- ifelse(present %in% le_set,
                          paste0(present, "*"),
                          present)
        if (length(labels) > highlight_max_listed) {
            labels <- c(head(labels, highlight_max_listed), "...")
        }
        paste(labels, collapse = ", ")
    }

    # ---- assemble grobs as a grid ----
    n_path <- length(pathway_names)

    # Header row (5 columns: pathway, gene-ranks, NES, padj, hits)
    header_grobs <- list(
        grid::textGrob("Pathway", x = 1, hjust = 1, gp = grid::gpar(fontface = "bold")),
        grid::textGrob("Gene ranks",     gp = grid::gpar(fontface = "bold")),
        grid::textGrob("NES",            gp = grid::gpar(fontface = "bold")),
        grid::textGrob("padj",           gp = grid::gpar(fontface = "bold")),
        grid::textGrob("GOI", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold"))
    )

    # Per-pathway rows
    row_grobs <- vector("list", n_path)
    for (i in seq_len(n_path)) {
        pname <- pathway_names[i]
        pgenes <- pathways[[pname]]
        fidx <- fgsea_idx[i]

        # ---- significance test (drives the name-only greying) ----
        # A pathway is "significant" if padj is present and strictly below the
        # threshold. Only the pathway NAME text is de-emphasised (drawn in
        # nonsig_color) when non-significant (padj >= threshold or NA) and
        # nonsig_color is not NULL. All other columns keep their normal colours.
        # >>> To grey the OTHER side (the significant pathways) instead, change
        #     `< padj_threshold` to `>= padj_threshold` on the next line. <<<
        padj_i <- fgseaRes$padj[fidx]
        is_sig <- !is.na(padj_i) && (padj_i < padj_threshold)
        name_col <- if (!is_sig && !is.null(nonsig_color)) nonsig_color else "black"

        row_grobs[[i]] <- list(
            grid::textGrob(pname, x = 1, hjust = 1,
                           gp = grid::gpar(fontsize = 10, col = name_col)),
            make_ticks_grob(pname, pgenes),
            grid::textGrob(fmt_nes(fgseaRes$NES[fidx]),
                           gp = grid::gpar(fontsize = 10)),
            grid::textGrob(fmt_padj(padj_i),
                           gp = grid::gpar(fontsize = 10)),
            grid::textGrob(highlight_listing(pname, pgenes),
                           x = 0, hjust = 0,
                           gp = grid::gpar(fontsize = 9, col = highlight_color))
        )
    }

    # Flatten into a single list in row-major order (header first, then rows)
    all_grobs <- c(header_grobs, unlist(row_grobs, recursive = FALSE))

    g <- gridExtra::arrangeGrob(
        grobs = all_grobs,
        ncol  = 5,
        nrow  = n_path + 1,
        widths = grid::unit(colwidths, "null"),
        heights = grid::unit(c(1.5, rep(1, n_path)), "null")
    )

    if (render) {
        grid::grid.newpage()
        grid::grid.draw(g)
    }
    invisible(g)
}