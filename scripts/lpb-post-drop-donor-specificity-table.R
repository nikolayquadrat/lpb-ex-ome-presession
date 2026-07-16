#' Drop the MSigDB collection prefix and underscores from pathway names
#'
#' Display helper for `plot_donor_specificity_table(label_fn = )`. The
#' collection prefix (REACTOME_, GOBP_, ...) costs 3-9 characters per label and
#' carries no information the figure needs -- state the collection in the
#' caption instead. Converting underscores to spaces also lets `label_wrap`
#' break lines (strwrap needs whitespace).
#'
#' @param x character vector of pathway names.
#' @return character vector of display labels.
strip_msigdb_prefix <- function(x) {
    x <- sub("^(REACTOME|GOBP|GOCC|GOMF|GOBP|KEGG_MEDICUS|KEGG|HALLMARK|WP|HP|BIOCARTA|PID|MODULE)_",
             "", x)
    gsub("_", " ", x)
}

#' Table-style plot of donor-specificity for SZ07-significant GSEA pathways
#'
#' One row per pathway, laid out like fgsea's plotGseaTable:
#'   Pathway | padj | NES | donors p<SZ07 | strip of every donor's NES
#' The strip places each donor's NES as a vertical tick on a scale SHARED by all
#' rows, so rows are directly comparable. SZ07 is highlighted; the other donors
#' are coloured by diagnosis group (e.g. patients vs. no psychiatric diagnosis).
#'
#' ROW ORDER. If `pathway_groups` is supplied, rows are blocked by functional
#' group first and ordered by specificity within each group; groups themselves
#' are ordered by `group_order` (manual) or by their best member (automatic).
#' Group blocks carry a header row. Without `pathway_groups`, rows are ordered
#' by specificity alone.
#'
#' SPECIFICITY. The default statistic is the RUNNER-UP GAP: how far SZ07's NES
#' sits beyond the most extreme other donor, in SZ07's own direction (so a
#' down-regulated pathway is judged on how much MORE negative SZ07 is). This is
#' the magnitude a plain rank throws away: rank 1 with a large gap is a
#' qualitatively distinct signal, rank 1 with a hair's-breadth gap is not.
#'
#' IMPORTANT -- donor identity is POSITIONAL. `other_donors_NES` is a
#' comma-separated string carrying no donor names, so `donor_ids` /
#' `donor_groups` must be given in EXACTLY the column order used to build that
#' string. A wrong order silently mislabels every donor and every group colour,
#' with no error. Verify the order against the code that produced the workbook.
#'
#' @param res         data.frame of the third-pass fgsea results (one row per
#'                    SZ07-significant pathway).
#' @param donor_ids   Character vector naming the other donors, in the column
#'                    order of `others_nes_col`. Length must match.
#' @param donor_groups Character vector (same length/order as `donor_ids`)
#'                    giving each donor's group, e.g. "SZ" / "HC". Values must
#'                    be names of `group_colors`.
#' @param pathway_groups Optional manual functional grouping: a NAMED character
#'                    vector mapping pathway name -> group label, e.g.
#'                      c(REACTOME_NEUREXINS_AND_NEUROLIGINS = "synaptic",
#'                        GOBP_MATURATION_OF_LSU_RRNA        = "RNA/ribosome")
#'                    Pathways absent from the vector fall into
#'                    `ungrouped_label`. NULL = no grouping (rank order only).
#' @param group_order Optional character vector fixing the top-to-bottom order
#'                    of groups. NULL = groups ordered by their best (most
#'                    specific) member. Groups listed here but absent from the
#'                    selected rows are skipped; groups present but unlisted are
#'                    appended in automatic order.
#' @param top_n       Total rows to show when `top_n_per_group` is NULL
#'                    (default 12).
#' @param top_n_per_group If set, take pathways from EACH group instead of
#'                    `top_n` overall. Guarantees every group is represented --
#'                    useful when a group's message is that it is NOT specific
#'                    (e.g. showing the immune block in context). Accepts:
#'                      scalar        -- same N for every group, e.g. 5
#'                      NAMED vector  -- per-group N, matched by label and so
#'                                       independent of ordering (recommended):
#'                                       c(synaptic = 5, immune = 3)
#'                      unnamed vector-- per-group N aligned POSITIONALLY to
#'                                       `group_order`, which must then be given
#'                                       and of equal length: c(5, 5, 4, 3)
#'                    Groups left uncovered are DROPPED, so per-group counts
#'                    also choose which groups appear (this is how to exclude an
#'                    auto-appended "other" block). A count of 0 drops a group.
#' @param rank_by     Specificity statistic: "gap" (runner-up gap, default),
#'                    "z" (SZ07's NES in SDs of the other donors' NES), or
#'                    "padj" (SZ07's own significance -- NOT a specificity
#'                    measure).
#' @param ungrouped_label Group label for pathways missing from
#'                    `pathway_groups` (default "other").
#' @param sz07_label  Label for the index donor (default "SZ07").
#' @param sz07_color,group_colors Colours. `group_colors` is a NAMED vector
#'                    keyed by the values in `donor_groups`.
#' @param show_zero   Draw a reference line at NES = 0 (default TRUE).
#' @param warn_duplicate_donors Check for donor columns identical across all
#'                    pathways (a duplicated donor inflates n and narrows the
#'                    null spread). Default TRUE.
#' @param pathway_col,padj_col,nes_col,sz07_nes_col,others_nes_col,n_smaller_col
#'                    Column names in `res`.
#' @param label_fn    Optional function applied to pathway names for DISPLAY
#'                    only -- join keys and the `pathway_groups` lookup still use
#'                    the raw names, so this cannot desynchronise them. Use
#'                    `strip_msigdb_prefix` to drop the REACTOME_/GOBP_/... 
#'                    prefix and underscores. NULL = show names verbatim.
#' @param label_wrap  Wrap display labels onto multiple lines at this width in
#'                    characters (NULL = no wrapping). Applied AFTER `label_fn`,
#'                    and it needs whitespace to break on -- so pair it with
#'                    `label_fn = strip_msigdb_prefix`, or underscore-joined
#'                    names will not wrap. When set, `shorten_names` is ignored
#'                    (wrapping keeps the whole label instead of cutting it).
#' @param label_size  Font size of the pathway labels (default 7.5). Lower it
#'                    to fit more text rather than widening the column.
#' @param shorten_names Truncate pathway labels to this many characters
#'                    (NULL = no truncation). This -- not `colwidths` -- is what
#'                    produces the trailing ellipsis; widening the column alone
#'                    will not reveal more text. Ignored when `label_wrap` is set.
#' @param colwidths   Relative widths: c(pathway, padj, NES, donors, strip).
#'                    Raise the first element to give names more room; the strip
#'                    is the natural donor (its message survives compression).
#' @param group_header_height Height of a group header row, relative to a data
#'                    row (default 0.85).
#' @param render      Draw on the current device (default TRUE).
#'
#' @return Invisibly, a list with $grob (the gtable) and $data (the plotted
#'   rows incl. functional group, gap, z, rank, n_more_extreme and the
#'   donors-p<SZ07 fraction), so the numbers behind the figure can be reported.
plot_donor_specificity_table <- function(
        res,
        donor_ids,
        donor_groups,
        pathway_groups  = NULL,
        group_order     = NULL,
        top_n           = 12,
        top_n_per_group = NULL,
        rank_by         = c("gap", "z", "padj"),
        ungrouped_label = "other",
        sz07_label      = "SZ07",
        sz07_color      = "#d62728",
        group_colors    = c(SZ = "#e69f00", HC = "#0072b2"),
        show_zero       = TRUE,
        warn_duplicate_donors = TRUE,
        pathway_col     = "pathway",
        padj_col        = "padj",
        nes_col         = "NES",
        sz07_nes_col    = "SZ07_simple_NES",
        others_nes_col  = "other_donors_NES",
        n_smaller_col   = "n_donors_psmallerSZ07",
        label_fn        = NULL,
        label_wrap      = NULL,
        label_size      = 7.5,
        shorten_names   = 44,
        colwidths       = c(8, 1.2, 1.0, 1.3, 5.5),
        group_header_height = 0.85,
        render          = TRUE) {

    for (p in c("ggplot2", "gridExtra", "grid")) {
        if (!requireNamespace(p, quietly = TRUE)) stop(p, " is required")
    }
    rank_by <- match.arg(rank_by)
    res <- as.data.frame(res, stringsAsFactors = FALSE)

    need <- c(pathway_col, padj_col, nes_col, sz07_nes_col, others_nes_col)
    miss <- need[!need %in% names(res)]
    if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))

    # ---- parse the positional donor-NES strings ---------------------------
    parse_row <- function(s) {
        if (is.na(s) || !nzchar(s)) return(numeric(0))
        as.numeric(strsplit(gsub("\\s", "", s), ",")[[1]])
    }
    other_list <- lapply(res[[others_nes_col]], parse_row)
    n_other <- unique(lengths(other_list))
    if (length(n_other) != 1L)
        stop("rows of ", others_nes_col, " have differing donor counts: ",
             paste(sort(unique(lengths(other_list))), collapse = ", "))
    if (length(donor_ids) != n_other || length(donor_groups) != n_other)
        stop(sprintf(paste0("donor_ids/donor_groups length (%d/%d) must equal the ",
                            "number of donors in %s (%d). Donor identity is POSITIONAL."),
                     length(donor_ids), length(donor_groups), others_nes_col, n_other))
    bad_grp <- setdiff(unique(donor_groups), names(group_colors))
    if (length(bad_grp))
        stop("donor_groups value(s) without a colour in group_colors: ",
             paste(bad_grp, collapse = ", "))

    M <- do.call(rbind, other_list)          # pathways x donors
    colnames(M) <- donor_ids

    # ---- duplicated-donor check -------------------------------------------
    if (warn_duplicate_donors && ncol(M) > 1) {
        dup <- character(0)
        for (i in seq_len(ncol(M) - 1)) for (j in (i + 1):ncol(M)) {
            if (isTRUE(all.equal(M[, i], M[, j], tolerance = 1e-10)))
                dup <- c(dup, sprintf("%s==%s", donor_ids[i], donor_ids[j]))
        }
        if (length(dup))
            warning("donor column(s) identical across ALL pathways -- a duplicated ",
                    "donor inflates n and narrows the null spread: ",
                    paste(dup, collapse = ", "))
    }

    # ---- specificity statistics (direction-aware) -------------------------
    sz  <- as.numeric(res[[sz07_nes_col]])
    sgn <- sign(sz); sgn[sgn == 0] <- 1
    Ms  <- M * sgn                      # orient so "more extreme" is larger
    szs <- sz * sgn
    res$.gap    <- szs - apply(Ms, 1, max, na.rm = TRUE)
    res$.z      <- (szs - rowMeans(Ms, na.rm = TRUE)) /
                    apply(Ms, 1, stats::sd, na.rm = TRUE)
    res$.n_more <- rowSums(Ms >= szs, na.rm = TRUE)
    res$.rank   <- res$.n_more + 1L
    res$.stat   <- switch(rank_by,
                          gap  =  res$.gap,
                          z    =  res$.z,
                          padj = -res[[padj_col]])   # negate: larger = better

    # ---- functional grouping ----------------------------------------------
    if (!is.null(pathway_groups)) {
        if (is.null(names(pathway_groups)))
            stop("pathway_groups must be a NAMED vector: pathway -> group label")
        res$.fgroup <- unname(pathway_groups[res[[pathway_col]]])
        res$.fgroup[is.na(res$.fgroup)] <- ungrouped_label
    } else {
        res$.fgroup <- NA_character_
    }

    # ---- resolve how many rows to take from each group ---------------------
    # `top_n_per_group` may be:
    #   scalar          -> same N for every group
    #   NAMED vector    -> matched by group label (order-independent, safest)
    #   unnamed vector  -> aligned POSITIONALLY to `group_order`, which must
    #                      then be supplied and of equal length
    # Groups left uncovered are DROPPED: specifying per-group counts also
    # selects which groups appear (this is how you exclude an auto-appended
    # "other" block). A count of 0 drops that group too.
    resolve_per_group <- function(tpg, groups_in_data) {
        if (length(tpg) == 1L && is.null(names(tpg)))
            return(setNames(rep(as.integer(tpg), length(groups_in_data)),
                            groups_in_data))
        if (!is.null(names(tpg)) && all(nzchar(names(tpg)))) {
            unknown <- setdiff(names(tpg), groups_in_data)
            if (length(unknown))
                warning("top_n_per_group name(s) absent from pathway_groups ",
                        "(ignored): ", paste(unknown, collapse = ", "))
            v <- tpg[names(tpg) %in% groups_in_data]
            if (!length(v)) stop("no top_n_per_group name matches any group")
            return(setNames(as.integer(v), names(v)))
        }
        if (is.null(group_order))
            stop("unnamed top_n_per_group of length ", length(tpg),
                 " needs group_order to say which count belongs to which group. ",
                 "Either set group_order, or pass a NAMED vector, e.g. ",
                 "c(synaptic = 5, immune = 3).")
        if (length(tpg) != length(group_order))
            stop(sprintf(paste0("top_n_per_group has %d value(s) but group_order has ",
                                "%d group(s); lengths must match for positional alignment."),
                         length(tpg), length(group_order)))
        setNames(as.integer(tpg), group_order)
    }

    # ---- select rows -------------------------------------------------------
    if (!is.null(top_n_per_group) && !is.null(pathway_groups)) {
        groups_in_data <- unique(res$.fgroup)
        npg <- resolve_per_group(top_n_per_group, groups_in_data)
        npg <- npg[npg > 0]
        dropped <- setdiff(groups_in_data, names(npg))
        if (length(dropped))
            message("[donor-specificity] group(s) not covered by top_n_per_group, ",
                    "omitted from the figure: ", paste(dropped, collapse = ", "))
        sel <- do.call(rbind, lapply(names(npg), function(gn) {
            d <- res[res$.fgroup == gn, , drop = FALSE]
            if (!nrow(d)) return(NULL)
            d <- d[order(-d$.stat), , drop = FALSE]
            head(d, min(npg[[gn]], nrow(d)))
        }))
        if (is.null(sel) || !nrow(sel))
            stop("no rows selected: check top_n_per_group against the group labels")
    } else {
        if (!is.null(top_n_per_group))
            warning("top_n_per_group ignored: needs pathway_groups; using top_n.")
        sel <- res[order(-res$.stat), , drop = FALSE]
        sel <- head(sel, min(top_n, nrow(sel)))
    }
    if (!nrow(sel)) stop("no rows selected")

    # ---- order groups, then rows within group ------------------------------
    if (!is.null(pathway_groups)) {
        auto <- names(sort(tapply(sel$.stat, sel$.fgroup, max), decreasing = TRUE))
        if (is.null(group_order)) {
            glev <- auto
        } else {
            glev <- c(intersect(group_order, auto), setdiff(auto, group_order))
        }
        sel$.fgroup <- factor(sel$.fgroup, levels = glev)
        sel <- sel[order(sel$.fgroup, -sel$.stat), , drop = FALSE]
    } else {
        sel <- sel[order(-sel$.stat), , drop = FALSE]
    }

    # rows of M matching sel, by pathway name
    Msel  <- M[match(sel[[pathway_col]], res[[pathway_col]]), , drop = FALSE]
    szsel <- sz[match(sel[[pathway_col]], res[[pathway_col]])]

    # ---- shared NES scale across every row --------------------------------
    xr <- range(c(as.vector(Msel), szsel, 0), na.rm = TRUE)
    xr <- xr + c(-1, 1) * 0.04 * diff(xr)

    # ---- per-row strip: one vertical tick per donor -----------------------
    make_strip <- function(i) {
        d_o <- data.frame(nes = as.numeric(Msel[i, ]),
                          grp = as.character(donor_groups),
                          stringsAsFactors = FALSE)
        d_s <- data.frame(nes = szsel[i])
        p <- ggplot2::ggplot() +
            ggplot2::scale_x_continuous(limits = xr, expand = c(0, 0)) +
            ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.margin = grid::unit(c(1, 2, 1, 2), "pt"))
        if (show_zero)
            p <- p + ggplot2::geom_vline(xintercept = 0, colour = "grey80",
                                         linewidth = 0.3)
        p <- p + ggplot2::geom_segment(
                data = d_o,
                ggplot2::aes(x = nes, xend = nes, y = 0.12, yend = 0.88, colour = grp),
                linewidth = 0.7, show.legend = FALSE) +
            ggplot2::geom_segment(
                data = d_s,
                ggplot2::aes(x = nes, xend = nes, y = 0.02, yend = 0.98),
                colour = sz07_color, linewidth = 1.5) +
            ggplot2::scale_colour_manual(values = group_colors)
        ggplot2::ggplotGrob(p)
    }

    # ---- formatters --------------------------------------------------------
    fmt_p <- function(p) if (is.na(p)) "NA" else if (p == 0) "<1e-300" else sprintf("%.1e", p)
    fmt_n <- function(x) if (is.na(x)) "NA" else sprintf("%.2f", x)
    # n / N -- the fraction conveys "how many of the cohort match SZ07"
    fmt_frac <- function(n, N) if (is.na(n)) "NA" else sprintf("%d/%d", as.integer(n), N)
    # display-label pipeline: label_fn -> wrap (or shorten). Applied to the
    # rendered text only; `pathway_col` stays the join key throughout.
    make_label <- function(s) {
        if (!is.null(label_fn)) {
            if (!is.function(label_fn)) stop("label_fn must be a function")
            s <- label_fn(s)
        }
        if (!is.null(label_wrap)) {
            # keep the whole label, break it over lines
            return(vapply(s, function(x)
                paste(strwrap(x, width = label_wrap), collapse = "\n"),
                character(1), USE.NAMES = FALSE))
        }
        if (!is.null(shorten_names))
            s <- ifelse(nchar(s) > shorten_names,
                        paste0(substr(s, 1, shorten_names - 1), "\u2026"), s)
        s
    }

    has_nsm <- n_smaller_col %in% names(sel)

    header <- list(
        grid::textGrob("Pathway", x = 1, hjust = 1, gp = grid::gpar(fontface = "bold", fontsize = 9)),
        grid::textGrob("padj",    gp = grid::gpar(fontface = "bold", fontsize = 9)),
        grid::textGrob("NES",     gp = grid::gpar(fontface = "bold", fontsize = 9)),
        grid::textGrob("donors\np<SZ07", gp = grid::gpar(fontface = "bold", fontsize = 8)),
        grid::textGrob("donor NES (shared scale)", gp = grid::gpar(fontface = "bold", fontsize = 9))
    )

    make_data_row <- function(i) {
        nsm <- if (has_nsm) sel[[n_smaller_col]][i] else sel$.n_more[i]
        list(
            grid::textGrob(make_label(sel[[pathway_col]][i]), x = 1, hjust = 1,
                           gp = grid::gpar(fontsize = label_size)),
            grid::textGrob(fmt_p(sel[[padj_col]][i]), gp = grid::gpar(fontsize = 7.5)),
            grid::textGrob(fmt_n(sel[[nes_col]][i]),  gp = grid::gpar(fontsize = 7.5)),
            grid::textGrob(fmt_frac(nsm, n_other),    gp = grid::gpar(fontsize = 7.5)),
            make_strip(i)
        )
    }
    make_group_row <- function(label) list(
        grid::textGrob(label, x = 1, hjust = 1,
                       gp = grid::gpar(fontface = "bold.italic", fontsize = 8.5,
                                       col = "grey25")),
        grid::nullGrob(), grid::nullGrob(), grid::nullGrob(),
        grid::linesGrob(x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(0.5, 0.5), "npc"),
                        gp = grid::gpar(col = "grey85", lwd = 0.8))
    )

    # ---- assemble rows (with group headers) --------------------------------
    row_grobs <- list(); row_h <- numeric(0)
    if (is.null(pathway_groups)) {
        for (i in seq_len(nrow(sel))) {
            row_grobs[[length(row_grobs) + 1L]] <- make_data_row(i)
            row_h <- c(row_h, 1)
        }
    } else {
        for (gname in levels(droplevels(sel$.fgroup))) {
            idx <- which(as.character(sel$.fgroup) == gname)
            if (!length(idx)) next
            row_grobs[[length(row_grobs) + 1L]] <- make_group_row(gname)
            row_h <- c(row_h, group_header_height)
            for (i in idx) {
                row_grobs[[length(row_grobs) + 1L]] <- make_data_row(i)
                row_h <- c(row_h, 1)
            }
        }
    }

    # ---- x-axis row under the strips --------------------------------------
    axis_p <- ggplot2::ggplot(data.frame(x = xr, y = c(0, 0))) +
        ggplot2::geom_blank(ggplot2::aes(x, y)) +
        ggplot2::scale_x_continuous(limits = xr, expand = c(0, 0), name = "NES") +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(
            panel.grid   = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y  = ggplot2::element_blank(),
            axis.line.x  = ggplot2::element_line(colour = "grey40"),
            axis.ticks.x = ggplot2::element_line(colour = "grey40"),
            axis.title.x = ggplot2::element_text(size = 8),
            plot.margin  = grid::unit(c(8, 2, 0, 2), "pt"))
    axis_row <- list(grid::nullGrob(), grid::nullGrob(), grid::nullGrob(),
                     grid::nullGrob(), ggplot2::ggplotGrob(axis_p))

    # ---- legend row --------------------------------------------------------
    leg_lab <- c(sz07_label, names(group_colors))
    leg_col <- c(sz07_color, unname(group_colors))
    leg_df  <- data.frame(x = seq_along(leg_lab), y = 1, lab = leg_lab,
                          stringsAsFactors = FALSE)
    leg_p <- ggplot2::ggplot(leg_df, ggplot2::aes(x, y)) +
        ggplot2::geom_segment(ggplot2::aes(x = x - 0.12, xend = x - 0.12,
                                           y = 0.6, yend = 1.4),
                              colour = leg_col, linewidth = 1.2) +
        ggplot2::geom_text(ggplot2::aes(x = x - 0.05, label = lab), hjust = 0, size = 2.7) +
        ggplot2::scale_x_continuous(limits = c(0.5, length(leg_lab) + 1.2)) +
        ggplot2::scale_y_continuous(limits = c(0.4, 1.6)) +
        ggplot2::theme_void()
    legend_row <- list(grid::nullGrob(), grid::nullGrob(), grid::nullGrob(),
                       grid::nullGrob(), ggplot2::ggplotGrob(leg_p))

    all_grobs <- c(header, unlist(row_grobs, recursive = FALSE), axis_row, legend_row)
    g <- gridExtra::arrangeGrob(
        grobs   = all_grobs,
        ncol    = 5,
        nrow    = length(row_h) + 3L,
        widths  = grid::unit(colwidths, "null"),
        heights = grid::unit(c(1.1, row_h, 1.1, 0.7), "null"))

    if (render) { grid::grid.newpage(); grid::grid.draw(g) }

    out <- sel[, c(pathway_col, padj_col, nes_col, sz07_nes_col), drop = FALSE]
    if (!is.null(pathway_groups)) out$functional_group <- as.character(sel$.fgroup)
    out$gap_runnerup   <- sel$.gap
    out$z_vs_donors    <- sel$.z
    out$sz07_rank      <- sel$.rank
    out$n_more_extreme <- sel$.n_more
    out$n_donors       <- n_other
    if (has_nsm)
        out$donors_p_lt_SZ07 <- sprintf("%d/%d", as.integer(sel[[n_smaller_col]]), n_other)
    rownames(out) <- NULL
    invisible(list(grob = g, data = out))
}