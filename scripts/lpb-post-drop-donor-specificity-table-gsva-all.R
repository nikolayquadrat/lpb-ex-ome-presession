#' Multi-region donor-specificity table for GSVA scores (SZ07 across regions).
#'
#' Shows the SAME pathways across several brain regions side by side, to test
#' whether a per-region GSVA signal replicates. Two stacked panels:
#'   Panel 1 ("replication")  -- a FIXED pathway set (e.g. the OUTRIDER/fgsea
#'                               significant pathways), grouped by theme, shown
#'                               across all regions.
#'   Panel 2 ("data-driven")  -- top-N up / top-N down pathways selected by
#'                               robust_z within chosen driver regions (blocks),
#'                               also shown across all regions.
#'
#' Layout per row:
#'   Pathway | neuronal spec | [region1: gap | <SZ07 | strip] | [region2: ...] | ...
#' Each region's strip is SELF-SCALED (GSVA scores are not comparable across
#' separately-run regions), and donor identity is BY COLUMN NAME (SZ07 found via
#' `sz07_pattern`, SZ vs HC via `sz_pattern`) -- no positional assumptions.
#'
#' Small-n regions (e.g. BA4) are shown but visually de-emphasised
#' (`deemph_regions`): thinner ticks + greyed header, so few-donor strips are
#' not over-read. Missing pathways in a region render as an empty strip.
#'
#' @param scores_list Named list of per-region score data.frames/matrices; each
#'   has pathway names in the first column (or row names) and one numeric column
#'   per donor, named "<DONOR>_<REGION>". Names of the list are the region labels
#'   and set column order.
#' @param panel1_pathways Character vector of the FIXED replication pathways
#'   (raw pathway ids). Shown in Panel 1.
#' @param panel1_groups Named vector pathway -> theme, for grouping Panel 1.
#' @param panel1_group_order Theme display order for Panel 1.
#' @param spec_all Named vector pathway -> whole-pathway neuronal specificity
#'   (mean delta_log2). Shown as the specificity column (spans both panels).
#' @param panel2_blocks Named list defining Panel 2 blocks; each element is
#'   list(region=<name in scores_list>, dir="up"|"down", n=<int>). The block name
#'   is its row-group label. Selection is by robust_z within that region.
#' @param panel1_title,panel2_title Panel header strings.
#' @param sz07_pattern,sz_pattern Regex to find the SZ07 column and to classify
#'   SZ (vs HC) donors within each region's columns.
#' @param deemph_regions Region names to de-emphasise (small n).
#' @param group_colors,sz07_color Strip colours.
#' @param label_fn,label_truncate Label display (see plot_gsva_specificity_table).
#' @param gap_digits,spec_digits Decimals.
#' @param colwidths Optional explicit column widths; if NULL, computed as
#'   c(pathway, spec, then per region: gap,frac,strip).
#' @param render draw to device.
#' @return invisibly list(grob, data) where data is the per-row, per-region long
#'   table (pathway, panel, block/group, region, sz07_score, gap, rank, n_more,
#'   n_donors, spec).
plot_gsva_multiregion_table <- function(
        scores_list,
        panel1_pathways,
        panel1_groups,
        panel1_group_order,
        spec_all,
        panel2_blocks,
        panel1_title    = "Replication of OUTRIDER pathways",
        panel2_title    = "Region-specific extremes (top by robust z)",
        sz07_pattern    = "^SZ07",
        sz_pattern      = "^SZ",
        deemph_regions  = character(0),
        group_colors    = c(SZ = "#e69f00", HC = "#0072b2"),
        sz07_color      = "#d62728",
        show_zero       = TRUE,
        label_fn        = NULL,
        label_truncate  = 46,
        label_size      = 7,
        gap_digits      = 2,
        spec_digits     = 2,
        spec_header     = "neuronal\nspecificity",
        colwidths       = NULL,
        header_lineheight = 0.9,
        header_height   = 1.6,
        panel_header_height = 0.8,
        group_header_height = 0.7,
        render          = TRUE) {

    for (p in c("ggplot2", "gridExtra", "grid"))
        if (!requireNamespace(p, quietly = TRUE)) stop(p, " is required")

    regions <- names(scores_list)
    if (is.null(regions) || !length(regions)) stop("scores_list must be named by region")

    # ---- normalise each region's matrix: pathway rownames + numeric donor mat -
    reg <- lapply(scores_list, function(d) {
        d <- as.data.frame(d, stringsAsFactors = FALSE, check.names = FALSE)
        if (!is.numeric(d[[1]])) { rn <- as.character(d[[1]]); d <- d[, -1, drop = FALSE]; rownames(d) <- rn }
        m <- as.matrix(d)
        storage.mode(m) <- "double"
        sz_col <- grep(sz07_pattern, colnames(m), value = TRUE)
        if (length(sz_col) != 1)
            stop("region needs exactly one SZ07 column matching '", sz07_pattern,
                 "'; found: ", paste(sz_col, collapse = ", "))
        others <- setdiff(colnames(m), sz_col)
        grp    <- ifelse(grepl(sz_pattern, others), "SZ", "HC")
        list(m = m, sz_col = sz_col, others = others, grp = grp)
    })
    names(reg) <- regions

    # ---- per-region per-pathway SZ07 gap/rank (direction-aware) --------------
    stat_one <- function(r, pw) {
        if (!(pw %in% rownames(r$m)))
            return(list(present = FALSE))
        sz <- r$m[pw, r$sz_col]
        M  <- r$m[pw, r$others]
        s  <- sign(sz); if (s == 0) s <- 1
        Ms <- M * s; szs <- sz * s
        list(present = TRUE, sz = sz,
             gap = szs - max(Ms, na.rm = TRUE),
             n_more = sum(Ms >= szs, na.rm = TRUE),
             n = length(M), M = M, grp = r$grp)
    }

    # ---- resolve Panel 2 selection (robust_z within the block's region) ------
    robz <- function(r, pw) {
        if (!(pw %in% rownames(r$m))) return(NA_real_)
        sz <- r$m[pw, r$sz_col]; M <- r$m[pw, r$others]
        s <- sign(sz); if (s == 0) s <- 1
        (sz * s - stats::median(M, na.rm = TRUE) * s) / stats::mad(M, na.rm = TRUE)
        # NB direction-aware robust z; ties fgsea 'z' convention
    }
    panel2_sel <- lapply(names(panel2_blocks), function(bn) {
        b <- panel2_blocks[[bn]]
        rr <- reg[[b$region]]
        if (is.null(rr)) stop("panel2 block '", bn, "' region '", b$region,
                              "' not in scores_list")
        pws <- rownames(rr$m)
        z   <- vapply(pws, function(pw) {
            sz <- rr$m[pw, rr$sz_col]; M <- rr$m[pw, rr$others]
            (sz - stats::median(M, na.rm = TRUE)) / stats::mad(M, na.rm = TRUE)
        }, numeric(1))
        ord <- if (b$dir == "up") order(-z) else order(z)
        data.frame(block = bn, pathway = pws[ord][seq_len(b$n)],
                   stringsAsFactors = FALSE)
    })
    panel2_df <- do.call(rbind, panel2_sel)

    # ---- assemble the row plan (panel, group/block, pathway) -----------------
    p1 <- data.frame(panel = "P1",
                     group = unname(panel1_groups[panel1_pathways]),
                     pathway = panel1_pathways, stringsAsFactors = FALSE)
    p1$group[is.na(p1$group)] <- "Other"
    glev <- c(intersect(panel1_group_order, unique(p1$group)),
              setdiff(unique(p1$group), panel1_group_order))
    p1$group <- factor(p1$group, levels = glev)
    p1 <- p1[order(p1$group), , drop = FALSE]

    p2 <- data.frame(panel = "P2", group = panel2_df$block,
                     pathway = panel2_df$pathway, stringsAsFactors = FALSE)
    p2$group <- factor(p2$group, levels = names(panel2_blocks))

    # ---- shared x-range PER REGION (self-scaled) over all shown pathways -----
    shown <- unique(c(p1$pathway, p2$pathway))
    xr <- lapply(reg, function(r) {
        v <- as.vector(r$m[intersect(shown, rownames(r$m)), , drop = FALSE])
        rng <- range(c(v, 0), na.rm = TRUE); rng + c(-1, 1) * 0.04 * diff(rng)
    })

    # ---- grobs -------------------------------------------------------------
    fmt <- function(x, d) if (is.na(x)) "" else sprintf(paste0("%+.", d, "f"), x)
    fmt_frac <- function(n, N) if (is.na(n)) "" else sprintf("%d/%d", as.integer(n), N)
    mk_label <- function(s) {
        if (!is.null(label_fn)) s <- label_fn(s)
        if (!is.null(label_truncate)) s <- truncate_words(s, width = label_truncate)
        s
    }
    strip_grob <- function(r, pw, region) {
        st <- stat_one(r, pw)
        thin <- region %in% deemph_regions
        p <- ggplot2::ggplot() +
            ggplot2::scale_x_continuous(limits = xr[[region]], expand = c(0, 0)) +
            ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.margin = grid::unit(c(1, 2, 1, 2), "pt"))
        if (show_zero)
            p <- p + ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3)
        if (isTRUE(st$present)) {
            d_o <- data.frame(v = as.numeric(st$M), grp = st$grp)
            p <- p +
                ggplot2::geom_segment(data = d_o,
                    ggplot2::aes(x = v, xend = v, y = 0.15, yend = 0.85, colour = grp),
                    linewidth = if (thin) 0.45 else 0.7, show.legend = FALSE) +
                ggplot2::geom_segment(
                    data = data.frame(v = st$sz),
                    ggplot2::aes(x = v, xend = v, y = 0.02, yend = 0.98),
                    colour = sz07_color, linewidth = if (thin) 1.0 else 1.5) +
                ggplot2::scale_colour_manual(values = group_colors)
        }
        ggplot2::ggplotGrob(p)
    }

    # column plan: [pathway][spec] then per region [gap][frac][strip]
    per_reg_cols <- 3L
    ncol_tab <- 2L + per_reg_cols * length(regions)
    if (is.null(colwidths)) {
        base <- c(8.5, 1.2)
        eachr <- c(0.9, 0.8, 3.2)
        colwidths <- c(base, rep(eachr, length(regions)))
    }
    stopifnot(length(colwidths) == ncol_tab)

    txt <- function(s, sz = 7.5, face = "plain", col = "black", x = 0.5, hjust = 0.5)
        grid::textGrob(s, x = x, hjust = hjust,
                       gp = grid::gpar(fontsize = sz, fontface = face, col = col,
                                       lineheight = header_lineheight))

    # header (2 rows: region banner + column subheaders)
    region_banner <- c(list(grid::nullGrob(), grid::nullGrob()),
        unlist(lapply(regions, function(rg) {
            col <- if (rg %in% deemph_regions) "grey55" else "black"
            list(grid::nullGrob(),
                 txt(paste0(rg, if (rg %in% deemph_regions) "\n(small n)" else ""),
                     8.5, "bold", col), grid::nullGrob())
        }), recursive = FALSE))
    subhead <- c(list(txt("Pathway", 9, "bold", x = 1, hjust = 1), txt(spec_header, 8, "bold")),
        unlist(lapply(regions, function(rg) list(
            txt("gap", 8, "bold"), txt("<SZ07", 7.5, "bold"),
            txt("donor GSVA (self-scaled)", 7.5, "bold"))), recursive = FALSE))

    data_row <- function(pw) {
        cells <- list(txt(mk_label(pw), label_size, x = 1, hjust = 1),
                      txt(fmt(unname(spec_all[pw]), spec_digits), 7.5))
        for (rg in regions) {
            r <- reg[[rg]]; st <- stat_one(r, pw)
            gap <- if (isTRUE(st$present)) st$gap else NA_real_
            nm  <- if (isTRUE(st$present)) st$n_more else NA_real_
            nn  <- if (isTRUE(st$present)) st$n else length(r$others)
            cells <- c(cells, list(
                txt(fmt(gap, gap_digits), 7.5),
                txt(fmt_frac(nm, nn), 7.5),
                strip_grob(r, pw, rg)))
        }
        cells
    }
    banner_row <- function(label, height_face = "bold.italic", col = "grey25") {
        c(list(txt(label, 8.5, height_face, col, x = 1, hjust = 1), grid::nullGrob()),
          unlist(lapply(regions, function(rg) list(grid::nullGrob(), grid::nullGrob(),
              grid::linesGrob(x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(.5, .5), "npc"),
                              gp = grid::gpar(col = "grey88", lwd = 0.8)))), recursive = FALSE))
    }
    panel_row <- function(label)
        c(list(txt(label, 9.5, "bold", "grey10", x = 0, hjust = 0), grid::nullGrob()),
          unlist(lapply(regions, function(rg) list(grid::nullGrob(), grid::nullGrob(), grid::nullGrob())),
                 recursive = FALSE))

    # ---- build rows --------------------------------------------------------
    rows <- list(); rh <- numeric(0)
    add <- function(g, h) { rows[[length(rows) + 1L]] <<- g; rh <<- c(rh, h) }

    add(panel_row(panel1_title), panel_header_height)
    for (gname in levels(droplevels(p1$group))) {
        idx <- which(as.character(p1$group) == gname); if (!length(idx)) next
        add(banner_row(gname), group_header_height)
        for (i in idx) add(data_row(p1$pathway[i]), 1)
    }
    add(panel_row(panel2_title), panel_header_height)
    for (bname in levels(droplevels(p2$group))) {
        idx <- which(as.character(p2$group) == bname); if (!length(idx)) next
        add(banner_row(bname), group_header_height)
        for (i in idx) add(data_row(p2$pathway[i]), 1)
    }

    # axis row (one self-scaled axis per region)
    axis_cells <- c(list(grid::nullGrob(), grid::nullGrob()))
    for (rg in regions) {
        ap <- ggplot2::ggplot(data.frame(x = xr[[rg]], y = c(0, 0))) +
            ggplot2::geom_blank(ggplot2::aes(x, y)) +
            ggplot2::scale_x_continuous(limits = xr[[rg]], expand = c(0, 0)) +
            ggplot2::theme_minimal(base_size = 7) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           axis.title = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.line.x = ggplot2::element_line(colour = "grey45"),
                           axis.ticks.x = ggplot2::element_line(colour = "grey45"),
                           plot.margin = grid::unit(c(6, 2, 0, 2), "pt"))
        axis_cells <- c(axis_cells, list(grid::nullGrob(), grid::nullGrob(),
                                         ggplot2::ggplotGrob(ap)))
    }
    add(axis_cells, 0.9)

    # legend row (spans strips)
    leg_df <- data.frame(x = 1:3, lab = c("SZ07", "SZ", "HC"),
                         col = c(sz07_color, unname(group_colors["SZ"]), unname(group_colors["HC"])))
    legp <- ggplot2::ggplot(leg_df, ggplot2::aes(x, 1)) +
        ggplot2::geom_segment(ggplot2::aes(x = x - .12, xend = x - .12, y = .6, yend = 1.4),
                              colour = leg_df$col, linewidth = 1.2) +
        ggplot2::geom_text(ggplot2::aes(x = x - .05, label = lab), hjust = 0, size = 2.6) +
        ggplot2::scale_x_continuous(limits = c(.5, 4.5)) +
        ggplot2::scale_y_continuous(limits = c(.4, 1.6)) + ggplot2::theme_void()
    leg_cells <- c(list(grid::nullGrob(), grid::nullGrob()),
                   list(grid::nullGrob(), grid::nullGrob(), ggplot2::ggplotGrob(legp)),
                   unlist(lapply(regions[-1], function(rg) list(grid::nullGrob(), grid::nullGrob(), grid::nullGrob())),
                          recursive = FALSE))
    add(leg_cells, 0.7)

    all_grobs <- c(region_banner, subhead, unlist(rows, recursive = FALSE))
    g <- gridExtra::arrangeGrob(
        grobs = all_grobs, ncol = ncol_tab, nrow = length(rh) + 2L,
        widths = grid::unit(colwidths, "null"),
        heights = grid::unit(c(1.0, header_height, rh), "null"))

    if (render) { grid::grid.newpage(); grid::grid.draw(g) }

    # tidy data out (long: pathway x region)
    out <- do.call(rbind, lapply(c(split(p1, seq_len(nrow(p1))), split(p2, seq_len(nrow(p2)))),
        function(rw) do.call(rbind, lapply(regions, function(rg) {
            st <- stat_one(reg[[rg]], rw$pathway)
            data.frame(panel = rw$panel, group = as.character(rw$group), pathway = rw$pathway,
                       region = rg,
                       sz07_score = if (isTRUE(st$present)) st$sz else NA_real_,
                       gap = if (isTRUE(st$present)) st$gap else NA_real_,
                       n_more = if (isTRUE(st$present)) st$n_more else NA_integer_,
                       n_donors = if (isTRUE(st$present)) st$n else NA_integer_,
                       spec = unname(spec_all[rw$pathway]), stringsAsFactors = FALSE)
        }))))
    rownames(out) <- NULL
    invisible(list(grob = g, data = out))
}
