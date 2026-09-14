#' Drop the MSigDB collection prefix and underscores from pathway names.
#'
#' Display helper for `plot_gsva_specificity_table(label_fn = )`. Removes the
#' collection prefix (REACTOME_, GOBP_, ...) and turns underscores into spaces
#' so labels read as text (and so `truncate_words`/`strwrap` can split on
#' whitespace).
#'
#' @param x character vector of pathway names.
#' @return character vector of display labels.
strip_msigdb_prefix <- function(x) {
    x <- sub("^(REACTOME|GOBP|GOCC|GOMF|KEGG_MEDICUS|KEGG|HALLMARK|WP|HP|BIOCARTA|PID|MODULE)_",
             "", x)
    gsub("_", " ", x)
}

#' Compress a label to <= `width` chars by truncating its LONGEST words first.
#'
#' Each over-long word is shortened to its initial letters + "." , starting with
#' the longest, until the whole label fits `width`. Unlike an end-cut, this keeps
#' the END of the name visible -- useful when the tail carries the meaning
#' (e.g. "... IN RESPONSE TO STRESS"). Labels already within `width` pass through
#' unchanged. Words are not shortened below `min_keep` letters (+ the dot), so
#' they stay recognisable ("REGULATION" -> "REG.", never "R.").
#'
#' Apply AFTER `strip_msigdb_prefix` (it needs spaces between words). Intended as
#' the `label_fn`, or via the `label_truncate` argument which calls it for you.
#'
#' @param x character vector of (already prefix-stripped) labels.
#' @param width target maximum characters (default 80).
#' @param min_keep minimum letters kept from a shortened word, before the dot
#'   (default 3).
#' @return character vector of compressed labels.
truncate_words <- function(x, width = 80, min_keep = 3) {
    vapply(x, function(s) {
        if (is.na(s) || nchar(s) <= width) return(s)
        w <- strsplit(s, " ", fixed = TRUE)[[1]]
        # Fully compress the longest word to its stem + "." , then the next
        # longest, until the label fits. Compressing one word all the way (rather
        # than nibbling every word by one char) fits the target with the FEWEST
        # words touched, so most of the name -- including its tail -- stays intact.
        repeat {
            if (nchar(paste(w, collapse = " ")) <= width) break
            len  <- nchar(w)
            elig <- which(len > (min_keep + 1) & !grepl("\\.$", w))  # not-yet-stemmed long words
            if (!length(elig)) break                                  # nothing left to shrink
            j <- elig[which.max(len[elig])]                           # longest remaining word
            w[j] <- paste0(substr(w[j], 1, min_keep), ".")            # stem it fully: KEEP+"."
        }
        paste(w, collapse = " ")
    }, character(1), USE.NAMES = FALSE)
}


#' GSVA analogue of `plot_donor_specificity_table()`. Same layout --
#'   Pathway | score | [neuronal specificity] | donors <SZ07 | strip of donor scores
#' -- but driven by a GSVA pathways x donors SCORE MATRIX rather than fgsea
#' output. The differences from the fgsea version are only in the INPUTS, because
#' GSVA has no NES, no per-donor enrichment p-value, and no leading edge:
#'   * "NES"            -> SZ07's GSVA score          (the enrichment column)
#'   * "padj"           -> DROPPED                    (GSVA has no pathway p-value)
#'   * donor-NES strip  -> donor GSVA-score strip     (identical mechanic; GSVA
#'                                                     scores are directly
#'                                                     comparable on a shared
#'                                                     scale, so it is if anything
#'                                                     cleaner than for NES)
#'   * leading-edge spec-> DROPPED                    (no leading edge); only the
#'                                                     WHOLE-PATHWAY neuronal
#'                                                     specificity (spec_all) is
#'                                                     shown -- which is exactly
#'                                                     the column that reveals
#'                                                     whether a GSVA signal is
#'                                                     composition-driven.
#'
#' Specificity/gap/z/rank are computed the SAME direction-aware way as the fgsea
#' version: SZ07's score is oriented by its own sign, so a pathway where SZ07 is
#' strongly NEGATIVE is judged on how much MORE negative it is than other donors.
#'
#' DONOR IDENTITY IS BY MATRIX COLUMN NAME here, not positional -- `scores`
#' columns carry donor names, so SZ07 is found by `sz07_id` and the other donors
#' by `setdiff`. `donor_groups` must be NAMED by donor id (or, if unnamed, given
#' in the column order of the non-SZ07 columns). This removes the positional
#' foot-gun of the fgsea version.
#'
#' @param scores      Numeric matrix, pathways (rows) x donors (columns), the
#'                    GSVA output. Row names = pathway names (raw, matching
#'                    `pathway_groups` / `spec_all` keys). Column names = donor
#'                    ids, including `sz07_id`.
#' @param sz07_id     Column name of the index donor (default "SZ07").
#' @param donor_groups Group per NON-SZ07 donor. Either a NAMED vector keyed by
#'                    donor id (recommended, order-independent), or an unnamed
#'                    vector aligned to the non-SZ07 columns of `scores` in order.
#'                    Values must be names of `group_colors`.
#' @param spec_all    Per-pathway WHOLE-PATHWAY neuronal specificity (mean
#'                    delta_log2 over the gene set). Named by pathway (or unnamed,
#'                    aligned to `rownames(scores)`). This is the key column for
#'                    GSVA: negative = non-neuronal-leaning (composition-driven),
#'                    positive = neuronal, ~0 = composition-independent. Column is
#'                    shown only if supplied.
#' @param rank_by     Specificity statistic to rank rows by: "gap" (runner-up
#'                    gap in SZ07's direction, default), "z" (SZ07 score in SDs of
#'                    the other donors), or "score" (SZ07's raw signed score).
#' @param score_header,strip_header Column/strip titles (default "GSVA" and
#'                    "donor GSVA score (shared scale)").
#' @param score_digits Decimals for the score and gap columns (default 2).
#' @param resid_z,resid_z_header,resid_z_digits,resid_z_width Optional SZ07
#'                    composition-residual z column (from
#'                    pathway_composition_correlation()), same meaning as in the
#'                    fgsea version: large |z| = NOT explained by composition.
#'                    Named by pathway (or aligned to rownames(scores)).
#' All other arguments (pathway_groups, group_order, top_n, top_n_per_group,
#' ungrouped_label, colours, label_fn, label_wrap, label_truncate, label_size,
#' shorten_names, colwidths, heights, render) behave as in
#' `plot_donor_specificity_table()`. `label_truncate` is GSVA-version-specific:
#' set it to a character width (e.g. 80) to compress over-long labels by
#' shortening their longest words to "WORD." (via `truncate_words`), keeping the
#' name's END visible; it takes precedence over `label_wrap` and `shorten_names`.
#'
#' @return Invisibly, list($grob, $data) as in the fgsea version, with SZ07's
#'   score, gap, z, rank, donors-below-SZ07 fraction and specificity.
plot_gsva_specificity_table <- function(
        scores,
        sz07_id         = "SZ07",
        donor_groups,
        pathway_groups  = NULL,
        group_order     = NULL,
        top_n           = 12,
        top_n_per_group = NULL,
        rank_by         = c("gap", "z", "score"),
        ungrouped_label = "other",
        sz07_label      = "SZ07",
        sz07_color      = "#d62728",
        group_colors    = c(SZ = "#e69f00", HC = "#0072b2"),
        show_zero       = TRUE,
        # score (enrichment) column
        score_header    = "GSVA",
        strip_header    = "donor GSVA score (shared scale)",
        score_digits    = 2,
        # whole-pathway neuronal specificity column
        spec_all        = NULL,
        spec_header     = "neuronal\nspecificity",
        spec_digits     = 2,
        spec_width      = 1.4,
        # optional SZ07 composition-residual column
        resid_z         = NULL,
        resid_z_header  = "SZ07\nresid-z",
        resid_z_digits  = 1,
        resid_z_width   = 1.1,
        header_lineheight = 0.9,
        header_height     = 1.5,
        label_fn        = NULL,
        label_wrap      = NULL,
        label_truncate  = NULL,
        label_size      = 7.5,
        shorten_names   = 44,
        colwidths       = NULL,     # base 4: pathway, score, donors, strip
        group_header_height = 0.85,
        render          = TRUE) {

    for (p in c("ggplot2", "gridExtra", "grid")) {
        if (!requireNamespace(p, quietly = TRUE)) stop(p, " is required")
    }
    rank_by <- match.arg(rank_by)

    # ---- inputs: matrix -> SZ07 vector + other-donor matrix ---------------
    scores <- as.matrix(scores)
    if (is.null(rownames(scores))) stop("`scores` needs pathway row names")
    if (is.null(colnames(scores))) stop("`scores` needs donor column names")
    if (!sz07_id %in% colnames(scores))
        stop("sz07_id '", sz07_id, "' is not a column of `scores`")

    other_cols <- setdiff(colnames(scores), sz07_id)
    if (!length(other_cols)) stop("no non-SZ07 donor columns in `scores`")

    # donor_groups: named by donor id, or unnamed aligned to other_cols
    if (!is.null(names(donor_groups)) && any(nzchar(names(donor_groups)))) {
        miss <- setdiff(other_cols, names(donor_groups))
        if (length(miss))
            stop("donor_groups missing entries for donor(s): ",
                 paste(miss, collapse = ", "))
        grp_vec <- unname(donor_groups[other_cols])         # align to other_cols
    } else {
        if (length(donor_groups) != length(other_cols))
            stop(sprintf(paste0("unnamed donor_groups length %d must equal the number of ",
                 "non-SZ07 donors (%d); or name it by donor id."),
                 length(donor_groups), length(other_cols)))
        grp_vec <- as.character(donor_groups)
    }
    bad_grp <- setdiff(unique(grp_vec), names(group_colors))
    if (length(bad_grp))
        stop("donor_groups value(s) without a colour in group_colors: ",
             paste(bad_grp, collapse = ", "))

    pw   <- rownames(scores)
    sz   <- as.numeric(scores[, sz07_id])                    # SZ07 score per pathway
    M    <- scores[, other_cols, drop = FALSE]               # pathways x other donors
    n_other <- length(other_cols)

    res <- data.frame(pathway = pw, sz07_score = sz,
                      stringsAsFactors = FALSE)

    # ---- specificity statistics (direction-aware, same as fgsea version) --
    sgn <- sign(sz); sgn[sgn == 0] <- 1
    Ms  <- M * sgn
    szs <- sz * sgn
    res$.gap    <- szs - apply(Ms, 1, max, na.rm = TRUE)
    res$.z      <- (szs - rowMeans(Ms, na.rm = TRUE)) /
                    apply(Ms, 1, stats::sd, na.rm = TRUE)
    res$.n_more <- rowSums(Ms >= szs, na.rm = TRUE)          # donors >= SZ07 (in its dir)
    res$.rank   <- res$.n_more + 1L
    res$.stat   <- switch(rank_by,
                          gap   = res$.gap,
                          z     = res$.z,
                          score = szs)

    # ---- functional grouping ----------------------------------------------
    if (!is.null(pathway_groups)) {
        if (is.null(names(pathway_groups)))
            stop("pathway_groups must be a NAMED vector: pathway -> group label")
        res$.fgroup <- unname(pathway_groups[res$pathway])
        res$.fgroup[is.na(res$.fgroup)] <- ungrouped_label
    } else {
        res$.fgroup <- NA_character_
    }

    # ---- per-group row-count resolution (same semantics as fgsea version) --
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
                 " needs group_order to map counts to groups, or pass a NAMED vector.")
        if (length(tpg) != length(group_order))
            stop(sprintf("top_n_per_group has %d value(s) but group_order has %d.",
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
            message("[gsva-specificity] group(s) not covered by top_n_per_group, ",
                    "omitted: ", paste(dropped, collapse = ", "))
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

    # ---- order groups then rows -------------------------------------------
    if (!is.null(pathway_groups)) {
        auto <- names(sort(tapply(sel$.stat, sel$.fgroup, max), decreasing = TRUE))
        glev <- if (is.null(group_order)) auto
                else c(intersect(group_order, auto), setdiff(auto, group_order))
        sel$.fgroup <- factor(sel$.fgroup, levels = glev)
        sel <- sel[order(sel$.fgroup, -sel$.stat), , drop = FALSE]
    } else {
        sel <- sel[order(-sel$.stat), , drop = FALSE]
    }

    # match M rows to selected pathways
    Msel  <- M[match(sel$pathway, pw), , drop = FALSE]
    szsel <- sz[match(sel$pathway, pw)]

    # ---- resolve optional per-pathway vectors to selected-row order --------
    show_spec    <- !is.null(spec_all)
    show_resid_z <- !is.null(resid_z)
    resolve_vec <- function(v, what) {
        if (is.null(v)) return(rep(NA_real_, nrow(sel)))
        v <- as.numeric(v)
        if (!is.null(names(v)) && any(nzchar(names(v)))) {
            unname(v[sel$pathway])
        } else {
            if (length(v) != length(pw))
                stop(sprintf(paste0("unnamed `%s` has length %d but must equal nrow(scores) ",
                     "(%d), or be NAMED by pathway."), what, length(v), length(pw)))
            v[match(sel$pathway, pw)]
        }
    }
    spec_all_sel <- resolve_vec(spec_all, "spec_all")
    resid_z_sel  <- resolve_vec(resid_z,  "resid_z")

    # ---- column layout: base 4 (pathway, score, donors, strip) + middles ---
    n_middle <- as.integer(show_spec) + as.integer(show_resid_z)
    ncol_tab <- 4L + n_middle
    if (is.null(colwidths)) colwidths <- c(9, 1.0, 1.3, 5.5)   # base 4
    extra_w <- c(if (show_spec) spec_width, if (show_resid_z) resid_z_width)
    if (length(colwidths) == 4L && length(extra_w) > 0)
        colwidths <- append(colwidths, extra_w, after = 2L)    # insert after score
    if (length(colwidths) != ncol_tab)
        stop(sprintf(paste0("colwidths has length %d but the table has %d columns. Pass 4 ",
             "base widths (pathway, score, donors, strip) -- extra widths are inserted ",
             "automatically -- or exactly %d."), length(colwidths), ncol_tab, ncol_tab))

    # ---- shared score scale across rows -----------------------------------
    xr <- range(c(as.vector(Msel), szsel, 0), na.rm = TRUE)
    xr <- xr + c(-1, 1) * 0.04 * diff(xr)

    make_strip <- function(i) {
        d_o <- data.frame(v = as.numeric(Msel[i, ]), grp = grp_vec,
                          stringsAsFactors = FALSE)
        d_s <- data.frame(v = szsel[i])
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
                ggplot2::aes(x = v, xend = v, y = 0.12, yend = 0.88, colour = grp),
                linewidth = 0.7, show.legend = FALSE) +
            ggplot2::geom_segment(
                data = d_s,
                ggplot2::aes(x = v, xend = v, y = 0.02, yend = 0.98),
                colour = sz07_color, linewidth = 1.5) +
            ggplot2::scale_colour_manual(values = group_colors)
        ggplot2::ggplotGrob(p)
    }

    # ---- formatters --------------------------------------------------------
    fmt_score <- function(x) if (is.na(x)) "NA"
                             else sprintf(paste0("%.", score_digits, "f"), x)
    fmt_frac  <- function(n, N) if (is.na(n)) "NA" else sprintf("%d/%d", as.integer(n), N)
    fmt_spec  <- function(x) if (is.na(x)) "NA"
                             else sprintf(paste0("%+.", spec_digits, "f"), x)
    fmt_z     <- function(z) if (is.na(z)) "NA"
                             else sprintf(paste0("%+.", resid_z_digits, "f"), z)

    make_label <- function(s) {
        if (!is.null(label_fn)) {
            if (!is.function(label_fn)) stop("label_fn must be a function")
            s <- label_fn(s)
        }
        # word-truncation (keeps the name's END visible) takes precedence over
        # wrapping / end-shortening; the three are mutually-exclusive strategies.
        if (!is.null(label_truncate))
            return(truncate_words(s, width = label_truncate))
        if (!is.null(label_wrap))
            return(vapply(s, function(x)
                paste(strwrap(x, width = label_wrap), collapse = "\n"),
                character(1), USE.NAMES = FALSE))
        if (!is.null(shorten_names))
            s <- ifelse(nchar(s) > shorten_names,
                        paste0(substr(s, 1, shorten_names - 1), "\u2026"), s)
        s
    }

    # ---- header ------------------------------------------------------------
    header <- c(list(
        grid::textGrob("Pathway", x = 1, hjust = 1,
                       gp = grid::gpar(fontface = "bold", fontsize = 9, lineheight = header_lineheight)),
        grid::textGrob(score_header,
                       gp = grid::gpar(fontface = "bold", fontsize = 9, lineheight = header_lineheight))),
        if (show_spec) list(grid::textGrob(spec_header,
            gp = grid::gpar(fontface = "bold", fontsize = 8, lineheight = header_lineheight))),
        if (show_resid_z) list(grid::textGrob(resid_z_header,
            gp = grid::gpar(fontface = "bold", fontsize = 8, lineheight = header_lineheight))),
        list(
        grid::textGrob("donors\n<SZ07", gp = grid::gpar(fontface = "bold", fontsize = 8, lineheight = header_lineheight)),
        grid::textGrob(strip_header, gp = grid::gpar(fontface = "bold", fontsize = 9, lineheight = header_lineheight))
    ))

    make_data_row <- function(i) {
        c(list(
            grid::textGrob(make_label(sel$pathway[i]), x = 1, hjust = 1,
                           gp = grid::gpar(fontsize = label_size)),
            grid::textGrob(fmt_score(sel$sz07_score[i]), gp = grid::gpar(fontsize = 7.5))),
            if (show_spec) list(grid::textGrob(fmt_spec(spec_all_sel[i]),
                gp = grid::gpar(fontsize = 7.5))),
            if (show_resid_z) list(grid::textGrob(fmt_z(resid_z_sel[i]),
                gp = grid::gpar(fontsize = 7.5))),
            list(
            grid::textGrob(fmt_frac(sel$.n_more[i], n_other), gp = grid::gpar(fontsize = 7.5)),
            make_strip(i)
        ))
    }
    make_group_row <- function(label) c(list(
        grid::textGrob(label, x = 1, hjust = 1,
                       gp = grid::gpar(fontface = "bold.italic", fontsize = 8.5, col = "grey25")),
        grid::nullGrob()),
        if (show_spec) list(grid::nullGrob()),
        if (show_resid_z) list(grid::nullGrob()),
        list(grid::nullGrob(),
        grid::linesGrob(x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(0.5, 0.5), "npc"),
                        gp = grid::gpar(col = "grey85", lwd = 0.8))
    ))

    # ---- assemble rows -----------------------------------------------------
    row_grobs <- list(); row_h <- numeric(0)
    if (is.null(pathway_groups)) {
        for (i in seq_len(nrow(sel))) {
            row_grobs[[length(row_grobs) + 1L]] <- make_data_row(i); row_h <- c(row_h, 1)
        }
    } else {
        for (gname in levels(droplevels(sel$.fgroup))) {
            idx <- which(as.character(sel$.fgroup) == gname)
            if (!length(idx)) next
            row_grobs[[length(row_grobs) + 1L]] <- make_group_row(gname)
            row_h <- c(row_h, group_header_height)
            for (i in idx) {
                row_grobs[[length(row_grobs) + 1L]] <- make_data_row(i); row_h <- c(row_h, 1)
            }
        }
    }

    # ---- axis + legend rows -----------------------------------------------
    axis_p <- ggplot2::ggplot(data.frame(x = xr, y = c(0, 0))) +
        ggplot2::geom_blank(ggplot2::aes(x, y)) +
        ggplot2::scale_x_continuous(limits = xr, expand = c(0, 0)) +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(
            panel.grid   = ggplot2::element_blank(),
            axis.title   = ggplot2::element_blank(),
            axis.text.y  = ggplot2::element_blank(),
            axis.line.x  = ggplot2::element_line(colour = "grey40"),
            axis.ticks.x = ggplot2::element_line(colour = "grey40"),
            plot.margin  = grid::unit(c(8, 2, 0, 2), "pt"))
    axis_row <- c(list(grid::nullGrob(), grid::nullGrob()),
                  if (show_spec) list(grid::nullGrob()),
                  if (show_resid_z) list(grid::nullGrob()),
                  list(grid::nullGrob(), ggplot2::ggplotGrob(axis_p)))

    leg_lab <- c(sz07_label, names(group_colors))
    leg_col <- c(sz07_color, unname(group_colors))
    leg_df  <- data.frame(x = seq_along(leg_lab), y = 1, lab = leg_lab, stringsAsFactors = FALSE)
    leg_p <- ggplot2::ggplot(leg_df, ggplot2::aes(x, y)) +
        ggplot2::geom_segment(ggplot2::aes(x = x - 0.12, xend = x - 0.12, y = 0.6, yend = 1.4),
                              colour = leg_col, linewidth = 1.2) +
        ggplot2::geom_text(ggplot2::aes(x = x - 0.05, label = lab), hjust = 0, size = 2.7) +
        ggplot2::scale_x_continuous(limits = c(0.5, length(leg_lab) + 1.2)) +
        ggplot2::scale_y_continuous(limits = c(0.4, 1.6)) +
        ggplot2::theme_void()
    legend_row <- c(list(grid::nullGrob(), grid::nullGrob()),
                    if (show_spec) list(grid::nullGrob()),
                    if (show_resid_z) list(grid::nullGrob()),
                    list(grid::nullGrob(), ggplot2::ggplotGrob(leg_p)))

    all_grobs <- c(header, unlist(row_grobs, recursive = FALSE), axis_row, legend_row)
    g <- gridExtra::arrangeGrob(
        grobs   = all_grobs,
        ncol    = ncol_tab,
        nrow    = length(row_h) + 3L,
        widths  = grid::unit(colwidths, "null"),
        heights = grid::unit(c(header_height, row_h, 1.1, 0.7), "null"))

    if (render) { grid::grid.newpage(); grid::grid.draw(g) }

    out <- sel[, c("pathway", "sz07_score"), drop = FALSE]
    if (!is.null(pathway_groups)) out$functional_group <- as.character(sel$.fgroup)
    out$gap_runnerup   <- sel$.gap
    out$z_vs_donors    <- sel$.z
    out$sz07_rank      <- sel$.rank
    out$n_more_extreme <- sel$.n_more
    out$n_donors       <- n_other
    out$donors_lt_SZ07 <- sprintf("%d/%d", as.integer(sel$.n_more), n_other)
    if (show_spec)    out$spec_neuronal_all <- spec_all_sel
    if (show_resid_z) out$sz07_residual_z   <- resid_z_sel
    rownames(out) <- NULL
    invisible(list(grob = g, data = out))
}
