#' Prioritise candidate genes by min-p enrichment overlap, matched on pathway
#' SIZE PROFILE -- run for BOTH membership bases (set + leading edge).
#'
#' Per-gene statistic: min-p = the smallest donor enrichment p-value across the
#' pathways the gene "belongs to". Belonging is defined TWO ways, each scored as
#' its own test and returned together:
#'   test = "set"          gene is a MEMBER of the pathway gene set (as given in
#'                         `pathway_list`) -- the original behaviour.
#'   test = "leading_edge" gene is in the pathway's GSEA LEADING EDGE (from
#'                         `fgsea_res[[le_col]]`) -- among the genes that
#'                         actually drive the enrichment. A strict subset of
#'                         membership, and the more stringent question.
#'
#' Both tests share all machinery and differ ONLY in the gene->pathway map.
#' Crucially, for each test the NULL POOL is scored under the SAME belonging
#' rule as the candidate. Leading-edge membership is a strict subset of set
#' membership, so a gene's min-p is drawn from fewer pathways and is expected to
#' be LARGER under the leading-edge rule; scoring the candidate on leading-edge
#' pathways while scoring the null on set membership would inflate the empirical
#' p. Both must use the same basis, so the matched null pool AND the matching
#' strata are rebuilt separately for each test (the size-profile distribution
#' differs between bases).
#'
#' MATCHING. The null pool is matched on the gene's pathway SIZE PROFILE -- the
#' counts of its pathways in each size class (small/specific vs large/broad) --
#' rather than on total pathway count (hubness). Total hubness saturates for
#' heavily-annotated genes; the size profile subsumes hubness (its bins sum to
#' the total count) while separating "has specific-pathway opportunity" from
#' "broad-only", the axis that carries the signal. Under the leading_edge test
#' the size class is still the FULL pathway size (a pathway's identity, hence
#' its p-value, is unchanged); the profile just counts pathways the gene is in
#' the leading edge OF, by size class -- the correct "opportunity" covariate for
#' that statistic. Bins keep tied counts together.
#'
#' CIRCULARITY / leave-one-out. min-p uses the FIXED donor pathway p-values in
#' `fgsea_res`. For a candidate that drives one of its pathways (is in that
#' pathway's leading edge), its min-p is self-vouching; recompute enrichment
#' with the candidate removed from the ranking and pass its leave-one-out run
#' via `fgsea_res_g` to score that candidate's min-p on de-biased p-values.
#' NOTE for the leading_edge test: the candidate's leading-edge MEMBERSHIP is
#' taken from the ORIGINAL run, never the leave-one-out run -- removing the gene
#' would by construction drop it from every leading edge and empty the
#' statistic. Only the pathway P-VALUES are overlaid from the leave-one-out run.
#' This is the same candidate-p / null-p asymmetry as the set test, now also
#' spanning membership. Membership-only candidates barely move their pathways'
#' p-values, so the shared `fgsea_res` is fine for them.
#'
#' @param candidates  Character vector of candidate gene symbols.
#' @param fgsea_res   data.frame/data.table with a pathway-name column, a
#'                    p-value column, and (for the leading_edge test) a
#'                    leading-edge column. Every contributing pathway appears.
#' @param pathway_list Named list: list("pathway" = c(genes...)). The tested
#'                    collection (uncollapsed).
#' @param universe    Character vector the null pool is drawn from -- the
#'                    exome-testable genes (NOT the expressed set). Candidates
#'                    are excluded from the null pool automatically.
#' @param size_breaks Numeric breakpoints defining pathway size classes, passed
#'                    to cut() on pathway gene-set size. Default c(0,100,Inf).
#' @param profile_bins Max bins per size-class count axis when forming matching
#'                    strata (default 3). Value-based, ties kept together.
#' @param pathway_col,p_col Column names in fgsea_res (default "pathway","pval").
#'                    Use raw "pval" not "padj".
#' @param size_col   Column in `fgsea_res` giving each pathway's EFFECTIVE size
#'                    -- the number of set members present in the ranked list,
#'                    post-intersection, which is what fgsea's minSize/maxSize
#'                    filtered on (default "size", as fgsea reports it). Using it
#'                    keeps path_size, size_breaks, and the fgsea size filter on
#'                    one scale. If absent, the function falls back to the full
#'                    annotation gene-set length and warns -- that scale differs
#'                    from fgsea's effective size and can misplace or silently
#'                    drop pathways in the size profile (e.g. a 312-gene set that
#'                    fgsea tested at effective size 91).
#' @param le_col     Column in `fgsea_res` holding each pathway's leading-edge
#'                    genes, used by the leading_edge test (default
#'                    "leadingEdge"). Accepts a fgsea list-column (a character
#'                    vector per pathway) OR a delimited string (",", ";", "|",
#'                    or whitespace between symbols). If the test is requested
#'                    but this column is absent, that test is skipped (warning).
#' @param tests      Which test(s) to run: any of "set","leading_edge"
#'                    (default both).
#' @param le_sig_col Column in `fgsea_res` used to gate which pathways may
#'                    contribute LEADING-EDGE membership, so that only ENRICHED
#'                    pathways count (default "padj"). A leading edge from a
#'                    non-enriched pathway (p ~ 1) is meaningless; because the
#'                    FULL fgsea table is passed for effective sizes, it also
#'                    holds such pathways, and without this gate a candidate
#'                    would be scored as "driving" them. Pathways with
#'                    le_sig_col <= le_sig_threshold contribute leading edges;
#'                    others are dropped from the leading-edge map only (they
#'                    remain in the size profile). Set NULL to disable the gate
#'                    (not recommended; warns). Affects the leading_edge test
#'                    only; the set test is unchanged.
#' @param le_sig_threshold Cutoff applied to le_sig_col (default 0.05).
#' @param le_direction Restrict leading-edge pathways by enrichment sign, using
#'                    the NES (or ES) column: "any" (default), "pos" (only
#'                    positively-enriched), or "neg". Use when your ranking is
#'                    signed and you only want drivers of one tail.
#' @param fgsea_res_g Optional per-candidate leave-one-out enrichment, used ONLY
#'                    to score the candidate's own min-p (never the null pool).
#'                    Single data.frame (applied to every candidate) or a named
#'                    list gene -> data.frame. FULL re-run or PARTIAL patch; the
#'                    p-values are overlaid on `fgsea_res`. Only p-values are
#'                    used -- leading-edge membership always comes from the
#'                    original run (see CIRCULARITY note). NULL = no correction.
#' @param min_stratum Warn if a candidate's matched null pool is smaller than
#'                    this. Default 20.
#' @param seed        Unused by the core test; kept for API parity.
#'
#' @return data.frame in LONG form, one row per (candidate x test), with a
#'   leading `test` column then: gene, min_p, best_pathway, best_pathway_size,
#'   n_pathways, profile, stratum, n_null, n_null_le, emp_p, emp_p_bh (BH
#'   computed WITHIN each test as its own family). Sorted by test, then emp_p.
#'   Split downstream with res[res$test == "leading_edge", ] etc.
minp_size_matched_test <- function(candidates,
                                   fgsea_res,
                                   pathway_list,
                                   universe,
                                   size_breaks  = c(0, 100, Inf),
                                   profile_bins = 3,
                                   pathway_col  = "pathway",
                                   p_col        = "pval",
                                   le_col       = "leadingEdge",
                                   le_sig_col       = "padj",
                                   le_sig_threshold = 0.05,
                                   le_direction     = "any",
                                   size_col     = "size",
                                   tests        = c("set", "leading_edge"),
                                   fgsea_res_g  = NULL,
                                   min_stratum  = 20,
                                   seed         = NULL) {

    fr <- as.data.frame(fgsea_res)
    stopifnot(pathway_col %in% names(fr), p_col %in% names(fr))
    candidates <- unique(as.character(candidates))
    universe   <- unique(as.character(universe))
    tests <- match.arg(tests, c("set", "leading_edge"), several.ok = TRUE)

    # ---- pathway-level lookups: p-value and size class (shared) -----------
    sets <- lapply(pathway_list, function(g) unique(as.character(g)))
    path_pval  <- setNames(as.numeric(fr[[p_col]]), fr[[pathway_col]])

    # Pathway SIZE must be the EFFECTIVE size fgsea used -- the count of set
    # members present in the ranked list (post-intersection), which is what
    # minSize/maxSize filtered on -- NOT the full annotation size. fgsea reports
    # it in `size_col`; use that so path_size, size_breaks, and the fgsea filter
    # all refer to the same quantity. Fall back to annotation length only if the
    # column is absent (warns: this reintroduces the annotation-vs-effective
    # mismatch and can misplace/drop pathways in the size profile).
    if (size_col %in% names(fr)) {
        path_size <- setNames(as.integer(fr[[size_col]]), fr[[pathway_col]])
        # pathways in pathway_list but not in fr (untested) have no effective
        # size; leave NA so they are excluded from size classing consistently.
        miss <- setdiff(names(sets), names(path_size))
        if (length(miss)) path_size[miss] <- NA_integer_
        path_size <- path_size[names(sets)]
    } else {
        warning("size_col '", size_col, "' not in fgsea_res; falling back to ",
                "annotation gene-set length. This can differ from fgsea's ",
                "effective (post-intersection) size and misplace pathways in ",
                "the size profile -- pass fgsea's `size` column for consistency.")
        path_size <- vapply(sets, length, integer(1))
    }
    size_class <- as.integer(cut(path_size, breaks = size_breaks,
                                 include.lowest = TRUE, labels = FALSE))
    names(size_class) <- names(sets)
    n_class <- length(size_breaks) - 1L
    class_labels <- paste0("nsize", seq_len(n_class))

    keep_genes <- union(universe, candidates)

    # ---- invert a membership definition to gene -> pathways ---------------
    # (restricted to universe+candidates). Used for both bases.
    build_g2p <- function(memb_sets) {
        glen <- lengths(memb_sets)
        long_gene <- unlist(memb_sets, use.names = FALSE)
        long_path <- rep(names(memb_sets), times = glen)
        inuniv <- long_gene %in% keep_genes
        split(long_path[inuniv], long_gene[inuniv])
    }

    g2p_set <- build_g2p(sets)

    # leading-edge membership map (parsed from le_col), if the test is on
    g2p_le <- NULL
    if ("leading_edge" %in% tests) {
        if (!(le_col %in% names(fr))) {
            warning(sprintf("leading_edge test requested but column '%s' is not "
                            , le_col),
                    "in fgsea_res; skipping the leading_edge test.")
            tests <- setdiff(tests, "leading_edge")
        } else {
            # A leading edge is only meaningful for an ENRICHED pathway. When the
            # full fgsea table is passed (needed for effective sizes), it also
            # contains non-enriched pathways (p ~ 1) whose "leading edge" is
            # noise -- a candidate must NOT be counted as driving those. Gate the
            # pathways that contribute leading-edge membership to the enriched
            # set: significance (le_sig_col <= le_sig_threshold) and, optionally,
            # NES direction. Pathways failing the gate are dropped from g2p_le
            # (but remain in the size profile via g2p_set / path_size).
            enr <- rep(TRUE, nrow(fr))
            if (!is.null(le_sig_col)) {
                if (le_sig_col %in% names(fr)) {
                    enr <- enr & !is.na(fr[[le_sig_col]]) &
                        as.numeric(fr[[le_sig_col]]) <= le_sig_threshold
                } else {
                    warning("le_sig_col '", le_sig_col, "' not in fgsea_res; ",
                            "leading-edge membership will NOT be significance-",
                            "filtered -- non-enriched pathways may contribute ",
                            "spurious leading edges. Pass the significant-",
                            "pathway p/padj column, or set le_sig_col=NULL to ",
                            "silence this.")
                }
            }
            if (!is.null(le_direction) && le_direction != "any") {
                nes_col <- if ("NES" %in% names(fr)) "NES" else
                           if ("ES" %in% names(fr)) "ES" else NA_character_
                if (is.na(nes_col)) {
                    warning("le_direction='", le_direction, "' requested but no ",
                            "NES/ES column in fgsea_res; ignoring direction.")
                } else {
                    v <- as.numeric(fr[[nes_col]])
                    enr <- enr & !is.na(v) &
                        (if (le_direction == "pos") v > 0 else v < 0)
                }
            }

            le_raw <- fr[[le_col]][enr]
            le_paths <- fr[[pathway_col]][enr]
            if (is.list(le_raw)) {
                le <- lapply(le_raw, function(g) unique(as.character(g)))
            } else {
                le <- lapply(strsplit(as.character(le_raw), "[,;|[:space:]]+"),
                             function(g) unique(g[nzchar(g)]))
            }
            names(le) <- le_paths
            le <- le[names(le) %in% names(sets)]     # size known for these only
            g2p_le <- build_g2p(le)
            message(sprintf("[minp] leading_edge test: %d of %d pathways pass the "
                            , sum(enr), nrow(fr)),
                    sprintf("enrichment gate (%s<=%.3g%s) and contribute leading edges.",
                            if (is.null(le_sig_col)) "none" else le_sig_col,
                            le_sig_threshold,
                            if (!is.null(le_direction) && le_direction != "any")
                                paste0(", NES ", le_direction) else ""))
        }
    }
    if (length(tests) == 0L)
        stop("no tests to run (leading_edge requested but ", le_col,
             " absent, and 'set' not selected).")

    # ---- per-gene statistic under a given gene->pathway map ---------------
    per_gene <- function(gene, g2p_x, pvsource) {
        paths <- g2p_x[[gene]]
        if (is.null(paths) || length(paths) == 0) {
            prof <- setNames(integer(n_class), class_labels)
            return(list(min_p = NA_real_, best = NA_character_,
                        best_size = NA_integer_, npath = 0L, prof = prof))
        }
        pv <- pvsource[paths]
        sc <- size_class[paths]
        prof <- setNames(tabulate(sc, nbins = n_class), class_labels)
        ok <- !is.na(pv)
        if (!any(ok)) {
            return(list(min_p = NA_real_, best = NA_character_,
                        best_size = NA_integer_, npath = length(paths), prof = prof))
        }
        wmin <- which.min(pv)
        list(min_p = unname(pv[wmin]),
             best  = paths[wmin],
             best_size = unname(path_size[paths[wmin]]),
             npath = length(paths),
             prof  = prof)
    }

    # p-value source for a candidate: original path_pval with its leave-one-out
    # run overlaid where provided (full re-run or partial patch). p-values only.
    get_cand_pvsource <- function(gene) {
        if (is.null(fgsea_res_g)) return(path_pval)
        tab <- if (is.data.frame(fgsea_res_g)) {
                   fgsea_res_g
               } else if (is.list(fgsea_res_g) && gene %in% names(fgsea_res_g)) {
                   fgsea_res_g[[gene]]
               } else NULL
        if (is.null(tab)) return(path_pval)
        tab <- as.data.frame(tab)
        gpv <- setNames(as.numeric(tab[[p_col]]), tab[[pathway_col]])
        src <- path_pval
        common <- intersect(names(gpv), names(src))
        src[common] <- gpv[common]
        src
    }

    # ---- score all candidates under one membership basis ------------------
    score_membership <- function(g2p_x, label) {
        # null-pool per-gene stats -- ALWAYS the original run
        univ_genes <- names(g2p_x)
        univ_genes <- univ_genes[univ_genes %in% universe]
        up <- lapply(univ_genes, per_gene, g2p_x = g2p_x, pvsource = path_pval)
        univ <- data.frame(
            gene  = univ_genes,
            min_p = vapply(up, `[[`, numeric(1), "min_p"),
            npath = vapply(up, `[[`, integer(1), "npath"),
            stringsAsFactors = FALSE
        )
        prof_mat <- do.call(rbind, lapply(up, `[[`, "prof"))
        colnames(prof_mat) <- class_labels
        univ <- cbind(univ, prof_mat)

        # matching strata from size-class counts (ties together), THIS basis
        make_axis_bin <- function(x, b) {
            br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = b + 1),
                                         na.rm = TRUE))
            if (length(br) < 2) return(rep(1L, length(x)))
            as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
        }
        axis_bins <- lapply(class_labels, function(cl) make_axis_bin(univ[[cl]], profile_bins))
        names(axis_bins) <- class_labels
        axis_breaks <- lapply(class_labels, function(cl) {
            br <- unique(stats::quantile(univ[[cl]], probs = seq(0, 1, length.out = profile_bins + 1),
                                         na.rm = TRUE))
            if (length(br) < 2) NULL else br
        })
        names(axis_breaks) <- class_labels
        assign_axis_bin <- function(value, cl) {
            br <- axis_breaks[[cl]]
            if (is.null(br)) return(1L)
            b <- as.integer(cut(value, breaks = br, include.lowest = TRUE, labels = FALSE))
            if (is.na(b)) b <- if (value <= br[1]) 1L else (length(br) - 1L)
            b
        }
        univ$stratum <- do.call(paste, c(lapply(class_labels, function(cl) axis_bins[[cl]]), sep = "|"))

        null_pool <- univ[!(univ$gene %in% candidates) & !is.na(univ$min_p), , drop = FALSE]

        out <- vector("list", length(candidates))
        for (i in seq_along(candidates)) {
            g <- candidates[i]
            pg <- per_gene(g, g2p_x = g2p_x, pvsource = get_cand_pvsource(g))

            cand_bins <- vapply(class_labels, function(cl)
                assign_axis_bin(pg$prof[[cl]], cl), integer(1))
            cand_stratum <- paste(cand_bins, collapse = "|")

            pool <- null_pool[null_pool$stratum == cand_stratum, , drop = FALSE]
            n_null <- nrow(pool)
            prof_str <- paste(sprintf("%s=%d", class_labels, pg$prof), collapse = ", ")

            if (is.na(pg$min_p)) {
                emp_p <- NA_real_; n_le <- NA_integer_
            } else {
                n_le  <- sum(pool$min_p <= pg$min_p)
                emp_p <- (1 + n_le) / (1 + n_null)
            }

            row <- data.frame(
                test = label, gene = g, min_p = pg$min_p, best_pathway = pg$best,
                best_pathway_size = pg$best_size, n_pathways = pg$npath,
                profile = prof_str, stratum = cand_stratum,
                n_null = n_null, n_null_le = n_le, emp_p = emp_p,
                stringsAsFactors = FALSE
            )
            if (!is.na(emp_p) && n_null < min_stratum) {
                warning(sprintf("[%s] candidate %s: matched null pool = %d (< min_stratum=%d); "
                                , label, g, n_null, min_stratum),
                        "empirical p may be unstable -- consider coarser size_breaks/profile_bins.")
            }
            out[[i]] <- row
        }
        res <- do.call(rbind, out)
        res$emp_p_bh <- p.adjust(res$emp_p, method = "BH")   # BH within this test
        res
    }

    # ---- run the requested tests and stack them ---------------------------
    parts <- list()
    if ("set" %in% tests)          parts[["set"]]          <- score_membership(g2p_set, "set")
    if ("leading_edge" %in% tests) parts[["leading_edge"]] <- score_membership(g2p_le, "leading_edge")

    res <- do.call(rbind, parts)
    res <- res[order(res$test, res$emp_p, na.last = TRUE), , drop = FALSE]
    rownames(res) <- NULL
    res
}
