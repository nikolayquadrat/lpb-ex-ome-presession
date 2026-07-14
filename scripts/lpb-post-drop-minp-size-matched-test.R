#' Prioritise candidate genes by min-p enrichment overlap, matched on pathway
#' SIZE PROFILE (not total hubness)
#'
#' Per-gene statistic: min-p = the smallest donor enrichment p-value across the
#' pathways the gene belongs to (its single most strongly-enriched pathway).
#' A candidate is "interesting" if its min-p is smaller (more extreme) than that
#' of matched null genes.
#'
#' MATCHING. The null pool is matched on the gene's pathway SIZE PROFILE -- the
#' counts of its pathways in each size class (e.g. small/specific vs large/broad)
#' -- rather than on total pathway count (hubness). Total hubness is a poor
#' covariate for min-p: it saturates for heavily-annotated genes (any gene in
#' ~all pathways has some low-p pathway by breadth alone), so hubness-matched
#' min-p cannot resolve candidates that are all near-maximal hubs. The size
#' profile subsumes hubness (its bins sum to the total count) while separating
#' "has specific-pathway opportunity" from "broad-only", which is the axis that
#' carries the signal. Matching bins keep tied counts together (genes with the
#' same profile are matchable), unlike an equal-frequency split.
#'
#' CIRCULARITY / leave-one-out. min-p uses the FIXED donor pathway p-values in
#' `fgsea_res`. For a candidate that itself drives one of its pathways (i.e. is
#' in that pathway's leading edge), its min-p is self-vouching. Handle that
#' upstream: recompute enrichment with the candidate removed from the ranking
#' and pass its leave-one-out enrichment run via `fgsea_res_g` for that gene.
#' This matters
#' only for leading-edge candidates; membership-only candidates barely move
#' their pathways' p-values, so the shared `fgsea_res` is fine for them.
#'
#' @param candidates  Character vector of candidate gene symbols.
#' @param fgsea_res   data.frame/data.table with a pathway-name column and a
#'                    p-value column (the donor enrichment result). Every pathway
#'                    in `pathway_list` that could contribute should appear here.
#' @param pathway_list Named list: list("pathway" = c(genes...)). The tested
#'                    collection (uncollapsed).
#' @param universe    Character vector of genes the null pool is drawn from --
#'                    the exome-testable genes (NOT the expressed set). Candidates
#'                    are excluded from the null pool automatically.
#' @param size_breaks Numeric breakpoints defining pathway size classes, passed
#'                    to cut() on pathway gene-set size. Default c(0,100,Inf) ->
#'                    two classes: "small" (<=100) and "large" (>100).
#' @param profile_bins Max bins per size-class count axis when forming matching
#'                    strata (default 3). Bins are value-based and KEEP TIES
#'                    TOGETHER, so equal counts share a bin; the realised number
#'                    of bins may be fewer if counts are highly tied (fine here).
#' @param pathway_col,p_col Column names in fgsea_res (default "pathway","pval").
#'                    Use raw "pval" not "padj": padj flooring creates ties at
#'                    the top that destroy min-p resolution.
#' @param fgsea_res_g Optional per-candidate leave-one-out enrichment, used ONLY
#'                    to score the candidate's own min-p (never the null pool,
#'                    which always uses `fgsea_res`). Accepts either:
#'                      - a single data.frame (same columns as fgsea_res): the
#'                        candidate's gene-removed run, applied to every candidate
#'                        in this call (natural when calling with one candidate);
#'                      - a named list mapping candidate gene -> its data.frame:
#'                        each candidate scored on its own run; candidates absent
#'                        from the list fall back to `fgsea_res` (no correction).
#'                    The table may be a FULL re-run (all pathways) or a PARTIAL
#'                    patch (only the pathways whose p changed under leave-one-
#'                    out); values are overlaid on `fgsea_res`, so both work.
#'                    The asymmetry is deliberate and essential: candidate min-p
#'                    on its own gene-removed run, null pool on the original run,
#'                    so empirical p's stay comparable across candidates. NULL
#'                    (default) = no leave-one-out (candidate scored on fgsea_res
#'                    like the null pool).
#' @param min_stratum Warn if a candidate's matched null pool is smaller than
#'                    this (unstable percentile). Default 20.
#' @param seed        Unused by the core test (matching is deterministic); kept
#'                    for API parity. Set for reproducibility if you later add
#'                    stochastic tie-breaking.
#'
#' @return data.frame, one row per candidate: gene, min_p, best_pathway,
#'   best_pathway_size, n_pathways, size profile counts, stratum, n_null,
#'   n_null_le (nulls at least as extreme), emp_p, emp_p_bh. Sorted by emp_p.
minp_size_matched_test <- function(candidates,
                                   fgsea_res,
                                   pathway_list,
                                   universe,
                                   size_breaks  = c(0, 100, Inf),
                                   profile_bins = 3,
                                   pathway_col  = "pathway",
                                   p_col        = "pval",
                                   fgsea_res_g = NULL,
                                   min_stratum  = 20,
                                   seed         = NULL) {

    fr <- as.data.frame(fgsea_res)
    stopifnot(pathway_col %in% names(fr), p_col %in% names(fr))
    candidates <- unique(as.character(candidates))
    universe   <- unique(as.character(universe))

    # ---- pathway-level lookups: p-value and size class --------------------
    sets <- lapply(pathway_list, function(g) unique(as.character(g)))
    path_size  <- vapply(sets, length, integer(1))
    path_pval  <- setNames(as.numeric(fr[[p_col]]), fr[[pathway_col]])
    # size class per pathway (integer 1..K); K = length(size_breaks)-1
    size_class <- as.integer(cut(path_size, breaks = size_breaks,
                                 include.lowest = TRUE, labels = FALSE))
    names(size_class) <- names(sets)
    n_class <- length(size_breaks) - 1L
    class_labels <- paste0("nsize", seq_len(n_class))

    # ---- invert to gene -> pathways (restricted to the universe+candidates) -
    keep_genes <- union(universe, candidates)
    glen <- lengths(sets)
    long_gene <- unlist(sets, use.names = FALSE)
    long_path <- rep(names(sets), times = glen)
    inuniv <- long_gene %in% keep_genes
    g2p <- split(long_path[inuniv], long_gene[inuniv])

    # ---- per-gene: min_p, best pathway, and size-class profile ------------
    # per-gene stats. `pvsource` is a named pathway->pval vector: the ORIGINAL
    # path_pval for null-pool genes, or a candidate's leave-one-out-overlaid
    # vector for a candidate. Size profile is p-value-independent (from sizes).
    per_gene <- function(gene, pvsource) {
        paths <- g2p[[gene]]
        if (is.null(paths) || length(paths) == 0) {
            prof <- setNames(integer(n_class), class_labels)
            return(list(min_p = NA_real_, best = NA_character_,
                        best_size = NA_integer_, npath = 0L, prof = prof))
        }
        pv <- pvsource[paths]
        # size-class profile counts
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

    # p-value source for a candidate: original path_pval, with the candidate's
    # leave-one-out run overlaid where provided (full re-run or partial patch).
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

    # compute for the whole universe (null pool) -- always the ORIGINAL run
    univ_genes <- names(g2p)
    univ_genes <- univ_genes[univ_genes %in% universe]
    up <- lapply(univ_genes, per_gene, pvsource = path_pval)
    univ <- data.frame(
        gene  = univ_genes,
        min_p = vapply(up, `[[`, numeric(1), "min_p"),
        npath = vapply(up, `[[`, integer(1), "npath"),
        stringsAsFactors = FALSE
    )
    prof_mat <- do.call(rbind, lapply(up, `[[`, "prof"))
    colnames(prof_mat) <- class_labels
    univ <- cbind(univ, prof_mat)

    # ---- build matching strata from size-class counts (ties kept together) -
    # value-based bins per axis: equal counts always share a bin.
    make_axis_bin <- function(x, b) {
        br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = b + 1),
                                     na.rm = TRUE))
        if (length(br) < 2) return(rep(1L, length(x)))
        as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
    }
    axis_bins <- lapply(class_labels, function(cl) make_axis_bin(univ[[cl]], profile_bins))
    names(axis_bins) <- class_labels
    # also need bin assignment as a FUNCTION of a count, to place candidates in
    # the same scheme; recover per-axis breakpoints and reuse.
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
        if (is.na(b)) b <- if (value <= br[1]) 1L else (length(br) - 1L)  # clamp
        b
    }
    univ$stratum <- do.call(paste, c(lapply(class_labels, function(cl) axis_bins[[cl]]), sep = "|"))

    # exclude candidates from the null pool, and drop genes with undefined min_p
    null_pool <- univ[!(univ$gene %in% candidates) & !is.na(univ$min_p), , drop = FALSE]

    # ---- score each candidate ---------------------------------------------
    out <- vector("list", length(candidates))
    for (i in seq_along(candidates)) {
        g <- candidates[i]
        pg <- per_gene(g, pvsource = get_cand_pvsource(g))

        # candidate's matching stratum, using the same per-axis binning
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
            gene = g, min_p = pg$min_p, best_pathway = pg$best,
            best_pathway_size = pg$best_size, n_pathways = pg$npath,
            profile = prof_str, stratum = cand_stratum,
            n_null = n_null, n_null_le = n_le, emp_p = emp_p,
            stringsAsFactors = FALSE
        )
        if (!is.na(emp_p) && n_null < min_stratum) {
            warning(sprintf("candidate %s: matched null pool = %d (< min_stratum=%d); "
                            , g, n_null, min_stratum),
                    "empirical p may be unstable -- consider coarser size_breaks/profile_bins.")
        }
        out[[i]] <- row
    }
    res <- do.call(rbind, out)
    res$emp_p_bh <- p.adjust(res$emp_p, method = "BH")
    res <- res[order(res$emp_p, na.last = TRUE), , drop = FALSE]
    rownames(res) <- NULL
    res
}