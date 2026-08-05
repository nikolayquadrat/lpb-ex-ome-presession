#' Per-pathway correlation of donor NES with cell-type composition, computed
#' WITHIN groups (healthy controls; patients excluding SZ07), plus a test of
#' whether SZ07's own NES is explained by its composition.
#'
#' WHY WITHIN-GROUP, EXCLUDING SZ07. The question is whether a pathway's
#' enrichment tracks cell composition. SZ07 is the outlier being explained, so it
#' must not be in the correlation (it would drive it). We measure the
#' composition->NES relationship in the OTHER donors, separately for HC and for
#' non-SZ07 patients, then predict where SZ07 *should* sit given its fraction and
#' compare to where it actually sits. Interpretation:
#'   * pathway with POSITIVE neuronal specificity (e.g. synaptic): a real
#'     composition effect shows NES rising with neuronal fraction within groups.
#'     If SZ07 has LOW neuronal fraction but HIGH NES (large positive residual),
#'     its enrichment runs AGAINST composition -> genuine regulatory signal.
#'   * pathway with NEGATIVE neuronal specificity (e.g. immune/complement): NES
#'     falls with neuronal fraction within groups. If SZ07 (low neurons) sits ON
#'     that line (small residual), its enrichment is EXPLAINED by composition.
#'
#' Reuses the plot's data contract: the same `res`, `donor_ids`, `donor_groups`.
#' `other_donors_NES` already excludes SZ07, so the HC/SZ groups here are exactly
#' "controls" and "patients without SZ07".
#'
#' @param res            data.frame used by plot_donor_specificity_table.
#' @param donor_ids      other-donor names, in the column order of
#'                       `others_nes_col` (POSITIONAL, same as the plot).
#' @param donor_groups   group per other-donor (same length/order), e.g. HC/SZ.
#' @param donor_fraction NAMED numeric vector: cell-type fraction per donor
#'                       (e.g. mean CIBERSORTx neuronal fraction across that
#'                       donor's samples), names matching `donor_ids`.
#' @param sz07_fraction  scalar: SZ07's own cell-type fraction (mean across its
#'                       samples). If given, adds the SZ07-residual columns.
#' @param method         "spearman" (default; robust, monotonic) or "pearson".
#' @param predict_from   which group(s) to fit for the SZ07 prediction. Default
#'                       c("HC","SZ") = all non-SZ07 donors (general relationship).
#' @param min_n          minimum donors in a group to report a correlation.
#' @param pathway_col,others_nes_col,sz07_nes_col  column names in `res`.
#' @return data.frame, one row per pathway, with per-group rho/p/n, and (if
#'   `sz07_fraction` given) SZ07 predicted NES, residual, and residual z-score
#'   (residual in SD units of the fitted group's residuals).
pathway_composition_correlation <- function(
        res, donor_ids, donor_groups, donor_fraction,
        sz07_fraction   = NULL,
        method          = c("spearman", "pearson"),
        predict_from    = c("HC", "SZ"),
        min_n           = 4,
        pathway_col     = "pathway",
        others_nes_col  = "other_donors_NES",
        sz07_nes_col    = "SZ07_simple_NES") {

    method <- match.arg(method)
    res <- as.data.frame(res, stringsAsFactors = FALSE)

    # ---- parse positional per-donor NES (same as the plot) ----------------
    parse_row <- function(s) {
        if (is.na(s) || !nzchar(s)) return(numeric(0))
        as.numeric(strsplit(gsub("\\s", "", s), ",")[[1]])
    }
    other_list <- lapply(res[[others_nes_col]], parse_row)
    nD <- unique(lengths(other_list))
    if (length(nD) != 1L)
        stop("rows of ", others_nes_col, " have differing donor counts")
    if (length(donor_ids) != nD || length(donor_groups) != nD)
        stop(sprintf("donor_ids/donor_groups length (%d/%d) must equal donors in %s (%d)",
                     length(donor_ids), length(donor_groups), others_nes_col, nD))
    M <- do.call(rbind, other_list)             # pathways x donors
    colnames(M) <- donor_ids

    # ---- align the composition vector to donor order ----------------------
    if (is.null(names(donor_fraction)))
        stop("donor_fraction must be NAMED by donor id")
    miss <- setdiff(donor_ids, names(donor_fraction))
    if (length(miss))
        stop("donor_fraction missing donor(s): ", paste(miss, collapse = ", "))
    fr <- as.numeric(donor_fraction[donor_ids])
    grp <- as.character(donor_groups)

    groups <- sort(unique(grp))
    all_lab <- "all_others"                     # HC + SZ combined (no SZ07)

    corr_in <- function(idx, y) {
        x <- fr[idx]; yy <- y[idx]
        ok <- is.finite(x) & is.finite(yy)
        if (sum(ok) < min_n) return(c(rho = NA, p = NA, n = sum(ok)))
        ct <- suppressWarnings(stats::cor.test(x[ok], yy[ok], method = method))
        c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
    }

    fit_predict_sz07 <- function(y) {
        # linear fit of NES ~ fraction over predict_from donors; predict SZ07 at
        # its own fraction, report residual and residual z (SD of fit residuals).
        if (is.null(sz07_fraction)) return(c(pred = NA, resid = NA, z = NA))
        idx <- which(grp %in% predict_from)
        x <- fr[idx]; yy <- y[idx]
        ok <- is.finite(x) & is.finite(yy)
        if (sum(ok) < max(min_n, 3)) return(c(pred = NA, resid = NA, z = NA))
        fit <- stats::lm(yy[ok] ~ x[ok])
        b <- unname(coef(fit))
        pred <- b[1] + b[2] * sz07_fraction
        sdr  <- stats::sd(stats::residuals(fit))
        obs  <- NA_real_                        # filled by caller (SZ07 NES)
        c(pred = pred, resid = NA, z = NA, sd = sdr, b0 = b[1], b1 = b[2])
    }

    sz07_nes <- if (sz07_nes_col %in% names(res)) as.numeric(res[[sz07_nes_col]]) else NA

    rows <- lapply(seq_len(nrow(M)), function(i) {
        y <- as.numeric(M[i, ])
        rec <- list(pathway = res[[pathway_col]][i])
        # per-group correlations (HC, SZ, ...) and combined
        for (g in groups) {
            c3 <- corr_in(which(grp == g), y)
            rec[[paste0("rho_", g)]] <- unname(c3["rho"])
            rec[[paste0("p_",   g)]] <- unname(c3["p"])
            rec[[paste0("n_",   g)]] <- unname(c3["n"])
        }
        c3 <- corr_in(seq_along(y), y)          # all others combined
        rec[[paste0("rho_", all_lab)]] <- unname(c3["rho"])
        rec[[paste0("p_",   all_lab)]] <- unname(c3["p"])
        rec[[paste0("n_",   all_lab)]] <- unname(c3["n"])
        # SZ07 residual vs the predict_from fit
        if (!is.null(sz07_fraction)) {
            fp <- fit_predict_sz07(y)
            obs <- sz07_nes[i]
            resid <- if (is.finite(fp["pred"]) && is.finite(obs)) obs - fp["pred"] else NA
            z <- if (is.finite(resid) && is.finite(fp["sd"]) && fp["sd"] > 0)
                     resid / fp["sd"] else NA
            rec$sz07_fraction     <- sz07_fraction
            rec$sz07_NES_observed <- obs
            rec$sz07_NES_predicted <- unname(fp["pred"])
            rec$sz07_residual     <- unname(resid)
            rec$sz07_residual_z   <- unname(z)
        }
        rec
    })

    out <- do.call(rbind, lapply(rows, function(r)
        as.data.frame(r, stringsAsFactors = FALSE)))
    rownames(out) <- NULL
    attr(out, "method") <- method
    attr(out, "predict_from") <- predict_from
    out
}
