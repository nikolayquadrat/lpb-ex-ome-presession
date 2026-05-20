#!/usr/bin/env python3
"""
tier_candidates.py - tier-based VUS prioritization for DROP-directed analysis.

Replaces the original r14_filter_annotations R script.

Pipeline position: VEP TSV (cohort) + per-sample VCF  ->  three TSV tiers.

Tier definitions
----------------
Tier A  (DROP-testable: LoF + splice)
    * LoF consequences with LOFTEE high-confidence flag where present
      (stop_gained, frameshift, splice_donor, splice_acceptor, start_lost).
      If LOFTEE is unavailable, falls back to consequence-only LoF detection.
    * SpliceAI max delta >= 0.2 in any of AG/AL/DG/DL.
    Rationale: predictable mRNA signature (NMD-driven expression loss or
    aberrant splicing) that DROP can empirically validate.

Tier B  (Missense, rank-only)
    * Missense variants with AlphaMissense pathogenicity class
      'likely_pathogenic' (score >= 0.564) in constrained genes
      (gnomAD LOEUF < 0.6), OR REVEL >= 0.75 in the same gene set.
    Rationale: DROP usually cannot confirm missense (normal mRNA, normal
    splicing); this tier feeds downstream hand-curation.

Tier C  (Gene-panel long-shots)
    * ANY rare protein-altering variant in SCHEMA FDR<0.25 / BipEx / ASC /
      DDG2P gene sets, regardless of predicted impact.
    Rationale: prior-probability boost for known psychiatric / NDD genes.

Common filters applied to all tiers:
    * gnomAD v4 popmax AF < 0.001 (ultra-rare; rare-large-effect hypothesis)
    * Sample genotype is heterozygous or homozygous-alt (ALT present)
    * Canonical / MANE-Select transcript only (collapsed upstream by VEP
      --pick_allele_gene)

Inputs handled by this version
------------------------------
The script auto-detects column naming variants from real-world annotation
sources, since file formats and column names drift between releases:

  * VEP TSV column for AlphaMissense:    am_pathogenicity / am_class
                                          (also accepts AlphaMissense_score
                                           and AlphaMissense_class)
  * gnomAD AF columns may appear as:     gnomADv4_AF, AF, AF_popmax, etc.
  * SCHEMA gene column:                  hgnc_symbol or gene_name
  * SCHEMA FDR column:                   "Q meta" or FDR
  * DDG2P gene column:                   gene_symbol / hgnc_symbol / "gene symbol"
  * gnomAD constraint LOEUF column:      lof.oe_ci.upper / oe_lof_upper / loeuf
  * gnomAD constraint pLI column:        lof.pLI / pLI

Gene-level gnomAD constraint (pLI + LOEUF) is joined directly from the
gnomAD constraint TSV (the same file as --loeuf_table) and lands in the
output as gnomADv4_pLI / gnomADv4_LOEUF. dbNSFP 5.3.1a does NOT include
these gene-level metrics, so they cannot come from the dbNSFP plugin.

Missing optional gene-set files (e.g., empty BipEx/ASC placeholder files)
are tolerated; load_gene_set returns an empty set with a warning, the
overall pipeline continues, and Tier C uses the panels it could load.

Outputs
-------
    --out_tier_a / _b / _c  Tier-specific TSVs (one variant per row)
    --out_master            All sample variants with tier column (TSV)
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from io import StringIO

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
AF_MAX               = 0.001    # gnomAD v4 popmax cutoff (ultra-rare)
SPLICEAI_THRESHOLD   = 0.20     # min of max(AG,AL,DG,DL) to flag splice
ALPHAMIS_LIKELY_PATH = 0.564    # AlphaMissense likely_pathogenic threshold
REVEL_THRESHOLD      = 0.773    # ClinGen-endorsed thresholds (0.773 for moderate evidence, 0.932 for strong evidence under ACMG PP3 criteria).
LOEUF_THRESHOLD      = 0.60     # gnomAD constraint: constrained genes (from gnomAD v4.0+ (2024))
SCHEMA_FDR_MAX       = 0.25     # include sub-significant SCHEMA genes
ASC_FDR_MAX          = 0.25     # include sub-significant ASC genes

LOF_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation",
}

MISSENSE_CONSEQUENCE = "missense_variant"


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def parse_vep_tsv(path):
    """
    VEP --tab output has a commented header block; column-name line starts
    with '#Uploaded_variation'. Read that line for columns, then load the
    rest as a DataFrame. VEP uses '-' as the missing-value marker.
    """
    open_fn = gzip.open if path.endswith(".gz") else open
    header = None
    data_lines = []
    with open_fn(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#Uploaded_variation"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            data_lines.append(line)
    if header is None:
        sys.exit("ERROR: could not find VEP TSV header line")
    df = pd.read_csv(StringIO("".join(data_lines)), sep="\t",
                     names=header, na_values=["-", "."], low_memory=False)
    return df


def extract_genotypes_for_sample(vcf_path, sample_id):
    """
    Pull GT/DP/AD/GQ for the sample. Tries bcftools first; falls back to a
    minimal pure-Python VCF parser if bcftools is unavailable in the
    environment (e.g., when this script runs in a pandas-only container).
    """
    cmd = [
        "bcftools", "query",
        "-s", sample_id,
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%DP]\t[%AD]\t[%GQ]\n",
        vcf_path,
    ]
    try:
        out = subprocess.check_output(cmd, text=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as err:
        print(f"WARNING: bcftools unavailable or failed ({err}); "
              f"falling back to pure-Python VCF parser",
              file=sys.stderr)
        return _parse_vcf_fallback(vcf_path, sample_id)

    if not out.strip():
        return pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT",
                                     "GT", "DP", "AD", "GQ"])

    gt = pd.read_csv(StringIO(out), sep="\t", header=None,
                     names=["CHROM", "POS", "REF", "ALT", "GT", "DP", "AD", "GQ"],
                     na_values=["."])
    gt["POS"] = gt["POS"].astype(int)
    return gt


def _parse_vcf_fallback(vcf_path, sample_id):
    """Minimal VCF parser - last resort when bcftools is unavailable."""
    rows = []
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    sample_idx = None
    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if sample_id not in cols:
                    sys.exit(f"ERROR: sample {sample_id} not in VCF header. "
                             f"Available: {cols[9:]}")
                sample_idx = cols.index(sample_id)
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = parts[:5]
            fmt = parts[8].split(":")
            vals = parts[sample_idx].split(":")
            d = dict(zip(fmt, vals))
            rows.append({
                "CHROM": chrom, "POS": int(pos), "REF": ref, "ALT": alt,
                "GT": d.get("GT"), "DP": d.get("DP"),
                "AD": d.get("AD"), "GQ": d.get("GQ"),
            })
    return pd.DataFrame(rows)


def split_location(df):
    """VEP 'Location' is 'chr1:12345' or 'chr1:12345-12346'. Split."""
    loc = df["Location"].astype(str).str.split(":", n=1, expand=True)
    df["CHROM"] = loc[0]
    pos_part = loc[1].str.split("-", n=1, expand=True)[0]
    df["POS"] = pd.to_numeric(pos_part, errors="coerce").astype("Int64")
    return df


def normalize_uploaded_variation(df):
    """
    VEP's '#Uploaded_variation' is CHROM_POS_REF/ALT after normalization.
    Parse REF and ALT so we can join cleanly with per-sample VCF rows.
    """
    pattern = re.compile(r"^(?P<chrom>[^_]+)_(?P<pos>\d+)_(?P<ref>[^/]+)/(?P<alt>.+)$")
    parsed = df["Uploaded_variation"].astype(str).str.extract(pattern)
    df["REF"] = parsed["ref"]
    df["ALT"] = parsed["alt"]
    return df


# ---------------------------------------------------------------------------
# Per-tier classification
# ---------------------------------------------------------------------------
def max_spliceai(row):
    """Max of SpliceAI delta scores across AG/AL/DG/DL (ignores NaN)."""
    vals = []
    for k in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
              "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"):
        v = row.get(k)
        try:
            v = float(v)
            if not np.isnan(v):
                vals.append(v)
        except (TypeError, ValueError):
            pass
    return max(vals) if vals else np.nan


def is_lof(row, require_hc):
    """
    LoF detection. If `require_hc` is True, additionally require LOFTEE's
    LoF=HC flag. If False (LOFTEE not present in this VEP run),
    fall back to consequence-only detection.
    """
    consequences = str(row.get("Consequence", "")).split(",")
    if not any(c in LOF_CONSEQUENCES for c in consequences):
        return False
    if not require_hc:
        return True
    lof_flag = str(row.get("LoF", "")).strip()
    return lof_flag == "HC"


def is_splice_hit(row):
    val = row.get("SpliceAI_max", np.nan)
    if val is None or pd.isna(val):
        return False
    return val >= SPLICEAI_THRESHOLD


def is_damaging_missense(row, constrained_genes,
                         am_score_col, am_class_col, revel_col):
    """
    Damaging missense in a constrained gene per AlphaMissense or REVEL.
    Column names are passed in because they vary across VEP plugin versions.
    """
    consequences = str(row.get("Consequence", "")).split(",")
    if MISSENSE_CONSEQUENCE not in consequences:
        return False
    if row.get("SYMBOL") not in constrained_genes:
        return False

    am_class = str(row.get(am_class_col, "")).lower() if am_class_col else ""
    am_score_raw = row.get(am_score_col) if am_score_col else None
    revel_raw    = row.get(revel_col) if revel_col else None

    try:
        am_score = float(am_score_raw)
    except (TypeError, ValueError):
        am_score = np.nan
    try:
        revel = float(revel_raw)
    except (TypeError, ValueError):
        revel = np.nan

    am_ok    = (am_class == "likely_pathogenic"
                or (not np.isnan(am_score) and am_score >= ALPHAMIS_LIKELY_PATH))
    revel_ok = (not np.isnan(revel) and revel >= REVEL_THRESHOLD)
    return am_ok or revel_ok


def is_protein_altering(row):
    """Any non-synonymous coding or splice-region consequence (Tier C gate)."""
    consequences = str(row.get("Consequence", "")).split(",")
    altering = {
        "missense_variant", "inframe_insertion", "inframe_deletion",
        "splice_region_variant", "protein_altering_variant",
    } | LOF_CONSEQUENCES
    return any(c in altering for c in consequences)


# ---------------------------------------------------------------------------
# Column-name auto-detection (different VEP plugin versions emit different
# column names for the same data)
# ---------------------------------------------------------------------------
def detect_alphamissense_cols(columns):
    """Return (score_col, class_col) names present in the VEP output, or
    (None, None) if AlphaMissense wasn't run."""
    score = next((c for c in ("am_pathogenicity",
                              "AlphaMissense_score",
                              "AlphaMissense_pathogenicity") if c in columns), None)
    klass = next((c for c in ("am_class",
                              "AlphaMissense_class") if c in columns), None)
    return score, klass


def detect_revel_col(columns):
    return next((c for c in ("REVEL", "REVEL_score") if c in columns), None)


def detect_loftee_present(columns):
    """LOFTEE plugin emits a 'LoF' column. Its presence enables HC filtering."""
    return "LoF" in columns


def detect_af_cols(columns):
    """
    Find AF columns from VEP --custom (typical names: gnomADv4_AF,
    gnomADv4_AF_nfe, ...) or from bcftools-injected INFO (AF, AF_nfe, ...).
    Population-stratified columns are included so popmax can be computed.
    """
    af_cols = [c for c in columns if c.startswith("gnomADv4_AF")]
    if af_cols:
        return af_cols
    # Fallback: bare AF / AF_<pop> from upstream bcftools annotate
    af_cols = [c for c in columns if c == "AF" or
               (c.startswith("AF_") and not c.startswith("AF_popmax"))]
    if af_cols:
        return af_cols
    # Last fallback: a precomputed popmax column
    if "AF_popmax" in columns:
        return ["AF_popmax"]
    return []


# ---------------------------------------------------------------------------
# Gene-set loading
# ---------------------------------------------------------------------------
def load_gene_set(path, gene_col=None, fdr_col=None, fdr_max=None,
                  confidence_col=None, min_confidence=None):
    """
    Load a gene set from CSV/TSV. Tolerant of:
      - missing or empty files (returns empty set with warning)
      - varied column names (auto-detection)
      - optional FDR filter (e.g. SCHEMA Q meta)
      - optional confidence-level filter (e.g. PanelApp green)
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"WARNING: {path} not found or empty; gene set will be empty",
              file=sys.stderr)
        return set()

    sep = "," if path.endswith(".csv") else "\t"
    try:
        df = pd.read_csv(path, sep=sep, low_memory=False)
    except Exception as err:
        print(f"WARNING: could not read {path}: {err}; gene set will be empty",
              file=sys.stderr)
        return set()

    if df.empty:
        print(f"WARNING: {path} has no rows; gene set will be empty",
              file=sys.stderr)
        return set()

    # Resolve gene-symbol column. Try the explicit hint first, then a list of
    # common synonyms, then a last-resort case-insensitive search.
    chosen_gene_col = None
    if gene_col and gene_col in df.columns:
        chosen_gene_col = gene_col
    else:
        for cand in ("hgnc_symbol", "gene_symbol", "gene symbol",
                     "gene_name", "gene", "symbol", "Gene", "Symbol"):
            if cand in df.columns:
                chosen_gene_col = cand
                break
    if chosen_gene_col is None:
        # Case-insensitive fallback
        ci_match = next((c for c in df.columns
                         if c.lower() in ("gene", "symbol", "gene_symbol",
                                          "gene_name", "hgnc_symbol")), None)
        chosen_gene_col = ci_match
    if chosen_gene_col is None:
        print(f"WARNING: no recognizable gene-symbol column in {path}; "
              f"available columns: {list(df.columns)[:8]}...",
              file=sys.stderr)
        return set()

    # Apply FDR / confidence filters if requested
    if fdr_col and fdr_col in df.columns and fdr_max is not None:
        before = len(df)
        df = df[pd.to_numeric(df[fdr_col], errors="coerce") <= fdr_max]
        print(f"  {path}: FDR filter {fdr_col}<={fdr_max} kept "
              f"{len(df)}/{before} rows", file=sys.stderr)

    if confidence_col and confidence_col in df.columns and min_confidence is not None:
        before = len(df)
        df = df[pd.to_numeric(df[confidence_col], errors="coerce") >= min_confidence]
        print(f"  {path}: confidence filter {confidence_col}>={min_confidence} "
              f"kept {len(df)}/{before} rows", file=sys.stderr)

    genes = (df[chosen_gene_col]
             .astype(str)
             .str.strip()
             .dropna())
    genes = set(g for g in genes if g and g.lower() not in ("nan", "none", ""))
    print(f"  loaded {len(genes)} genes from {path} (col: {chosen_gene_col})",
          file=sys.stderr)
    return genes


def load_gene_annotations(path, prefix,
                          gene_col=None,
                          fdr_col=None, fdr_max=None,
                          confidence_col=None, min_confidence=None,
                          report_columns=None):
    """
    Like load_gene_set, but returns BOTH the filtered gene set AND a
    per-gene annotation table indexed on the gene symbol. The annotation
    columns (chosen via `report_columns`) are renamed with `prefix_` so
    they can be merged into the variant report without column-name collisions.

    Parameters
    ----------
    path : str
        Path to the gene-set TSV/CSV.
    prefix : str
        Short prefix for output column names, e.g. "SCHEMA" or "DDG2P".
    gene_col, fdr_col, fdr_max, confidence_col, min_confidence
        Same meaning as in load_gene_set.
    report_columns : list[str] or None
        Column names from the source file to carry through to the report.
        Names match against the source file's literal column names; missing
        columns are silently skipped. If None, no annotations are returned.

    Returns
    -------
    genes : set[str]
        The gene symbols that pass the FDR / confidence filter.
    annot : pandas.DataFrame
        DataFrame indexed on gene symbol with columns named `<prefix>_<col>`
        for every column in `report_columns` that was found in the source.
        Empty DataFrame if the file is missing or no report_columns matched.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"WARNING: {path} not found or empty; gene set will be empty",
              file=sys.stderr)
        return set(), pd.DataFrame()

    sep = "," if path.endswith(".csv") else "\t"
    try:
        df = pd.read_csv(path, sep=sep, low_memory=False)
    except Exception as err:
        print(f"WARNING: could not read {path}: {err}; gene set will be empty",
              file=sys.stderr)
        return set(), pd.DataFrame()

    if df.empty:
        print(f"WARNING: {path} has no rows; gene set will be empty",
              file=sys.stderr)
        return set(), pd.DataFrame()

    # Resolve gene-symbol column (same logic as load_gene_set)
    chosen_gene_col = None
    if gene_col and gene_col in df.columns:
        chosen_gene_col = gene_col
    else:
        for cand in ("hgnc_symbol", "gene_symbol", "gene symbol",
                     "gene_name", "gene", "symbol", "Gene", "Symbol"):
            if cand in df.columns:
                chosen_gene_col = cand
                break
    if chosen_gene_col is None:
        ci_match = next((c for c in df.columns
                         if c.lower() in ("gene", "symbol", "gene_symbol",
                                          "gene_name", "hgnc_symbol")), None)
        chosen_gene_col = ci_match
    if chosen_gene_col is None:
        print(f"WARNING: no recognizable gene-symbol column in {path}; "
              f"available columns: {list(df.columns)[:8]}...",
              file=sys.stderr)
        return set(), pd.DataFrame()

    # Apply FDR / confidence filters BEFORE extracting annotations -- so
    # the annotation table only carries genes that actually pass the filter.
    if fdr_col and fdr_col in df.columns and fdr_max is not None:
        before = len(df)
        df = df[pd.to_numeric(df[fdr_col], errors="coerce") <= fdr_max]
        print(f"  {path}: FDR filter {fdr_col}<={fdr_max} kept "
              f"{len(df)}/{before} rows", file=sys.stderr)

    if confidence_col and confidence_col in df.columns and min_confidence is not None:
        before = len(df)
        df = df[pd.to_numeric(df[confidence_col], errors="coerce") >= min_confidence]
        print(f"  {path}: confidence filter {confidence_col}>={min_confidence} "
              f"kept {len(df)}/{before} rows", file=sys.stderr)

    # Build the gene set
    genes_series = (df[chosen_gene_col]
                    .astype(str)
                    .str.strip())
    genes_series = genes_series[genes_series.str.lower().isin(
        ["nan", "none", ""]) == False]
    genes = set(genes_series.dropna())
    print(f"  loaded {len(genes)} genes from {path} (col: {chosen_gene_col})",
          file=sys.stderr)

    # Build the annotation table
    if not report_columns:
        return genes, pd.DataFrame()

    keep_cols = [c for c in report_columns if c in df.columns]
    missing = [c for c in report_columns if c not in df.columns]
    if missing:
        print(f"  {path}: requested annotation columns not found and "
              f"will be empty in the report: {missing}",
              file=sys.stderr)
    if not keep_cols:
        return genes, pd.DataFrame()

    annot = df[[chosen_gene_col] + keep_cols].copy()
    annot[chosen_gene_col] = annot[chosen_gene_col].astype(str).str.strip()
    # Drop duplicates by gene -- some panels list a gene multiple times
    # (e.g. one row per disease in DDG2P). Keep the first occurrence;
    # for users who want full multiplicity, see the "aggregate=" hook below.
    annot = annot.drop_duplicates(subset=[chosen_gene_col], keep="first")
    annot = annot.set_index(chosen_gene_col)

    # Rename columns with the panel prefix
    annot.columns = [f"{prefix}_{c}" for c in annot.columns]
    return genes, annot


def aggregate_panel_rows(path, prefix, gene_col=None, agg_columns=None,
                         joiner="; "):
    """
    Some panels (notably DDG2P / PanelApp) list a gene many times -- once per
    associated disorder -- and useful information is the *concatenation* of
    all rows. This helper produces a per-gene table where the requested
    columns are concatenated with `joiner` across all rows for the gene.
    Used as an alternative to load_gene_annotations when first-row semantics
    aren't right (e.g., to keep all phenotype names per gene).

    Returns a DataFrame indexed on gene symbol with columns named
    `<prefix>_<col>`. Empty if file missing or columns absent.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    sep = "," if path.endswith(".csv") else "\t"
    try:
        df = pd.read_csv(path, sep=sep, low_memory=False)
    except Exception:
        return pd.DataFrame()
    if df.empty or not agg_columns:
        return pd.DataFrame()

    # Resolve gene column the same way as the other loaders
    chosen_gene_col = None
    if gene_col and gene_col in df.columns:
        chosen_gene_col = gene_col
    else:
        for cand in ("hgnc_symbol", "gene_symbol", "gene symbol",
                     "gene_name", "gene", "symbol", "Gene", "Symbol"):
            if cand in df.columns:
                chosen_gene_col = cand
                break
    if chosen_gene_col is None:
        return pd.DataFrame()

    keep = [c for c in agg_columns if c in df.columns]
    if not keep:
        return pd.DataFrame()

    df = df[[chosen_gene_col] + keep].copy()
    df[chosen_gene_col] = df[chosen_gene_col].astype(str).str.strip()
    df = df[df[chosen_gene_col].str.lower().isin(["", "nan", "none"]) == False]

    def _concat_unique(series):
        seen = []
        for v in series.astype(str):
            v = v.strip()
            if v and v.lower() not in ("nan", "none") and v not in seen:
                seen.append(v)
        return joiner.join(seen)

    aggregated = df.groupby(chosen_gene_col).agg({c: _concat_unique for c in keep})
    aggregated.columns = [f"{prefix}_{c}" for c in aggregated.columns]
    return aggregated


def load_constrained_genes(path, threshold=LOEUF_THRESHOLD):
    """Load constrained genes from gnomAD constraint table; filter by LOEUF."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"WARNING: {path} not found; no constraint filter applied",
              file=sys.stderr)
        return None  # sentinel: no filter
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as err:
        print(f"WARNING: could not read {path}: {err}", file=sys.stderr)
        return None

    loeuf_col = (next((c for c in df.columns if c == "lof.oe_ci.upper"), None)
                 or next((c for c in df.columns if "loeuf" in c.lower()), None)
                 or next((c for c in df.columns if "oe_ci.upper" in c.lower()), None))
    if loeuf_col is None:
        print(f"WARNING: no LOEUF column in {path}", file=sys.stderr)
        return None

    gene_col = next((c for c in df.columns
                     if c.lower() in ("gene", "symbol", "gene_symbol",
                                      "hgnc_symbol")), "gene")
    df[loeuf_col] = pd.to_numeric(df[loeuf_col], errors="coerce")
    constrained = set(df.loc[df[loeuf_col] < threshold, gene_col]
                        .astype(str).str.strip())
    constrained.discard("")
    constrained.discard("nan")
    print(f"  loaded {len(constrained)} constrained genes "
          f"(LOEUF < {threshold}) from {path}", file=sys.stderr)
    return constrained


def load_gnomad_constraint_per_gene(path):
    """Load per-gene pLI and LOEUF from the gnomAD constraint table.

    dbNSFP 5.3.1a does NOT include gene-level constraint metrics (pLI / LOEUF
    are gene-level, while dbNSFP is per-variant), so they have to come from
    gnomAD directly. Returns a DataFrame with columns:
        SYMBOL          -- gene symbol (HGNC), match key
        gnomADv4_pLI    -- probability of LoF intolerance
        gnomADv4_LOEUF  -- loss-of-function observed/expected upper-bound
    The "gnomADv4_" prefix mirrors the existing prefix used for gnomAD v4 AF
    columns (from --custom file=...,short_name=gnomADv4) so both groups of
    gnomAD-sourced columns share a single visually-obvious source label.

    If the gnomAD constraint file has multiple rows per gene (one per
    transcript), the row flagged as MANE-select or canonical is kept. If
    neither flag column is present, the first row per gene wins.

    The column-name detection accepts both the v2.1.1 convention ("pLI",
    "oe_lof_upper") and the v4.x convention ("lof.pLI", "lof.oe_ci.upper")
    so the function survives schema variation across gnomAD releases.
    Returns None if the file is missing or has neither pLI nor LOEUF columns
    (the caller then skips the merge with a warning).
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"WARNING: {path} not found; gnomADv4_pLI/LOEUF columns "
              "will be absent from output", file=sys.stderr)
        return None
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as err:
        print(f"WARNING: could not read {path}: {err}; gnomADv4_pLI/LOEUF "
              "columns will be absent from output", file=sys.stderr)
        return None

    # Detect column names: accept gnomAD v2.1.1 or v4.x conventions
    gene_col = next((c for c in df.columns
                     if c.lower() in ("gene", "symbol", "gene_symbol",
                                      "hgnc_symbol")), None)
    pli_col = (next((c for c in df.columns if c == "lof.pLI"), None)
               or next((c for c in df.columns if c.lower() == "pli"), None))
    loeuf_col = (next((c for c in df.columns if c == "lof.oe_ci.upper"), None)
                 or next((c for c in df.columns if c == "oe_lof_upper"), None)
                 or next((c for c in df.columns if "loeuf" in c.lower()), None))

    if gene_col is None:
        print(f"WARNING: no gene-symbol column in {path}; "
              "gnomADv4_pLI/LOEUF columns will be absent", file=sys.stderr)
        return None
    if pli_col is None and loeuf_col is None:
        print(f"WARNING: neither pLI nor LOEUF columns found in {path}; "
              "gnomADv4_pLI/LOEUF columns will be absent", file=sys.stderr)
        return None

    # Deduplicate rows per gene. Prefer mane_select == True, then canonical
    # == True, then keep the first row per gene as a fallback. gnomAD's
    # canonical / mane_select are boolean strings ("true"/"false") in the
    # TSV, so cast first.
    for flag in ("mane_select", "canonical"):
        if flag in df.columns:
            df[flag] = df[flag].astype(str).str.lower().eq("true")
    if "mane_select" in df.columns and df["mane_select"].any():
        df = df[df["mane_select"]].copy()
    elif "canonical" in df.columns and df["canonical"].any():
        df = df[df["canonical"]].copy()
    df = df.drop_duplicates(subset=gene_col, keep="first")

    # Build minimal output table
    keep = [gene_col] + [c for c in (pli_col, loeuf_col) if c is not None]
    out = df[keep].copy()
    rename = {gene_col: "SYMBOL"}
    if pli_col is not None:
        rename[pli_col]   = "gnomADv4_pLI"
        out[pli_col]      = pd.to_numeric(out[pli_col],   errors="coerce")
    if loeuf_col is not None:
        rename[loeuf_col] = "gnomADv4_LOEUF"
        out[loeuf_col]    = pd.to_numeric(out[loeuf_col], errors="coerce")
    out = out.rename(columns=rename)
    out["SYMBOL"] = out["SYMBOL"].astype(str).str.strip()
    out = out[out["SYMBOL"].ne("") & out["SYMBOL"].ne("nan")]

    print(f"  loaded gnomAD constraint per-gene table: {len(out)} genes "
          f"({'pLI ' if pli_col else ''}{'LOEUF' if loeuf_col else ''}) "
          f"from {path}", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--vep_tsv", required=True)
    ap.add_argument("--sample_vcf", required=True)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--schema_genes", required=True)
    ap.add_argument("--bipex_genes", required=True)
    ap.add_argument("--asc_genes", required=True)
    ap.add_argument("--ndd_genes", required=True)
    ap.add_argument("--loeuf_table", required=True)
    ap.add_argument("--artifact_sites", default=None,
                    help="TSV of (CHROM, POS, REF, ALT) tuples to exclude as "
                         "internal cohort artifacts (recurrent in unrelated "
                         "samples + absent/rare in gnomAD). Optional; if "
                         "omitted, no internal artifact filtering is applied.")
    ap.add_argument("--out_tier_a", required=True)
    ap.add_argument("--out_tier_b", required=True)
    ap.add_argument("--out_tier_c", required=True)
    ap.add_argument("--out_master", required=True)
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # 1. Load VEP annotations
    # ------------------------------------------------------------------
    print("Loading VEP TSV...", file=sys.stderr)
    vep = parse_vep_tsv(args.vep_tsv)
    print(f"  {len(vep)} VEP rows, {len(vep.columns)} columns", file=sys.stderr)
    vep = normalize_uploaded_variation(vep)
    vep = split_location(vep)
    vep["SpliceAI_max"] = vep.apply(max_spliceai, axis=1)

    # Detect the AlphaMissense, REVEL, LOFTEE column variants the VEP run
    # actually produced
    am_score_col, am_class_col = detect_alphamissense_cols(vep.columns)
    revel_col                  = detect_revel_col(vep.columns)
    loftee_present             = detect_loftee_present(vep.columns)
    af_cols                    = detect_af_cols(vep.columns)

    print(f"  AlphaMissense: score={am_score_col}, class={am_class_col}",
          file=sys.stderr)
    print(f"  REVEL: {revel_col}", file=sys.stderr)
    print(f"  LOFTEE present: {loftee_present}", file=sys.stderr)
    print(f"  AF columns ({len(af_cols)}): {af_cols[:5]}{'...' if len(af_cols)>5 else ''}",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # 2. Per-sample genotypes; keep only this sample's ALT-carrying calls
    # ------------------------------------------------------------------
    print(f"Extracting genotypes for {args.sample_id}...", file=sys.stderr)
    gt = extract_genotypes_for_sample(args.sample_vcf, args.sample_id)
    if gt.empty:
        print(f"WARNING: no genotype rows from {args.sample_vcf}",
              file=sys.stderr)
    else:
        gt = gt[gt["GT"].isin(["0/1", "1/0", "1/1",
                               "0|1", "1|0", "1|1"])].copy()
    print(f"  {len(gt)} ALT-carrying genotype rows", file=sys.stderr)

    merged = vep.merge(gt, on=["CHROM", "POS", "REF", "ALT"], how="inner")
    print(f"  {len(merged)} variants merged with VEP annotations",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # 2b. Compute BAF (B-allele frequency) = AD[ALT] / DP
    # ------------------------------------------------------------------
    # AD comes from GATK as "REF_count,ALT_count" (or longer for multi-allelic
    # sites; since the upstream pipeline normalizes/splits multi-allelics in
    # r09, every row here is biallelic and AD is exactly two integers).
    # BAF is the standard QC / allele-balance metric for a heterozygous call:
    #   - ~0.5 for a true heterozygous diploid germline variant
    #   - skewed toward 0 or 1 for mosaicism, allelic dropout, mapping bias
    # NaN for missing AD/DP (rows from samples with no coverage at the site).
    def _baf(ad_str, dp_val):
        try:
            parts = str(ad_str).split(",")
            if len(parts) < 2:
                return np.nan
            alt_count = float(parts[1])
            dp = float(dp_val)
            return alt_count / dp if dp > 0 else np.nan
        except (TypeError, ValueError):
            return np.nan

    merged["BAF"] = [
        _baf(ad, dp) for ad, dp in zip(merged.get("AD", []),
                                       merged.get("DP", []))
    ]

    # ------------------------------------------------------------------
    # 3. gnomAD AF filter
    # ------------------------------------------------------------------
    if af_cols:
        merged["gnomAD_popmax"] = merged[af_cols].apply(
            lambda r: pd.to_numeric(r, errors="coerce").max(), axis=1)
    else:
        print("WARNING: no gnomAD AF columns; AF filter skipped",
              file=sys.stderr)
        merged["gnomAD_popmax"] = np.nan

    # NaN AF means "absent from gnomAD" -> rare by default
    rare = merged[merged["gnomAD_popmax"].isna() |
                  (merged["gnomAD_popmax"] < AF_MAX)].copy()
    print(f"  {len(rare)} rare variants (popmax AF < {AF_MAX} or absent)",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # 3b. Internal artifact filter (optional; relatedness-aware)
    # ------------------------------------------------------------------
    # The pipeline's r09e_artifact_blacklist rule produces a list of
    # (CHROM, POS, REF, ALT) tuples that recurred across unrelated
    # representatives but are absent or ultra-rare in gnomAD - i.e.
    # batch / mapping / sample-prep artifacts. Drop them.
    if args.artifact_sites and os.path.exists(args.artifact_sites):
        try:
            artifacts = pd.read_csv(args.artifact_sites, sep="\t",
                                    header=None,
                                    names=["CHROM", "POS", "REF", "ALT"],
                                    dtype={"CHROM": str, "POS": int,
                                           "REF": str, "ALT": str})
        except pd.errors.EmptyDataError:
            artifacts = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT"])

        if len(artifacts) > 0:
            artifact_keys = set(zip(artifacts["CHROM"],
                                    artifacts["POS"].astype(int),
                                    artifacts["REF"],
                                    artifacts["ALT"]))
            before = len(rare)
            rare_pos_int = rare["POS"].astype("Int64").astype(int, errors="ignore")
            mask_artifact = [
                (c, int(p) if pd.notna(p) else -1, r, a) in artifact_keys
                for c, p, r, a in zip(rare["CHROM"], rare["POS"],
                                      rare["REF"], rare["ALT"])
            ]
            rare = rare[~pd.Series(mask_artifact, index=rare.index)].copy()
            print(f"  internal artifact filter: dropped "
                  f"{before - len(rare)} of {before} rare variants "
                  f"({len(artifacts)} blacklist entries)",
                  file=sys.stderr)
        else:
            print(f"  internal artifact filter: blacklist is empty, skipping",
                  file=sys.stderr)
    else:
        print(f"  internal artifact filter: not requested (--artifact_sites unset)",
              file=sys.stderr)

    # ------------------------------------------------------------------
    # 4. Gene sets (with per-panel annotations attached to each variant)
    # ------------------------------------------------------------------
    # Each panel contributes:
    #   - a gene set used for Tier C membership
    #   - a small block of per-gene columns added to the variant report,
    #     prefixed with the panel name to avoid column collisions.
    # Panels with files that don't exist (or are empty placeholders) just
    # contribute empty annotations -- the script keeps running.
    print("Loading gene sets...", file=sys.stderr)

    # SCHEMA: meta-analysis statistics. Q meta is the FDR; OR (PTV) and
    # OR (Class I) are the protein-truncating and damaging-missense odds
    # ratios respectively. P meta is redundant with Q meta but is included
    # because it's often the value reported in publications.
    schema_set, schema_annot = load_gene_annotations(
        args.schema_genes, prefix="SCHEMA",
        gene_col="hgnc_symbol",
        fdr_col="Q meta", fdr_max=SCHEMA_FDR_MAX,
        report_columns=["Q meta", "P meta", "OR (PTV)", "OR (Class I)",
                        "Case PTV", "Ctrl PTV"],
    )

    # BipEx and ASC: the public exome-browser TSVs use slightly different
    # column conventions across releases. We list several common candidates
    # and let load_gene_annotations skip the ones that don't exist in the
    # actual file. Empty placeholder files yield empty annotations.
    # BipEx and ASC TSVs from the Broad exome browsers are keyed only on
    # Ensembl gene_id, so they have been pre-joined with HGNC symbols via
    # join_bipex_asc_with_hgnc.sh. The pre-join also collapses multiple-
    # group rows to one row per gene (BipEx: "Bipolar Disorder" group;
    # ASC: "All" group) and surfaces the columns we actually use.
    
    # BipEx has no broad published FDR gene panel in this browser table.
    # Keep the full table for annotations, but do not use all tested genes
    # as Tier C prior genes.
    bipex_set, bipex_annot = load_gene_annotations(
        args.bipex_genes, prefix="BipEx",
        gene_col="hgnc_symbol",
        report_columns=[
            "ptv_case_count", "ptv_control_count",
            "ptv_fisher_gnom_non_psych_pval", "ptv_fisher_gnom_non_psych_OR",
            "damaging_missense_case_count", "damaging_missense_control_count",
            "damaging_missense_fisher_gnom_non_psych_pval",
            "damaging_missense_fisher_gnom_non_psych_OR",
        ],
    )

    # Empty set means BipEx does not contribute to Tier C.
    bipex_set = set()

    asc_set, asc_annot = load_gene_annotations(
        args.asc_genes, prefix="ASC",
        gene_col="hgnc_symbol",
        fdr_col="qval", fdr_max=ASC_FDR_MAX, # filter fdr
        report_columns=[
            "xcase_dn_ptv", "xcont_dn_ptv",
            "xcase_dn_misb", "xcont_dn_misb",  # damaging missense de novo
            "xcase_dn_misa", "xcont_dn_misa", # other missense de novo
            "qval",
        ],
    )

    # DDG2P (PanelApp): the most clinically meaningful columns are the
    # confidence level and the mode of inheritance. Phenotypes are listed
    # for context. Multiple-disease genes are aggregated across rows so a
    # gene with three different developmental disorders shows all three.
    ndd_set, _ndd_first = load_gene_annotations(
        args.ndd_genes, prefix="DDG2P",
        gene_col="gene_symbol",
        confidence_col="confidence_level",
        min_confidence=2,  # green + amber
        report_columns=["confidence_level", "mode_of_inheritance"],
    )
    # Phenotypes need concatenation (one gene -> many disorders), not
    # first-row semantics
    ndd_phen = aggregate_panel_rows(
        args.ndd_genes, prefix="DDG2P",
        gene_col="gene_symbol",
        agg_columns=["phenotypes"],
    )
    # Combine: confidence/MOI from first row, phenotypes aggregated
    if not _ndd_first.empty and not ndd_phen.empty:
        ndd_annot = _ndd_first.join(ndd_phen, how="outer")
    elif not _ndd_first.empty:
        ndd_annot = _ndd_first
    elif not ndd_phen.empty:
        ndd_annot = ndd_phen
    else:
        ndd_annot = pd.DataFrame()

    panel_set = schema_set | bipex_set | asc_set | ndd_set
    print(f"  combined panel gene set: {len(panel_set)} genes",
          file=sys.stderr)

    constrained_set = load_constrained_genes(args.loeuf_table)
    if constrained_set is None:
        # Fall back to "any panel gene" if no constraint table available
        print("WARNING: using panel set as constrained-gene fallback",
              file=sys.stderr)
        constrained_set = panel_set if panel_set else set()

    # ------------------------------------------------------------------
    # 5. Tier classification
    # ------------------------------------------------------------------
    print("Classifying variants...", file=sys.stderr)
    rare["is_lof"] = rare.apply(is_lof, axis=1, require_hc=loftee_present)
    rare["is_splice"] = rare.apply(is_splice_hit, axis=1)
    rare["is_damaging_mis"] = rare.apply(
        is_damaging_missense, axis=1,
        constrained_genes=constrained_set,
        am_score_col=am_score_col,
        am_class_col=am_class_col,
        revel_col=revel_col,
    )
    rare["is_panel_gene"] = rare["SYMBOL"].isin(panel_set)
    rare["is_protein_altering"] = rare.apply(is_protein_altering, axis=1)

    tier_a_mask = rare["is_lof"] | rare["is_splice"]
    tier_b_mask = rare["is_damaging_mis"] & ~tier_a_mask
    tier_c_mask = (rare["is_panel_gene"] & rare["is_protein_altering"]
                   & ~tier_a_mask & ~tier_b_mask)

    rare["tier"] = np.select(
        [tier_a_mask, tier_b_mask, tier_c_mask],
        ["A", "B", "C"],
        default="",
    )

    # Counts (computed from masks; per-tier DataFrames are sliced after the
    # panel annotations are merged and the columns are reordered, see below)
    n_tier_a = int(tier_a_mask.sum())
    n_tier_b = int(tier_b_mask.sum())
    n_tier_c = int(tier_c_mask.sum())

    # Per-panel boolean membership flags -- handy for filtering reports
    rare["in_SCHEMA"] = rare["SYMBOL"].isin(schema_set)

    # Annotation-only flag: gene is present in the BipEx browser/result table.
    rare["has_BipEx_annotation"] = rare["SYMBOL"].isin(bipex_set)

    rare["in_ASC"]   = rare["SYMBOL"].isin(asc_set)
    rare["in_DDG2P"] = rare["SYMBOL"].isin(ndd_set)

    # Attach per-panel annotation columns (panel-prefixed) by joining each
    # annotation table on SYMBOL. Genes not in a given panel get NaN for
    # that panel's columns, which is the correct semantics.
    panel_annotations = []
    for name, annot in (("SCHEMA", schema_annot),
                        ("BipEx",  bipex_annot),
                        ("ASC",    asc_annot),
                        ("DDG2P",  ndd_annot)):
        if not annot.empty:
            rare = rare.merge(annot, how="left", left_on="SYMBOL",
                              right_index=True)
            panel_annotations.extend(annot.columns.tolist())
            print(f"  attached {len(annot.columns)} {name} annotation column(s)",
                  file=sys.stderr)

    print(f"  Tier A (LoF + splice):       {n_tier_a} variants",
          file=sys.stderr)
    print(f"  Tier B (damaging missense):  {n_tier_b} variants",
          file=sys.stderr)
    print(f"  Tier C (panel gene):         {n_tier_c} variants",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # 5b. Join gene-level gnomAD constraint (pLI + LOEUF)
    # ------------------------------------------------------------------
    # dbNSFP 5.3.1a does NOT include gene-level constraint metrics (these
    # are gene-level, while dbNSFP is per-variant), so we add them now by
    # joining the dedicated gnomAD constraint TSV by HGNC gene symbol.
    # The constraint table was already opened once by load_constrained_genes
    # above; here we re-read it as a per-gene table with the two columns we
    # want.
    constraint_per_gene = load_gnomad_constraint_per_gene(args.loeuf_table)
    if constraint_per_gene is not None:
        before_cols = set(rare.columns)
        rare = rare.merge(constraint_per_gene, on="SYMBOL", how="left")
        added = sorted(set(rare.columns) - before_cols)
        if added:
            n_pli   = (rare["gnomADv4_pLI"].notna().sum()
                       if "gnomADv4_pLI"   in rare.columns else 0)
            n_loeuf = (rare["gnomADv4_LOEUF"].notna().sum()
                       if "gnomADv4_LOEUF" in rare.columns else 0)
            print(f"  added gnomAD constraint columns ({added}): "
                  f"pLI populated in {n_pli}/{len(rare)} rows, "
                  f"LOEUF populated in {n_loeuf}/{len(rare)} rows",
                  file=sys.stderr)

    # ------------------------------------------------------------------
    # 6. Source-attributed column renaming + 5-group output ordering
    # ------------------------------------------------------------------
    # Two operations, in this order:
    #
    # (a) Rename plugin-sourced columns with their source as a prefix
    #     (LOFTEE_*, AlphaMissense_*, SpliceAI_*, dbNSFP_*, VEP_*), so the
    #     output table makes clear which tool produced each annotation.
    #     This matters for license attribution, version tracking, and for
    #     anyone reading the table without access to the VEP invocation.
    #
    # (b) Group columns into 5 ordered sections in the final report:
    #         1. General info (variant ID + sample genotype data)
    #         2. Tier A criteria (LoF + splice -- DROP-testable)
    #         3. Tier B criteria (missense pathogenicity + gene constraint)
    #         4. Tier C criteria (gene panels + ClinVar)
    #         5. Others (population frequencies, misc VEP flags, tier label)
    #
    # Rename happens AFTER classification (lines 849-870) -- the classifier
    # functions use the original VEP plugin column names; we only rename
    # for the final output table.

    # Build the rename map from columns that actually exist in the table.
    # Missing columns are silently skipped, so this works whether or not
    # each plugin was run.
    rename_map = {}

    # LOFTEE plugin
    for old, new in [("LoF",        "LOFTEE_LoF"),
                     ("LoF_filter", "LOFTEE_LoF_filter"),
                     ("LoF_flags",  "LOFTEE_LoF_flags"),
                     ("LoF_info",   "LOFTEE_LoF_info")]:
        if old in rare.columns:
            rename_map[old] = new

    # AlphaMissense plugin (handles either of the two naming conventions
    # that the VEP plugin uses depending on version)
    for old, new in [("am_pathogenicity",          "AlphaMissense_pathogenicity"),
                     ("am_class",                  "AlphaMissense_class"),
                     ("AlphaMissense_pathogenicity", "AlphaMissense_pathogenicity"),
                     ("AlphaMissense_class",        "AlphaMissense_class"),
                     ("AlphaMissense_score",        "AlphaMissense_pathogenicity")]:
        if old in rare.columns and old != new:
            rename_map[old] = new

    # REVEL plugin
    if "REVEL" in rare.columns:
        rename_map["REVEL"] = "REVEL_score"

    # SpliceAI plugin -- with split_output=1 the four delta-scores and four
    # delta-positions are emitted as separate columns. Drop the "_pred_"
    # infix for compactness: SpliceAI_pred_DS_AG -> SpliceAI_DS_AG.
    for kind in ("DS", "DP"):
        for direction in ("AG", "AL", "DG", "DL"):
            old = f"SpliceAI_pred_{kind}_{direction}"
            new = f"SpliceAI_{kind}_{direction}"
            if old in rare.columns:
                rename_map[old] = new
    # The SYMBOL column SpliceAI emits is the gene it scored against; rename
    # under the same convention so it groups with the other SpliceAI columns.
    if "SpliceAI_pred_SYMBOL" in rare.columns:
        rename_map["SpliceAI_pred_SYMBOL"] = "SpliceAI_SYMBOL"

    # VEP-cache built-in predictors (SIFT / PolyPhen-2 come from --sift b
    # --polyphen b, NOT from dbNSFP, so prefix them with VEP_ to make the
    # source obvious. Without the prefix, readers commonly assume these
    # are dbNSFP-sourced and treat them as independent evidence -- they
    # are not, since dbNSFP's other predictors already incorporate them.)
    for old, new in [("SIFT",     "VEP_SIFT"),
                     ("PolyPhen", "VEP_PolyPhen")]:
        if old in rare.columns:
            rename_map[old] = new

    # dbNSFP plugin predictors. Pandas converts hyphens in column names
    # to dots when the VEP TSV is read, so "M-CAP_pred" can become
    # "M.CAP_pred" -- handle both spellings.
    # NOTE: gnomAD_pLI / gnomAD_LOEUF used to be requested here too, but
    # dbNSFP 5.3.1a does NOT include those gene-level metrics (they are
    # gene-level, while dbNSFP is per-variant). They now come from the
    # gnomAD constraint TSV directly, merged in section 5b above, and
    # land in the output as gnomADv4_pLI / gnomADv4_LOEUF.
    dbnsfp_predictors = [
        "MutationTaster_pred", "PROVEAN_pred",
        "MetaLR_pred", "MetaRNN_pred",
        "M-CAP_pred", "M.CAP_pred",
        "PrimateAI_pred", "ClinPred_pred",
        "BayesDel_addAF_pred",
    ]
    for old in dbnsfp_predictors:
        if old in rare.columns:
            rename_map[old] = f"dbNSFP_{old}"

    rare = rare.rename(columns=rename_map)
    if rename_map:
        print(f"  renamed {len(rename_map)} source-attributed columns",
              file=sys.stderr)

    # ------------------------------------------------------------------
    # Column-ordering: 5 groups
    # ------------------------------------------------------------------
    # Each list is a candidate set -- we filter to columns that actually
    # exist in the table, then concatenate the groups in order. Anything
    # not assigned to a group lands at the end of group 5 ("others").

    # Group 1: General info -- variant identification, genotype data
    # SOMATIC and PHENO are VEP --check_existing output flags (whether the
    # variant has been found in somatic databases like COSMIC, and whether
    # it has any phenotype associations from VEP's existing-variation
    # lookup). They are variant-identification context rather than tier
    # criteria; group them with the other VEP identifier columns.
    group1 = [
        "SYMBOL", "Consequence", "IMPACT",
        "CHROM", "POS", "REF", "ALT",
        "Uploaded_variation", "Location", "Allele",
        "Existing_variation", "SOMATIC", "PHENO",
        "Gene", "Feature", "Feature_type", "BIOTYPE",
        "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SOURCE",
        "SYMBOL_SOURCE", "HGNC_ID",
        "STRAND", "FLAGS", "DISTANCE",
        "cDNA_position", "CDS_position", "Protein_position",
        "Amino_acids", "Codons",
        "HGVSc", "HGVSp", "HGVS_OFFSET", "EXON", "INTRON",
        "GT", "DP", "AD", "BAF", "GQ",
    ]

    # Group 2: Tier A criteria -- LoF + splice (DROP-testable)
    group2 = [
        "is_lof", "is_splice",
        # LOFTEE outputs (renamed)
        "LOFTEE_LoF", "LOFTEE_LoF_filter", "LOFTEE_LoF_flags", "LOFTEE_LoF_info",
        # SpliceAI: max aggregate, then individual delta-scores, then
        # delta-positions, then the SpliceAI-internal gene symbol used for
        # prediction (renamed from SpliceAI_pred_SYMBOL), then the legacy
        # packed-string column if present.
        "SpliceAI_max",
        "SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL",
        "SpliceAI_DP_AG", "SpliceAI_DP_AL", "SpliceAI_DP_DG", "SpliceAI_DP_DL",
        "SpliceAI_SYMBOL",
        "SpliceAI_pred",
    ]

    # Group 3: Tier B criteria -- missense pathogenicity + gene constraint
    group3 = [
        "is_damaging_mis", "is_protein_altering",
        # AlphaMissense
        "AlphaMissense_pathogenicity", "AlphaMissense_class",
        # Other meta-predictors
        "REVEL_score",
        # CADD
        "CADD_PHRED", "CADD_RAW",
        # VEP-cache predictors
        "VEP_SIFT", "VEP_PolyPhen",
        # dbNSFP predictors (renamed)
        "dbNSFP_MutationTaster_pred", "dbNSFP_PROVEAN_pred",
        "dbNSFP_MetaLR_pred", "dbNSFP_MetaRNN_pred",
        "dbNSFP_M-CAP_pred", "dbNSFP_M.CAP_pred",
        "dbNSFP_PrimateAI_pred", "dbNSFP_ClinPred_pred",
        "dbNSFP_BayesDel_addAF_pred",
        # Gene-level constraint from gnomAD v4.1 (joined directly from the
        # gnomAD constraint TSV in section 5b -- NOT from dbNSFP, which
        # does not include these gene-level metrics).
        "gnomADv4_pLI", "gnomADv4_LOEUF",
    ]

    # Group 4: Tier C criteria -- gene panels + ClinVar
    group4 = [
        "is_panel_gene",
        # Per-panel boolean membership flags
        "in_SCHEMA", "has_BipEx_annotation", "in_ASC", "in_DDG2P",
    ]
    # Per-panel annotation columns (already panel-prefixed). Add them in
    # panel order so columns from the same panel stay grouped together.
    for prefix in ("SCHEMA", "BipEx", "ASC", "DDG2P"):
        group4 += [c for c in rare.columns
                   if c.startswith(prefix + "_") and c not in group4]
    # ClinVar fields
    group4 += [
        "ClinVar",
        "ClinVar_CLNSIG", "ClinVar_CLNDN",
        "ClinVar_CLNREVSTAT", "ClinVar_CLNDISDB",
        # VEP's own existing-variation flag (different from ClinVar plugin)
        "CLIN_SIG",
    ]

    # Group 5: Others -- population frequencies, misc VEP flags, tier
    # (SOMATIC and PHENO have moved up to group1 alongside the other
    # VEP existing-variation columns.)
    group5 = ["gnomAD_popmax"]
    group5 += [c for c in rare.columns if c.startswith("gnomADv4")
               and c not in ("gnomADv4_pLI", "gnomADv4_LOEUF")]
    group5 += ["tier"]  # final tier label goes last

    # Concatenate groups, keeping only columns that exist; preserve order
    ordered = []
    seen = set()
    for group in (group1, group2, group3, group4, group5):
        for c in group:
            if c in rare.columns and c not in seen:
                ordered.append(c)
                seen.add(c)

    # Anything not explicitly placed in a group lands at the end (before
    # 'tier', so 'tier' stays last). Sorted for determinism.
    leftover = sorted(c for c in rare.columns if c not in seen)
    if leftover:
        print(f"  {len(leftover)} unassigned columns appended at the end: "
              f"{leftover[:5]}{'...' if len(leftover) > 5 else ''}",
              file=sys.stderr)
    # Pull 'tier' out of `ordered` (we know it's there) and re-append after
    # the leftovers, so it stays as the final column.
    if "tier" in ordered:
        ordered.remove("tier")
    ordered = ordered + leftover + (["tier"] if "tier" in rare.columns else [])

    rare = rare[ordered]

    # Slice per-tier views from the reordered, fully-annotated `rare` so the
    # tier-specific TSVs share the same columns and column order as the
    # master file.
    tier_a = rare[tier_a_mask].copy()
    tier_b = rare[tier_b_mask].copy()
    tier_c = rare[tier_c_mask].copy()

    # ------------------------------------------------------------------
    # 7. Write outputs
    # ------------------------------------------------------------------
    def _write(df, path, drop_tier=True):
        if drop_tier and "tier" in df.columns:
            df = df.drop(columns=["tier"])
        df.to_csv(path, sep="\t", index=False)
        print(f"  wrote {len(df):>6} rows -> {path}", file=sys.stderr)

    print("Writing outputs...", file=sys.stderr)
    _write(tier_a, args.out_tier_a)
    _write(tier_b, args.out_tier_b)
    _write(tier_c, args.out_tier_c)
    _write(rare,   args.out_master, drop_tier=False)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
