#!/usr/bin/env python3
"""Shared logic for the deconvolution-reference builder scripts.

Every reference script in scripts_to_make_deconv_reference/ imports this module
and provides only the parts that differ between references: the source URLs,
the obs filters, and a function that turns obs into a per-cell final CLASS label
(the string "DROP" excludes a cell). Everything else -- backed reading, class
balancing, pulling the selected rows, and writing the canonical hspe contract --
lives here so a fix lands once.

CANONICAL OUTPUT CONTRACT (identical for every reference; consumed by hspe):
  <outdir>/
    matrix.mtx.gz      genes x cells, raw integer counts (MatrixMarket)
    features.tsv.gz     gene symbol per row
    cells.tsv.gz        cell_id \t subject_id \t source_label \t class
    provenance.json     inputs, sha256, filters, seed, per-class counts

Design invariants (learned the hard way in this project):
  * read h5ad BACKED -- never densify (full atlas densified ~= 750 GB).
  * raw COUNTS only (hspe needs linear counts); warn if the layer is non-integer.
  * every source label must map to a class or be explicit DROP -- a silently
    dropped label is how a reference quietly stops matching the biology.
  * fail loud: unknown columns / unmapped labels / empty selection are errors.
"""
import gzip, hashlib, json, os, sys, csv, urllib.request
from datetime import datetime, timezone


def eprint(*a, **k):
    print(*a, file=sys.stderr, flush=True, **k)


# --------------------------------------------------------------------------- #
# download management (stdlib only; scripts manage their own sources)
# --------------------------------------------------------------------------- #
def fetch(url, target, expect_min_bytes=1 << 20):
    """Download url -> target if absent. Atomic (.partial + rename), skip if a
    non-empty file already sits at target. Streams, so multi-GB is fine."""
    if os.path.isfile(target) and os.path.getsize(target) > 0:
        eprint(f"[fetch] present, skip: {target} "
               f"({os.path.getsize(target)/1e9:.1f} GB)")
        return target
    if not url:
        sys.exit(f"[ERR] need to download {os.path.basename(target)} but no URL "
                 "given (set the appropriate *_URL env var; see script header).")
    os.makedirs(os.path.dirname(target), exist_ok=True)
    tmp = target + ".partial"
    eprint(f"[fetch] {url}\n        -> {target}")
    try:
        with urllib.request.urlopen(url) as r, open(tmp, "wb") as out:
            total = int(r.headers.get("Content-Length", 0))
            got = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
                got += len(chunk)
                if total:
                    eprint(f"\r[fetch]   {got/1e9:5.1f} / {total/1e9:5.1f} GB",
                           end="")
        eprint("")
    except Exception as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        sys.exit(f"[ERR] download failed: {e}")
    if os.path.getsize(tmp) < expect_min_bytes:
        os.remove(tmp)
        sys.exit(f"[ERR] {target} suspiciously small; refusing.")
    os.replace(tmp, target)
    return target


CXG_API = "https://api.cellxgene.cziscience.com/curation/v1"


def resolve_cxg_h5ad(collection_id, dataset_id):
    """dataset_id -> current CELLxGENE h5ad asset URL (URL is presigned and
    ephemeral; dataset_id is stable, so we resolve fresh each run)."""
    import json
    req = urllib.request.Request(f"{CXG_API}/collections/{collection_id}",
                                 headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        coll = json.load(r)
    for d in coll.get("datasets", []):
        if (d.get("dataset_id") or d.get("id")) == dataset_id:
            for a in (d.get("assets") or d.get("dataset_assets") or []):
                if str(a.get("filetype", "")).upper() == "H5AD" and a.get("url"):
                    return a["url"]
            sys.exit(f"[ERR] dataset {dataset_id} has no H5AD asset")
    sys.exit(f"[ERR] dataset {dataset_id} not in collection {collection_id}")


def ensure_source(target, url=None, cxg=None):
    """Local-first source resolution for the recipe scripts:
      1. target already staged -> use it, NO network (works behind YC allowlist);
      2. else direct `url` -> download;
      3. else `cxg`=(collection_id, dataset_id) -> resolve fresh URL + download.
    """
    if os.path.isfile(target) and os.path.getsize(target) > 0:
        eprint(f"[source] using staged file: {target} "
               f"({os.path.getsize(target)/1e9:.1f} GB)")
        return target
    if url:
        return fetch(url, target)
    if cxg:
        eprint(f"[source] resolving CELLxGENE dataset {cxg[1]} ...")
        return fetch(resolve_cxg_h5ad(cxg[0], cxg[1]), target)
    sys.exit(f"[ERR] {target} not present and no url/cxg to fetch it. "
             "Stage the file there, or provide a source.")


def sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(buf), b""):
            h.update(c)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# label construction helpers (used by the scripts to build the class Series)
# --------------------------------------------------------------------------- #
def load_mapping(path):
    """TSV with columns source_label, class. class 'DROP' (any case) excludes."""
    m = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if not rd.fieldnames or "source_label" not in rd.fieldnames \
                or "class" not in rd.fieldnames:
            sys.exit(f"[ERR] mapping {path} needs columns: source_label, class")
        for row in rd:
            m[row["source_label"].strip()] = row["class"].strip()
    if not m:
        sys.exit(f"[ERR] mapping {path} is empty")
    return m


def apply_flat_mapping(obs, col, mapping):
    """obs[col] -> class Series via mapping; HARD FAIL on any unmapped value."""
    vals = obs[col].astype(str)
    unmapped = sorted(set(vals) - set(mapping))
    if unmapped:
        sys.exit(f"[ERR] labels in '{col}' absent from mapping "
                 f"(add them or map to DROP): {unmapped}")
    return vals.map(mapping)


def compose_hierarchical_label(obs, coarse_col, fine_col, coarse_rules,
                               coarse_name="Neuron"):
    """Asymmetric annotation: resolve some branches fine, collapse others.

    coarse_rules maps each value of coarse_col to an action:
      'COARSE' -> emit the fixed `coarse_name` (collapse, e.g. all neurons)
      'FINE'   -> emit this cell's fine_col value (resolve, e.g. glial subtypes)
      'DROP'   -> exclude (e.g. non-neuroectoderm lineages)
    Every value of coarse_col must appear in coarse_rules (HARD FAIL otherwise).
    Returns a class Series ('DROP' string for excluded cells).
    """
    import pandas as pd
    for c in (coarse_col, fine_col):
        if c not in obs.columns:
            sys.exit(f"[ERR] obs has no column '{c}'. available: "
                     f"{list(obs.columns)[:30]}")
    cv = obs[coarse_col].astype(str)
    unruled = sorted(set(cv) - set(coarse_rules))
    if unruled:
        sys.exit(f"[ERR] values of '{coarse_col}' without a rule "
                 f"(add COARSE/FINE/DROP): {unruled}")
    fine = obs[fine_col].astype(str)
    out = []
    for c, f in zip(cv, fine):
        act = coarse_rules[c]
        out.append(coarse_name if act == "COARSE"
                   else ("DROP" if act == "DROP" else f))
    return pd.Series(out, index=obs.index)


def cap_by_group(df, col, cap, rng):
    """Index labels of df, capped to `cap` rows per value of col (seeded).
    Version-robust: uses .groupby().groups, avoids apply() semantics that
    changed across pandas versions."""
    import numpy as np
    keep = []
    for _, idx in df.groupby(col).groups.items():
        idx = np.asarray(idx)
        if len(idx) > cap:
            idx = idx[rng.permutation(len(idx))[:cap]]
        keep.append(idx)
    return np.concatenate(keep) if keep else np.array([], dtype=int)


# --------------------------------------------------------------------------- #
# enumerate-first support (for references whose obs schema is unknown)
# --------------------------------------------------------------------------- #
def enumerate_obs(h5ad_path, max_card=60, cross=None):
    """Print obs columns and value counts for low-cardinality columns, plus an
    optional cross-tab of two columns. Reads BACKED (obs only). Writes nothing.
    Used so a reference with an unknown annotation vocabulary can surface it
    before anyone pins label columns to it."""
    import anndata as ad
    A = ad.read_h5ad(h5ad_path, backed="r")
    obs = A.obs
    eprint(f"[enumerate] {h5ad_path}: {obs.shape[0]:,} cells, "
           f"{obs.shape[1]} obs columns\n")
    for c in obs.columns:
        try:
            vc = obs[c].value_counts()
        except Exception:
            continue
        if 1 < len(vc) <= max_card:
            eprint(f"--- {c}  ({len(vc)} values) ---")
            for v, n in vc.items():
                eprint(f"    {n:>10,}  {v}")
            eprint("")
        else:
            eprint(f"--- {c}  ({len(vc)} values; too many to list) ---\n")
    if cross and all(c in obs.columns for c in cross):
        import pandas as pd
        eprint(f"=== cross-tab {cross[0]} (rows) x {cross[1]} (cols) ===")
        ct = pd.crosstab(obs[cross[0]], obs[cross[1]])
        eprint(ct.to_string())


# --------------------------------------------------------------------------- #
# the core builder
# --------------------------------------------------------------------------- #
def build_reference(h5ad_paths, outdir, label_fn, subject_col,
                    tissue_col=None, tissue_terms=None,
                    dissection_col=None, dissection_terms=None,
                    counts_layer="raw",
                    cells_per_class=500, min_cells_per_class=100,
                    balance_col=None, cap_per_balance=2000,
                    seed=42, provenance_extra=None):
    """Two-pass backed build -> canonical contract under outdir.

    label_fn(obs) -> class Series (index-aligned; 'DROP' excludes a cell).
    balance_col: if set, cap per that column BEFORE the per-class cap -- use it
      when several source labels collapse into one class and you want that class
      internally balanced (e.g. neuron/non-neuron: keeps the non-neuronal class
      from being all oligodendrocyte). Leave None when each class is already a
      single source label (e.g. glia resolved to superclusters/clusters).
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import scipy.io as sio

    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    tissue_terms = [t.strip() for t in (tissue_terms or []) if t.strip()]
    diss_terms = [t.strip().lower() for t in (dissection_terms or []) if t.strip()]

    # ---- pass 1: obs only, decide which cells to keep ----
    picks = []
    for fi, path in enumerate(h5ad_paths):
        eprint(f"[obs] backed read: {path}")
        A = ad.read_h5ad(path, backed="r")
        obs = A.obs
        if subject_col not in obs.columns:
            sys.exit(f"[ERR] {path}: no subject column '{subject_col}'. "
                     f"available: {list(obs.columns)[:30]}")
        keep = pd.Series(True, index=obs.index)
        if tissue_col and tissue_terms:
            if tissue_col not in obs.columns:
                sys.exit(f"[ERR] {path}: no tissue column '{tissue_col}'")
            keep &= obs[tissue_col].astype(str).isin(tissue_terms)
        if dissection_col and diss_terms:
            if dissection_col not in obs.columns:
                sys.exit(f"[ERR] {path}: no dissection column '{dissection_col}'")
            d = obs[dissection_col].astype(str).str.lower()
            keep &= d.apply(lambda x: any(t in x for t in diss_terms))
        sub = obs.loc[keep].copy()
        eprint(f"[obs]   kept {len(sub):,}/{len(obs):,} after tissue/dissection")

        cls = label_fn(sub).astype(str)
        sub = sub.assign(_class=cls.values)
        sub = sub[(sub["_class"].str.upper() != "DROP") & sub["_class"].notna()]
        sub["_row"] = np.array([obs.index.get_loc(i) for i in sub.index])

        if balance_col:
            if balance_col not in sub.columns:
                sys.exit(f"[ERR] balance_col '{balance_col}' not in obs")
            sub = sub.loc[cap_by_group(sub, balance_col, cap_per_balance, rng)]
        picks.append((fi, sub))

    allsub = pd.concat([s.assign(_file=fi) for fi, s in picks],
                       ignore_index=False)
    if allsub.empty:
        sys.exit("[ERR] no cells selected -- check filters / label rules")

    # ---- balance per class, floor-check ----
    take = allsub.loc[cap_by_group(allsub, "_class", cells_per_class, rng)]
    kept_counts = {}
    for cls, g in take.groupby("_class"):
        kept_counts[cls] = int(len(g))
        if len(g) < min_cells_per_class:
            eprint(f"[WARN] class '{cls}': only {len(g)} cells (< "
                   f"{min_cells_per_class}); consider merging it upward.")
    eprint(f"[balance] final class counts: {kept_counts}")

    final_rows = {}
    for fi, r in zip(take["_file"], take["_row"]):
        final_rows.setdefault(int(fi), []).append(int(r))

    # ---- pass 2: pull ONLY selected rows, extract counts ----
    mats, cell_ids, subj, srclab, clslab = [], [], [], [], []
    var_ref = None
    for fi, path in enumerate(h5ad_paths):
        rows = sorted(final_rows.get(fi, []))
        if not rows:
            continue
        eprint(f"[pull] {path}: {len(rows):,} nuclei")
        A = ad.read_h5ad(path, backed="r")
        sel = A[rows].to_memory()
        if counts_layer == "raw":
            if sel.raw is None:
                sys.exit(f"[ERR] {path}: --counts-layer raw but .raw is None. "
                         "On CELLxGENE, raw.X is counts and X is log-normalised.")
            X = sel.raw.X
            var = sel.raw.var if sel.raw.var is not None else sel.var
        elif counts_layer == "X":
            X, var = sel.X, sel.var
        else:
            X, var = sel.layers[counts_layer], sel.var
        X = sp.csr_matrix(X)
        d = X.data[:10000]
        if d.size and not np.allclose(d, np.round(d)):
            eprint(f"[WARN] {path}: counts layer '{counts_layer}' is non-integer "
                   "-- looks normalised, not raw counts. Check --counts-layer.")
        mats.append(X)
        o = sel.obs
        cls = label_fn(o).astype(str)
        cell_ids += [f"{fi}:{i}" for i in o.index.astype(str)]
        subj += o[subject_col].astype(str).tolist()
        srclab += cls.tolist()          # the final class doubles as source_label
        clslab += cls.tolist()
        if var_ref is None:
            var_ref = var

    counts = sp.vstack(mats).tocsc()
    sym_col = next((c for c in ("feature_name", "gene_symbol", "Gene", "symbol")
                    if c in var_ref.columns), None)
    genes = (var_ref[sym_col].astype(str).tolist() if sym_col
             else var_ref.index.astype(str).tolist())
    gxc = counts.transpose().tocsc()
    eprint(f"[write] {gxc.shape[0]} genes x {gxc.shape[1]} cells, "
           f"nnz={gxc.nnz:,}")

    # ---- write canonical outputs (atomic) ----
    def gz_lines(path, it):
        tmp = path + ".partial"
        with gzip.open(tmp, "wt") as fh:
            for r in it:
                fh.write(r)
        os.replace(tmp, path)

    mtx = os.path.join(outdir, "matrix.mtx")
    sio.mmwrite(mtx, gxc, field="integer")
    with open(mtx, "rb") as fi, gzip.open(mtx + ".gz", "wb") as fo:
        fo.writelines(fi)
    os.remove(mtx)
    gz_lines(os.path.join(outdir, "features.tsv.gz"), (g + "\n" for g in genes))
    gz_lines(os.path.join(outdir, "cells.tsv.gz"),
             ("\t".join(x) + "\n" for x in zip(cell_ids, subj, srclab, clslab)))

    # ---- per-gene cell-type specificity (from the SAME raw pseudobulk) --------
    # Complete, unfiltered values (unlike CIBERSORTx's signature GEP, which drops
    # sub-population markers like GFAP under a coarse partition). Binary refs get
    # a signed log2 ratio; multi-class refs get tau. Written to specificity.tsv.
    spec_meta = write_specificity(gxc, genes, clslab, outdir)

    prov = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [{"path": p, "sha256": sha256(p)} for p in h5ad_paths],
        "n_cells": int(gxc.shape[1]), "n_genes": int(gxc.shape[0]),
        "n_subjects": int(len(set(subj))), "class_counts": kept_counts,
        "counts_layer": counts_layer, "seed": seed,
        "filters": {"tissue_col": tissue_col, "tissue_terms": tissue_terms,
                    "dissection_col": dissection_col,
                    "dissection_terms": diss_terms},
    }
    if provenance_extra:
        prov.update(provenance_extra)
    prov["specificity"] = spec_meta
    with open(os.path.join(outdir, "provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    eprint(f"[done] {prov['n_cells']} cells, {prov['n_subjects']} subjects, "
           f"{len(kept_counts)} classes -> {outdir}")
    return prov


def write_specificity(gxc, genes, cell_classes, outdir, pseudocount=1.0):
    """Per-gene cell-type specificity from raw per-class pseudobulk.

    Aggregates counts per class (sparse, no densify), CPM-normalises each class
    profile, then:
      * 2 classes  -> signed log2 ratio delta = log2((A+eps)/(B+eps)), where A is
                      the alphabetically-FIRST class. Sign gives the compartment,
                      magnitude the degree. 'higher_in' names the class.
      * >2 classes -> tau = sum(1 - x_i/x_max)/(n-1) on linear CPM (Yanai index):
                      0 = uniform (composition-independent), 1 = single-class
                      (pure composition proxy). 'dominant_class' = argmax.
    Writes <outdir>/specificity.tsv (gene + per-class CPM + score columns).
    Returns a small metadata dict for provenance.
    """
    import numpy as np, scipy.sparse as sp, csv, os
    classes = sorted(set(cell_classes))
    cc = np.asarray(cell_classes)
    # sparse (cells x class) 0/1 indicator -> gxc(genes x cells) %*% ind = per-class sums
    col = {c: j for j, c in enumerate(classes)}
    j_idx = np.fromiter((col[c] for c in cc), dtype=int, count=len(cc))
    ind = sp.csc_matrix((np.ones(len(cc)), (np.arange(len(cc)), j_idx)),
                        shape=(len(cc), len(classes)))
    pb = np.asarray((gxc @ ind).todense())               # genes x classes (summed)
    # CPM per class (column-wise); sum==mean-per-cell in CPM so summing is fine
    colsum = pb.sum(axis=0, keepdims=True)
    colsum[colsum == 0] = 1.0
    cpm = pb / colsum * 1e6

    header = ["gene"] + [f"cpm_{c}" for c in classes]
    rows = []
    if len(classes) == 2:
        method = "log2_ratio"
        a, b = classes[0], classes[1]                    # a = first alphabetically
        A, B = cpm[:, 0], cpm[:, 1]
        delta = np.log2((A + pseudocount) / (B + pseudocount))
        higher = np.where(delta >= 0, a, b)
        header += ["delta_log2", "higher_in", "positive_class", "negative_class"]
        for i, g in enumerate(genes):
            rows.append([g] + [f"{v:.6g}" for v in cpm[i]]
                        + [f"{delta[i]:.6g}", higher[i], a, b])
        meta = {"method": method, "score_column": "delta_log2",
                "convention": f"positive = higher in '{a}', negative = '{b}'",
                "pseudocount": pseudocount, "classes": classes}
    else:
        method = "tau"
        xmax = cpm.max(axis=1)
        n = len(classes)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = cpm / xmax[:, None]                 # x_i / x_max
        tau = np.where(xmax > 0, (n - ratios.sum(axis=1)) / (n - 1), np.nan)
        dom = np.array(classes)[cpm.argmax(axis=1)]
        dom = np.where(xmax > 0, dom, "NA")
        header += ["tau", "dominant_class", "max_cpm"]
        for i, g in enumerate(genes):
            tstr = "NA" if np.isnan(tau[i]) else f"{tau[i]:.6g}"
            rows.append([g] + [f"{v:.6g}" for v in cpm[i]]
                        + [tstr, dom[i], f"{xmax[i]:.6g}"])
        meta = {"method": method, "score_column": "tau",
                "note": "tau in [0,1]; 0=uniform, 1=single-class; NA if gene "
                        "unexpressed in all classes", "classes": classes}

    out = os.path.join(outdir, "specificity.tsv")
    tmp = out + ".partial"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    os.replace(tmp, out)
    eprint(f"[spec] {method}: {len(genes)} genes x {len(classes)} classes "
           f"-> specificity.tsv")
    meta["file"] = "specificity.tsv"
    meta["n_genes"] = len(genes)
    return meta


def already_built(outdir):
    """Idempotency check the scripts use to no-op on re-run. Requires the
    specificity file too, so references built before specificity was added are
    rebuilt to generate it."""
    return os.path.isfile(os.path.join(outdir, "matrix.mtx.gz")) \
        and os.path.isfile(os.path.join(outdir, "provenance.json")) \
        and os.path.isfile(os.path.join(outdir, "specificity.tsv"))
