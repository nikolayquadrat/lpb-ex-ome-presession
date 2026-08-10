#!/usr/bin/env python3
"""
Per-cell-type expression of a gene (default IZUMO4) in a large snRNA-seq atlas.
MEMORY-EFFICIENT: never loads the full matrix.

Why the naive version OOMs: Siletti X is stored CSR (cells-major), so slicing a
single gene COLUMN with anndata (`adata[:, g].X`) touches the whole matrix. Here
we instead stream the gene's column straight from the HDF5 sparse arrays in
row-chunks, reading only the `indices` slice per chunk plus the handful of
matching `data` values -- peak memory is one chunk of indices, not the matrix.

For each cell-type group (supercluster/cluster) it reports mean raw counts,
% of nuclei expressing (dropout-aware), and n cells, ranked. Optional by region.

Usage (run with the deconv venv's python):
    /mnt/data/rnaseq/annotation/deconv/venv/bin/python3 izumo4_celltype_expression.py \
        --h5ad neurons.h5ad nonneurons.h5ad \
        --groupby supercluster_term [--region-col tissue] \
        [--counts raw|X] [--chunk-cells 250000] [--min-cells 50] \
        [--out izumo4_by_celltype.tsv]

Deps: h5py + numpy + pandas (matrix); anndata (obs only, backed -- cheap).
      All present in the deconv venv. tqdm optional (nicer progress).
"""

import os
import math
import argparse
import numpy as np
import pandas as pd

DEFAULT_ENSEMBL = "ENSG00000099840"   # IZUMO4 (used when neither gene flag is given)
DEFAULT_SYMBOL = "IZUMO4"


# ---------------------------------------------------------------------------
def _progress(iterable, total=None, desc=""):
    try:
        from tqdm import tqdm
        return tqdm(iterable, total=total, desc=desc, unit="chunk")
    except ImportError:
        def gen():
            import sys
            n = str(total) if total is not None else "?"
            i = 0
            for x in iterable:
                i += 1
                sys.stderr.write(f"\r  {desc}: {i}/{n}   ")
                sys.stderr.flush()
                yield x
            sys.stderr.write("\n")
        return gen()


def _decode(arr):
    return np.array([x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
                     for x in arr])


def _find_gene_idx(f, varkey, gene_ensembl=None, gene_symbol=None):
    """Locate the gene's column in var: Ensembl in _index, then symbol in _index,
    then symbol in feature_name. Pass either identifier (or both)."""
    import h5py
    vi = _decode(f[f"{varkey}/_index"][:])
    if gene_ensembl:
        hit = np.where(vi == gene_ensembl)[0]
        if hit.size:
            return int(hit[0]), f"ensembl _index ({gene_ensembl})"
    if gene_symbol:
        hit = np.where(vi == gene_symbol)[0]          # some h5ads key var by symbol
        if hit.size:
            return int(hit[0]), f"symbol _index ({gene_symbol})"
        fnk = f"{varkey}/feature_name"
        if fnk in f:
            node = f[fnk]
            if isinstance(node, h5py.Group):          # categorical: codes + categories
                cats = _decode(node["categories"][:])
                names = cats[node["codes"][:]]
            else:
                names = _decode(node[:])
            hit = np.where(names == gene_symbol)[0]
            if hit.size:
                return int(hit[0]), f"symbol feature_name ({gene_symbol})"
    return None, None


def read_gene_vector(path, gene_ensembl=None, gene_symbol=None,
                     counts="raw", chunk=250000, progress=True):
    """Return (gene_vec float32 [n_cells], how, source_desc) reading only what's
    needed. CSC -> direct column slice; CSR -> chunked row scan."""
    import h5py
    with h5py.File(path, "r") as f:
        use_raw = (counts == "raw" and "raw" in f and "X" in f["raw"])
        xkey = "raw/X" if use_raw else "X"
        varkey = "raw/var" if use_raw else "var"
        gene_idx, how = _find_gene_idx(f, varkey, gene_ensembl, gene_symbol)
        if gene_idx is None:
            return None, None, None

        X = f[xkey]
        enc = X.attrs.get("encoding-type", "")
        if isinstance(enc, (bytes, bytearray)):
            enc = enc.decode()
        shape = tuple(int(s) for s in X.attrs["shape"])
        n_cells = shape[0]
        gene_vec = np.zeros(n_cells, dtype=np.float32)
        indptr = X["indptr"][:]

        if "csc" in enc:
            # gene is one contiguous column -> read only its nonzeros (cheap)
            s, e = int(indptr[gene_idx]), int(indptr[gene_idx + 1])
            if e > s:
                gene_vec[X["indices"][s:e]] = X["data"][s:e]
            return gene_vec, how, f"{xkey} [CSC direct col {gene_idx}]"

        # CSR: scan indices in row-chunks; read data only at matches
        idx_ds, dat_ds = X["indices"], X["data"]
        n_chunks = math.ceil(n_cells / chunk)
        it = range(n_chunks)
        if progress:
            it = _progress(it, total=n_chunks, desc=f"scan {os.path.basename(path)}")
        for c in it:
            i0, i1 = c * chunk, min(n_cells, (c + 1) * chunk)
            s, e = int(indptr[i0]), int(indptr[i1])
            if e <= s:
                continue
            idx = idx_ds[s:e]                      # the one big read per chunk
            hits = np.flatnonzero(idx == gene_idx)
            if hits.size:
                local = indptr[i0:i1 + 1] - s      # row boundaries within chunk
                rows = np.searchsorted(local, hits, side="right") - 1
                vals = dat_ds[(s + hits).tolist()]  # sparse fancy read (few values)
                gene_vec[i0 + rows] = vals
        return gene_vec, how, f"{xkey} [CSR chunked scan, chunk={chunk}]"


def read_obs_column(path, col, region_col=None):
    """Obs is small; read it backed (X stays on disk). Returns (groups, regions)."""
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    if col not in a.obs.columns:
        raise SystemExit(f"--groupby '{col}' not in obs of {os.path.basename(path)}. "
                         f"available: {list(a.obs.columns)}")
    groups = a.obs[col].astype(str).values
    regions = None
    if region_col:
        if region_col not in a.obs.columns:
            raise SystemExit(f"--region-col '{region_col}' not in obs")
        regions = a.obs[region_col].astype(str).values
    return groups, regions


def summarise(counts, groups, min_cells):
    df = pd.DataFrame({"raw": counts, "grp": np.asarray(groups)})
    g = df.groupby("grp", observed=True)["raw"]
    out = g.agg(
        n_cells="size",
        mean_raw="mean",
        pct_expressing=lambda x: 100.0 * (x > 0).mean(),
        n_expressing=lambda x: int((x > 0).sum()),
    ).reset_index()
    out = out[out["n_cells"] >= min_cells]
    return out.sort_values(["pct_expressing", "mean_raw"], ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", nargs="+", required=True)
    ap.add_argument("--gene-ensembl", default=None,
                    help="Ensembl gene id (e.g. ENSG00000099840). Default: IZUMO4 "
                         "if neither --gene-ensembl nor --gene-symbol is given.")
    ap.add_argument("--gene-symbol", default=None,
                    help="gene symbol (e.g. IZUMO4); matched in var _index or feature_name")
    ap.add_argument("--groupby", default="supercluster_term")
    ap.add_argument("--region-col", default=None)
    ap.add_argument("--counts", choices=["raw", "X"], default="raw",
                    help="raw counts (raw/X, default) or main X")
    ap.add_argument("--chunk-cells", type=int, default=250000,
                    help="CSR scan chunk size; smaller = less RAM, more passes")
    ap.add_argument("--min-cells", type=int, default=50)
    ap.add_argument("--out", default=None,
                    help="output tsv (default: <symbol_or_ensembl>_by_celltype.tsv)")
    a = ap.parse_args()

    # resolve gene: use whatever was given; fall back to IZUMO4 only if BOTH absent
    gene_ensembl, gene_symbol = a.gene_ensembl, a.gene_symbol
    if gene_ensembl is None and gene_symbol is None:
        gene_ensembl, gene_symbol = DEFAULT_ENSEMBL, DEFAULT_SYMBOL
    gene_label = gene_symbol or gene_ensembl
    out_path = a.out or f"{gene_label}_by_celltype.tsv"
    print(f"[gene] {gene_label}  (ensembl={gene_ensembl}, symbol={gene_symbol})")

    all_counts, all_groups, all_regions = [], [], []
    for path in a.h5ad:
        print(f"[read] {path}")
        gene_vec, how, src = read_gene_vector(path, gene_ensembl, gene_symbol,
                                              counts=a.counts, chunk=a.chunk_cells)
        if gene_vec is None:
            print(f"   [warn] {gene_label} not found in var; skipping")
            continue
        print(f"   gene via {how}; matrix {src}; {len(gene_vec):,} cells")
        groups, regions = read_obs_column(path, a.groupby, a.region_col)
        if len(groups) != len(gene_vec):
            raise SystemExit(f"obs ({len(groups)}) vs matrix ({len(gene_vec)}) "
                             f"cell-count mismatch in {path}")
        all_counts.append(gene_vec)
        all_groups.append(groups)
        if regions is not None:
            all_regions.append(regions)

    if not all_counts:
        raise SystemExit(f"{gene_label} not found in any input")

    counts = np.concatenate(all_counts)
    groups = np.concatenate(all_groups)
    print(f"\n[data] {len(counts):,} nuclei; {len(pd.unique(groups))} groups")
    print(f"[data] overall pct-expressing: {100*(counts>0).mean():.3f}%  "
          f"mean raw: {counts.mean():.4f}")

    res = summarise(counts, groups, a.min_cells)
    res.to_csv(out_path, sep="\t", index=False)
    print(f"\n[top cell types by % expressing]  (>= {a.min_cells} cells)")
    print(res.head(20).to_string(index=False))
    print(f"\n[written] {out_path}")

    if all_regions:
        regions = np.concatenate(all_regions)
        combo = pd.Series(groups).str.cat(pd.Series(regions), sep=" | ").values
        res_r = summarise(counts, combo, a.min_cells)
        rp = out_path.replace(".tsv", "_by_region.tsv")
        res_r.to_csv(rp, sep="\t", index=False)
        print(f"[written] {rp}  (cell type x region)")


if __name__ == "__main__":
    main()
