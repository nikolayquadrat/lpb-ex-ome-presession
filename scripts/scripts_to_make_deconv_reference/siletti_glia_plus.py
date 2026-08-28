#!/usr/bin/env python3
"""Build siletti_glia_plus: non-neuronal reference with GLIA *and* the other
non-neuronal lineages RESOLVED, neurons pooled, from Siletti 2023, ALL regions.

This is the postmortem-tissue counterpart to siletti_glia. siletti_glia was built
for the iPSC methodological study and therefore DROPPED the lineages a forebrain
iPSC culture cannot make (microglia, vascular, fibroblast, choroid plexus,
Bergmann glia). A real postmortem brain (e.g. SZ07) contains all of them, and two
of them are central to this cohort's biology -- microglia carry the innate/TLR
signature and the vasculature indexes the ischaemic pathology -- so here we ADD
them back as resolved deconvolution classes instead of dropping them.

Class set (relative to siletti_glia, the additions are marked +):
  glia (neuroectoderm):  Astrocyte, Oligodendrocyte, Oligodendrocyte precursor,
                         Committed oligodendrocyte precursor, Ependymal
  + Bergmann glia        (cerebellar radial astrocyte; real in postmortem tissue)
  + Microglia            (innate-immune; drives the TLR/NF-kB signature)
  + Vascular             (endothelial/mural; indexes ischaemia)
  + Fibroblast           (meningeal/perivascular)
  + Choroid plexus       (CSF-secreting epithelium)
  Neuron (nuisance):     ALL neuronal superclusters pooled to one class

Two outputs, same rationale as siletti_glia (resolution is bounded by how
distinct the reference profiles are):
  <outdir>_supercluster/   non-neuronal classes at supercluster level -- the
                           robust baseline; always trustworthy.
  <outdir>_cluster/        the glial-neuroectoderm classes resolved to CLUSTER
                           level; the *added* non-neuronal lineages are kept at
                           supercluster level even here (several are single-
                           supercluster / low-n -- Bergmann glia is 1 cluster,
                           Choroid plexus / Ependymal are small -- so cluster
                           splitting them buys nothing and only adds noise).
                           Treat cluster level as provisional; trust only where
                           it aggregates to the supercluster result.

Neurons are pooled at both levels (the question is non-neuronal composition).

Sources: the same two Siletti h5ad as siletti_cortex / siletti_glia (shared,
staged once; not re-downloaded if a sibling already fetched them).
  export SILETTI_NEURONS_URL / SILETTI_NONNEURONS_URL   (see siletti_cortex.py)

The CLUSTER column name is not hard-known; pass --cluster-col to match your
h5ad (the builder hard-errors listing available columns if it is absent).
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _deconv_common as C

# Same stable Siletti sources as siletti_cortex / siletti_glia (shared, staged once).
CXG_COLLECTION = "283d65eb-dd53-496d-adb7-7570c7caa443"
SILETTI_SOURCES = {
    "siletti_neurons.h5ad":     "8e10f1c4-8e98-41e5-b65f-8cd89a887122",
    "siletti_nonneurons.h5ad":  "b165f033-9dec-468a-9248-802fc6902a74",
}

# supercluster -> role. Real Siletti supercluster_term strings (31 total),
# grouped per the S3 'Class auto-annotation' column.
# NEURON superclusters (Class=NEUR) collapse to one "Neuron" nuisance class.
NEURON_SC = {
    "Upper-layer intratelencephalic", "Deep-layer intratelencephalic",
    "Deep-layer near-projecting", "Deep-layer corticothalamic and 6b",
    "MGE interneuron", "CGE interneuron", "LAMP5-LHX6 and Chandelier",
    "Amygdala excitatory", "Hippocampal CA1-3", "Hippocampal CA4",
    "Hippocampal dentate gyrus", "Medium spiny neuron",
    "Eccentric medium spiny neuron", "Splatter", "Cerebellar inhibitory",
    "Thalamic excitatory", "Midbrain-derived inhibitory", "Mammillary body",
    "Lower rhombic lip", "Upper rhombic lip",
    "Miscellaneous",                 # pooled into the Neuron nuisance class
}
# Glial (neuroectoderm) superclusters -- resolved FINE at cluster level.
GLIA_SC = {
    "Astrocyte", "Oligodendrocyte", "Oligodendrocyte precursor",
    "Committed oligodendrocyte precursor", "Ependymal",
}
# The ADDED non-neuronal lineages (vs siletti_glia, which DROPs these). Kept at
# supercluster level even in the cluster build: several are single-supercluster
# or low-n, so cluster-splitting them adds noise without resolution.
OTHER_NONNEURON_SC = {
    "Bergmann glia",     # ASTRO (cerebellar radial); 1 cluster / 8k cells
    "Microglia",         # MGL;  innate-immune, the TLR/NF-kB signature
    "Vascular",          # ENDO/mural; indexes ischaemia
    "Fibroblast",        # FIB;  meningeal/perivascular
    "Choroid plexus",    # CHRP; CSF epithelium
}
# Nothing is dropped in glia_plus (contrast siletti_glia's DROP_SC).
DROP_SC = set()

_ALL_KNOWN = NEURON_SC | GLIA_SC | OTHER_NONNEURON_SC | DROP_SC


def rules_for(level):
    """coarse_rules for compose_hierarchical_label at the chosen glia level.
    Actions: COARSE -> 'Neuron'; FINE -> cluster value; SELF -> supercluster
    value (label_fn maps SELF to FINE with fine_col = supercluster_col); DROP."""
    r = {}
    for sc in NEURON_SC:
        r[sc] = "COARSE"                       # -> "Neuron"
    for sc in DROP_SC:
        r[sc] = "DROP"
    for sc in GLIA_SC:
        r[sc] = "FINE" if level == "cluster" else "SELF"
    for sc in OTHER_NONNEURON_SC:
        r[sc] = "SELF"                         # always supercluster-level
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--levels", default="supercluster,cluster",
                    help="comma list of {supercluster,cluster}")
    ap.add_argument("--supercluster-col", default="supercluster_term")
    ap.add_argument("--cluster-col", default="cluster_id",
                    help="461-cluster label column in the h5ad obs")
    ap.add_argument("--subject-col", default="donor_id")
    # Raw integer counts are in X for this Siletti distribution (no .raw/layers,
    # verified integer). Builder warns if the chosen layer is non-integer.
    ap.add_argument("--counts-layer", default="X")
    # Per-class cap. Cluster level benefits from more cells per subtype because
    # abundant glial superclusters fan into many well-populated clusters
    # (Astrocyte 155k -> 13, Oligodendrocyte 490k -> 8). Rare classes cap out
    # where they are; min_cells_per_class only warns, never drops them.
    ap.add_argument("--cells-per-class", type=int, default=500,
                    help="supercluster-level per-class cap")
    ap.add_argument("--cells-per-class-cluster", type=int, default=1500,
                    help="cluster-level per-class cap (larger: abundant glial "
                         "clusters are well-populated)")
    # explicit source-file overrides (mirror siletti_glia / siletti_cortex):
    # flag > env > staged default > CxG resolve.
    ap.add_argument("--neurons-h5ad",
                    default=os.environ.get("SILETTI_NEURONS_H5AD", ""))
    ap.add_argument("--nonneurons-h5ad",
                    default=os.environ.get("SILETTI_NONNEURONS_H5AD", ""))
    a = ap.parse_args()

    neu = C.ensure_source(
        a.neurons_h5ad or os.path.join(a.source_dir, "siletti_neurons.h5ad"),
        cxg=(CXG_COLLECTION, SILETTI_SOURCES["siletti_neurons.h5ad"]))
    non = C.ensure_source(
        a.nonneurons_h5ad or os.path.join(a.source_dir, "siletti_nonneurons.h5ad"),
        cxg=(CXG_COLLECTION, SILETTI_SOURCES["siletti_nonneurons.h5ad"]))

    for level in [x.strip() for x in a.levels.split(",") if x.strip()]:
        # Flat sibling output per level: --outdir .../siletti_glia_plus gives
        # .../siletti_glia_plus_supercluster and .../siletti_glia_plus_cluster,
        # each a complete independent canonical reference.
        out = f"{a.outdir}_{level}"
        if C.already_built(out):
            C.eprint(f"[skip] already built: {out}")
            continue
        rules = rules_for(level)

        def label_fn(obs, _rules=rules, _level=level):
            # "SELF" means emit the supercluster value; implement by composing
            # with fine_col = supercluster_col for those, and cluster_col for the
            # FINE glial classes. Because SELF and FINE need DIFFERENT fine
            # columns, resolve per-cell in two passes and merge.
            import pandas as pd
            sc_col, cl_col = a.supercluster_col, a.cluster_col
            # pass 1: SELF/DROP/COARSE resolved against supercluster as fine_col
            r_self = {k: ("FINE" if v == "SELF" else v) for k, v in _rules.items()}
            lab_self = C.compose_hierarchical_label(
                obs, sc_col, sc_col, r_self, coarse_name="Neuron")
            if _level != "cluster":
                return lab_self
            # pass 2 (cluster level): the GLIA_SC classes should instead be their
            # cluster value. Recompute with cluster as fine_col and overwrite only
            # those cells whose supercluster is a FINE-glia one.
            r_fine = {k: ("FINE" if v == "SELF" else v) for k, v in _rules.items()}
            lab_fine = C.compose_hierarchical_label(
                obs, sc_col, cl_col, r_fine, coarse_name="Neuron")
            scv = obs[sc_col].astype(str)
            use_cluster = scv.isin(GLIA_SC)          # only glia go to cluster res
            return lab_self.where(~use_cluster, lab_fine)

        cpc = (a.cells_per_class_cluster if level == "cluster"
               else a.cells_per_class)
        C.build_reference(
            h5ad_paths=[neu, non], outdir=out,
            label_fn=label_fn, subject_col=a.subject_col,
            tissue_col=None, tissue_terms=None,          # ALL regions
            dissection_col=None, dissection_terms=None,
            counts_layer=a.counts_layer,
            cells_per_class=cpc,
            balance_col=None,          # each non-neuronal class is one label already
            seed=a.seed,
            provenance_extra={"reference": f"siletti_glia_plus_{level}",
                              "purpose": "postmortem non-neuronal deconvolution "
                                         "(glia + microglia/vascular/fibroblast/"
                                         "choroid plexus/Bergmann)",
                              "glia_level": level,
                              "cells_per_class": cpc,
                              "neuroectoderm_only": False,
                              "added_vs_siletti_glia": sorted(OTHER_NONNEURON_SC),
                              "dropped_superclusters": sorted(DROP_SC)})


if __name__ == "__main__":
    main()
