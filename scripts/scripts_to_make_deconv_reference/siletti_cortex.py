#!/usr/bin/env python3
"""Build siletti_cortex: neuronal/non-neuronal reference from Siletti 2023,
restricted to NEOCORTEX matching the bulk regions (BA9/BA4/BA22p). For the
SZ07 composition study (is the synaptic signal a neuron:glia shift?).

Sources (two h5ad, shared with siletti_glia; downloaded to <source-dir>):
  export SILETTI_NEURONS_URL=...      # CELLxGENE Siletti "Neurons" h5ad
  export SILETTI_NONNEURONS_URL=...   # CELLxGENE Siletti "Nonneurons" h5ad
No stable programmatic URL exists for these; get the permanent links from the
CELLxGENE collection linked at github.com/linnarsson-lab/adult-human-brain.

Output: <outdir>/{matrix.mtx.gz, features.tsv.gz, cells.tsv.gz, provenance.json}
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _deconv_common as C

# Stable CELLxGENE identifiers (verified by cell count: 2,480,956 / 888,263).
# Collection "Human Brain Cell Atlas v1.0". Files are resolved local-first:
# if already staged at <source-dir>, used as-is with no network (works behind
# the YC allowlist); otherwise a fresh asset URL is resolved from these IDs.
CXG_COLLECTION = "283d65eb-dd53-496d-adb7-7570c7caa443"
SILETTI_SOURCES = {
    "siletti_neurons.h5ad":     "8e10f1c4-8e98-41e5-b65f-8cd89a887122",  # All neurons
    "siletti_nonneurons.h5ad":  "b165f033-9dec-468a-9248-802fc6902a74",  # All non-neuronal cells
}

# supercluster -> class. Real Siletti supercluster_term strings (31 total: 21
# neuronal + 10 non-neuronal). Two classes; balance handled by
# cap-per-supercluster so the non_neuronal class is not all oligodendrocyte.
# NB: the "All neurons" file is WHOLE-BRAIN; the neocortex dissection filter
# (below) is what restricts this reference to BA9/BA4/BA22p territory.
MAPPING = {
    # --- 21 neuronal superclusters ---
    "Upper-layer intratelencephalic": "neuronal",
    "Deep-layer intratelencephalic": "neuronal",
    "Deep-layer near-projecting": "neuronal",
    "Deep-layer corticothalamic and 6b": "neuronal",
    "MGE interneuron": "neuronal", "CGE interneuron": "neuronal",
    "LAMP5-LHX6 and Chandelier": "neuronal", "Amygdala excitatory": "neuronal",
    "Hippocampal CA1-3": "neuronal", "Hippocampal CA4": "neuronal",
    "Hippocampal dentate gyrus": "neuronal", "Medium spiny neuron": "neuronal",
    "Eccentric medium spiny neuron": "neuronal", "Splatter": "neuronal",
    "Cerebellar inhibitory": "neuronal", "Thalamic excitatory": "neuronal",
    "Midbrain-derived inhibitory": "neuronal", "Mammillary body": "neuronal",
    "Lower rhombic lip": "neuronal", "Upper rhombic lip": "neuronal",
    "Miscellaneous": "neuronal",
    # --- 10 non-neuronal superclusters ---
    "Astrocyte": "non_neuronal", "Oligodendrocyte": "non_neuronal",
    "Oligodendrocyte precursor": "non_neuronal",
    "Committed oligodendrocyte precursor": "non_neuronal",
    "Microglia": "non_neuronal", "Vascular": "non_neuronal",
    "Fibroblast": "non_neuronal", "Ependymal": "non_neuronal",
    "Bergmann glia": "non_neuronal", "Choroid plexus": "non_neuronal",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--label-col", default="supercluster_term")
    ap.add_argument("--subject-col", default="donor_id")
    ap.add_argument("--tissue-col", default="")           # region handled below
    ap.add_argument("--tissue-terms", default="")
    # Neocortex filter: ROIGroupCoarse == "Cerebral cortex" is a single clean
    # term (1,156,201 neurons) that excludes hippocampus/allocortex (listed
    # separately) and all subcortical structures. The neocortex is NOT lobe-
    # subdivided in ROIGroup* (frontal/temporal would need the 105-value 'roi'
    # field), so we take all neocortex -- correct for a neuron:glia ratio
    # question, where more reference cells stabilise the glial classes and the
    # ratio, not the exact areal match, is what matters.
    ap.add_argument("--dissection-col", default="ROIGroupCoarse")
    ap.add_argument("--dissection-terms", default="Cerebral cortex")
    # This CELLxGENE Siletti distribution stores raw integer counts directly in
    # X (no .raw, no layers) -- verified integer-valued. Override if a different
    # asset puts counts elsewhere; the builder warns if the chosen layer is
    # non-integer (i.e. accidentally log/normalised).
    ap.add_argument("--counts-layer", default="X")
    ap.add_argument("--cells-per-class", type=int, default=500)
    ap.add_argument("--cap-per-supercluster", type=int, default=2000)
    # explicit source-file overrides (mirror hnoca.py): use these when the h5ad
    # live somewhere other than <source-dir>/siletti_{neurons,nonneurons}.h5ad
    # or under different names. Precedence: flag > env > staged default > CxG.
    ap.add_argument("--neurons-h5ad",
                    default=os.environ.get("SILETTI_NEURONS_H5AD", ""))
    ap.add_argument("--nonneurons-h5ad",
                    default=os.environ.get("SILETTI_NONNEURONS_H5AD", ""))
    a = ap.parse_args()

    if C.already_built(a.outdir):
        C.eprint(f"[skip] already built: {a.outdir}")
        return

    neu = C.ensure_source(
        a.neurons_h5ad or os.path.join(a.source_dir, "siletti_neurons.h5ad"),
        cxg=(CXG_COLLECTION, SILETTI_SOURCES["siletti_neurons.h5ad"]))
    non = C.ensure_source(
        a.nonneurons_h5ad or os.path.join(a.source_dir, "siletti_nonneurons.h5ad"),
        cxg=(CXG_COLLECTION, SILETTI_SOURCES["siletti_nonneurons.h5ad"]))

    C.build_reference(
        h5ad_paths=[neu, non], outdir=a.outdir,
        label_fn=lambda obs: C.apply_flat_mapping(obs, a.label_col, MAPPING),
        subject_col=a.subject_col,
        tissue_col=a.tissue_col, tissue_terms=a.tissue_terms.split(","),
        dissection_col=a.dissection_col,
        dissection_terms=a.dissection_terms.split(","),
        counts_layer=a.counts_layer,
        cells_per_class=a.cells_per_class,
        balance_col=a.label_col, cap_per_balance=a.cap_per_supercluster,
        seed=a.seed,
        provenance_extra={"reference": "siletti_cortex",
                          "purpose": "SZ07 neuron/non-neuron composition",
                          "mapping": MAPPING})


if __name__ == "__main__":
    main()
