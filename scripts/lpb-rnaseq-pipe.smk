# =============================================================================
# RNA-seq alignment pipeline, GTEx V11-compatible
# =============================================================================
# Purpose: produce STAR-aligned BAM files for the 8-donor / ~36-sample brain
# RNA-seq cohort, using exactly the alignment parameters and reference files
# of the GTEx V11 release (https://gtexportal.org/home/downloads/adult-gtex/
# bulk_tissue_expression). This makes the BAMs and gene-level read counts
# directly mergeable with the GTEx V11 expression matrices for DROP/OUTRIDER
# reference panel expansion.
#
# Reference setup (GTEx V11):
#   - genome:    Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta (TOPMed ref)
#   - annotation: GENCODE v47 (gencode.v47.GRCh38.annotation.gtf)
#   - STAR:      v2.7.11b (use this exact version for full reproducibility)
#   - RSEM:      v1.3.3
#   - sjdbOverhang: depends on YOUR read length (NOT 75, unless your reads are
#                   2x76bp like GTEx). Set READ_LENGTH below to match your data.
#
# Outputs per sample:
#   - 03_bam_star/{sample}/Aligned.sortedByCoord.out.patched.md.bam
#       (the canonical GTEx-style output: STAR-aligned, coordinate-sorted,
#        duplicate-marked; this is the one to feed to DROP)
#   - 03_bam_star/{sample}/Aligned.toTranscriptome.out.bam
#       (transcriptome-coordinate BAM, for downstream RSEM if desired)
#   - 03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab
#       (STAR's gene-level counts, mergeable with GTEx V11 count matrix)
# =============================================================================

import pandas as pd
import os
import sys

# -----------------------------------------------------------------------------
# User-tunable parameters
# -----------------------------------------------------------------------------

# Read length in your FASTQs. CRITICAL: STAR's --sjdbOverhang must be
# (read_length - 1). GTEx itself uses 75 because their reads are 2x76bp.
# To check: zcat your_R1.fq.gz | head -2 | tail -1 | wc -c   (subtract 1 for newline)
READ_LENGTH = int(os.environ.get("READ_LENGTH", "150")) # 150 bp is the default
SJDB_OVERHANG = READ_LENGTH - 1

# Sample table
# columns: name, path, dataset1, dataset2, extension, include_in_analysis
SAMPLE_FILE = "/tmp/data/00_additional_files/sample_data.txt"
sample_data = pd.read_csv(SAMPLE_FILE, sep="\t")
sample_data = sample_data[sample_data["include_in_analysis"] == 1]

# Filter out samples whose forward FASTQ is not present on disk.
# This makes the pipeline robust to a master sample table that lists more
# samples than are currently available (e.g., when sequencing is rolling in
# in batches, or some samples have been moved/archived). The warning lists
# every dropped sample by name so a typo'd filename still surfaces clearly.
def _fwd_fastq_path(row):
    return f"/tmp/data/{row['path']}/{row['dataset1']}.{row['extension']}"

_present_mask = sample_data.apply(
    lambda r: os.path.exists(_fwd_fastq_path(r)), axis=1)

if (~_present_mask).any():
    missing = sample_data.loc[~_present_mask, "name"].tolist()
    print("=" * 72)
    print(f"WARNING: {len(missing)} sample(s) listed in {SAMPLE_FILE}")
    print("with include_in_analysis=1 but their forward FASTQ file is missing")
    print("on disk. These samples will be SILENTLY SKIPPED in this run.")
    print("If this is unintended (e.g., typo in 'path' or 'dataset1'),")
    print("fix the sample table or restore the missing files before re-running.")
    print()
    print("Skipped samples (and their expected R1 paths):")
    for name in missing:
        row = sample_data[sample_data["name"] == name].iloc[0]
        print(f"  - {name}: {_fwd_fastq_path(row)}")
    print("=" * 72)
    sample_data = sample_data.loc[_present_mask].reset_index(drop=True)

samples = sample_data["name"].tolist()
print(f"Aligning {len(samples)} samples")

# Reference paths.
# IMPORTANT: download these files to the indicated locations before running.
# See "Reference download instructions" at the bottom of this file.
REF_DIR     = "/tmp/data/00_additional_files/gtex_v11_refs"
GENOME_FA   = REF_DIR + "/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
GTF         = REF_DIR + "/gencode.v47.GRCh38.annotation.gtf"
STAR_INDEX  = REF_DIR + f"/star_index_oh{SJDB_OVERHANG}"
RSEM_INDEX  = REF_DIR + "/rsem_reference"

if not 25 <= READ_LENGTH <= 500:
    sys.exit(f"[ERR] READ_LENGTH={READ_LENGTH} is implausible for Illumina "
             "RNA-seq; expected ~50-250 (per mate, not the sum of both mates)")
print(f"[config] READ_LENGTH={READ_LENGTH}  -> sjdbOverhang={SJDB_OVERHANG}, "
      f"index={STAR_INDEX}", flush=True)

# Container versions (matching GTEx V11 release where possible)
STAR_CONTAINER   = "docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_0"
SAMTOOLS_CONTAINER = "docker://biocontainers/samtools:v1.9-4-deb_cv1"
PICARD_CONTAINER = "docker://broadinstitute/picard:2.27.5"
TRIM_CONTAINER   = "docker://mskaccess/trim_galore" 
# bwa-mem2 + samtools 1.13 mulled image. Both binaries are on PATH inside,
# so r05b can stream from samtools fastq -> bwa-mem2 mem -> samtools view
# in a single container invocation.
BWA_CONTAINER = (
    "docker://quay.io/biocontainers/"
    "mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:"
    "c5b8c4b7735290369693e2b63cfc1ea0732fde07-0"
)

# BBTools (bbduk.sh in particular). Used by r05c_entropy_filter to drop
# low-complexity / low-entropy reads (poly-A tails, simple repeats,
# AT-rich runs) before they're aligned to the contamination reference.
# Removing these reads at the FASTQ level is the most effective single
# intervention for reducing cross-mapping false positives, since
# low-complexity reads align cleanly to many references regardless of
# bwa-mem2 stringency settings.
BBMAP_CONTAINER = "docker://quay.io/biocontainers/bbmap:39.06--h92535d8_0"

# bedtools, used by r05a0_mask_conserved_regions to hard-mask conserved
# rRNA/tRNA regions in the contamination reference. These regions are
# extremely sequence-conserved across species (bacterial 16S/23S has
# stretches 100% identical to human rRNA, mitochondrial rRNA reads
# cross-map to apicomplexan plastid rRNA, etc.), so masking them at the
# reference level prevents host-derived rRNA reads from generating
# spurious contamination hits regardless of alignment stringency.
BEDTOOLS_CONTAINER = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"

# NCBI datasets CLI container, used by r05_prep_contamination_refs
# DATASETS_CONTAINER = "docker://quay.io/biocontainers/ncbi-datasets-cli:14.26.0"
DATASETS_CONTAINER = "docker://biocontainers/ncbi-datasets-cli:16.22.1_cv1"

# Env-driven, ON by default. Disable per run with -e USE_TRIMMING=0 at container
# start (same convention as ARCASHLA_ENABLED / CONTAMINATION_ENABLED / DECONV_ENABLED).
USE_TRIMMING = os.environ.get("USE_TRIMMING", "1").strip().lower() \
    in ("1", "true", "yes", "on")

# -----------------------------------------------------------------------------
# Optional contamination screen
# -----------------------------------------------------------------------------
# Three modes, chosen by which files are present in /tmp/data/00_additional_files/contamination/:
#
#   (A) PREBUILT: if contamination.fa + contamination.gtf + gtf_seqnames.tsv
#       are all present and non-empty, they're used as-is.
#
#   (B) ON-DEMAND: if (A)'s files are not all present but species.txt IS
#       present, the pipeline downloads each species' genome+GTF from NCBI
#       using the `datasets` CLI tool, concatenates the FASTAs into
#       contamination.fa, the GTFs into contamination.gtf, and emits a
#       gtf_seqnames.tsv that maps each FASTA seqname (NCBI accession) to
#       a species_id of the form <name>_<assembly_accession>.
#
#   (C) DISABLED: neither (A) nor (B) is satisfied. r05* rules are skipped
#       and the QC summary omits contamination columns.
#
# species.txt format (tab-separated, no header):
#   species_name<TAB>assembly_accession
# e.g.:
#   mesomycoplasma_hyorhinis<TAB>GCF_001705605.1
#
# The contig naming inside the auto-built contamination.fa preserves the
# native NCBI seqnames (e.g. ">NC_010163.1 Acholeplasma laidlawii PG-8A").
# r05b's awk uses the gtf_seqnames.tsv to attribute reads to species.
#
# Alignment parameters (r05b): bwa-mem2 mem is run with conservative
# settings to reduce false positives from rRNA homologs etc:
#   -T 60  : minimum alignment score (default 30) -- requires ~60bp of
#            near-perfect identity for a read to be reported
#   -k 25  : minimum seed length (default 19)    -- longer seeds reject
#            short noisy matches
#   -r 2.0 : re-seed threshold (default 1.5)     -- less aggressive
#            re-seeding; fewer marginal alignments
# Combined with the MAPQ>=20 post-filter, the result is much closer to
# "reads that genuinely come from this organism" than the defaults.

CONTAM_DIR          = "/tmp/data/00_additional_files/contamination"
CONTAM_FA           = f"{CONTAM_DIR}/contamination.fa"
CONTAM_GTF          = f"{CONTAM_DIR}/contamination.gtf"
CONTAM_SEQNAMES     = f"{CONTAM_DIR}/gtf_seqnames.tsv"
CONTAM_SPECIES_LIST = f"{CONTAM_DIR}/species.txt"

_prebuilt = all(
    os.path.exists(p) and os.path.getsize(p) > 0
    for p in (CONTAM_FA, CONTAM_GTF, CONTAM_SEQNAMES)
)
_have_species_list = (
    os.path.exists(CONTAM_SPECIES_LIST)
    and os.path.getsize(CONTAM_SPECIES_LIST) > 0
)
def _read_species_order(species_txt_path, seqname_map_path):
    """
    Return a list of species_id values in the canonical order they appear
    in species.txt. If species.txt is not present (i.e., prebuilt mode
    without the source list), fall back to deriving the ordering from the
    seqname map -- which guarantees consistent column ordering even when
    species.txt is unavailable, but reflects whatever order the seqname
    map was emitted in.

    species.txt parsing rules (matching r05_prep_contamination_reference):
      - Tab-separated: <species_name><TAB><accession>
      - Blank lines skipped
      - Lines starting with '#' skipped (after optional leading whitespace)
      - species_id = "<species_name>_<accession>", with whitespace in the
        species name replaced by underscores (matching how the seqname-map
        first column is constructed)
      - First occurrence wins on duplicates
    """
    order = []
    seen = set()

    # Path A: species.txt is present -- use it as the canonical source
    if species_txt_path and os.path.exists(species_txt_path) and os.path.getsize(species_txt_path) > 0:
        with open(species_txt_path) as f:
            for raw in f:
                line = raw.rstrip("\r\n")
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#"):
                    continue
                # Tab-split; tolerate stray whitespace around fields
                parts = [p.strip() for p in line.split("\t")]
                if len(parts) < 2:
                    # Malformed row (no tab) -- skip silently; the autobuild
                    # rule already errors on these and won't write the FASTA
                    continue
                species_name, accession = parts[0], parts[1]
                if not species_name or not accession:
                    continue
                # Same sanitisation as the autobuild rule:
                # collapse internal whitespace, '/' and '\\' to '_'
                safe = species_name
                for ch in (" ", "\t", "/", "\\"):
                    safe = safe.replace(ch, "_")
                # Collapse multiple underscores
                while "__" in safe:
                    safe = safe.replace("__", "_")
                species_id = f"{safe}_{accession}"
                if species_id not in seen:
                    seen.add(species_id)
                    order.append(species_id)
        return order

    # Path B: fallback -- read from the seqname map. Each species appears
    # on multiple lines (one per contig); we keep the first occurrence to
    # respect whatever order the autobuild rule emitted.
    if seqname_map_path and os.path.exists(seqname_map_path):
        with open(seqname_map_path) as f:
            header_skipped = False
            for raw in f:
                if not header_skipped:
                    header_skipped = True
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) == 2 and parts[0] and parts[0] not in seen:
                    seen.add(parts[0])
                    order.append(parts[0])
    return order

# CONTAMINATION_ENABLED can be forced on/off via the CONTAMINATION_ENABLED
# environment variable. Set to "0"/"false"/"no"/"off" to force-disable (skip
# r05* contamination rules regardless of which input files are present), or
# "1"/"true"/"yes"/"on" to force-enable. Anything else, including the default
# "auto", uses filesystem auto-detection (the original behaviour).
_contamination_override = os.environ.get("CONTAMINATION_ENABLED", "auto").lower()
if _contamination_override in ("0", "false", "no", "off"):
    CONTAMINATION_ENABLED = False
elif _contamination_override in ("1", "true", "yes", "on"):
    CONTAMINATION_ENABLED = True
elif _contamination_override == "auto":
    CONTAMINATION_ENABLED = _prebuilt or _have_species_list
else:
    print(f"[contamination screen] unrecognised "
          f"CONTAMINATION_ENABLED='{_contamination_override}', "
          f"falling back to auto-detection", file=sys.stderr)
    CONTAMINATION_ENABLED = _prebuilt or _have_species_list

# AUTOBUILD only makes sense when no prebuilt files exist AND we have a
# species list to download from. Independent of the override.
CONTAMINATION_AUTOBUILD = (not _prebuilt) and _have_species_list

# Warn the user if a forced-on flag will produce broken downstream rules
# because there is nothing for r05* to consume.
if CONTAMINATION_ENABLED and not _prebuilt and not _have_species_list:
    print(f"[contamination screen] WARNING: CONTAMINATION_ENABLED is forced ON "
          f"but neither prebuilt files nor {CONTAM_SPECIES_LIST} are present; "
          f"contamination rules will fail with missing-input errors.",
          file=sys.stderr)

if CONTAMINATION_ENABLED:
    if CONTAMINATION_AUTOBUILD:
        print(f"[contamination screen] enabled (auto-build mode)",
              file=sys.stderr)
        print(f"[contamination screen] species list: {CONTAM_SPECIES_LIST}",
              file=sys.stderr)
        print(f"[contamination screen] will download genomes from NCBI and "
              f"build {CONTAM_FA} / {CONTAM_GTF} / {CONTAM_SEQNAMES}",
              file=sys.stderr)
    else:
        print(f"[contamination screen] enabled (pre-built files)",
              file=sys.stderr)
        print(f"[contamination screen] FASTA: {CONTAM_FA}", file=sys.stderr)
        print(f"[contamination screen] seqname map: {CONTAM_SEQNAMES}",
              file=sys.stderr)
else:
    print(f"[contamination screen] disabled (need either pre-built "
          f"contamination.fa+gtf+gtf_seqnames.tsv or a species.txt in "
          f"{CONTAM_DIR}, or unset CONTAMINATION_ENABLED if it is forced off)",
          file=sys.stderr)

# Typical aligned read length (bp), used by r05f_contamination_uniformity
# in the Lander-Waterman expected-breadth formula. RNA-seq read lengths
# typically range 75-150 bp; the uniformity ratio is robust to +/-20%
# error so 75 is a reasonable default. Override here if your data is
# very different (e.g., long-read sequencing).
CONTAM_READ_LENGTH = 150

# -----------------------------------------------------------------------------
# Optional HLA class I typing with arcasHLA
# -----------------------------------------------------------------------------
# arcasHLA types HLA class I (and II) directly from RNA-seq reads. It is run
# here OPTIONALLY: the rules only activate if a fully built arcasHLA reference
# is present at the expected (post-bind) location. The reference is large and
# must be built once on the host with the companion installer script
# (install-arcashla-ref.sh), which handles all the bioconda-biocontainer
# quirks (the container ships only a partial reference). The host build dir is:
#
#   /mnt/data/rnaseq/rnaseq-drop/00_additional_files/arcashla_ref
#
# At run time that host directory is bind-mounted so the pipeline sees it at:
#
#   /tmp/data/00_additional_files/arcashla_ref
#
# and the dat/ subdirectory inside it must additionally be bound OVER the
# biocontainer's own (incomplete) dat/ directory. See the snakemake invocation
# note in the installer script and at the bottom of this file: the arcasHLA
# rules expect both
#   (a) /tmp:/tmp                                              (outer bind)
#   (b) arcashla_ref/dat : /usr/local/share/arcas-hla-0.6.0-2/dat   (inner bind)
#
# If the reference is not present (e.g., the installer was not run, or this is
# a fresh machine), ARCASHLA_ENABLED is False and the HLA rules + the HLA
# summary are simply omitted from the DAG. The rest of the pipeline is
# unaffected.

ARCASHLA_REF_DIR = "/tmp/data/00_additional_files/arcashla_ref"
ARCASHLA_DAT     = f"{ARCASHLA_REF_DIR}/dat"
ARCASHLA_CONTAINER = "docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2"

# Container-internal mount point for the dat/ directory. The host
# arcashla_ref/dat must be bound over this path at run time (see above).
ARCASHLA_CONTAINER_DAT_PATH = "/usr/local/share/arcas-hla-0.6.0-2/dat"

# Key files that must all exist (non-empty) for the reference to count as
# fully built. These mirror what install-arcashla-ref.sh verifies after it
# runs `arcasHLA reference --rebuild`. If any is missing, arcasHLA cannot run.
_arcashla_key_files = (
    f"{ARCASHLA_DAT}/IMGTHLA/hla.dat",
    f"{ARCASHLA_DAT}/IMGTHLA/wmda/hla_nom_p.txt",
    f"{ARCASHLA_DAT}/ref/hla.p.json",
    f"{ARCASHLA_DAT}/ref/hla.convert.json",
    f"{ARCASHLA_DAT}/ref/hla.idx",
    f"{ARCASHLA_DAT}/ref/hla_partial.idx",
)

_arcashla_override = os.environ.get("ARCASHLA_ENABLED", "auto").lower()
if _arcashla_override == "auto":
    ARCASHLA_ENABLED = all(
        os.path.exists(p) and os.path.getsize(p) > 0
        for p in _arcashla_key_files
    )
elif _arcashla_override in ("0", "false", "no", "off"):
    ARCASHLA_ENABLED = False
elif _arcashla_override in ("1", "true", "yes", "on"):
    ARCASHLA_ENABLED = True
else:
    print(f"[arcasHLA] unrecognised ARCASHLA_ENABLED='{_arcashla_override}', "
          f"falling back to auto-detection", file=sys.stderr)
    ARCASHLA_ENABLED = all(
        os.path.exists(p) and os.path.getsize(p) > 0
        for p in _arcashla_key_files
    )

if ARCASHLA_ENABLED:
    _missing = [p for p in _arcashla_key_files
                if not (os.path.exists(p) and os.path.getsize(p) > 0)]
    if _missing:
        print(f"[arcasHLA] WARNING: ARCASHLA_ENABLED is forced ON but the "
              f"reference is incomplete; arcasHLA rules will fail at run time",
              file=sys.stderr)
        print(f"[arcasHLA] missing/empty: {_missing[:3]}"
              f"{'...' if len(_missing) > 3 else ''}", file=sys.stderr)
    else:
        print(f"[arcasHLA] enabled -- reference found at {ARCASHLA_DAT}",
              file=sys.stderr)
        print(f"[arcasHLA] REMINDER: the snakemake run must bind the host "
              f"arcashla_ref/dat over {ARCASHLA_CONTAINER_DAT_PATH} "
              f"(see installer script for the exact --apptainer-args / "
              f"--singularity-args line)", file=sys.stderr)
else:
    _missing = [p for p in _arcashla_key_files
                if not (os.path.exists(p) and os.path.getsize(p) > 0)]
    print(f"[arcasHLA] disabled -- reference not fully built at "
          f"{ARCASHLA_DAT}", file=sys.stderr)
    print(f"[arcasHLA] missing/empty: {_missing[:3]}"
          f"{'...' if len(_missing) > 3 else ''}", file=sys.stderr)
    print(f"[arcasHLA] run install-arcashla-ref.sh on the host to build it, "
          f"then ensure arcashla_ref/ is bind-mounted into the container",
          file=sys.stderr)



def get_fastq_path(wildcards):
    """Same logic as existing get_fastq_path(): paired or single-end."""
    row = sample_data[sample_data.name == wildcards.sample].iloc[0]
    fwd = f"/tmp/data/{row['path']}/{row['dataset1']}.{row['extension']}"
    if str(row['dataset2']).strip().lower() == "nan":
        return [fwd, fwd]   # double up for SE
    rev = f"/tmp/data/{row['path']}/{row['dataset2']}.{row['extension']}"
    return [fwd, rev]

# -----------------------------------------------------------------------------
# Sex inference (RNA-seq, expression-based)
# -----------------------------------------------------------------------------
# RNA-seq sex inference uses MARKER-GENE EXPRESSION:
#   - XIST: long non-coding RNA that coats and silences the inactive X. It is
#     highly expressed only in cells with >=2 X chromosomes (i.e. XX, and also
#     XXY), and near-silent in XY. -> female / 2nd-X marker.
#   - Y-linked genes: expressed only when a Y chromosome is present. -> male
#     marker. We use a panel (sum) rather than one gene for robustness; several
#     of these are expressed in brain specifically (NLGN4Y, especially).
# The two signals together also distinguish sex-chromosome aneuploidy:
#   XX  -> XIST high, Y low ;  XY -> XIST low, Y high ;
#   XXY -> XIST high AND Y high ; X0 -> both low.

SEX_MARKER_GENES = [
    "XIST",                                   # female / inactive-X marker
    "RPS4Y1", "DDX3Y", "UTY", "USP9Y",        # Y panel (X-degenerate, broadly expressed)
    "KDM5D", "EIF1AY", "ZFY", "TXLNGY",
    "NLGN4Y",                                 # Y, notably expressed in brain
]
SEX_FEMALE_GENES = {"XIST"}                   # the rest are treated as the Y panel

# CPM thresholds for the sex call. Chosen to sit in the empirical gap between
# each marker's real signal and its opposite-sex cross-mapping baseline:
#   - true males carry XIST up to ~3 CPM (X/Y-adjacent noise); true females
#     start at ~17 CPM -> XIST cutoff 10 sits cleanly between.
#   - true females carry Y-panel up to ~7 CPM (X paralog cross-mapping);
#     true males start at ~64 CPM -> Y cutoff 20 sits cleanly between.
# A genuine XXY shows BOTH markers above threshold simultaneously.
SEX_XIST_CPM_MIN = 10.0
SEX_Y_CPM_MIN    = 20.0

# -----------------------------------------------------------------------------
# Cell-composition marker scoring (RNA-seq) / r08 rules
# -----------------------------------------------------------------------------
# Purpose: a cheap, reference-free check on whether a sample's neuron:glia
# balance is an outlier. This is the first-pass answer to "is an apparent
# synaptic/neuronal signature actually cell composition (grey-matter content,
# dissection) rather than biology?". It is NOT deconvolution: the output is a
# relative score across samples, not a cell-type proportion. Reference-based
# proportions (hspe / Bisque) are a future (?) separate, optional block.
#
# METHOD. Per gene: CPM from unique reads over the gene span -> log2(CPM+1) ->
# z-score ACROSS SAMPLES WITHIN A REGION. Per cell type: mean of its markers' z.
#
# NB: SOX10 is a canonical OPC/oligodendrocyte marker but is ALSO a VUS
# candidate gene in this study. It is deliberately omitted so the
# composition check stays independent of candidate evaluation.

# delete the 04_qc/cellcomp folder after applying any changes

CELL_MARKER_GENES = {
    # neuronal markers:
    "neuron": ["RBFOX3", "MAP2", "NEFL", "NEFM", "NEFH", "INA", "ELAVL3", "ELAVL4"], # non-synaptic

    # non-neuronal markers: 
    "astrocyte":       ["GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "GLUL"],
    "oligodendroglial": ["OLIG1", "OLIG2", "PLP1", "CNP", "MBP", "MOG", "MAG", "MOBP", "CLDN11"], # mature + OPC
    "endothelial":     ["CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "PTPRB"],
    "mural":           ["RGS5", "PDGFRB", "KCNJ8", "ABCC9", "HIGD1B", "NOTCH3", "MYH11"], # pericytes, VSMC
    "fibroblast":      ["DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL15A1"], # brain/perivascular fibroblasts
    "microglia":       ["P2RY12", "TMEM119", "SALL1", "GPR34", "HEXB"],
    "ependymal":       ["FOXJ1","PIFO", "TPPP3", "RSPH1", "FAM183A", "CCDC153", "SNTN"],
    "bam_macrophage":  ["MRC1", "CD163", "LYVE1", "F13A1", "MS4A7", "STAB1"],
    "choroid_plexus":  ["TTR", "FOLR1", "OTX2", "CLIC6", "KRT18", "AQP1"],

    # additional markers:
    "ruzicka24_down": ["BSN", "SHANK2", "RASGRF1", "NEURL1", "PSEN1", "BCR"],
    "fromer16_up":    ["TACR3", "HGF", "C3orf52", "MC4R", "ADCYAP1", "MST1L", "BEND4", "SCN9A", "TMEM26", "BIRC3"],
    "fromer16_down":  ["IGF2", "RERGL", "HPSE2", "SLCO2A1", "CLEC3B", "RAMP2", "ITIH2", "COL5A3", "SELPLG", "ALDH1A1"],
    "bowen19_type2_up": ["BAG3", "MT1X", "ANGPTL4", "F3", "PDK4", "EMP1", "ADM", "YBX3",
        "IFITM3", "EFEMP1", "CEBPD", "APOLD1", "IFITM2", "DDIT4", "CD44", "FGF2", "SOX9",
        "PLPP3", "SLC16A9", "MT2A", "BMPR1B", "NTRK2", "PLSCR4", "RANBP3L" ],
    "bowen19_type2_down": [ "NPY", "RERGL", "RTN1", "NEUROD6", "ATP6V1G2", "NSF",
        "ENSA", "AP1S1", "DCLK1", "G3BP2", "MAL2",  "BEND6", "CRYM", "TAC1", "ARMCX3",
        "MAP1B", "AP1S1", "VGF", "TRIM23", "ACAT2", "PGAM1", "OPN3"],
    "lanz19_dlpfc_neuronal_down": ["GAD1","PVALB","SST","CXCL12","ATP6V1A","RGS4","CNTNAP2"],
    "cytokine_response_genes": ["IL6", "STAT3", "SOCS3", "IL1RL1","SHC1", "AKT2"],
    "microglia_activated": ["SPP1", "GPNMB", "LPL", "CST7", "LGALS3", "ITGAX", "CLEC7A", "AXL", "CD68", "CTSB", "CTSD"],
    "sz07_vascular_depleted": ["VWF", "CD34", "CLEC14A", "ESAM", "NOS3", "ABCB1", "ABCG2", "HSPA12B", "GATA2", "ARHGEF15"],
    "sz07_neuronal_enriched": ["RTN4R", "FABP3"],
    "synaptic_readout": [
        "SYP", "SYN1", "SNAP25", "STX1A", "VAMP2", "SYT1", "SV2A", "BSN", "PCLO", "RIMS1", # presynaptic
        "DLG4", "HOMER1", "SHANK1", "SHANK2" # postsynaptic
    ],
    "fibroblast_core": ["DCN", "LUM", "APOD", "COL15A1", "FBLN1"],
    "endothelial_core": ["CLDN5", "PECAM1", "CDH5", "PTPRB", "FLT1"],
    "endothelial_activation": ["ICAM1","VCAM1","SELE","SELP","ANGPT2","MCAM"],
    "bbb_transport": ["ABCB1", "ABCG2", "MFSD2A","SLC2A1", "SLC7A5"],
    "fibrotic_ecm": ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1", "COL6A2"],
    "endmt": ["FN1","VIM","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2"],
    "fibromyocyte": ["CCL19", "IGFBP5", "KCNT2"],
}

# =============================================================================
# r09* -- hspe cell-type deconvolution of the bulk RNA-seq against the canonical
#         single-cell references built by scripts_to_make_deconv_reference/.
#
# Consumes, per reference <ref> at DECONV_REF_DIR/<ref>/:
#     matrix.mtx.gz  features.tsv.gz  cells.tsv.gz  provenance.json
# Produces:
#     09_deconv/{ref}/{sample}/proportions.tsv     per sample x reference
#     09_deconv/{ref}.tsv                          summary (samples x classes)
#
# Flow:
#   r09a_mixture     bulk STAR gene counts -> symbol-space CPM (once, shared)
#   r09b_reference   canonical sc ref -> pseudobulk profiles + markers (per ref)
#   r09c_run_hspe    deconvolve one sample vs one reference (per ref x sample)
#   r09d_summary     collate all samples for a reference (per ref)
#
# Assumes the including Snakefile defines: `samples` (list) and the STAR output
# "/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab", GTF, and a
# scripts dir. Paths below mirror the main pipeline's /tmp/data layout.
# =============================================================================
import os

# ---- config / paths ----------------------------------------------------------
DECONV_REF_DIR = "/tmp/data/00_additional_files/deconv/reference_canonical" 
DECONV_OUT     = "/tmp/data/09_deconv"
R09_SCRIPTS    = "/tmp/data/00_additional_files/deconv/scripts"
# GTF for gene_id -> symbol (same annotation as alignment; see main pipe)
GTF = REF_DIR + "/gencode.v47.GRCh38.annotation.gtf" if "REF_DIR" in globals() \
      else "/tmp/data/00_additional_files/gtex_v11_refs/gencode.v47.GRCh38.annotation.gtf"

# STAR ReadsPerGene strand column, as a DATA index (0=unstranded,1=fwd,2=reverse)
# GTEx / Illumina TruSeq stranded mRNA is reverse-stranded -> 2. Override with
# -e DECONV_STRAND=... at container start.
DECONV_STRAND    = int(os.environ.get("DECONV_STRAND", "2"))
DECONV_N_MARKERS = int(os.environ.get("DECONV_N_MARKERS", "25"))

# an R image carrying hspe + Matrix + data.table (override -e DECONV_R_CONTAINER=...)
# R_HSPE_CONTAINER = os.environ.get("DECONV_R_CONTAINER",
    # "docker://bioconductor/bioconductor_docker:RELEASE_3_19")
# R_HSPE_CONTAINER = "docker://nikolayquadrat/toolbox:hspe" # R container with hspe
R_HSPE_CONTAINER = "/tmp/data/00_additional_files/deconv/hspe_container/hspe.sif" 

# ---- enable switch --------------------------------------------------------
# Deconvolution is OPTIONAL and OFF by default. It is toggled by the
# DECONV_ENABLED environment variable, set with -e at container start, exactly
# like ARCASHLA_ENABLED / CONTAMINATION_ENABLED:
#     docker ... -e DECONV_ENABLED=1 ... snakemake/snakemake:v9.16.3
# Accepts 1/true/yes/on (any case). Even when enabled, the stage stays inert
# unless DECONV_REF_DIR exists and holds at least one canonical reference
# (a subfolder with matrix.mtx.gz).
DECONV_ENABLED = os.environ.get("DECONV_ENABLED", "0").strip().lower() \
    in ("1", "true", "yes", "on")

# ---- discover references (folders that contain a canonical matrix) -----------
def _discover_references():
    if not DECONV_ENABLED:
        return []
    if not os.path.isdir(DECONV_REF_DIR):
        return []
    return sorted(d for d in os.listdir(DECONV_REF_DIR)
                  if os.path.isfile(os.path.join(DECONV_REF_DIR, d, "matrix.mtx.gz")))

DECONV_REFERENCES = _discover_references()
DECONV_ACTIVE = DECONV_ENABLED and len(DECONV_REFERENCES) > 0

# constrain wildcards so {ref}/{sample} never eat slashes or cross each other
wildcard_constraints:
    ref    = "|".join(re.escape(r) for r in DECONV_REFERENCES) or "NO_REFERENCE",
    sample = "|".join(re.escape(s) for s in samples) or "NO_SAMPLE",

# Optionally also emit CIBERSORTx-ready input files (does NOT run CIBERSORTx --
# run it yourself manually from the emitted folder). Off by default.
DECONV_CIBERSORTX = os.environ.get("DECONV_CIBERSORTX", "0").strip().lower() \
    in ("1", "true", "yes", "on")
DECONV_CSX_MAXCELLS = int(os.environ.get("DECONV_CSX_MAXCELLS", "500"))

# ---- convenience target: everything deconvolution ----------------------------
# When the stage is inactive (disabled, or no references found) this target has
# NO inputs, so `snakemake r09_deconv_all` is a clean no-op rather than an error.
if not DECONV_ACTIVE:
    if DECONV_ENABLED:
        _reason = (f"no canonical references under {DECONV_REF_DIR} "
                   "(expected subfolders each with matrix.mtx.gz)")
    else:
        _reason = "disabled (set -e DECONV_ENABLED=1 at container start to turn on)"
    print(f"[r09] deconvolution stage inactive: {_reason}")

# =============================================================================
# 11. RNA-seq genotype extraction, chip concordance and kinship  (r11a .. r11n)
# =============================================================================
# Produce a PLINK fileset of CONFIDENT genotypes called from the RNA-seq BAMs,
# in the SAME BUILD and with the SAME VARIANT IDs as an existing SNP-array
# (chip) fileset, then use it for (a) identity / swap checking and (b) kinship.
#
# -----------------------------------------------------------------------------
# THE CENTRAL INVARIANT
# -----------------------------------------------------------------------------
#   The output fileset is always in the CHIP's coordinate system, carrying the
#   CHIP's variant IDs (rsIDs where the manifest has them).
#
# The BAMs are hg38 (GENOME_FA is the TOPMed hg38 reference, GENCODE v47), but
# the chip is typically hg19. So coordinates travel in both directions:
#
#   chip .bim (hg19)  --liftOver-->  hg38 positions  -->  pileup / call in hg38
#        ^                                                        |
#        +------------- crosswalk maps calls back ----------------+
#
# Only the SITE LIST is lifted over -- never a VCF, never genotypes. Because we
# genotype at known array sites only, the liftover is a 1 bp-per-row BED job,
# and every mapped site keeps its hg19 coordinates and rsID alongside, which is
# exactly the crosswalk needed to write the result back out in hg19.
#
# Liftover artefacts are handled rather than ignored:
#   * sites that fail to lift are dropped and counted
#   * sites that lift to more than one hg38 location are dropped (ambiguous)
#   * sites that lift to an inverted region come back strand-flipped. These are
#     detected in r11e (neither chip allele matches the hg38 reference base, but
#     the complement does), counted, and the alleles are complemented back to
#     the chip's orientation in r11i. This is only safe because A/T and C/G
#     sites were already excluded, where a flip is indistinguishable from a
#     genotype swap.
#
# -----------------------------------------------------------------------------
# DESIGN DECISION -- targeted, not de novo
# -----------------------------------------------------------------------------
# We do NOT discover variants in the RNA. We genotype ONLY at the sites the chip
# already assays, with bcftools call -C alleles constraining each site to the
# chip's own REF/ALT pair. Reasons:
#   1. Concordance only ever uses the intersection; calling elsewhere is waste.
#   2. Identical alleles at identical positions means the filesets merge with no
#      post-hoc allele reconciliation.
#   3. It sidesteps de novo RNA calling problems entirely: splice-junction
#      mapping artefacts, ADAR A-to-I editing read as A>G, allele-specific
#      expression. The residual issue is LOSS of hets, which degrades
#      concordance gracefully -- it never manufactures false identity.
#
# -----------------------------------------------------------------------------
# WHAT "CONFIDENT" MEANS (all tunable, see ENV VARS)
# -----------------------------------------------------------------------------
#   site:  single-letter ACGT alleles only (indel/CNV/multi-base probes dropped);
#          biallelic; NOT strand-ambiguous (A/T and C/G excluded); lifts
#          uniquely to hg38; falls inside a union-exon interval (transcribed)
#   reads: MAPQ >= 255 (STAR unique -- see NOTE), base quality >= 20,
#          unmapped/secondary/QC-fail/duplicate excluded, BAQ off, indels off
#   call:  DP >= 10 and GQ >= 20, else set MISSING (not dropped)
#   site kept if non-missing in >= 50% of samples
#
# NOTE ON MAPQ 255 -- this pipeline's "patched.md.bam" is patched in NAME ONLY
# (see the comment above r03f_markduplicates: bamsync/MAPQ patching is skipped
# because we start from FASTQ). STAR's native convention is therefore intact:
# 255 = unique, 3 = 2 loci, 1 = 3-4 loci, 0 = >=5. So -q 255 is the correct
# unique-only filter. If a GTEx-style MAPQ 255->60 patch is ever added upstream,
# set GENOCHECK_MIN_MQ=60 or this stage silently returns zero coverage.
#
# -----------------------------------------------------------------------------
# INPUTS  (/tmp/data/00_additional_files/genetics/)   -- basenames are free-form
# -----------------------------------------------------------------------------
#   <anything>.bim              REQUIRED. The .bim is the anchor: it is found by
#                               glob, so the fileset can be called whatever you
#                               like. Several .bim files in the directory means
#                               GENOCHECK_CHIP_PREFIX must say which to use.
#   <same>.bed / <same>.fam     OPTIONAL. Present -> the IDENTITY tier runs
#                               (merge + kinship + similarity report). Absent ->
#                               SITES tier: the .bim still supplies the site
#                               list, so RNA genotypes, the PLINK deliverable
#                               and RNA-vs-RNA kinship are produced, but no
#                               identity check is possible.
#   hg19ToHg38.over.chain.gz    REQUIRED when the array is hg19. Found by that
#                               name, or as the single *.chain.gz in the dir.
#   sample_map.tsv  OPTIONAL    2 cols, no header: rna_sample  array_id
#                               A HYPOTHESIS TO TEST, never an assumption --
#                               matching is always done from the genotypes, and
#                               r11q reports agreement or disagreement with it.
#   donor_map.tsv   OPTIONAL    2 cols, no header: rna_sample  donor_id
#                               Lets r11l flag expected same-donor RNA pairs.
#
# -----------------------------------------------------------------------------
# OUTPUTS  (/tmp/data/11_genocheck/)
# -----------------------------------------------------------------------------
#   rnaseq_genotypes.{bed,bim,fam}   <- the deliverable, in ARRAY build
#   rnaseq_genotypes.hg38.vcf.gz     the raw calls, hg38, for reference
#   00_genocheck_site_stats.tsv      per-sample call rate / depth / het
#   00_genocheck_crosswalk.tsv.gz    hg19 <-> hg38 site table actually used
#   00_kinship_rna.kin0              KING, RNA samples only
#   00_kinship_rna_report.tsv        RNA pairs, labelled (donor structure)
#   -- IDENTITY tier only --
#   00_kinship_merged.kin0           KING, RNA + array together
#   00_similarity_pairs.tsv          every pair: concordance, kinship, IBS0
#   00_similarity_report.tsv         one row per RNA sample: best array match,
#                                    margin to runner-up, verdict
#
# ENV VARS (same convention as CONTAMINATION_ENABLED / DECONV_ENABLED)
#   GENOCHECK_ENABLED      auto (default) | 1 | 0   auto = on if a .bim is found
#   GENOCHECK_IDENTITY     auto (default) | 0        0 forces the SITES tier
#   GENOCHECK_DIR          default /tmp/data/00_additional_files/genetics
#   GENOCHECK_CHIP_PREFIX  fileset path without extension; only needed when the
#                          directory holds more than one .bim
#   GENOCHECK_CHIP_BUILD   hg19 (default) | hg38    build of the array fileset
#   GENOCHECK_CHAIN        chain file path (else auto-discovered in the dir)
#   GENOCHECK_MIN_DP       default 10     per-genotype depth
#   GENOCHECK_MIN_GQ       default 20     per-genotype quality
#   GENOCHECK_MIN_MQ       default 255    read mapping quality (STAR unique)
#   GENOCHECK_MIN_BQ       default 20     base quality
#   GENOCHECK_MAX_DEPTH    default 10000  mpileup per-file depth cap
#   GENOCHECK_MAX_MISSING  default 0.5    drop site if more missing than this
#   GENOCHECK_KEEP_DUPS    default 0      1 = do NOT filter duplicate reads
#   GENOCHECK_INCLUDE_X    default 0      1 = also genotype chrX
#   GENOCHECK_LD_PRUNE     default 0      1 = LD-prune before kinship
#   GENOCHECK_DUP_KINSHIP  default 0.35   kinship above this = same individual
#   GENOCHECK_MATCH_CONC   default 0.85   concordance for a `match` verdict
#   GENOCHECK_MAX_IBS0     default 0.005  max IBS0 rate for a `match` verdict
#   GENOCHECK_MIN_MARGIN   default 0.05   min concordance gap to the runner-up
# =============================================================================

# -----------------------------------------------------------------------------
# Genotype fileset discovery
# -----------------------------------------------------------------------------
# The array fileset lives in GENETICS_DIR under ANY basename; the .bim is the
# anchor. Two tiers of capability, decided by what is actually present:
#
#   .bim only                 -> SITES tier. The .bim supplies the site list, so
#                                RNA genotypes, the PLINK deliverable and the
#                                RNA-vs-RNA kinship (donor structure) are all
#                                produced. No identity check is possible -- there
#                                are no array genotypes to compare against.
#   .bim + .bed + .fam        -> IDENTITY tier. Everything above, plus the merge,
#                                the cross-dataset kinship and the similarity
#                                report.
#
# Ambiguity is never resolved silently: several .bim files in the directory and
# no GENOCHECK_CHIP_PREFIX disables the stage with the candidate list printed
# (or errors, if the stage was force-enabled).
# -----------------------------------------------------------------------------
GENETICS_DIR = os.environ.get("GENOCHECK_DIR",
                              "/tmp/data/00_additional_files/genetics")
GENOCHECK_OUT = "/tmp/data/11_genocheck"
GENOCHECK_SAMPLE_MAP = GENETICS_DIR + "/sample_map.tsv"
GENOCHECK_DONOR_MAP  = GENETICS_DIR + "/donor_map.tsv"


def _genocheck_find_prefix(d):
    """(prefix, note, candidates). Basename is free-form; *.bim is the anchor."""
    import glob as _glob
    env = os.environ.get("GENOCHECK_CHIP_PREFIX", "").strip()
    if env:
        pre = env[:-4] if env.endswith(".bim") else env
        return pre, "set via GENOCHECK_CHIP_PREFIX", []
    bims = sorted(_glob.glob(os.path.join(d, "*.bim")))
    if len(bims) == 1:
        return bims[0][:-4], f"auto-discovered in {d}", bims
    if len(bims) > 1:
        return "", f"{len(bims)} .bim files in {d} -- ambiguous", bims
    return "", f"no *.bim in {d}", []


def _genocheck_find_chain(d):
    """Chain file for hg19->hg38: env > conventional name > a single *.chain*."""
    import glob as _glob
    env = os.environ.get("GENOCHECK_CHAIN", "").strip()
    if env:
        return env
    conventional = os.path.join(d, "hg19ToHg38.over.chain.gz")
    if os.path.isfile(conventional):
        return conventional
    hits = sorted(_glob.glob(os.path.join(d, "*.chain.gz"))
                  + _glob.glob(os.path.join(d, "*.chain")))
    return hits[0] if len(hits) == 1 else conventional


GENOCHECK_CHIP_PREFIX, _chip_note, _chip_candidates = \
    _genocheck_find_prefix(GENETICS_DIR)
GENOCHECK_CHAIN = _genocheck_find_chain(GENETICS_DIR)


def _nonempty(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


_have_bim = bool(GENOCHECK_CHIP_PREFIX) and _nonempty(GENOCHECK_CHIP_PREFIX + ".bim")
_have_full = _have_bim and _nonempty(GENOCHECK_CHIP_PREFIX + ".bed") \
    and _nonempty(GENOCHECK_CHIP_PREFIX + ".fam")

_genocheck_override = os.environ.get("GENOCHECK_ENABLED", "auto").lower()
if _genocheck_override in ("0", "false", "no", "off"):
    GENOCHECK_ENABLED = False
elif _genocheck_override in ("1", "true", "yes", "on"):
    GENOCHECK_ENABLED = True
elif _genocheck_override == "auto":
    GENOCHECK_ENABLED = _have_bim
else:
    print(f"[genocheck] unrecognised GENOCHECK_ENABLED='{_genocheck_override}', "
          f"falling back to auto-detection", file=sys.stderr)
    GENOCHECK_ENABLED = _have_bim

# Identity tier needs real array genotypes, not just the site list.
_identity_override = os.environ.get("GENOCHECK_IDENTITY", "auto").lower()
if _identity_override in ("0", "false", "no", "off"):
    GENOCHECK_IDENTITY = False
else:
    GENOCHECK_IDENTITY = GENOCHECK_ENABLED and _have_full

GENOCHECK_CHIP_BUILD  = os.environ.get("GENOCHECK_CHIP_BUILD", "hg19").strip().lower()
GENOCHECK_MIN_DP      = int(os.environ.get("GENOCHECK_MIN_DP", "10"))
GENOCHECK_MIN_GQ      = int(os.environ.get("GENOCHECK_MIN_GQ", "20"))
GENOCHECK_MIN_MQ      = int(os.environ.get("GENOCHECK_MIN_MQ", "255"))
GENOCHECK_MIN_BQ      = int(os.environ.get("GENOCHECK_MIN_BQ", "20"))
GENOCHECK_MAX_DEPTH   = int(os.environ.get("GENOCHECK_MAX_DEPTH", "10000"))
GENOCHECK_MAX_MISSING = float(os.environ.get("GENOCHECK_MAX_MISSING", "0.5"))
GENOCHECK_KEEP_DUPS   = os.environ.get("GENOCHECK_KEEP_DUPS", "0").strip().lower() \
    in ("1", "true", "yes", "on")
GENOCHECK_INCLUDE_X   = os.environ.get("GENOCHECK_INCLUDE_X", "0").strip().lower() \
    in ("1", "true", "yes", "on")
GENOCHECK_LD_PRUNE    = os.environ.get("GENOCHECK_LD_PRUNE", "0").strip().lower() \
    in ("1", "true", "yes", "on")
GENOCHECK_DUP_KINSHIP = float(os.environ.get("GENOCHECK_DUP_KINSHIP", "0.35"))
# Similarity thresholds for the r11q verdict. Heuristics, not hard science --
# see the interpretation note in the KINSHIP header.
GENOCHECK_MATCH_CONC  = float(os.environ.get("GENOCHECK_MATCH_CONC", "0.85"))
GENOCHECK_MAX_IBS0    = float(os.environ.get("GENOCHECK_MAX_IBS0", "0.005"))
GENOCHECK_MIN_MARGIN  = float(os.environ.get("GENOCHECK_MIN_MARGIN", "0.05"))

# hg38 contigs we genotype on (the pileup space)
GENOCHECK_CHROMS = [f"chr{i}" for i in range(1, 23)]
if GENOCHECK_INCLUDE_X:
    GENOCHECK_CHROMS.append("chrX")

GENOCHECK_NEEDS_LIFT = GENOCHECK_CHIP_BUILD in ("hg19", "grch37", "b37", "37")

# Duplicates excluded by default. For RNA a "duplicate" at a highly expressed
# gene is often a genuinely independent molecule, so dropping them costs depth,
# but keeping them lets one over-amplified original molecule masquerade as a
# high-confidence homozygote. For identity work, losing depth is the cheaper error.
GENOCHECK_FF = "UNMAP,SECONDARY,QCFAIL" if GENOCHECK_KEEP_DUPS \
    else "UNMAP,SECONDARY,QCFAIL,DUP"

# Container tags: VERIFY these resolve before the first run (`singularity pull`);
# biocontainers build-hash suffixes drift. If your nikolayquadrat/toolbox image
# carries bcftools / plink / plink2 / CrossMap, point these at it instead.
BCFTOOLS_CONTAINER = os.environ.get(
    "BCFTOOLS_CONTAINER", "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0")
PLINK2_CONTAINER = os.environ.get(
    "PLINK2_CONTAINER", "docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0")
PLINK1_CONTAINER = os.environ.get(
    "PLINK1_CONTAINER", "docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1")
CROSSMAP_CONTAINER = os.environ.get(
    "CROSSMAP_CONTAINER", "docker://quay.io/biocontainers/crossmap:0.7.0--pyhdfd78af_0")

if GENOCHECK_ENABLED and not _have_bim:
    if _chip_candidates:
        sys.exit("[genocheck] " + _chip_note + ":\n"
                 + "\n".join("              " + b for b in _chip_candidates)
                 + "\n            Set GENOCHECK_CHIP_PREFIX to the one to use "
                   "(path without the .bim suffix).")
    sys.exit(f"[genocheck] enabled but no usable .bim found ({_chip_note})")

if GENOCHECK_ENABLED:
    if GENOCHECK_CHIP_BUILD not in ("hg19", "grch37", "b37", "37",
                                    "hg38", "grch38", "b38", "38"):
        sys.exit(f"[genocheck] GENOCHECK_CHIP_BUILD='{GENOCHECK_CHIP_BUILD}' "
                 "not understood; use hg19 or hg38")
    if GENOCHECK_NEEDS_LIFT and not os.path.isfile(GENOCHECK_CHAIN):
        sys.exit(
            f"[genocheck] chip build is {GENOCHECK_CHIP_BUILD}, so the site list "
            f"must be lifted to hg38, but no chain file at:\n    {GENOCHECK_CHAIN}\n"
            "            Download hg19ToHg38.over.chain.gz from UCSC "
            "(hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) and stage it in "
            f"{GENETICS_DIR}, or set GENOCHECK_CHAIN.")
    print(f"[genocheck] enabled; fileset={GENOCHECK_CHIP_PREFIX} ({_chip_note})",
          file=sys.stderr)
    print(f"[genocheck] tier: "
          + ("IDENTITY (.bed/.bim/.fam present -- merge, kinship and "
             "similarity report will run)" if GENOCHECK_IDENTITY else
             "SITES ONLY (.bim without .bed/.fam -- RNA genotypes and RNA-vs-RNA "
             "kinship only, no identity check)"), file=sys.stderr)
    print(f"[genocheck] build={GENOCHECK_CHIP_BUILD} "
          f"lift={'hg19->hg38' if GENOCHECK_NEEDS_LIFT else 'none'}; output "
          f"fileset will be in {GENOCHECK_CHIP_BUILD} coords with chip variant IDs",
          file=sys.stderr)
    print(f"[genocheck] filters: MQ>={GENOCHECK_MIN_MQ} BQ>={GENOCHECK_MIN_BQ} "
          f"DP>={GENOCHECK_MIN_DP} GQ>={GENOCHECK_MIN_GQ} "
          f"max_missing={GENOCHECK_MAX_MISSING} ff={GENOCHECK_FF}",
          file=sys.stderr)
    for _lbl, _p in (("sample map (hypothesis to TEST, not assumed)",
                      GENOCHECK_SAMPLE_MAP),
                     ("donor map", GENOCHECK_DONOR_MAP)):
        if _nonempty(_p):
            print(f"[genocheck] {_lbl}: {_p}", file=sys.stderr)
else:
    print(f"[genocheck] disabled ({_chip_note}, or forced off)", file=sys.stderr)


# -----------------------------------------------------------------------------
# Small shared helpers (module level so several rules can use them)
# -----------------------------------------------------------------------------

def _genocheck_read_map(path):
    """Two-column whitespace file -> dict. Missing/empty file -> {}."""
    out = {}
    if not (os.path.isfile(path) and os.path.getsize(path) > 0):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and not f[0].startswith("#"):
                out[f[0]] = f[1]
    return out


def _genocheck_read_kin0(path):
    """
    Read a plink2 --make-king-table .kin0, normalising the column names, which
    differ between plink2 builds (#FID1/#IID1 vs ID1/ID2, KINSHIP vs Kinship).
    The resolved names are returned in df.attrs: 'kinship', 'ibs0', 'nsnp'
    ('ibs0'/'nsnp' may be None if that build does not emit them).
    """
    import pandas as pd
    df = pd.read_csv(path, sep=r"\s+")
    ren = {"#FID1": "FID1", "#IID1": "IID1", "#ID1": "IID1",
           "ID1": "IID1", "ID2": "IID2"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "IID1" not in df.columns or "IID2" not in df.columns:
        raise ValueError(f"{path}: cannot find the two sample-ID columns; "
                         f"found {list(df.columns)}")

    def pick(*names):
        for n in names:
            if n in df.columns:
                return n
        return None

    kinship = pick("KINSHIP", "Kinship", "kinship")
    if kinship is None:
        raise ValueError(f"{path}: no KINSHIP column; found {list(df.columns)}")
    df.attrs["kinship"] = kinship
    ibs0 = pick("IBS0", "Ibs0", "ibs0")
    nsnp = pick("NSNP", "N_SNP", "Nsnp", "nsnp")
    df.attrs["ibs0"] = ibs0
    df.attrs["nsnp"] = nsnp

    # IBS0 is a PROPORTION in KING-style .kin0 output (plink2 mimics KING), but
    # some builds/tools emit a count. Decide from the data rather than assuming:
    # an IBS0 count over tens of thousands of SNPs exceeds 1 for any unrelated
    # pair, while a proportion never does. If every pair is a duplicate both
    # readings give ~0, so the ambiguity is harmless.
    if ibs0 is not None:
        mx = pd.to_numeric(df[ibs0], errors="coerce").max()
        if nsnp is not None and pd.notna(mx) and mx > 1.0:
            df["ibs0_rate"] = pd.to_numeric(df[ibs0], errors="coerce") / \
                pd.to_numeric(df[nsnp], errors="coerce").replace(0, pd.NA)
            df.attrs["ibs0_scale"] = "count/NSNP"
        else:
            df["ibs0_rate"] = pd.to_numeric(df[ibs0], errors="coerce")
            df.attrs["ibs0_scale"] = "proportion"
    return df





# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

rule all:
    input:
        # STAR BAMs ready for DROP, plus splice-junction files
        expand("/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam", sample=samples),
        # expand("/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai", sample=samples),
        expand("/tmp/data/03_bam_star/{sample}/{sample}.SJ.out.tab", sample=samples),
        expand("/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab", sample=samples),
        # Mapping stat summary
        "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
        # QC: per-sample Picard CollectRnaSeqMetrics + cohort summary table
        # expand("/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt", sample=samples),
        "/tmp/data/04_qc/00_qc_summary.tsv",
        # Optional standalone contamination summary (wide TSV, one row per
        # sample, one column per species in species.txt order). Independent
        # of the QC summary -- same data but no percentages, flags, or
        # Picard merging. Gated on CONTAMINATION_ENABLED.
        *(["/tmp/data/05_contamination/00_contamination_summary.tsv"]
          if CONTAMINATION_ENABLED else []),
        # Optional per-species coverage uniformity (breadth-of-coverage,
        # fraction 0-1 of genome covered by >=1 read). Distinguishes real
        # presence (high read count + high uniformity) from off-target
        # alignment artifacts (high read count + low uniformity).
        *(["/tmp/data/05_contamination/00_contamination_uniformity.tsv"]
          if CONTAMINATION_ENABLED else []),
        # Optional per-species read fractions (of the post-entropy-filter
        # read pool fed to bwa-mem2). Same wide format as the summary.
        *(["/tmp/data/05_contamination/00_contamination_fractions.tsv"]
          if CONTAMINATION_ENABLED else []),
        # Optional HLA class I typing (only if the arcasHLA reference is present;
        # see ARCASHLA_ENABLED gating above).
        *(["/tmp/data/06_hla/00_hla_summary.tsv"] if ARCASHLA_ENABLED else []),
        # Optional HLA class II typing (DRB1/DQA1/DQB1/DPA1/DPB1). Same gating;
        # class II calls are expected to be sparse in brain RNA-seq (expressed
        # mainly on microglia), so this table will have empty cells for
        # low-expression samples.
        *(["/tmp/data/06_hla/00_hla_summary_classII.tsv"] if ARCASHLA_ENABLED else []),
        "/tmp/data/04_qc/00_inferred_sex.tsv",
        "/tmp/data/04_qc/00_cell_marker_expression.tsv",
        "/tmp/data/04_qc/00_expression_summary.xlsx",
        ((expand(f"{DECONV_OUT}/{{ref}}.tsv", ref=DECONV_REFERENCES)
          + expand(f"{DECONV_OUT}/{{ref}}/{{sample}}/proportions.tsv",
                   ref=DECONV_REFERENCES, sample=samples)) if DECONV_ACTIVE else [])
        + ((expand(f"{DECONV_OUT}/_cibersortx/{{ref}}/refsample.txt",
                   ref=DECONV_REFERENCES)
            + expand(f"{DECONV_OUT}/_cibersortx/{{ref}}/mixture.txt",
                     ref=DECONV_REFERENCES))
           if (DECONV_ACTIVE and DECONV_CIBERSORTX) else []),
        *([GENOCHECK_OUT + "/rnaseq_genotypes.bed",
            GENOCHECK_OUT + "/00_genocheck_site_stats.tsv",
            GENOCHECK_OUT + "/00_kinship_rna_report.tsv"] if GENOCHECK_ENABLED else []),
        *([GENOCHECK_OUT + "/00_similarity_report.tsv",
            GENOCHECK_OUT + "/00_similarity_pairs.tsv"] if GENOCHECK_IDENTITY else []),
    shell: "echo 'GTEx-V11-compatible alignment + QC complete.'"

# -----------------------------------------------------------------------------
# Trimming
# -----------------------------------------------------------------------------
# Note: GTEx itself does NOT trim. They feed raw FASTQ to STAR. Trimming is
# generally OK for STAR but produces slightly different alignments than GTEx.
# If you want maximum GTEx fidelity, set USE_TRIMMING = False below and edit
# rule r03d_STAR_mapping_GTEx to read directly from the un-trimmed FASTQ.
# For DROP use, trimmed-vs-untrimmed does not significantly affect outlier
# detection, so trimming is left enabled by default.

rule r02a_trim_galore:
    input:
        get_fastq_path,
    output:
        R1 = temp("/tmp/data/02_trimmed_fastq/{sample}_R1_val_1.fq.gz"),
    singularity: TRIM_CONTAINER
    shell:
        r"""
        mkdir -p /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
        n=$(perl -e 'if (qq|{input}| =~ /^(\S+)\s\1$/) {{print 1}} else {{print 0}}')
        if (($n==1)); then
            trim_galore {input[0]} -o /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
            sleep 5
            file_to_move=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                             -type f -name "*trimmed\.fq\.gz*")
            mv "$file_to_move" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R1_val_1.fq.gz
        else
            trim_galore --paired {input[0]} {input[1]} \
                -o /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
            sleep 5
            file_to_move1=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                              -type f -name "*_val_1\.fq\.gz*")
            mv "$file_to_move1" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R1_val_1.fq.gz
            file_to_move2=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                              -type f -name "*_val_2\.fq\.gz*")
            mv "$file_to_move2" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R2_val_2.fq.gz
        fi
        rm -rf /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
        """


# -----------------------------------------------------------------------------
# STAR index, GTEx V11 parameters
# -----------------------------------------------------------------------------

rule r03a_STAR_index:
    """
    Build the STAR index using GENCODE v47 annotation and sjdbOverhang
    matching YOUR read length (read_length - 1).
    GTEx itself uses sjdbOverhang=75 because their reads are 2x76bp.
    """
    input:
        genome = GENOME_FA,
        gtf    = GTF,
    output:
        checkpoint = STAR_INDEX + "/SA",
    params:
        index_dir = STAR_INDEX,
        sjdb_overhang = SJDB_OVERHANG,
    threads: 16
    singularity: STAR_CONTAINER
    shell:
        """
        mkdir -p {params.index_dir}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --runThreadN {threads}
        """


# -----------------------------------------------------------------------------
# STAR alignment with GTEx V11 parameters
# -----------------------------------------------------------------------------

rule r03d_STAR_mapping_GTEx:
    """
    STAR alignment using exactly the parameter set from the GTEx V11 pipeline
    (Broad's run_STAR.py defaults). The flags below come verbatim from
    https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/src/run_STAR.py
    Together with --twopassMode Basic, these are what GTEx V11 used to
    produce the count matrices on the GTEx portal. Reproducing them exactly
    means our gene-level counts can be merged into a unified expression
    matrix with GTEx samples for DROP/OUTRIDER reference panel expansion.

    Notable parameters:
      --twopassMode Basic              : two-pass alignment for better novel-junction detection
      --outFilterMultimapNmax 20       : up to 20 multi-mapping locations per read
      --outFilterMismatchNoverLmax 0.1 : <=10% mismatches per read length
      --outFilterType BySJout          : only retain junctions supported in 2nd pass
      --quantMode TranscriptomeSAM GeneCounts : output transcriptome BAM + STAR gene counts
      --chimSegmentMin 15              : enable chimeric detection (TOPMed/GTEx convention)
      --outSAMtype BAM SortedByCoordinate : GTEx writes a coordinate-sorted BAM directly
      --outSAMunmapped Within          : keep unmapped reads in the BAM (GTEx convention)
    """
    input:
        R1 = "/tmp/data/02_trimmed_fastq/{sample}_R1_val_1.fq.gz" if USE_TRIMMING
              else lambda wc: get_fastq_path(wc)[0],
        index_checkpoint = ancient(STAR_INDEX + "/SA"),
    output:
        bam_sorted = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bam_transcripts = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
        gene_counts = "/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab",
        sj_tab = "/tmp/data/03_bam_star/{sample}/{sample}.SJ.out.tab",
        chimeric = "/tmp/data/03_bam_star/{sample}/{sample}.Chimeric.out.junction",
        log_final = "/tmp/data/03_bam_star/{sample}/{sample}.Log.final.out",
    params:
        out_prefix = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/{wc.sample}.",
        out_dir    = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/",
        tmp_dir    = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/_STAR_tmp",
        star_index = STAR_INDEX,
        # Read group line for downstream tools. Replace SM tag with sample id;
        # ID is set to the sample id as well for simplicity. PL=ILLUMINA, LB=lib1.
        rg_line = lambda wc: f"ID:{wc.sample} SM:{wc.sample} PL:ILLUMINA LB:lib1",
        r2_path = lambda wc: f"/tmp/data/02_trimmed_fastq/{wc.sample}_R2_val_2.fq.gz",
    threads: 8
    singularity: STAR_CONTAINER
    shell:
        r"""
        # STAR refuses to write into a non-empty tmp dir; clean before running
        rm -rf {params.tmp_dir} || true
        mkdir -p {params.out_dir}

        # Decide PE vs SE based on whether the R2 file exists
        if [ -f "{params.r2_path}" ]; then
            R2_ARG="{params.r2_path}"
        else
            R2_ARG=""
        fi

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.star_index} \
            --twopassMode Basic \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterType BySJout \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --limitSjdbInsertNsj 1200000 \
            --readFilesIn {input.R1} $R2_ARG \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --outTmpDir {params.tmp_dir} \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs None \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --genomeLoad NoSharedMemory \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType Junctions WithinBAM SoftClip \
            --chimMainSegmentMultNmax 1 \
            --outSAMattributes NH HI AS nM NM ch \
            --outSAMattrRGline {params.rg_line}

        rm -rf {params.tmp_dir} || true
        """


# -----------------------------------------------------------------------------
# Index BAM (samtools)
# -----------------------------------------------------------------------------

rule r03e_index_bam:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    threads: 4
    singularity: SAMTOOLS_CONTAINER
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"


# -----------------------------------------------------------------------------
# Picard MarkDuplicates
# -----------------------------------------------------------------------------
# GTEx applies MarkDuplicates AFTER coordinate-sorting from STAR. The output
# is the canonical "<sample>.Aligned.sortedByCoord.out.patched.md.bam" -
# this is the file that gets fed to RNA-SeQC and RSEM in the GTEx pipeline,
# and it is also the file DROP MAE / OUTRIDER / FRASER expect as input.
#
# Note: bamsync (which carries QC flags from the original GTEx CRAM) is NOT
# applied here because we are starting from FASTQ, not from existing BAM.
# The "patched" suffix is kept in the filename only for naming consistency
# with the GTEx convention.

rule r03f_markduplicates:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        metrics = "/tmp/data/03_bam_star/{sample}/{sample}.markdup_metrics.txt",
    threads: 4
    singularity: PICARD_CONTAINER
    shell:
        """
        java -jar /usr/picard/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            ASSUME_SORT_ORDER=coordinate \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            VALIDATION_STRINGENCY=LENIENT
        """


rule r03g_index_md_bam:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
    output:
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
    threads: 4
    singularity: SAMTOOLS_CONTAINER
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"


# -----------------------------------------------------------------------------
# Aggregate mapping stats
# -----------------------------------------------------------------------------

rule r03h_mapping_stat:
    input:
        log_files = expand("/tmp/data/03_bam_star/{sample}/{sample}.Log.final.out",
                           sample=samples),
    output:
        "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
    run:
        with open(output[0], "w") as fout:
            fout.write("Sample\tInput_reads\tUniquely_mapped\tUniquely_mapped_pct\t"
                       "Multi_mapped\tMulti_mapped_pct\tToo_many_loci\tUnmapped\n")
            for log_file in input.log_files:
                sample = os.path.basename(log_file).replace(".Log.final.out", "")
                metrics = {"input_reads": "", "unique": "", "unique_pct": "",
                           "multi": "", "multi_pct": "", "too_many": "", "unmapped": ""}
                with open(log_file) as f:
                    for line in f:
                        line = line.strip()
                        if "Number of input reads" in line:
                            metrics["input_reads"] = line.split("|")[1].strip()
                        elif "Uniquely mapped reads number" in line:
                            metrics["unique"] = line.split("|")[1].strip()
                        elif "Uniquely mapped reads %" in line:
                            metrics["unique_pct"] = line.split("|")[1].strip()
                        elif "Number of reads mapped to multiple loci" in line:
                            metrics["multi"] = line.split("|")[1].strip()
                        elif "% of reads mapped to multiple loci" in line:
                            metrics["multi_pct"] = line.split("|")[1].strip()
                        elif "% of reads mapped to too many loci" in line:
                            metrics["too_many"] = line.split("|")[1].strip()
                        elif "% of reads unmapped: too short" in line:
                            metrics["unmapped"] = line.split("|")[1].strip()
                fout.write(f"{sample}\t{metrics['input_reads']}\t"
                           f"{metrics['unique']}\t{metrics['unique_pct']}\t"
                           f"{metrics['multi']}\t{metrics['multi_pct']}\t"
                           f"{metrics['too_many']}\t{metrics['unmapped']}\n")


# -----------------------------------------------------------------------------
# Picard CollectRnaSeqMetrics: RNA-seq QC including 3' bias
# -----------------------------------------------------------------------------
# Picard CollectRnaSeqMetrics reports several metrics critical for catching
# RNA degradation and library-prep issues:
#
#   - MEDIAN_3PRIME_BIAS / MEDIAN_5PRIME_BIAS / MEDIAN_5PRIME_TO_3PRIME_BIAS:
#     coverage uniformity along transcripts. A high 3' bias (>0.5 typically,
#     or much higher than the cohort median) indicates RNA degradation -
#     polyA-selected libraries lose 5' coverage as RNA degrades because
#     fragments without a polyA tail are lost preferentially.
#   - PCT_RIBOSOMAL_BASES: rRNA contamination. Should be <5% in a well-
#     prepared polyA-selected library; anything >10% indicates poor depletion.
#   - PCT_CODING_BASES / PCT_UTR_BASES / PCT_INTRONIC_BASES / PCT_INTERGENIC_BASES:
#     where reads come from. Low coding+UTR fraction (<60% combined) suggests
#     genomic DNA contamination or unspliced pre-mRNA dominance, both of which
#     can confound expression analysis.
#   - PCT_MRNA_BASES: percent of reads in coding+UTR (the "useful" fraction
#     for differential expression).
#
# Why this matters for cell-type composition: RNA degradation can selectively
# affect long transcripts (e.g., GFAP at ~10 kb) more than short ones,
# producing apparent "depletion" of cell types whose markers happen to be
# long transcripts. A high 3' bias on a sample with apparently missing
# astrocyte signal (GFAP, AQP4 are both long) would explain the pattern
# without needing a real biological loss-of-cell-type story.
# -----------------------------------------------------------------------------

rule r04a_make_refflat:
    """
    Convert GENCODE v47 GTF to UCSC refFlat format, which is what Picard
    CollectRnaSeqMetrics requires. We use ucsc-gtfToGenePred (with
    -genePredExt to retain gene_id) and then reformat columns to refFlat.

    refFlat columns: geneName, name (transcript), chrom, strand, txStart,
    txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds.
    genePredExt has these in different order; the awk fixes that.
    """
    input:
        gtf = GTF,
    output:
        refflat = REF_DIR + "/gencode.v47.GRCh38.refFlat.txt",
    singularity: "docker://quay.io/biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0"
    shell:
        r"""
        # gtfToGenePred -genePredExt produces 15-column genePredExt; refFlat
        # wants 11 columns with geneName as the FIRST column. We swap and trim.
        gtfToGenePred -genePredExt -geneNameAsName2 {input.gtf} /dev/stdout \
            | awk -v OFS='\t' '{{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \
            > {output.refflat}
        """


rule r04b_make_rrna_intervals:
    """
    Build a Picard interval-list of ribosomal RNA loci from the GENCODE GTF.
    Picard CollectRnaSeqMetrics uses this to compute PCT_RIBOSOMAL_BASES.

    Selects gene_type "rRNA" and "Mt_rRNA" entries. The interval-list format
    is a SAM-like header (from the BAM) plus one BED-like line per region.
    We extract the header from any sample BAM (they all share the same
    reference, so any one works).
    """
    input:
        gtf = GTF,
        # Borrow the header from the first sample's STAR BAM. The header
        # contains @SQ lines that identify the reference contigs - required
        # for a valid interval-list file.
        ref_bam = expand(
            "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            sample=samples[:1]),
    output:
        rrna_list = REF_DIR + "/gencode.v47.GRCh38.rRNA.interval_list",
    singularity: SAMTOOLS_CONTAINER
    shell:
        r"""
        # 1. Header from the BAM (only @HD and @SQ lines; @PG lines from STAR
        #    confuse Picard's interval-list parser)
        samtools view -H {input.ref_bam[0]} \
            | grep -E '^@HD|^@SQ' > {output.rrna_list}

        # 2. rRNA gene loci: from the GTF, gene_type "rRNA" or "Mt_rRNA",
        #    feature == "gene", emit chrom, start, end, strand, gene_id.
        #    Picard interval-list uses 1-based coordinates so no -1 on start.
        awk 'BEGIN{{OFS="\t"}} \
            $3=="gene" && /gene_type "rRNA"|gene_type "Mt_rRNA"/ {{ \
                match($0, /gene_id "[^"]+"/); \
                gid=substr($0, RSTART+9, RLENGTH-10); \
                print $1, $4, $5, $7, gid \
            }}' {input.gtf} >> {output.rrna_list}
        """


rule r04c_picard_rnaseq_metrics:
    """
    Per-sample CollectRnaSeqMetrics. Strand-specificity could be set to NONE
    update this to SECOND_READ_TRANSCRIPTION_STRAND (TruSeq dUTP/Stranded mRNA) or
    FIRST_READ_TRANSCRIPTION_STRAND. Wrong strand setting inflates the
    PCT_CORRECT_STRAND_READS metric but does NOT affect 3'/5' bias - so
    NONE is safe for primary QC purpose.
    """
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        refflat = REF_DIR + "/gencode.v47.GRCh38.refFlat.txt",
        rrna_list = REF_DIR + "/gencode.v47.GRCh38.rRNA.interval_list",
    output:
        metrics = "/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt",
    threads: 2
    singularity: PICARD_CONTAINER
    shell:
        """
        mkdir -p /tmp/data/04_qc/{wildcards.sample}
        java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
            I={input.bam} \
            O={output.metrics} \
            REF_FLAT={input.refflat} \
            RIBOSOMAL_INTERVALS={input.rrna_list} \
            STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORTED=true
        """

# =============================================================================
# Optional contamination screen
# =============================================================================
# Rules:
#   r05_prep_contamination_refs : (only when CONTAMINATION_AUTOBUILD)
#       Downloads each species' genome+GTF from NCBI using `datasets` and
#       concatenates them into contamination.fa / contamination.gtf /
#       gtf_seqnames.tsv.
#   r05a0_mask_conserved_regions: hard-mask rRNA/tRNA regions in
#                                 contamination.fa using bedtools
#                                 maskfasta. These regions cross-map to
#                                 host rRNA/mitochondrial RNA at 100%
#                                 identity regardless of alignment
#                                 stringency, so masking them is the only
#                                 way to eliminate that noise source.
#   r05a_contamination_index    : builds the bwa-mem2 index (run once)
#   r05b_extract_unmapped       : per-sample BAM -> unmapped.fq.gz
#   r05c_entropy_filter         : per-sample bbduk entropy filter to drop
#                                 low-complexity reads (poly-A, simple
#                                 repeats, etc.) before alignment. This is
#                                 the strongest single intervention
#                                 against the cross-mapping false-positive
#                                 floor.
#   r05d_contamination_screen   : per-sample bwa-mem2 alignment of the
#                                 entropy-filtered reads to the
#                                 contamination index, plus per-species
#                                 read attribution via the seqname map.
#
# r04d_qc_summary depends on r05d's per-sample counts via the conditional
# `unpack(...)` input, and folds the contamination metrics into the
# single cohort QC TSV (no separate summary file).
#
# Contig naming: contamination.fa uses the native NCBI seqnames
# (e.g. ">NC_010163.1 Acholeplasma laidlawii PG-8A"). The seqname->species
# mapping comes from gtf_seqnames.tsv (header: species_id<TAB>seqname);
# r05_prep_contamination_refs generates this file in auto-build mode by
# reading FASTA headers of each downloaded genome.

if CONTAMINATION_ENABLED:

    if CONTAMINATION_AUTOBUILD:

        rule r05_prep_contamination_refs:
            """
            On-demand contamination-reference builder. Reads species.txt
            (tab-separated: species_name<TAB>assembly_accession), downloads
            each genome+GTF from NCBI using the `datasets` CLI, then
            concatenates them into the three files r05a/r05b/r04d expect.

            Some NCBI assemblies, especially fungi/environmental genomes, provide
            genome FASTA but no GTF. For contamination screening we only need a
            syntactically valid annotation file, so generate one mock whole-contig
            gene/transcript/exon per contig. Do not interpret these as real genes.

            Internet access is required at run time (NCBI download). The
            species.txt file should be small (a few dozen species at
            most); the per-species download is ~3-50 MB and serialized,
            so the total wall-clock for a 10-species panel is ~2-5 min.

            Idempotent re-runs: if all three output files exist and are
            non-empty, Snakemake will not re-fire this rule. To force a
            rebuild after editing species.txt, delete the three output
            files.
            """
            input:
                species_list = CONTAM_SPECIES_LIST,
            output:
                fa       = CONTAM_FA,
                gtf      = CONTAM_GTF,
                seqnames = CONTAM_SEQNAMES,
            threads: 2
            singularity: DATASETS_CONTAINER
            shell:
                r"""
                set -euo pipefail
                CONTAM_DIR=$(dirname {output.fa})
                WORK="$CONTAM_DIR/_dl"
                mkdir -p "$WORK"

                # Initialize outputs to empty so any failure mid-way leaves
                # zero-byte files (which CONTAMINATION_ENABLED treats as
                # absent on next run).
                : > {output.fa}
                : > {output.gtf}
                printf "species_id\tseqname\n" > {output.seqnames}

                # Iterate species.txt. Format mirrors download_species.sh:
                #   species_name<TAB>assembly_accession   (no header)
                while IFS=$'\t' read -r species accession || \
                    [[ -n "${{species:-}}${{accession:-}}" ]]; do
                    [[ -z "${{species// }}" ]] && continue
                    [[ "${{species:0:1}}" == "#" ]] && continue

                    # Strip CRLF and leading/trailing whitespace
                    species=$(printf '%s' "$species" | sed 's/\r$//; s/^[[:space:]]*//; s/[[:space:]]*$//')
                    accession=$(printf '%s' "$accession" | sed 's/\r$//; s/^[[:space:]]*//; s/[[:space:]]*$//')

                    [[ -z "$species" ]] && continue
                    [[ "${{species:0:1}}" == "#" ]] && continue
                    [[ -z "$accession" ]] && continue

                    if [[ ! "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
                        echo "[r05_prep] ERROR: malformed accession for $species: [$accession]" >&2
                        printf '[r05_prep] shell-escaped accession: %q\n' "$accession" >&2
                        exit 1
                    fi

                    safe_species=$(printf '%s' "$species" \
                        | sed 's/[[:space:]]\+/_/g; s#[/\\]#_#g')
                    species_id="${{safe_species}}_${{accession}}"

                    zipfile="$WORK/${{species_id}}.zip"
                    destdir="$WORK/${{species_id}}"

                    printf '[r05_prep] accession shell-escaped: %q\n' "$accession" >&2
                    datasets version >&2 || true

                    echo "[r05_prep] Downloading: $species ($accession)" >&2
                    rm -f "$zipfile"

                    if ! datasets download genome accession "$accession" \
                            --include genome,gtf \
                            --filename "$zipfile" \
                            --no-progressbar; then
                        echo "[r05_prep] WARNING: genome,gtf download failed for $species ($accession); trying genome,gff3" >&2
                        rm -f "$zipfile"
                        datasets download genome accession "$accession" \
                            --include genome,gff3 \
                            --filename "$zipfile" \
                            --no-progressbar
                    fi

                    mkdir -p "$destdir"
                    unzip -oq "$zipfile" -d "$destdir"
                    rm -f "$zipfile"

                    # Locate the genome FASTA and GTF inside the unzipped
                    # tree (datasets puts them under
                    # ncbi_dataset/data/<accession>/).
                    asm_dir="$destdir/ncbi_dataset/data/${{accession}}"
                    if [[ ! -d "$asm_dir" ]]; then
                        echo "[r05_prep] ERROR: assembly dir missing for $species" >&2
                        exit 1
                    fi

                    # FASTA: usually *_genomic.fna (no fixed prefix).
                    # Use the first .fna found that is not the genomic_gaps
                    # variant.
                    asm_fa=$(find "$asm_dir" -maxdepth 2 -name "*.fna" \
                            ! -name "*_genomic_gaps.fna" | head -1)
                    asm_gtf=$(find "$asm_dir" -maxdepth 2 -name "*.gtf" | head -1)

                    if [[ -z "$asm_fa" ]]; then
                        echo "[r05_prep] ERROR: missing genome FASTA for $species" >&2
                        echo "[r05_prep]   asm_fa=$asm_fa" >&2
                        exit 1
                    fi

                    if [[ -z "$asm_gtf" ]]; then
                        echo "[r05_prep] WARNING: missing GTF for $species ($accession); generating mock whole-contig GTF" >&2

                        asm_gtf="$destdir/${{species_id}}.mock.gtf"

                        awk -v sp="$species_id" '
                            BEGIN {{
                                OFS = "\t"
                                seq = ""
                                len = 0
                            }}

                            /^>/ {{
                                if (seq != "" && len > 0) {{
                                    gene_id = sp "__" seq
                                    tx_id   = gene_id ".t1"

                                    print seq, "mock_gtf", "gene",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "transcript", 1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "exon",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; exon_number \"1\";"
                                }}

                                seq = substr($1, 2)
                                len = 0
                                next
                            }}

                            {{
                                gsub(/[[:space:]]/, "", $0)
                                len += length($0)
                            }}

                            END {{
                                if (seq != "" && len > 0) {{
                                    gene_id = sp "__" seq
                                    tx_id   = gene_id ".t1"

                                    print seq, "mock_gtf", "gene",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "transcript", 1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "exon",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; exon_number \"1\";"
                                }}
                            }}
                        ' "$asm_fa" > "$asm_gtf"
                    fi

                    # Concatenate FASTA + GTF onto the cohort files (we do
                    # NOT rename seqnames; native NCBI accessions are used
                    # throughout, with gtf_seqnames.tsv carrying the
                    # accession-to-species mapping).
                    cat "$asm_fa"  >> {output.fa}
                    cat "$asm_gtf" >> {output.gtf}

                    # Extract seqnames from the FASTA headers and emit
                    # one row per seqname in gtf_seqnames.tsv.
                    grep '^>' "$asm_fa" | awk -v sp="$species_id" '
                        BEGIN {{ OFS="\t" }}
                        {{
                            seqname = substr($1, 2)
                            print sp, seqname
                        }}
                    ' >> {output.seqnames}

                    rm -rf "$destdir"
                done < {input.species_list}

                rm -rf "$WORK"

                # Sanity: ensure each output is non-empty
                for f in {output.fa} {output.gtf} {output.seqnames}; do
                    if [[ ! -s "$f" ]]; then
                        echo "[r05_prep] ERROR: output is empty: $f" >&2
                        exit 1
                    fi
                done

                n_species=$(awk 'NR>1 {{print $1}}' {output.seqnames} \
                            | sort -u | wc -l)
                n_seqs=$(awk 'NR>1' {output.seqnames} | wc -l)
                echo "[r05_prep] Built $n_species species, $n_seqs contigs" >&2
                """

    rule r05a0_mask_conserved_regions:
        """
        Hard-mask conserved rRNA/tRNA regions in the contamination FASTA
        so that host-derived rRNA reads don't generate cross-mapping
        contamination hits.

        Mechanism: bacterial 16S/23S rRNA has stretches 100% identical
        to human cytosolic rRNA; mitochondrial rRNA reads (very abundant
        in any RNA-seq library that didn't deplete mt-RNA) cross-map to
        apicomplexan apicoplast rRNA and to bacterial rRNA broadly;
        fungal small-subunit rRNA shares conserved cores with all of the
        above. No alignment-parameter tightening can fix this -- the
        reads ARE identical sequences. The only effective fix is to
        prevent bwa-mem2 from seeing those regions in the reference at
        all, by hard-masking them with N's. bwa-mem2 won't seed into N
        runs, so reads hitting masked regions fail to align entirely.

        Feature types masked (extracted from contamination.gtf):
          - rRNA      : the primary problem; all bacterial / eukaryotic
                        species have annotated rRNA features
          - tRNA      : also evolutionarily conserved; secondary noise
                        source
          - ncRNA     : misc conserved ncRNAs (riboswitches, etc.)
          - misc_RNA  : catch-all for additional conserved RNA features

        Predicted effect: drops Plasmodium/Toxoplasma/Penicillium/
        Cladosporium/Cutibacterium "contamination" counts to near-zero
        in samples that don't actually carry those organisms, since
        their current signal is dominated by rRNA cross-mapping. Real
        contamination (which has reads across the whole genome, not
        just at rRNA loci) is preserved.
        """
        input:
            fa = CONTAM_FA,
            gtf = CONTAM_GTF,
        output:
            masked_fa = "/tmp/data/05_contamination/index/contamination.masked.fa",
            bed = "/tmp/data/05_contamination/index/conserved_regions.bed",
        threads: 2
        singularity: BEDTOOLS_CONTAINER
        shell:
            r"""
            set -euo pipefail
            OUTDIR=$(dirname {output.masked_fa})
            mkdir -p "$OUTDIR"

            # Extract conserved-RNA features from the GTF/GFF. Handles
            # both formats because GTF column 3 (feature type) is the
            # same in GFF3; we just look for feature-type values
            # matching rRNA / tRNA / ncRNA / misc_RNA. Note: column 4 is
            # 1-based inclusive (GTF/GFF), but BED is 0-based half-open,
            # so we subtract 1 from start. We also drop features with
            # malformed coordinates (start > end), which can happen for
            # malformed annotations.
            awk -F'\t' 'BEGIN {{ OFS="\t" }}
                /^#/ {{ next }}
                NF >= 8 && ($3 == "rRNA" || $3 == "tRNA" ||
                            $3 == "ncRNA" || $3 == "misc_RNA") {{
                    if ($4 + 0 > 0 && $5 + 0 >= $4 + 0) {{
                        print $1, $4 - 1, $5
                    }}
                }}
            ' {input.gtf} | sort -k1,1 -k2,2n -u > {output.bed}

            N_FEATURES=$(wc -l < {output.bed})
            echo "[r05a0_mask] extracted $N_FEATURES conserved-RNA features to mask" >&2

            # Even if no features were found (e.g., user supplied a
            # FASTA-only reference or all-virus species.txt), continue
            # gracefully: just copy the FASTA through unchanged.
            if [[ "$N_FEATURES" -eq 0 ]]; then
                echo "[r05a0_mask] WARNING: no rRNA/tRNA features in GTF; passing FASTA through unmasked" >&2
                cp {input.fa} {output.masked_fa}
                exit 0
            fi

            # Hard-mask with N's. bwa-mem2 cannot seed into N runs, so
            # reads that would have mapped into these regions simply
            # fail to align in r05d.
            bedtools maskfasta \
                -fi {input.fa} \
                -bed {output.bed} \
                -fo {output.masked_fa} \
                -mc N

            # Sanity check: report how many bases were masked
            N_BASES_ORIG=$(grep -v '^>' {input.fa}        | tr -d '\n' | tr -cd 'Nn' | wc -c)
            N_BASES_MASK=$(grep -v '^>' {output.masked_fa} | tr -d '\n' | tr -cd 'Nn' | wc -c)
            N_BASES_ADDED=$((N_BASES_MASK - N_BASES_ORIG))
            echo "[r05a0_mask] masked $N_BASES_ADDED bp of conserved-RNA across $N_FEATURES features" >&2
            """

    rule r05a_contamination_index:
        """
        Build a bwa-mem2 index for the contamination reference. Run once,
        shared across all samples. Consumes the masked FASTA from
        r05a0_mask_conserved_regions so that conserved rRNA/tRNA
        regions are N-masked before indexing. The seqname->species
        mapping lives in gtf_seqnames.tsv (either supplied by the user
        or auto-generated by r05_prep_contamination_refs); we do NOT
        re-derive it from FASTA headers because the contamination FASTA
        uses bare NCBI accessions (e.g. NC_010163.1) which don't carry
        species info on their own.

        The index is built under a name that points at the masked
        FASTA, so r05b/r05d (which derive the FASTA path from the
        index_done file's directory) automatically pick up the masked
        version too.
        """
        input:
            fa = "/tmp/data/05_contamination/index/contamination.masked.fa",
        output:
            done = "/tmp/data/05_contamination/index/index.done",
        threads: 4
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            INDEX_DIR=$(dirname {output.done})
            mkdir -p "$INDEX_DIR"
            # r05d's $INDEX_FA expects the file named contamination.fa
            # in the index dir. Since the masked FASTA already lives
            # there, just symlink the conventional name.
            ln -sf $(basename {input.fa}) "$INDEX_DIR/contamination.fa"
            cd "$INDEX_DIR"
            bwa-mem2 index contamination.fa
            touch {output.done}
            """

    rule r05b_extract_unmapped:
        """
        Extract a sample's unmapped reads from the STAR BAM as a single
        gzipped interleaved FASTQ. Separated from the contamination
        screen so the entropy-filter step (r05c) can operate on the
        FASTQ before alignment.

        Filter flags:
          -f 4    : keep only unmapped reads
          -F 256  : drop secondary alignments
          -N      : retain /1, /2 read-pair suffix in QNAME so the
                    downstream interleaved-stream consumers can pair
                    reads correctly
        """
        input:
            bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        output:
            fq = temp("/tmp/data/05_contamination/{sample}/{sample}.unmapped.fq.gz"),
        threads: 2
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.fq})
            samtools fastq -f 4 -F 256 -N {input.bam} 2>/dev/null \
                | gzip > {output.fq}
            """

    rule r05c_entropy_filter:
        """
        Drop low-complexity / low-entropy reads using bbduk.sh before the
        contamination alignment. This is the strongest single intervention
        against the cross-mapping false-positive floor: reads like
        poly-A tails, ATx50 simple repeats, and other low-complexity
        sequences align cleanly to many references (especially AT-rich
        ones like Plasmodium falciparum or AT-rich bacterial regions)
        regardless of bwa-mem2 score / seed-length stringency settings.

        Parameters:
          entropy=0.7        Shannon entropy threshold (1.0 = max entropy).
                             Typical RNA-seq reads sit at 0.85-0.95; poly-A
                             at ~0; simple repeats at ~0.5. 0.7 drops the
                             low-complexity tail without eating into real
                             reads.
          entropywindow=20   Window size for entropy calculation (default 50;
                             20 catches shorter low-complexity stretches).
          entropyk=5         k-mer size for entropy calculation (default 5).
          int=t              Treat stream as interleaved paired-end (matches
                             samtools fastq -N output for paired data; safe
                             for single-end since singletons are emitted to
                             the same stream and processed independently).
        """
        input:
            fq = "/tmp/data/05_contamination/{sample}/{sample}.unmapped.fq.gz",
        output:
            fq = temp("/tmp/data/05_contamination/{sample}/{sample}.unmapped.filtered.fq.gz"),
        threads: 2
        singularity: BBMAP_CONTAINER
        shell:
            r"""
            set -euo pipefail
            bbduk.sh \
                in={input.fq} \
                out={output.fq} \
                entropy=0.7 \
                entropywindow=20 \
                entropyk=5 \
                int=t \
                threads={threads} \
                overwrite=t 2>&1 | grep -E "Input|Result|entropy|Time" || true
            """

    rule r05d_contamination_screen:
        """
        Align the entropy-filtered unmapped reads to the contamination
        reference with conservative bwa-mem2 parameters, then attribute
        per-species alignment counts via the seqname->species map.

        bwa-mem2 parameters:
          -T 95   : minimum alignment score (default 30). A 95bp read
                    needs near-perfect identity to be reported. Pairs
                    with the entropy filter to push the false-positive
                    rate as low as practical.
          -k 35   : minimum seed length (default 19). Longer seeds
                    reject short noisy matches before extension.
          -r 2.0  : re-seed threshold multiplier (default 1.5). Less
                    aggressive re-seeding, fewer marginal alignments.
          -h 50   : report up to 50 tied-best alternatives in the XA
                    tag (default 5). Ensures all same-score alignments
                    to closely related references are captured.

        Multi-mapping handling: a read aligning equally well to N
        references is COUNTED ONCE PER REFERENCE it hits (primary +
        each XA-tag alternative). This means per-species counts
        represent "alignments supporting this species" rather than
        strict per-read counts; the sum across species can EXCEED
        total_mapped_to_contamination. Rationale: the strict alignment
        parameters (-T 95 -k 35) together with the upstream rRNA /
        conserved-region masking (r05a0_mask_conserved_regions) make
        per-genome attribution of multi-mappers a sensitivity recovery
        rather than a specificity risk. Without this change, closely
        related references (e.g., multiple Streptomyces species)
        cancel each other out at MAPQ filtering and ALL get
        undercounted.

        Per-species attribution uses the user-supplied gtf_seqnames.tsv
        (header: species_id<TAB>seqname). Reads whose accession is
        absent from the map are bucketed into 'unknown' rather than
        dropped, so additions to the FASTA without matching seqname-map
        entries are still counted and visibly attributed.
        """
        input:
            fq = "/tmp/data/05_contamination/{sample}/{sample}.unmapped.filtered.fq.gz",
            bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
            index_done = "/tmp/data/05_contamination/index/index.done",
            species_map = CONTAM_SEQNAMES,
        output:
            counts = "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
            breadth = "/tmp/data/05_contamination/{sample}/{sample}.contamination_breadth.tsv",
            bam = temp("/tmp/data/05_contamination/{sample}/{sample}.contamination.bam"),
        threads: 4
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            OUTDIR=$(dirname {output.counts})
            mkdir -p "$OUTDIR"
            echo "$OUTDIR"

            INDEX_FA=$(dirname {input.index_done})/contamination.fa

            # Count post-entropy-filter reads (the read pool that actually
            # entered bwa-mem2). This is the right denominator for "fraction
            # of attempted reads attributed to each species" in r05g. We
            # count the input FASTQ directly because samtools view -c on the
            # output BAM counts ALIGNMENT records, not reads (a read with
            # supplementary alignments would be counted more than once).
            # 4 lines per FASTQ record -> divide by 4. zcat is fine since
            # the entropy filter always emits gzip-compressed output.

            TOTAL_READS_ALIGNED=$(zcat {input.fq} | awk 'END {{print int(NR/4)}}')

            # bwa-mem2 with -p reads interleaved paired-end from stdin
            # (matches bbduk's int=t output and samtools fastq -N).
            # Works correctly for SE input too since singletons stream
            # through the same path. -h 50 raises the XA-tag alternative
            # cap from 5 to 50 so all tied-best multi-mappers are kept.

            zcat {input.fq} \
                | bwa-mem2 mem -t {threads} -p -T 95 -k 35 -r 2.0 -B 6 -O 8 -L 10 -h 50 \
                      "$INDEX_FA" - 2>/dev/null \
                | samtools view -b -F 4 - > {output.bam}

            # Totals: total unmapped reads in the host BAM (denominator
            # context only) and total READS confidently mapped to any
            # contaminant genome (no MAPQ filter; the strict alignment
            # parameters already filter out marginal alignments). Note
            # that the sum of per-species counts below can exceed this
            # because multi-mappers are counted once per sequence.

            TOTAL_UNMAPPED=$(samtools view -c -f 4 -F 256 {input.bam})
            TOTAL_MAPPED=$(samtools view -c -F 260 {output.bam})

            # Per-contig coverage breadth via samtools coverage. The
            # 'coverage' column is the percentage of contig bases with
            # >=1 read, i.e. the breadth of coverage. r05f aggregates
            # this across contigs of the same species (weighted by
            # contig length) into a per-species uniformity score.
            # Note: samtools coverage needs a sorted+indexed BAM. We
            # sort to a side-temp BAM rather than the rule's main BAM
            # (the main BAM is consumed unsorted by the awk pipeline
            # above; sorting it in place would invalidate the layout
            # other tools may expect from a contamination BAM).

            TMP_SORT_DIR=$(mktemp -d "$OUTDIR/sorted.XXXXXX")
            SORTED_BAM="$TMP_SORT_DIR/sorted.bam"
            trap 'rm -rf -- "${{TMP_SORT_DIR:-}}"' EXIT INT TERM
            samtools sort -@ {threads} -o "$SORTED_BAM" {output.bam:q}
            samtools index "$SORTED_BAM"

            # Header line begins with '#rname'. samtools coverage emits
            # one row per reference contig.

            samtools coverage "$SORTED_BAM" > {output.breadth}

            # Per-species counts via samtools view + awk join on the
            # gtf_seqnames.tsv map. The awk parses both the primary
            # alignment reference (column 3) AND the XA tag (when
            # present) so reads aligning equally well to multiple
            # references are counted once per reference.

            samtools view -F 260 {output.bam} \
                | awk -v map="{input.species_map}" \
                      -v sample="{wildcards.sample}" \
                      -v total_unmapped="$TOTAL_UNMAPPED" \
                      -v total_mapped="$TOTAL_MAPPED" \
                      -v total_reads_aligned="$TOTAL_READS_ALIGNED" \
                      -v out_counts="{output.counts}" '
                    BEGIN {{
                        OFS="\t"
                        # load seqname -> species_id; build species set
                        while ((getline line < map) > 0) {{
                            n = split(line, a, "\t")
                            if (n == 2 && a[1] != "species_id") {{
                                # a[1] = species_id, a[2] = seqname
                                seqname2species[a[2]] = a[1]
                                all_species[a[1]] = 1
                            }}
                        }}
                        close(map)
                    }}
                    function attribute(seqname,    species) {{
                        species = (seqname in seqname2species) ? seqname2species[seqname] : "unknown"
                        counts[species]++
                        all_species[species] = 1
                    }}
                    {{
                        # Always count the primary alignment.
                        attribute($3)
                        # Scan optional fields (columns 12..NF) for an XA
                        # tag and credit each comma-delimited alternative.
                        # XA format: XA:Z:chr,pos,CIGAR,NM;chr,pos,CIGAR,NM;
                        # The accession is the substring before the first
                        # comma in each ";"-delimited entry.

                        for (i = 12; i <= NF; i++) {{
                            if (substr($i, 1, 5) == "XA:Z:") {{
                                xa_payload = substr($i, 6)
                                ne = split(xa_payload, xa_entries, ";")
                                for (e = 1; e <= ne; e++) {{
                                    if (xa_entries[e] == "") continue
                                    split(xa_entries[e], xa_fields, ",")
                                    if (xa_fields[1] != "" && xa_fields[1] != $3) {{
                                        attribute(xa_fields[1])
                                    }}
                                }}
                                break
                            }}
                        }}
                    }}
                    END {{
                        print "metric", "value" > out_counts
                        print "sample", sample > out_counts
                        print "total_unmapped_in_star_bam", total_unmapped > out_counts
                        print "total_reads_aligned", total_reads_aligned > out_counts
                        print "total_mapped_to_contamination", total_mapped > out_counts
                        # emit a deterministic order of species rows
                        n = 0
                        for (s in all_species) sp_arr[++n] = s
                        # naive insertion sort
                        for (i = 2; i <= n; i++) {{
                            key = sp_arr[i]; j = i - 1
                            while (j > 0 && sp_arr[j] > key) {{
                                sp_arr[j + 1] = sp_arr[j]; j--
                            }}
                            sp_arr[j + 1] = key
                        }}
                        for (i = 1; i <= n; i++) {{
                            print "species:" sp_arr[i], (sp_arr[i] in counts ? counts[sp_arr[i]] : 0) > out_counts
                        }}
                    }}
                '
            """

    rule r05e_contamination_summary:
            """
            Aggregate per-sample contamination counts into a single wide TSV
            (one row per sample, one column per species). Independent of
            r04d_qc_summary: no Input_reads percentages, no flags, no Picard
            merging -- just the raw per-species alignment counts from r05d
            for inspection or downstream analysis.

            Column order matches the species.txt source file (with comments
            and blank lines stripped), preserving the user's intended
            ordering. When species.txt is not present (prebuilt-only mode),
            falls back to first-occurrence order in the seqname map.

            Columns:
            - sample
            - total_unmapped_in_star_bam        (reads unmapped in host BAM)
            - total_mapped_to_contamination     (reads aligned to contam ref;
                                                note multi-mappers are
                                                counted once per matched
                                                sequence in the per-species
                                                columns below, so
                                                sum(species_*) >= this)
            - one column per species_id, in species.txt order, holding the
                "alignments supporting this species" count from r05d (see
                r05d docstring for multi-mapping semantics)
            - "unknown" column appended at the end IF any sample has reads
                attributed to the unknown bucket (contig not in seqname map)
            """
            input:
                counts = expand(
                    "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
                    sample=samples,
                ),
                species_map = CONTAM_SEQNAMES,
            output:
                tsv = "/tmp/data/05_contamination/00_contamination_summary.tsv",
            run:
                # Canonical species ordering from species.txt (or seqname map
                # fallback). Helper handles blank lines, comments, whitespace
                # sanitisation -- see _read_species_order() definition near
                # the top of this file.
                species_order = _read_species_order(
                    CONTAM_SPECIES_LIST,
                    input.species_map,
                )
                # Catch species in the seqname map that aren't in species.txt
                # (post-build edits to the map). Append at end, alphabetically,
                # so they're visible rather than silently dropped.
                in_order = set(species_order)
                extras = set()
                with open(input.species_map) as f:
                    next(f)
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) == 2 and parts[0] and parts[0] not in in_order:
                            extras.add(parts[0])
                if extras:
                    print(f"[contam summary] WARN: {len(extras)} species_ids in "
                        f"seqname map but not in species.txt; appending at end: "
                        f"{sorted(extras)[:5]}...", file=sys.stderr)
                    species_order = species_order + sorted(extras)

                base_cols = [
                    "sample",
                    "total_unmapped_in_star_bam",
                    "total_mapped_to_contamination",
                ]

                # Two-pass: scan all files to (a) parse counts and (b) detect
                # whether any sample has an "unknown" bucket, so we only emit
                # that column when relevant.
                parsed = {}
                saw_unknown = False
                for path in input.counts:
                    data = {}
                    with open(path) as fh:
                        next(fh)  # header (metric<TAB>value)
                        for line in fh:
                            parts = line.rstrip("\n").split("\t", 1)
                            if len(parts) == 2:
                                data[parts[0]] = parts[1]
                    sample = data.get("sample", "")
                    if not sample:
                        print(f"[contam summary] WARN: no sample id in {path}, skipping",
                            file=sys.stderr)
                        continue
                    parsed[sample] = data
                    try:
                        if int(data.get("species:unknown", "0") or 0) > 0:
                            saw_unknown = True
                    except ValueError:
                        pass

                with open(output.tsv, "w") as fout:
                    header = base_cols + species_order + (["unknown"] if saw_unknown else [])
                    fout.write("\t".join(header) + "\n")

                    n_rows_written = 0
                    # Emit rows in cohort sample-list order (consistent with
                    # other summary tables) rather than file-discovery order.
                    for sample in samples:
                        if sample not in parsed:
                            print(f"[contam summary] WARN: no contamination counts "
                                f"for sample {sample}; row will be blank",
                                file=sys.stderr)
                            row = [sample] + ["" for _ in header[1:]]
                            fout.write("\t".join(row) + "\n")
                            continue
                        data = parsed[sample]
                        row = [
                            sample,
                            data.get("total_unmapped_in_star_bam", "0"),
                            data.get("total_mapped_to_contamination", "0"),
                        ]
                        for sp in species_order:
                            row.append(data.get(f"species:{sp}", "0"))
                        if saw_unknown:
                            row.append(data.get("species:unknown", "0"))
                        fout.write("\t".join(row) + "\n")
                        n_rows_written += 1

                print(f"[contam summary] wrote {output.tsv} -- "
                    f"{n_rows_written} samples x {len(species_order)} species"
                    f"{' (+ unknown bucket)' if saw_unknown else ''}",
                    flush=True)

    rule r05f_contamination_uniformity:
        """
        Aggregate per-sample coverage data into a per-species uniformity
        ratio. One row per sample, one column per species in species.txt
        order.

        SCORE: Lander-Waterman uniformity ratio
        ---------------------------------------
        For each (sample, species), the score is:

            rho_s = B_observed / E[B]_uniform

        where
          - B_observed  = sum(covbases) across the species's contigs
                          (bases of the species's genome covered by >=1 read)
          - E[B]_uniform = L * (1 - exp(-N*r/L))
            is the expected covered bases under random uniform placement
            of N reads of length r on L bases of genome
          - N            = total read count across the species's contigs
          - L            = sum of contig lengths for the species
          - r            = typical read length (CONTAM_READ_LENGTH, default 75)

        Interpretation
        --------------
          - rho ~ 1.0  : reads are placed about as uniformly as random
                         placement would predict (consistent with reads
                         genuinely originating from many loci across the
                         organism's genome)
          - rho ~ 0    : reads are heavily clustered relative to uniform
                         expectation (consistent with off-target alignment
                         to a few conserved loci, OR with biology that
                         legitimately concentrates reads, e.g., bradyzoite
                         stage of latent Toxoplasma where only a subset
                         of genes is expressed)
          - rho > 1    : observed breadth exceeds uniform expectation;
                         can happen due to read-length under-estimation
                         or under-dispersion of reads. Values slightly
                         above 1 are normal; values much above 1 suggest
                         the CONTAM_READ_LENGTH constant is too low for
                         this dataset
          - empty cell : no reads aligned to this species in this sample

        Key advantage over raw breadth (covbases/L)
        -------------------------------------------
        Raw breadth is bounded above by N*r/L, which is tiny for large
        genomes at low read counts. So real and spurious alignments give
        similar tiny breadth values when N*r << L (e.g., 10000 reads on
        a 65 Mb Toxoplasma genome gives raw breadth ~0.01 even at perfect
        uniformity). The ratio rho factors out N and L, making the score
        comparable across species of vastly different genome sizes and
        across samples with vastly different read counts.

        Caveats
        -------
          - r is approximated by CONTAM_READ_LENGTH = 75 (configurable at
            the top of this file). Robust to ~20% error.
          - At very low N (say N < 10), the expectation is small and
            noisy; rho can swing wildly between samples. The "raw
            breadth" was also unreliable in this regime. Read this
            score alongside the read-count from r05e.
          - rho measures spatial uniformity, not biological reality. A
            real organism in a developmental stage with restricted
            expression (e.g., bradyzoite tissue cysts) will give LOW rho
            even when present. Combine with the fraction-vs-controls
            analysis (r05g + cell-line controls) for a complete picture.

        Reference: Lander & Waterman, Genomics 1988 (original derivation
        of expected coverage for random shotgun sequencing).
        """
        input:
            breadth = expand(
                "/tmp/data/05_contamination/{sample}/{sample}.contamination_breadth.tsv",
                sample=samples,
            ),
            species_map = CONTAM_SEQNAMES,
        output:
            tsv = "/tmp/data/05_contamination/00_contamination_uniformity.tsv",
        run:
            import math

            # Canonical species ordering from species.txt (or seqname map
            # fallback). Same helper r05e uses.
            species_order = _read_species_order(
                CONTAM_SPECIES_LIST,
                input.species_map,
            )
            # Append any species_ids in the seqname map that aren't in
            # species.txt (post-build edits) at end, alphabetically.
            in_order = set(species_order)
            extras = set()
            seqname_to_species = {}
            with open(input.species_map) as f:
                next(f)  # header
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2 and parts[0]:
                        seqname_to_species[parts[1]] = parts[0]
                        if parts[0] not in in_order:
                            extras.add(parts[0])
            if extras:
                print(f"[contam uniformity] WARN: {len(extras)} species_ids "
                      f"in seqname map but not in species.txt; appending at end: "
                      f"{sorted(extras)[:5]}...", file=sys.stderr)
                species_order = species_order + sorted(extras)

            # Parse each per-sample breadth file (samtools coverage output).
            # Columns from samtools coverage:
            #   1: #rname (contig)
            #   2: startpos
            #   3: endpos
            #   4: numreads
            #   5: covbases  (bases covered by >=1 read)
            #   6: coverage  (covbases / contig_length * 100, percentage)
            #   7-9: meandepth, meanbaseq, meanmapq
            #
            # We aggregate per (sample, species):
            #   N_s = sum(numreads), L_s = sum(contig_length), B_s = sum(covbases)
            r = CONTAM_READ_LENGTH

            def lw_ratio(B_obs, N, L, r=r):
                """
                Lander-Waterman uniformity ratio.
                E[B]_uniform = L * (1 - exp(-N*r/L))
                For very small expected B (i.e., N*r/L approx 0), returns
                None to signal "undefined" rather than a noisy ratio.
                """
                if N <= 0 or L <= 0:
                    return None
                lam = (N * r) / L  # mean coverage in reads-per-base
                E_B = L * (1.0 - math.exp(-lam))
                if E_B <= 0:
                    return None
                return B_obs / E_B

            parsed = {}  # sample -> {species_id -> (B_obs, N, L)}
            for path in input.breadth:
                sample = os.path.basename(path).replace(".contamination_breadth.tsv", "")
                per_species = {}  # species_id -> [B_obs, N, L]
                with open(path) as fh:
                    header = fh.readline()
                    if not header.startswith("#"):
                        parsed[sample] = {}
                        continue
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 6:
                            continue
                        rname = parts[0]
                        try:
                            startpos_i = int(parts[1])
                            endpos_i = int(parts[2])
                            numreads = int(parts[3])
                            covbases = int(parts[4])
                        except ValueError:
                            continue
                        contig_len = endpos_i - startpos_i + 1
                        species = seqname_to_species.get(rname, "unknown")
                        if species not in per_species:
                            per_species[species] = [0, 0, 0]  # B, N, L
                        per_species[species][0] += covbases
                        per_species[species][1] += numreads
                        per_species[species][2] += contig_len
                parsed[sample] = per_species

            # Decide whether to emit an "unknown" column
            saw_unknown = any(
                "unknown" in d and d["unknown"][1] > 0
                for d in parsed.values()
            )

            with open(output.tsv, "w") as fout:
                header_cols = ["sample"] + species_order + (["unknown"] if saw_unknown else [])
                fout.write("\t".join(header_cols) + "\n")
                for sample in samples:
                    if sample not in parsed:
                        print(f"[contam uniformity] WARN: no breadth file "
                              f"for sample {sample}; row will be blank",
                              file=sys.stderr)
                        row = [sample] + ["" for _ in header_cols[1:]]
                        fout.write("\t".join(row) + "\n")
                        continue
                    d = parsed[sample]
                    row = [sample]
                    for sp in species_order:
                        if sp in d:
                            B, N, L = d[sp]
                            rho = lw_ratio(B, N, L)
                            row.append(f"{rho:.4f}" if rho is not None else "")
                        else:
                            row.append("")
                    if saw_unknown:
                        if "unknown" in d:
                            B, N, L = d["unknown"]
                            rho = lw_ratio(B, N, L)
                            row.append(f"{rho:.4f}" if rho is not None else "")
                        else:
                            row.append("")
                    fout.write("\t".join(row) + "\n")

            print(f"[contam uniformity] wrote {output.tsv} -- "
                  f"{len(samples)} samples x {len(species_order)} species "
                  f"(read length assumed {r} bp)"
                  f"{' (+ unknown bucket)' if saw_unknown else ''}",
                  flush=True)


    rule r05g_contamination_fractions:
        """
        Aggregate per-sample contamination counts as FRACTIONS of the
        post-entropy-filter read pool (i.e., the reads that actually
        entered bwa-mem2). Wide TSV with the same species column order
        as r05e/r05f.

        Denominator: total_reads_aligned (the post-entropy-filter
        FASTQ read count, tracked per sample by r05d). This is the
        correct denominator for "of the reads we tried to attribute,
        what fraction landed in each species" -- it accounts for
        entropy-filter dropouts (low-complexity reads removed before
        alignment) so the fractions sum to <= the total fraction of
        reads attributed to ANY contaminant.

        Because multi-mapping reads are counted once per matched
        species (r05d's XA-tag-aware counting), the sum of per-species
        fractions can exceed total_mapped_to_contamination / total_reads_aligned.
        That's by design -- a read with same-best alignments to N species
        contributes 1 to each species's numerator.

        Values are emitted as plain decimal fractions (e.g., 0.0123 for
        1.23 percent). Multiply by 100 downstream for percentage display.
        """
        input:
            counts = expand(
                "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
                sample=samples,
            ),
            species_map = CONTAM_SEQNAMES,
        output:
            tsv = "/tmp/data/05_contamination/00_contamination_fractions.tsv",
        run:
            species_order = _read_species_order(
                CONTAM_SPECIES_LIST,
                input.species_map,
            )
            in_order = set(species_order)
            extras = set()
            with open(input.species_map) as f:
                next(f)
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2 and parts[0] and parts[0] not in in_order:
                        extras.add(parts[0])
            if extras:
                print(f"[contam fractions] WARN: {len(extras)} species_ids "
                        f"in seqname map but not in species.txt; appending at end: "
                        f"{sorted(extras)[:5]}...", file=sys.stderr)
                species_order = species_order + sorted(extras)

            parsed = {}
            saw_unknown = False
            for path in input.counts:
                data = {}
                with open(path) as fh:
                    next(fh)  # header
                    for line in fh:
                        parts = line.rstrip("\n").split("\t", 1)
                        if len(parts) == 2:
                            data[parts[0]] = parts[1]
                sample = data.get("sample", "")
                if not sample:
                    print(f"[contam fractions] WARN: no sample id in {path}, skipping",
                            file=sys.stderr)
                    continue
                parsed[sample] = data
                try:
                    if int(data.get("species:unknown", "0") or 0) > 0:
                        saw_unknown = True
                except ValueError:
                    pass

            with open(output.tsv, "w") as fout:
                header_cols = (
                    ["sample", "total_reads_aligned"]
                    + species_order
                    + (["unknown"] if saw_unknown else [])
                )
                fout.write("\t".join(header_cols) + "\n")
                n_rows_written = 0
                for sample in samples:
                    if sample not in parsed:
                        print(f"[contam fractions] WARN: no contamination counts "
                                f"for sample {sample}; row will be blank",
                                file=sys.stderr)
                        row = [sample] + ["" for _ in header_cols[1:]]
                        fout.write("\t".join(row) + "\n")
                        continue
                    data = parsed[sample]
                    try:
                        denom = int(data.get("total_reads_aligned", "0") or 0)
                    except ValueError:
                        denom = 0
                    if denom <= 0:
                        # Avoid division by zero; emit blank fractions.
                        print(f"[contam fractions] WARN: sample {sample} has "
                                f"total_reads_aligned={denom}; fractions emitted as blank",
                                file=sys.stderr)
                        row = [sample, str(denom)] + ["" for _ in species_order]
                        if saw_unknown:
                            row.append("")
                        fout.write("\t".join(row) + "\n")
                        continue
                    row = [sample, str(denom)]
                    for sp in species_order:
                        try:
                            cnt = int(data.get(f"species:{sp}", "0") or 0)
                        except ValueError:
                            cnt = 0
                        row.append(f"{cnt / denom:.6f}")
                    if saw_unknown:
                        try:
                            cnt = int(data.get("species:unknown", "0") or 0)
                        except ValueError:
                            cnt = 0
                        row.append(f"{cnt / denom:.6f}")
                    fout.write("\t".join(row) + "\n")
                    n_rows_written += 1

            print(f"[contam fractions] wrote {output.tsv} -- "
                    f"{n_rows_written} samples x {len(species_order)} species"
                    f"{' (+ unknown bucket)' if saw_unknown else ''}",
                    flush=True)

rule r04d_qc_summary:
    """
    Aggregate per-sample CollectRnaSeqMetrics outputs (Picard) and STAR
    mapping statistics into one table, flagging outliers. Picard's metrics
    file has two main sections:
      - METRICS section: one line of values per sample
      - HISTOGRAM section: normalized coverage by percentile of transcript
        position (0=5'-end, 100=3'-end); we don't aggregate this here but
        leave the per-sample file available for plotting.

    The STAR Log.final.out values come from r03h_mapping_stat's aggregated
    output (mapping_stat.txt) so we don't re-parse the per-sample logs here.
    The uniquely-mapped read count drives the library-complexity flag.

    If the optional contamination screen (r05a/r05b) is enabled, the
    per-sample contamination counts are folded into the same cohort table.
    For each species in the contamination reference we add two data
    columns (<species>_reads and <species>_rpm) plus a per-species flag.
    A combined flag_contamination fires if total contaminant reads exceed
    100 reads per million unmapped (the Olarerin-George & Hogenesch 2015
    threshold).

    Thresholds chosen for DROP / OUTRIDER analysis:
      - 3' bias high: cohort-relative outlier (median + 3*MAD on MEDIAN_3PRIME_BIAS)
      - rRNA high:    > 10%  (PCT_RIBOSOMAL_BASES > 0.10)
      - mRNA low:     < 60%  (PCT_MRNA_BASES < 0.60)
      - low complexity: < 10 M uniquely-mapped reads (OUTRIDER/FRASER lose
                        statistical power below this; ENCODE's full-quality
                        bar is 30 M, but 10 M is the working floor)
      - contamination: combined contaminant reads > 0.01% of total input reads
      - per-species:   single-species contaminant reads > 0.001% of total input reads

    Contamination percentages use TOTAL input reads (STAR's "Number of
    input reads") as denominator, so percentages are directly comparable
    across samples regardless of host-alignment rate. The 0.01% threshold
    corresponds to 100 contaminant reads per million total reads -- the
    same Olarerin-George & Hogenesch (2015, NAR) bar as the original
    rpm-of-unmapped formulation, just rescaled to a denominator that
    doesn't change with sample alignment quality.
    """
    input:
        unpack(lambda wc: {
            "metrics_files": expand(
                "/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt",
                sample=samples),
            "mapping_stat":
                "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
            **(
                {
                    "contam_counts": expand(
                        "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
                        sample=samples),
                    "contam_species_map": CONTAM_SEQNAMES,
                }
                if CONTAMINATION_ENABLED else {}
            ),
        }),
    output:
        summary = "/tmp/data/04_qc/00_qc_summary.tsv",
    run:
        import os, statistics

        # Columns of interest from Picard's METRICS section
        cols_of_interest = [
            "PF_BASES", "PF_ALIGNED_BASES",
            "PCT_RIBOSOMAL_BASES",
            "PCT_CODING_BASES", "PCT_UTR_BASES",
            "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES",
            "PCT_MRNA_BASES", "PCT_USABLE_BASES",
            "MEDIAN_CV_COVERAGE",
            "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS",
            "MEDIAN_5PRIME_TO_3PRIME_BIAS",
            "PCT_CORRECT_STRAND_READS",
        ]

        # Extra columns we surface from STAR's mapping_stat.txt. Input_reads
        # is needed as the denominator for contamination percentages so that
        # the value is comparable across samples (independent of host-alignment
        # rate, library size, etc).
        mapping_cols = ["Input_reads", "Uniquely_mapped", "Uniquely_mapped_pct"]

        # Outlier flags from RNA-quality / complexity metrics
        flag_cols = [
            "flag_3p_bias_high",
            "flag_rrna_high",
            "flag_mrna_low",
            "flag_low_complexity",
        ]

        # -- Parse Picard metrics files --
        rows = {}
        for f in input.metrics_files:
            sample = os.path.basename(f).replace(".rnaseq_metrics.txt", "")
            with open(f) as fh:
                lines = fh.readlines()
            # Find the METRICS header line (starts with "PF_BASES")
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith("PF_BASES"):
                    header_idx = i
                    break
            if header_idx is None:
                continue
            header = lines[header_idx].strip().split("\t")
            values = lines[header_idx + 1].strip().split("\t")
            row = dict(zip(header, values))
            rows[sample] = {c: row.get(c, "") for c in cols_of_interest}

        # -- Parse STAR mapping_stat.txt for read-count columns --
        # Format (from r03h_mapping_stat):
        #   Sample  Input_reads  Uniquely_mapped  Uniquely_mapped_pct
        #           Multi_mapped Multi_mapped_pct Too_many_loci Unmapped
        with open(input.mapping_stat) as f:
            ms_lines = f.readlines()
        if ms_lines:
            ms_header = ms_lines[0].rstrip("\n").split("\t")
            for line in ms_lines[1:]:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0]:
                    continue
                sample = parts[0]
                if sample not in rows:
                    # STAR-only sample; will appear as Picard-NaN row below
                    rows[sample] = {c: "" for c in cols_of_interest}
                ms_row = dict(zip(ms_header, parts))
                for c in mapping_cols:
                    rows[sample][c] = ms_row.get(c, "")

        # -- Parse contamination counts (when CONTAMINATION_ENABLED) --
        # Format per-sample (long-form key-value, from r05b):
        #   metric                              value
        #   sample                              S001
        #   total_unmapped_in_star_bam          12345678
        #   total_mapped_to_contamination       42
        #   species:mycoplasma_hyorhinis_GCF... 37
        #   species:candida_albicans_GCF...     5
        #   ...
        #
        # Percentages use the sample's TOTAL input reads (STAR's
        # "Number of input reads", surfaced as Input_reads) as denominator,
        # so they are directly comparable across samples regardless of how
        # much of the library aligned to the host genome.
        contam_cols = []           # Data columns appended to row table
        contam_flag_cols = []      # Flag columns appended after the RNA flags
        species_sorted = []        # Stable ordering for column emission
        if CONTAMINATION_ENABLED:
            # Column order matches the species.txt source file (with comments
            # and blank lines stripped). When species.txt is not present
            # (prebuilt-only mode), falls back to first-occurrence order in
            # the seqname map. Both branches deduplicate, so a species listed
            # twice still occupies a single column at its first position.
            # See _read_species_order() near the top of this file for details.
            species_sorted = _read_species_order(
                CONTAM_SPECIES_LIST,
                input.contam_species_map,
            )
            # Defensive: any species_ids present in the seqname map but not
            # in species.txt are appended at the end in alphabetical order so
            # they don't silently disappear from the summary. This catches
            # post-build edits to the seqname map (rare but possible).
            in_order = set(species_sorted)
            extras = set()
            with open(input.contam_species_map) as f:
                next(f)
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2 and parts[0] and parts[0] not in in_order:
                        extras.add(parts[0])
            if extras:
                print(f"[qc summary] WARN: {len(extras)} species_ids in "
                      f"seqname map are not in species.txt; appending at end: "
                      f"{sorted(extras)[:5]}...", file=sys.stderr)
                species_sorted = species_sorted + sorted(extras)

            contam_cols = (
                ["total_mapped_to_contamination",
                 "contamination_pct_of_input"]
                + [f"{sp}_reads" for sp in species_sorted]
                + [f"{sp}_pct"   for sp in species_sorted]
            )
            contam_flag_cols = (
                ["flag_contamination"]
                + [f"flag_{sp}" for sp in species_sorted]
            )

            for f_path in input.contam_counts:
                data = {}
                with open(f_path) as fh:
                    next(fh)  # skip header
                    for line in fh:
                        parts = line.rstrip("\n").split("\t", 1)
                        if len(parts) == 2:
                            data[parts[0]] = parts[1]

                sample = data.get("sample", "")
                if not sample:
                    continue
                if sample not in rows:
                    rows[sample] = {c: "" for c in cols_of_interest}

                # Denominator: TOTAL input reads (from mapping_stat, already
                # in the row dict via the mapping_cols pass above).
                try:
                    total_input = int(rows[sample].get("Input_reads", 0) or 0)
                except (ValueError, TypeError):
                    total_input = 0
                total_contam = int(data.get("total_mapped_to_contamination", 0) or 0)
                pct = (total_contam / total_input * 100) if total_input > 0 else 0.0

                rows[sample]["total_mapped_to_contamination"] = total_contam
                rows[sample]["contamination_pct_of_input"] = f"{pct:.6f}"

                for sp in species_sorted:
                    sp_reads = int(data.get(f"species:{sp}", 0) or 0)
                    sp_pct = (sp_reads / total_input * 100) if total_input > 0 else 0.0
                    rows[sample][f"{sp}_reads"] = sp_reads
                    rows[sample][f"{sp}_pct"]   = f"{sp_pct:.6f}"

        # -- Compute cohort statistics for outlier detection --
        def floats(col):
            return [float(rows[s][col]) for s in rows
                    if rows[s].get(col, "") not in ("", "?", "NA")]
        bias_3p = floats("MEDIAN_3PRIME_BIAS")
        med_3p_bias = statistics.median(bias_3p) if bias_3p else 0
        mad_3p_bias = statistics.median([abs(x - med_3p_bias) for x in bias_3p]) \
                      if bias_3p else 0
        bias_3p_threshold = med_3p_bias + 3 * mad_3p_bias

        # Library-complexity threshold: 10M uniquely-mapped reads. ENCODE's
        # gold-standard bar is 30M, but OUTRIDER/FRASER autoencoders can still
        # work down to ~10M with reduced power; below that, statistical
        # outlier detection becomes unreliable.
        complexity_threshold = 10_000_000

        # Contamination thresholds (only used when CONTAMINATION_ENABLED).
        # Both are expressed as percent of total input reads to match the
        # contamination column units. 100 reads per million total reads =
        # 0.01% combined; per-species 10 reads per million = 0.001%.
        contam_total_pct_threshold   = 0.01
        contam_species_pct_threshold = 0.001

        def is_outlier(sample):
            r = rows[sample]
            flags = {}
            # 3' bias: cohort-relative MAD outlier
            try:
                flags["flag_3p_bias_high"] = float(r["MEDIAN_3PRIME_BIAS"]) > bias_3p_threshold
            except (ValueError, KeyError):
                flags["flag_3p_bias_high"] = ""
            # rRNA: absolute threshold
            try:
                flags["flag_rrna_high"] = float(r["PCT_RIBOSOMAL_BASES"]) > 0.10
            except (ValueError, KeyError):
                flags["flag_rrna_high"] = ""
            # mRNA: absolute threshold
            try:
                flags["flag_mrna_low"] = float(r["PCT_MRNA_BASES"]) < 0.60
            except (ValueError, KeyError):
                flags["flag_mrna_low"] = ""
            # Library complexity: absolute threshold on uniquely-mapped reads
            try:
                flags["flag_low_complexity"] = int(r["Uniquely_mapped"]) < complexity_threshold
            except (ValueError, KeyError):
                flags["flag_low_complexity"] = ""

            # Contamination flags (only when CONTAMINATION_ENABLED)
            if CONTAMINATION_ENABLED:
                try:
                    flags["flag_contamination"] = (
                        float(r["contamination_pct_of_input"]) > contam_total_pct_threshold
                    )
                except (ValueError, KeyError):
                    flags["flag_contamination"] = ""
                for sp in species_sorted:
                    try:
                        flags[f"flag_{sp}"] = (
                            float(r[f"{sp}_pct"]) > contam_species_pct_threshold
                        )
                    except (ValueError, KeyError):
                        flags[f"flag_{sp}"] = ""
            return flags

        # -- Write summary --
        all_data_cols = cols_of_interest + mapping_cols + contam_cols
        all_flag_cols = flag_cols + contam_flag_cols
        with open(output.summary, "w") as fout:
            fout.write("sample\t" + "\t".join(all_data_cols)
                       + "\t" + "\t".join(all_flag_cols) + "\n")
            for sample in sorted(rows):
                vals = [str(rows[sample].get(c, "")) for c in all_data_cols]
                flags = is_outlier(sample)
                flag_vals = [str(flags.get(c, "")) for c in all_flag_cols]
                fout.write(sample + "\t" + "\t".join(vals)
                           + "\t" + "\t".join(flag_vals) + "\n")



# =============================================================================
# Optional HLA class I + class II typing with arcasHLA (RNA-seq)
# =============================================================================
# These rules are only added to the DAG when ARCASHLA_ENABLED is True (i.e.
# the reference at /tmp/data/00_additional_files/arcashla_ref/dat is fully
# built; see the gating block near the top of this file). arcasHLA is the
# correct tool for HLA typing from RNA-seq: it pseudo-aligns reads to a cDNA
# reference, which is exactly what RNA-seq provides.
#
# Workflow per sample:
#   r06a_arcashla_extract           : pull HLA-region reads out of the markdup
#                                     BAM into a pair of FASTQs (chr6 + decoys,
#                                     handled by arcasHLA internally).
#   r06b_arcashla_genotype          : class I alleles A, B, C. Emits
#                                     {sample}.genotype.json.
#   r06c_hla_summary                : aggregate class I JSONs into one TSV
#                                     (00_hla_summary.tsv).
#   r06d_arcashla_genotype_classII  : class II alleles DRB1, DQA1, DQB1, DPA1,
#                                     DPB1 from the SAME extracted FASTQs. Kept
#                                     separate because class II expression (and
#                                     thus typing reliability) is much lower and
#                                     patchier in brain RNA-seq. Emits
#                                     {sample}.genotype.classII.json.
#   r06e_hla_summary_classII        : aggregate class II JSONs into a separate
#                                     TSV (00_hla_summary_classII.tsv).
#
# IMPORTANT run-time bind requirement: the biocontainer ships an incomplete
# dat/ directory. The host arcashla_ref/dat must be bound OVER the container's
# internal dat/. With apptainer that means the snakemake invocation needs
# something like:
#
#   snakemake --use-singularity \
#     --singularity-args "--bind /mnt/data/rnaseq:/tmp/data ... \
#       --bind /mnt/data/rnaseq/rnaseq-drop/00_additional_files/arcashla_ref/dat:/usr/local/share/arcas-hla-0.6.0-2/dat"
#
# (Adjust the outer bind to whatever maps data root to /tmp/data.) The
# installer script prints the exact line for this machine.
#
# arcasHLA CLI quirks baked into these rules (learned the hard way):
#   - NO --paired flag; arcasHLA auto-detects paired vs single from the input.
#   - -g A,B,C (or the class II gene list) is comma-separated.
#   - --temp takes a writable scratch dir; we make one per job and clean it up
#     with a trap. JOB_TMP is initialised empty first so the trap is safe even
#     under set -u if mktemp has not run yet.
#   - genotype output is named from the input FASTQ basename, so the class II
#     rule runs in a private temp dir and renames its output to avoid
#     clobbering the class I genotype.json.
# =============================================================================

if ARCASHLA_ENABLED:

    rule r06a_arcashla_extract:
        """
        Extract HLA-region reads from the markdup STAR BAM into FASTQs for
        arcasHLA. arcasHLA extract handles the region selection internally
        (chr6 plus the IMGT/HLA decoy contigs), so we just hand it the BAM.
        """
        input:
            bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
            bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        output:
            fq1 = "/tmp/data/06_hla/{sample}/{sample}.extracted.1.fq.gz",
            fq2 = "/tmp/data/06_hla/{sample}/{sample}.extracted.2.fq.gz",
        params:
            outdir = "/tmp/data/06_hla/{sample}",
        threads: 4
        singularity: ARCASHLA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.outdir}

            # Initialise JOB_TMP before the trap so cleanup is safe under
            # set -u even if mktemp below fails.
            JOB_TMP=""
            trap 'rm -rf "${{JOB_TMP:-}}"' EXIT INT TERM
            JOB_TMP=$(mktemp -d -p {params.outdir} arcas_extract.XXXXXX)

            # arcasHLA names its outputs from the BAM basename. The markdup BAM
            # basename is long <wildcards.sample>.Aligned.sortedByCoord.out.patched.md),
            # which would yield FASTQs that don't match our expected output
            # names. Symlink the BAM (and index) to a canonical <sample>.bam in
            # the temp dir so the output FASTQs are named <sample>.extracted.*.

            ln -sf "$(realpath {input.bam})" "$JOB_TMP/{wildcards.sample}.bam"
            ln -sf "$(realpath {input.bai})" "$JOB_TMP/{wildcards.sample}.bam.bai"

            arcasHLA extract \
                -t {threads} \
                -o {params.outdir} \
                --temp "$JOB_TMP" \
                -v \
                "$JOB_TMP/{wildcards.sample}.bam"
            """

    rule r06b_arcashla_genotype:
        """
        arcasHLA class I genotyping (A, B, C) from the extracted FASTQs.
        Emits {sample}.genotype.json with the called alleles.
        """
        input:
            fq1 = "/tmp/data/06_hla/{sample}/{sample}.extracted.1.fq.gz",
            fq2 = "/tmp/data/06_hla/{sample}/{sample}.extracted.2.fq.gz",
        output:
            json = "/tmp/data/06_hla/{sample}/{sample}.genotype.json",
        params:
            outdir = "/tmp/data/06_hla/{sample}",
        threads: 4
        singularity: ARCASHLA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.outdir}

            JOB_TMP=""
            trap 'rm -rf "${{JOB_TMP:-}}"' EXIT INT TERM
            JOB_TMP=$(mktemp -d -p {params.outdir} arcas_genotype.XXXXXX)

            # No --paired flag (arcasHLA auto-detects). -g A,B,C is the
            # comma-separated gene list. --temp gives it scratch space.
            arcasHLA genotype \
                -g A,B,C \
                -t {threads} \
                -o {params.outdir} \
                --temp "$JOB_TMP" \
                -v \
                {input.fq1} {input.fq2}
            """

    rule r06c_hla_summary:
        """
        Aggregate per-sample arcasHLA genotype.json files into one TSV with
        two alleles per class I gene. Pairs are sorted within each gene for
        determinism (arcasHLA's output order is not guaranteed).
        """
        input:
            jsons = expand("/tmp/data/06_hla/{sample}/{sample}.genotype.json",
                           sample=samples),
        output:
            tsv = "/tmp/data/06_hla/00_hla_summary.tsv",
        run:
            import json

            def parse_arcashla(path):
                # genotype.json maps gene -> list of allele strings, e.g.
                #   {"A": ["A*02:01:01", "A*30:01:01"], "B": [...], "C": [...]}
                # We keep 2-field resolution (A*02:01) and sort each pair so
                # the column assignment is stable across samples.
                with open(path) as f:
                    data = json.load(f)
                out = {}
                for gene in ("A", "B", "C"):
                    alleles = data.get(gene, []) or []
                    twofield = []
                    for a in alleles:
                        fields = a.split(":")
                        twofield.append(":".join(fields[:2]) if len(fields) >= 2 else a)
                    twofield = sorted(twofield)
                    # Pad to two entries (homozygous calls may list one allele).
                    while len(twofield) < 2:
                        twofield.append(twofield[0] if twofield else "")
                    out[f"{gene}1"] = twofield[0]
                    out[f"{gene}2"] = twofield[1]
                return out

            cols = ["A1", "A2", "B1", "B2", "C1", "C2"]
            with open(output.tsv, "w") as fout:
                fout.write("sample\t" + "\t".join(f"arcashla_{c}" for c in cols) + "\n")
                for path in input.jsons:
                    sample = os.path.basename(path).replace(".genotype.json", "")
                    try:
                        calls = parse_arcashla(path)
                    except Exception as e:
                        print(f"[HLA summary] WARN: could not parse {path}: {e}",
                              file=sys.stderr)
                        calls = {c: "" for c in cols}
                    row = [calls.get(c, "") for c in cols]
                    fout.write(sample + "\t" + "\t".join(row) + "\n")
            n = len(input.jsons)
            print(f"[HLA summary] wrote {output.tsv} -- {n} samples", flush=True)

    rule r06d_arcashla_genotype_classII:
        """
        arcasHLA class II genotyping from the same extracted FASTQs as class I
        (class II loci sit in the same chr6 MHC region, so r06a_arcashla_extract
        already pulled the relevant reads).

        Kept SEPARATE from the class I genotype rule on purpose: class II HLA is
        expressed mainly on antigen-presenting cells (in brain, chiefly
        microglia), so in brain RNA-seq class II read depth is low and highly
        variable -- much patchier than class I. Isolating it keeps the reliable
        class I summary clean and makes the (less reliable) class II calls easy
        to evaluate on their own.

        Gene set: the clinically standard DRB1, DQA1, DQB1, DPA1, DPB1. DRB3/4/5
        are copy-number-variable (present only in some DRB1 haplotypes), so
        arcasHLA may legitimately report them absent; they are not requested
        here to keep the output columns fixed.

        IMPORTANT: arcasHLA names its output {basename}.genotype.json from the
        input FASTQ basename, which is the SAME basename the class I rule uses.
        Running genotype again in the sample dir would OVERWRITE the class I
        genotype.json. To avoid that, this rule runs entirely inside a private
        temp dir and then moves the result to a distinct class II filename.
        """
        input:
            fq1 = "/tmp/data/06_hla/{sample}/{sample}.extracted.1.fq.gz",
            fq2 = "/tmp/data/06_hla/{sample}/{sample}.extracted.2.fq.gz",
        output:
            json = "/tmp/data/06_hla/{sample}/{sample}.genotype.classII.json",
        params:
            outdir = "/tmp/data/06_hla/{sample}",
        threads: 4
        singularity: ARCASHLA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.outdir}

            JOB_TMP=""
            trap 'rm -rf "${{JOB_TMP:-}}"' EXIT INT TERM
            JOB_TMP=$(mktemp -d -p {params.outdir} arcas_genotype2.XXXXXX)

            # Run genotyping with output directed into the private temp dir so
            # the class I <sample>.genotype.json in <params.outdir> is not
            # clobbered. -g lists the class II genes; no --paired (auto-detect).

            arcasHLA genotype \
                -g DRB1,DQA1,DQB1,DPA1,DPB1 \
                -t {threads} \
                -o "$JOB_TMP" \
                --temp "$JOB_TMP" \
                -v \
                {input.fq1} {input.fq2}

            # arcasHLA writes <basename>.genotype.json into $JOB_TMP. The
            # basename is derived from the FASTQ name (<sample>.extracted ->
            # arcasHLA strips known suffixes to <sample>). Find whatever
            # genotype.json it produced and move it to the class II output name.

            produced=$(find "$JOB_TMP" -maxdepth 1 -name '*.genotype.json' | head -n1)
            if [ -z "$produced" ]; then
                echo "ERROR: arcasHLA produced no genotype.json in $JOB_TMP" >&2
                ls -la "$JOB_TMP" >&2
                exit 1
            fi
            mv "$produced" {output.json}
            """

    rule r06e_hla_summary_classII:
        """
        Aggregate per-sample arcasHLA CLASS II genotype JSONs into one TSV.
        Same 2-field, sorted-pair logic as the class I summary. Expect missing
        or sparse calls for low-expression samples -- empty cells are normal for
        class II in brain RNA-seq and indicate insufficient read support, not a
        pipeline error.
        """
        input:
            jsons = expand("/tmp/data/06_hla/{sample}/{sample}.genotype.classII.json",
                           sample=samples),
        output:
            tsv = "/tmp/data/06_hla/00_hla_summary_classII.tsv",
        run:
            import json

            classII_genes = ("DRB1", "DQA1", "DQB1", "DPA1", "DPB1")

            def parse_arcashla_classII(path):
                with open(path) as f:
                    data = json.load(f)
                out = {}
                for gene in classII_genes:
                    alleles = data.get(gene, []) or []
                    twofield = []
                    for a in alleles:
                        fields = a.split(":")
                        twofield.append(":".join(fields[:2]) if len(fields) >= 2 else a)
                    twofield = sorted(twofield)
                    while len(twofield) < 2:
                        twofield.append(twofield[0] if twofield else "")
                    out[f"{gene}1"] = twofield[0]
                    out[f"{gene}2"] = twofield[1]
                return out

            cols = []
            for g in classII_genes:
                cols.extend([f"{g}1", f"{g}2"])

            with open(output.tsv, "w") as fout:
                fout.write("sample\t" + "\t".join(f"arcashla_{c}" for c in cols) + "\n")
                for path in input.jsons:
                    sample = os.path.basename(path).replace(".genotype.classII.json", "")
                    try:
                        calls = parse_arcashla_classII(path)
                    except Exception as e:
                        print(f"[HLA summary II] WARN: could not parse {path}: {e}",
                              file=sys.stderr)
                        calls = {c: "" for c in cols}
                    row = [calls.get(c, "") for c in cols]
                    fout.write(sample + "\t" + "\t".join(row) + "\n")
            n = len(input.jsons)
            print(f"[HLA summary II] wrote {output.tsv} -- {n} samples", flush=True)


rule r07a_sex_marker_regions:
    """
    Extract genomic spans (chr, start, end) for the sex-marker genes from the
    GENCODE v47 GTF, by gene_name. Written once and reused by every per-sample
    job. POSIX-awk parsing (no gawk-only match(,,arr)) for container safety.
    """
    input:
        gtf = GTF,
    output:
        regions = "/tmp/data/04_qc/sex/sex_marker_regions.tsv",
    run:
        wanted = set(SEX_MARKER_GENES)
        os.makedirs(os.path.dirname(output.regions), exist_ok=True)
        found = {}
        with open(input.gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attr = f[8]
                # Portable gene_name extraction
                key = 'gene_name "'
                i = attr.find(key)
                if i < 0:
                    continue
                gname = attr[i + len(key):].split('"', 1)[0]
                if gname in wanted and gname not in found:
                    # GTF is 1-based inclusive; samtools regions are 1-based too.
                    found[gname] = (f[0], f[3], f[4])
        missing = wanted - set(found)
        if missing:
            print(f"[sex markers] WARN: {len(missing)} marker genes not found "
                  f"in GTF (skipped): {sorted(missing)}", file=sys.stderr)
        with open(output.regions, "w") as out:
            # one row per gene: gene<TAB>chr<TAB>start<TAB>end
            for g in SEX_MARKER_GENES:
                if g in found:
                    chrom, start, end = found[g]
                    out.write(f"{g}\t{chrom}\t{start}\t{end}\n")
        print(f"[sex markers] wrote {len(found)} gene regions to {output.regions}",
              flush=True)

rule r07b_sex_counts_per_sample:
    """
    Per-sample read counts over the sex-marker gene spans, plus library size
    and chrX/chrY read fractions, from the markduplicates STAR BAM.

    Efficiency: uses only index-based access (samtools idxstats reads the .bai;
    samtools view -c <region> seeks to the gene span). No full-BAM scan, no
    per-base depth text. Each per-sample job is independent so they run in
    parallel across cores.

    Counting flags:
      -q 30  : STAR assigns MAPQ 255 to uniquely-mapped reads and 0-3 to
               multi-mappers, so -q 30 keeps ONLY unique reads. This is
               important for the Y panel: several Y genes (DDX3Y/USP9Y/UTY)
               have closely homologous X paralogs (DDX3X/USP9X/KDM6A), and
               unique-read filtering prevents X-derived reads from being
               miscounted as Y signal (and vice versa).
      -F 3332: exclude unmapped(4) + secondary(256) + duplicate(1024) +
               supplementary(2048) records, so each fragment is counted once.

    Library size (for CPM) is the total mapped read records from idxstats.
    It includes duplicates (idxstats can't filter), but that is a uniform
    per-sample scaling that does not affect the XIST-vs-Y comparison or the
    bimodal call.
    """
    input:
        bam     = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        bai     = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        regions = "/tmp/data/04_qc/sex/sex_marker_regions.tsv",
    output:
        tsv = "/tmp/data/04_qc/sex/{sample}.sex_counts.tsv",
    threads: 2
    singularity: SAMTOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})"

        IDX=$(mktemp)
        trap 'rm -f "${{IDX:-}}"' EXIT INT TERM

        # One idxstats pass (index only) -> library size + chrX/chrY records.
        samtools idxstats {input.bam} > "$IDX"
        LIB=$(awk '{{m += $3}} END {{print m + 0}}' "$IDX")
        CHRX=$(awk '$1 == "chrX" {{print $3 + 0}}' "$IDX"); [ -n "$CHRX" ] || CHRX=0
        CHRY=$(awk '$1 == "chrY" {{print $3 + 0}}' "$IDX"); [ -n "$CHRY" ] || CHRY=0

        XIST=0
        YSUM=0
        # Per-gene unique-read counts via indexed region seeks.
        while IFS=$'\t' read -r gene chr start end; do
            [ -z "$gene" ] && continue
            cnt=$(samtools view -c -q 30 -F 3332 {input.bam} "${{chr}}:${{start}}-${{end}}")
            if [ "$gene" = "XIST" ]; then
                XIST=$cnt
            else
                YSUM=$(( YSUM + cnt ))
            fi
        done < {input.regions}

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "{wildcards.sample}" "$LIB" "$CHRX" "$CHRY" "$XIST" "$YSUM" \
            > {output.tsv}
        """

rule r07c_infer_sex:
    """
    Aggregate per-sample marker counts into the cohort sex-inference table.
    Computes CPM-normalised XIST and Y-panel expression, calls genetic sex,
    and compares against the clinical 'gender' field from sample_data.txt.

    Calling logic (XIST >= SEX_XIST_CPM_MIN? / Y >= SEX_Y_CPM_MIN?):
      XIST hi, Y lo  -> "XX"
      XIST lo, Y hi  -> "XY"
      XIST hi, Y hi  -> "ambiguous_possible_XXY"  (inactive X + Y present)
      both low       -> "ambiguous_low_signal"

    Concordance vs. reported gender (male/female/other):
      - reported gender is normalised (lower-cased, trimmed); anything that
        is not exactly 'male' or 'female' (incl. 'unknown', 'NA', '' ) is
        treated as UNKNOWN and never counts as a mismatch.
      - concordance is one of:
          "concordant"      genetic XX==female or XY==male
          "DISCORDANT"      genetic XX==male or XY==female  <- identity flag
          "check_aneuploidy" genetic call is an ambiguous_* class
          "no_reported_sex" reported gender is unknown/blank
      A DISCORDANT row is a sample-swap / mislabel flag (e.g. KH36), NOT a
      biological conclusion.
    """
    input:
        counts = expand("/tmp/data/04_qc/sex/{sample}.sex_counts.tsv", sample=samples),
    output:
        tsv = "/tmp/data/04_qc/00_inferred_sex.tsv",
    run:
        # Build sample -> reported gender map from sample_data (already loaded
        # as `sample_data` with a 'name' column). Normalise defensively.
        def norm_gender(v):
            if v is None:
                return "unknown"
            s = str(v).strip().lower()
            if s in ("male", "m"):
                return "male"
            if s in ("female", "f"):
                return "female"
            return "unknown"   # 'unknown', 'na', 'nan', '', anything else

        gender_map = {}
        if "gender" in sample_data.columns:
            for _, row in sample_data.iterrows():
                gender_map[row["name"]] = norm_gender(row["gender"])
        else:
            print("[infer_sex] WARN: no 'gender' column in sample_data; "
                  "concordance will be 'no_reported_sex' for all samples",
                  file=sys.stderr)

        n_discordant = 0
        with open(output.tsv, "w") as out:
            out.write("sample\tlib_size\txist_cpm\ty_panel_cpm\t"
                      "chrX_frac\tchrY_frac\tinferred_sex\t"
                      "reported_gender\tconcordance\n")
            for path in input.counts:
                sample, lib, chrx, chry, xist, ysum = open(path).read().split()
                lib = float(lib); chrx = float(chrx); chry = float(chry)
                xist = float(xist); ysum = float(ysum)
                reported = gender_map.get(sample, "unknown")

                if lib <= 0:
                    out.write(f"{sample}\t0\tNA\tNA\tNA\tNA\t"
                              f"ambiguous_low_signal\t{reported}\t"
                              f"{'no_reported_sex' if reported=='unknown' else 'check_aneuploidy'}\n")
                    continue

                xist_cpm = xist / lib * 1e6
                y_cpm    = ysum / lib * 1e6
                chrx_frac = chrx / lib
                chry_frac = chry / lib
                xist_present = xist_cpm >= SEX_XIST_CPM_MIN
                y_present    = y_cpm    >= SEX_Y_CPM_MIN

                if xist_present and not y_present:
                    sex = "XX"
                elif y_present and not xist_present:
                    sex = "XY"
                elif xist_present and y_present:
                    sex = "ambiguous_possible_XXY"
                else:
                    sex = "ambiguous_low_signal"

                # Concordance vs. reported gender
                if reported == "unknown":
                    conc = "no_reported_sex"
                elif sex == "XX":
                    conc = "concordant" if reported == "female" else "DISCORDANT"
                elif sex == "XY":
                    conc = "concordant" if reported == "male" else "DISCORDANT"
                else:
                    conc = "check_aneuploidy"
                if conc == "DISCORDANT":
                    n_discordant += 1

                out.write(f"{sample}\t{int(lib)}\t{xist_cpm:.4f}\t{y_cpm:.4f}\t"
                          f"{chrx_frac:.6f}\t{chry_frac:.6f}\t{sex}\t"
                          f"{reported}\t{conc}\n")

        msg = f"[infer_sex] wrote {output.tsv} for {len(input.counts)} samples"
        if n_discordant:
            msg += f" -- {n_discordant} DISCORDANT sample(s) flagged for identity review"
        print(msg, flush=True)

rule r08a_cell_marker_regions:
    """
    Extract genomic spans for the cell-type marker genes from the GENCODE v47
    GTF, by gene_name. Mirrors r07a_sex_marker_regions. Written once, reused by
    every per-sample job. POSIX-safe attribute parsing (no gawk-only features).

    Emits cell_type as a column so the aggregation rule never has to re-derive
    the panel mapping, and so the panel actually used is auditable on disk.
    """
    input:
        gtf = GTF,
    output:
        regions = "/tmp/data/04_qc/cellcomp/cell_marker_regions.tsv",
    run:
        gene2type = {}
        for ct, genes in CELL_MARKER_GENES.items():
            for g in genes:
                gene2type.setdefault(g, ct)
        wanted = set(gene2type)
        os.makedirs(os.path.dirname(output.regions), exist_ok=True)

        found = {}
        with open(input.gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attr = f[8]
                key = 'gene_name "'
                i = attr.find(key)
                if i < 0:
                    continue
                gname = attr[i + len(key):].split('"', 1)[0]
                if gname in wanted and gname not in found:
                    # GTF is 1-based inclusive; samtools regions are 1-based too.
                    found[gname] = (f[0], f[3], f[4])

        missing = wanted - set(found)
        if missing:
            print(f"[cell markers] WARN: {len(missing)} marker genes not in GTF "
                  f"(skipped): {sorted(missing)}", file=sys.stderr)

        with open(output.regions, "w") as out:
            for ct, genes in CELL_MARKER_GENES.items():
                for g in genes:
                    if g in found:
                        chrom, start, end = found[g]
                        out.write(f"{g}\t{ct}\t{chrom}\t{start}\t{end}\n")
        print(f"[cell markers] wrote {len(found)} gene regions to "
              f"{output.regions}", flush=True)

rule r08b_cell_marker_counts_per_sample:
    """
    Per-sample unique-read counts over each marker gene span, plus library size.
    Mirrors r07b: index-based only (idxstats reads the .bai; view -c seeks to
    each span), so no full-BAM scan. ~35 indexed region counts per sample.

    -q 30  : STAR gives MAPQ 255 to uniquely-mapped reads, so this keeps only
             unique reads. Matters here for the paralogous marker families
             (NEFL/NEFM/NEFH; C1QA/C1QB; OLIG1/OLIG2).
    -F 3332: exclude unmapped + secondary + duplicate + supplementary.

    lib_size comes from idxstats and includes duplicates, which view -c excludes.
    That is a uniform per-sample scale factor and cancels in the per-gene z, so
    it does not affect the score -- same reasoning as r07b.
    """
    input:
        bam     = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        bai     = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        regions = "/tmp/data/04_qc/cellcomp/cell_marker_regions.tsv",
    output:
        tsv = "/tmp/data/04_qc/cellcomp/{sample}.marker_counts.tsv",
    threads: 2
    singularity: SAMTOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})"

        IDX=$(mktemp)
        trap 'rm -f "${{IDX:-}}"' EXIT INT TERM

        samtools idxstats {input.bam} > "$IDX"
        LIB=$(awk '{{m += $3}} END {{print m + 0}}' "$IDX")

        printf 'sample\tlib_size\tgene\tcell_type\tcount\n' > {output.tsv}
        while IFS=$'\t' read -r gene ctype chr start end; do
            [ -z "$gene" ] && continue
            cnt=$(samtools view -c -q 30 -F 3332 {input.bam} "${{chr}}:${{start}}-${{end}}")
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "{wildcards.sample}" "$LIB" "$gene" "$ctype" "$cnt" >> {output.tsv}
        done < {input.regions}
        """

rule r08c_cell_marker_expression:
    """
    Aggregate per-sample marker counts into one cohort table.

    SCOPE: this rule emits per-sample MEASUREMENTS only -- counts, library size,
    CPM, log2(CPM+1). Every value depends solely on its own sample, so adding or
    removing a sample never changes another sample's numbers, and the table is a
    general-purpose marker-expression artifact rather than an answer to one
    question.

    Deliberately NOT done here (all cohort-dependent or study-specific, so they
    belong in downstream analysis):
      - z-scoring across samples          (changes when the cohort changes)
      - stratifying by brain region       (parses study-specific sample names)
      - cell-type composite scores        (an analysis choice, e.g. which types
                                           count as glia)
      - outlier / donor comparisons
    Scoring recipe for the composition check is in the analysis scripts; see the
    note on r08a for why the neuron panel is non-synaptic and the glial panels
    are the independent evidence.

    Long format (one row per sample x marker) keeps it lossless and easy to
    extend with new panels; pivot downstream as needed.
    """
    input:
        counts = expand("/tmp/data/04_qc/cellcomp/{sample}.marker_counts.tsv",
                        sample=samples),
    output:
        tsv = "/tmp/data/04_qc/00_cell_marker_expression.tsv",
    run:
        import math

        order = {g: i for i, g in enumerate(
            g for ct in CELL_MARKER_GENES for g in CELL_MARKER_GENES[ct])}

        rows = []
        for path in input.counts:
            with open(path) as fh:
                next(fh)                                  # header
                for line in fh:
                    s, L, gene, ctype, c = line.rstrip("\n").split("\t")
                    lib, cnt = float(L), float(c)
                    if lib > 0:
                        cpm = cnt / lib * 1e6
                        rows.append((s, gene, ctype, int(cnt), int(lib),
                                     f"{cpm:.6f}", f"{math.log2(cpm + 1):.6f}"))
                    else:
                        # no mapped reads: emit the row, flag the metric
                        rows.append((s, gene, ctype, int(cnt), int(lib),
                                     "NA", "NA"))

        rows.sort(key=lambda r: (r[0], order.get(r[1], 10**6)))

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w") as out:
            out.write("sample\tgene\tcell_type\tcount\tlib_size\tcpm\tlog2cpm\n")
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")

        n_s = len({r[0] for r in rows})
        n_g = len({r[1] for r in rows})
        print(f"[cell markers] wrote {output.tsv}: {len(rows)} rows "
              f"({n_s} samples x {n_g} markers)", flush=True)


# ---- r09a: bulk mixture in symbol space (once; shared by all references) ------
rule r09a_mixture:
    input:
        star = expand("/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab",
                      sample=samples),
        gtf  = GTF,
    output:
        matrix = f"{DECONV_OUT}/_mixture/mixture_symbol_cpm.tsv.gz",
        genes  = f"{DECONV_OUT}/_mixture/mixture_genes.txt",
    params:
        strand = DECONV_STRAND,
        script = os.path.join(R09_SCRIPTS, "lpb-rnaseq-deconv-r09a_mixture.R"),
        pairs  = lambda wc: " ".join(
            f"{s} /tmp/data/03_bam_star/{s}/{s}.ReadsPerGene.out.tab"
            for s in samples),
    singularity: R_HSPE_CONTAINER
    shell:
        r"""
        Rscript {params.script} {input.gtf} {params.strand} \
            {output.matrix} {output.genes} {params.pairs}
        """


# ---- r09b: canonical sc reference -> hspe reference object (per ref; cache) ---
rule r09b_reference:
    input:
        mtx   = f"{DECONV_REF_DIR}/{{ref}}/matrix.mtx.gz",
        feat  = f"{DECONV_REF_DIR}/{{ref}}/features.tsv.gz",
        cells = f"{DECONV_REF_DIR}/{{ref}}/cells.tsv.gz",
        genes = rules.r09a_mixture.output.genes,
    output:
        rds = f"{DECONV_OUT}/{{ref}}/_reference/reference.rds",
        qc  = f"{DECONV_OUT}/{{ref}}/_reference/qc_markers.tsv",
    params:
        n_markers = DECONV_N_MARKERS,
        script    = os.path.join(R09_SCRIPTS, "lpb-rnaseq-deconv-r09b_reference.R"),
    singularity: R_HSPE_CONTAINER
    shell:
        r"""
        Rscript {params.script} {input.mtx} {input.feat} {input.cells} \
            {input.genes} {params.n_markers} {output.rds} {output.qc}
        """


# ---- r09c: deconvolve one sample vs one reference (per ref x sample) ----------
rule r09c_run_hspe:
    input:
        ref = rules.r09b_reference.output.rds,
        mix = rules.r09a_mixture.output.matrix,
    output:
        props = f"{DECONV_OUT}/{{ref}}/{{sample}}/proportions.tsv",
        rds   = f"{DECONV_OUT}/{{ref}}/{{sample}}/hspe.rds",
    params:
        script = os.path.join(R09_SCRIPTS, "lpb-rnaseq-deconv-r09c_run_hspe.R"),
    singularity: R_HSPE_CONTAINER
    shell:
        r"""
        Rscript {params.script} {input.ref} {input.mix} {wildcards.sample} \
            {output.props} {output.rds}
        """


# ---- r09d: summary TSV per reference (samples x classes) ---------------------
rule r09d_summary:
    input:
        props = lambda wc: expand(
            f"{DECONV_OUT}/{wc.ref}/{{sample}}/proportions.tsv", sample=samples),
    output:
        tsv = f"{DECONV_OUT}/{{ref}}.tsv",
    run:
        import csv, os
        rows, classes = {}, []
        for p in input.props:
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    rows.setdefault(r["sample"], {})[r["class"]] = r["proportion"]
                    if r["class"] not in classes:
                        classes.append(r["class"])
        classes = sorted(classes)
        os.makedirs(os.path.dirname(output.tsv) or ".", exist_ok=True)
        with open(output.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample"] + classes)
            for s in sorted(rows):
                w.writerow([s] + [rows[s].get(c, "0") for c in classes])
        print(f"[r09d] {len(rows)} samples x {len(classes)} classes -> {output.tsv}")

# ---- r09e: emit CIBERSORTx input files (per ref) -- prepare only, don't run ---
# Writes 09_deconv/_cibersortx/{ref}/{refsample.txt, mixture.txt}. Run
# CIBERSORTx yourself from that folder (see r09e_cibersortx_prep.R header). Gated
# by DECONV_CIBERSORTX=1; independent of whether hspe (r09c/d) runs.
rule r09e_cibersortx_prep:
    input:
        mtx   = f"{DECONV_REF_DIR}/{{ref}}/matrix.mtx.gz",
        feat  = f"{DECONV_REF_DIR}/{{ref}}/features.tsv.gz",
        cells = f"{DECONV_REF_DIR}/{{ref}}/cells.tsv.gz",
        mix   = rules.r09a_mixture.output.matrix,
    output:
        ref = f"{DECONV_OUT}/_cibersortx/{{ref}}/refsample.txt",
        mix = f"{DECONV_OUT}/_cibersortx/{{ref}}/mixture.txt",
    params:
        script   = os.path.join(R09_SCRIPTS, "lpb-rnaseq-deconv-r09e_cibersortx_prep.R"),
        maxcells = DECONV_CSX_MAXCELLS,
    singularity: R_HSPE_CONTAINER
    shell:
        """
        Rscript {params.script} {input.mtx} {input.feat} {input.cells} \
            {input.mix} {output.ref} {output.mix} {params.maxcells}
        """

# =============================================================================
# r10 -- cohort expression summary (counts + TPM)
#   r10a : compute the two gene x sample matrices, emit as gzip TSVs
#          (pure-Python run: block; needs nothing beyond pandas on the host)
#   r10b : render the two TSVs into a single two-sheet .xlsx workbook
#          (R + writexl, in the nikolayquadrat/toolbox:rnaseq3 container, so the
#           host Snakemake env needs no openpyxl/xlsxwriter)
#
# Both TSV/sheet matrices are genes x samples with two leading gene-annotation
# columns from the GTF:
#     gene_id   : Ensembl id WITH GENCODE version suffix (e.g. ENSG00000141510.17)
#     gene_name : HGNC symbol (human) / gene_name from the GTF for non-human
#                 annotations (common name); "" if the gene is unnamed.
#
# Strand column: STAR ReadsPerGene.out.tab has 4 columns -- col1 gene_id,
#   col2 unstranded, col3 forward, col4 reverse. We take DECONV_STRAND (default
#   2 = reverse, the GTEx / Illumina TruSeq stranded-mRNA convention), matching
#   the deconvolution stage. Override with -e DECONV_STRAND=... (0/1/2).
#
# TPM: gene length is the UNION-EXON length from the GTF (bp covered by >=1 exon,
#   overlapping/alternative exons merged, NOT summed). TPM_i = count_i / len_i,
#   normalised so each sample column sums to 1e6. Genes with no GTF exon length
#   are dropped from the TPM table (kept in counts); reported in the log.
# =============================================================================

# Container with R + writexl (for r10b). Same override convention as the
# other *_CONTAINER constants.
RNASEQ3_CONTAINER = os.environ.get(
    "RNASEQ3_CONTAINER", "docker://nikolayquadrat/toolbox:rnaseq3")


rule r10a_expression_tables:
    input:
        star = expand(
            "/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab",
            sample=samples),
        gtf = GTF,
    output:
        counts_tsv = "/tmp/data/04_qc/00_expression_counts.tsv.gz",
        tpm_tsv    = "/tmp/data/04_qc/00_expression_tpm.tsv.gz",
    params:
        strand = DECONV_STRAND,   # 0=unstranded, 1=fwd, 2=reverse (data index)
    run:
        import os, re, gzip
        import pandas as pd

        strand_col = 1 + int(params.strand)   # -> col2/3/4 of ReadsPerGene

        # -- 1. GTF: union-exon length + gene_name, keyed by full gene_id ------
        def _open(p):
            return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)

        gid_re   = re.compile(r'gene_id "([^"]+)"')
        gname_re = re.compile(r'gene_name "([^"]+)"')
        exons, names = {}, {}
        with _open(input.gtf) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("exon", "gene"):
                    continue
                m = gid_re.search(f[8])
                if not m:
                    continue
                gid = m.group(1)
                if gid not in names:
                    gn = gname_re.search(f[8])
                    if gn:
                        names[gid] = gn.group(1)
                if f[2] == "exon":
                    exons.setdefault(gid, []).append((int(f[3]), int(f[4])))

        lengths = {}
        for gid, iv in exons.items():
            iv.sort()
            merged = []
            for s_, e_ in iv:
                if merged and s_ <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e_))
                else:
                    merged.append((s_, e_))
            lengths[gid] = sum(e_ - s_ + 1 for s_, e_ in merged)
        print(f"[r10a] GTF: {len(lengths)} genes with exon length, "
              f"{len(names)} named", flush=True)

        # -- 2. per-sample STAR counts (chosen strand column) -----------------
        count_cols = {}
        for path in input.star:
            sample = os.path.basename(path).replace(".ReadsPerGene.out.tab", "")
            ids, vals = [], []
            with open(path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4 or p[0].startswith("N_"):
                        continue
                    ids.append(p[0])
                    vals.append(int(p[strand_col]))
            count_cols[sample] = pd.Series(vals, index=ids, dtype="int64")

        first = os.path.basename(input.star[0]).replace(".ReadsPerGene.out.tab", "")
        counts = pd.DataFrame(count_cols).reindex(count_cols[first].index)
        counts = counts.fillna(0).astype("int64")
        counts = counts[[s for s in samples if s in counts.columns]]
        print(f"[r10a] counts: {counts.shape[0]} genes x {counts.shape[1]} samples",
              flush=True)

        # -- 3. TPM (union-exon length; per-sample column sums to 1e6) ---------
        len_series = pd.Series(
            {g: lengths.get(g, lengths.get(g.split(".")[0])) for g in counts.index},
            dtype="float64")
        have_len = len_series.notna()
        n_missing = int((~have_len).sum())
        if n_missing:
            print(f"[r10a] WARN: {n_missing} genes have no GTF exon length; "
                  f"excluded from TPM table (present in counts table)",
                  file=sys.stderr)

        c_for_tpm = counts.loc[have_len]
        l_kb = (len_series.loc[have_len] / 1000.0)
        rpk = c_for_tpm.div(l_kb, axis=0)
        col_sums = rpk.sum(axis=0)
        tpm = rpk.div(col_sums.replace(0, pd.NA), axis=1) * 1e6
        tpm = tpm.fillna(0.0).round(4)

        # -- 4. annotation columns + write gzip TSVs (gene_id as first column) -
        def _annotate(df):
            out = df.copy()
            out.insert(0, "gene_name", [names.get(g, "") for g in out.index])
            out.index.name = "gene_id"
            return out.reset_index()          # gene_id becomes a real column

        os.makedirs(os.path.dirname(output.counts_tsv), exist_ok=True)
        _annotate(counts).to_csv(output.counts_tsv, sep="\t", index=False)
        _annotate(tpm).to_csv(output.tpm_tsv, sep="\t", index=False)

        print(f"[r10a] wrote:\n  {output.counts_tsv} "
              f"({counts.shape[0]}x{counts.shape[1]})\n  {output.tpm_tsv} "
              f"({tpm.shape[0]}x{tpm.shape[1]}) [strand col {strand_col}]",
              flush=True)

        # -- diagnostics: byte sizes + head preview, loud warning if empty -----
        if counts.shape[0] == 0:
            print("[r10a] WARN: counts table has 0 GENES -- the .xlsx will be "
                  "empty. Check that the ReadsPerGene.out.tab files contain "
                  "gene rows (not just the 4 N_* header lines) and that the "
                  "strand column index is valid.", file=sys.stderr)
        if counts.shape[1] == 0:
            print("[r10a] WARN: 0 SAMPLES matched -- check `samples` vs the "
                  "STAR output basenames.", file=sys.stderr)
        for _p in (output.counts_tsv, output.tpm_tsv):
            print(f"[r10a] {_p}: {os.path.getsize(_p)} bytes", flush=True)
        with gzip.open(output.counts_tsv, "rt") as _fh:
            _head = [next(_fh, "").rstrip("\n") for _ in range(3)]
        print("[r10a] counts head (first 3 lines):\n  " + "\n  ".join(_head),
              flush=True)


rule r10b_expression_xlsx:
    """
    Render the two r10a TSVs into a single .xlsx with sheets 'counts' and 'TPM',
    using R + writexl inside the rnaseq3 container (so the host env needs no
    Python Excel writer). writexl::write_xlsx takes a named list of data frames
    -> one sheet per name. read.delim keeps gene_id/gene_name as plain character
    columns; check.names=FALSE preserves sample column names verbatim.
    """
    input:
        counts_tsv = rules.r10a_expression_tables.output.counts_tsv,
        tpm_tsv    = rules.r10a_expression_tables.output.tpm_tsv,
    output:
        xlsx = "/tmp/data/04_qc/00_expression_summary.xlsx",
    singularity: RNASEQ3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # --- shell-side: prove the container sees non-empty input TSVs --------
        echo "[r10b] inputs on disk (inside container):"
        ls -l {input.counts_tsv} {input.tpm_tsv}
        # decompressed line counts (header + one row per gene). '|| true' guards
        # against SIGPIPE from head under pipefail.
        CN=$(gzip -dc {input.counts_tsv} | wc -l)
        TN=$(gzip -dc {input.tpm_tsv}    | wc -l)
        echo "[r10b] counts TSV: $CN lines (incl header); tpm TSV: $TN lines"
        echo "[r10b] counts TSV head:"
        gzip -dc {input.counts_tsv} | head -3 || true
        if [ "$CN" -le 1 ]; then
            echo "[r10b] ERROR: counts TSV has no data rows (<=1 line). The" >&2
            echo "       emptiness originates in r10a, not here. Re-run r10a" >&2
            echo "       and read its [r10a] diagnostics." >&2
            exit 1
        fi

        # --- R side: read, report dimensions, refuse to write an empty book ---
        Rscript -e '
            suppressPackageStartupMessages(library(writexl));
            args <- commandArgs(trailingOnly=TRUE);
            cat(sprintf("[r10b][R] received %d args: %s\n",
                        length(args), paste(args, collapse=" | ")));
            if (length(args) < 3) stop("[r10b][R] need 3 args: counts_tsv tpm_tsv out_xlsx");
            counts <- read.delim(gzfile(args[1]), header=TRUE, sep="\t",
                                 check.names=FALSE, stringsAsFactors=FALSE,
                                 quote="", comment.char="");
            tpm    <- read.delim(gzfile(args[2]), header=TRUE, sep="\t",
                                 check.names=FALSE, stringsAsFactors=FALSE,
                                 quote="", comment.char="");
            cat(sprintf("[r10b][R] counts parsed: %d rows x %d cols (%s)\n",
                        nrow(counts), ncol(counts),
                        paste(head(colnames(counts), 6), collapse=",")));
            cat(sprintf("[r10b][R] tpm    parsed: %d rows x %d cols\n",
                        nrow(tpm), ncol(tpm)));
            if (nrow(counts) == 0 || nrow(tpm) == 0)
                stop("[r10b][R] a parsed table has 0 rows -- refusing to write an empty workbook");
            if (ncol(counts) < 3 || ncol(tpm) < 3)
                stop("[r10b][R] <3 columns -- expected gene_id, gene_name, and >=1 sample");
            write_xlsx(list(counts=counts, TPM=tpm), path=args[3]);
            cat(sprintf("[r10b][R] wrote %s (%d bytes)\n",
                        args[3], file.info(args[3])$size));
        ' {input.counts_tsv} {input.tpm_tsv} {output.xlsx}

        echo "[r10b] output on disk:"
        ls -l {output.xlsx}
        """


# -----------------------------------------------------------------------------
# r11a: chip .bim -> filtered site list, in chip coordinates
# -----------------------------------------------------------------------------

rule r11a_genocheck_sites:
    """
    Keep only sites that can be both genotyped from RNA and unambiguously
    matched back to the array:

      * single-letter ACGT alleles on BOTH sides. This drops indels, CNV and
        multi-base probes, PLINK I/D coding, '-' and monomorphic '0' entries in
        one test -- nothing that is not a clean SNP survives.
      * biallelic (a1 != a2)
      * NOT strand-ambiguous: A/T and C/G pairs are excluded outright. For
        these, a strand flip is indistinguishable from a genotype swap, which
        is the exact error this QC exists to detect.
      * autosomes (+chrX if GENOCHECK_INCLUDE_X); PLINK 23/24/25/26 codes
        resolved; the ORIGINAL chip chromosome string is carried through
        unchanged so the final fileset merges with the chip without recoding.

    Transitions (A/G, C/T) are deliberately NOT dropped. They are ADAR-editing
    prone, but editing inflates het calls at a modest rate while costing 2/3 of
    all SNPs to exclude -- a bad trade when the surviving set only has to
    establish identity.
    """
    input:
        bim = GENOCHECK_CHIP_PREFIX + ".bim",
    output:
        sites = GENOCHECK_OUT + "/sites_chipbuild.tsv",
        bed   = GENOCHECK_OUT + "/sites_chipbuild.bed",
        stats = GENOCHECK_OUT + "/sites_chipbuild.log",
    run:
        import numpy as np
        import pandas as pd

        AMBIG = {frozenset(("A", "T")), frozenset(("C", "G"))}
        ACGT = list("ACGT")
        # hg38 contig names we will end up genotyping on, as bare numbers/letters
        want = {c.replace("chr", "") for c in GENOCHECK_CHROMS}
        counts = {}

        bim = pd.read_csv(input.bim, sep=r"\s+", header=None,
                          names=["chrom", "id", "cm", "bp", "a1", "a2"],
                          dtype={"chrom": str, "id": str, "bp": np.int64,
                                 "a1": str, "a2": str})
        counts["chip_variants"] = len(bim)

        # normalised contig (for liftover + selection); chrom_orig preserved
        bim["chrom_orig"] = bim["chrom"]
        c = bim["chrom"].str.replace(r"^chr", "", regex=True, case=False)
        c = c.replace({"23": "X", "24": "Y", "25": "X", "26": "MT"})
        bim["contig"] = c
        bim = bim[bim["contig"].isin(want)]
        counts["on_requested_chroms"] = len(bim)

        bim["a1"] = bim["a1"].str.upper().str.strip()
        bim["a2"] = bim["a2"].str.upper().str.strip()
        # single-letter ACGT on both sides, and biallelic
        ok = bim["a1"].isin(ACGT) & bim["a2"].isin(ACGT) & (bim["a1"] != bim["a2"])
        bim = bim[ok]
        counts["single_letter_biallelic_snp"] = len(bim)

        pairs = [frozenset((x, y)) for x, y in zip(bim["a1"], bim["a2"])]
        bim = bim[[p not in AMBIG for p in pairs]]
        counts["after_drop_AT_and_CG"] = len(bim)

        bim = bim.drop_duplicates(subset=["contig", "bp"], keep="first")
        bim = bim.drop_duplicates(subset=["id"], keep="first")
        counts["after_dedup"] = len(bim)

        bim = bim.sort_values(["contig", "bp"])
        bim[["chrom_orig", "contig", "bp", "a1", "a2", "id"]].to_csv(
            output.sites, sep="\t", header=False, index=False)

        # BED for liftover. Name field packs everything needed to rejoin, with
        # '|' as separator (never present in rsIDs or PLINK chrom codes).
        with open(output.bed, "w") as fh:
            for co, ct, bp, a1, a2, vid in bim[
                    ["chrom_orig", "contig", "bp", "a1", "a2", "id"]].itertuples(
                        index=False):
                fh.write(f"chr{ct}\t{bp - 1}\t{bp}\t"
                         f"{vid}|{co}|{bp}|{a1}|{a2}\n")

        with open(output.stats, "w") as fh:
            for k, v in counts.items():
                fh.write(f"{k}\t{v}\n")
                print(f"[genocheck][sites] {k}: {v}", file=sys.stderr)

        if len(bim) < 1000:
            raise ValueError(
                f"only {len(bim)} usable chip SNPs; expected tens of thousands. "
                "Check the .bim allele coding and chromosome naming.")


# -----------------------------------------------------------------------------
# r11b: lift the site list hg19 -> hg38  (no-op when the chip is already hg38)
# -----------------------------------------------------------------------------

rule r11b_genocheck_liftover:
    """
    CrossMap on a 1 bp-per-row BED. Only the SITE LIST is lifted; genotypes
    never are. CrossMap writes successes to the output BED and failures to
    <output>.unmap, which is counted in r11c.
    """
    input:
        bed = GENOCHECK_OUT + "/sites_chipbuild.bed",
        **({"chain": GENOCHECK_CHAIN} if GENOCHECK_NEEDS_LIFT else {}),
    output:
        bed   = GENOCHECK_OUT + "/sites_hg38_raw.bed",
        unmap = GENOCHECK_OUT + "/sites_hg38_raw.bed.unmap",
    params:
        needs_lift = GENOCHECK_NEEDS_LIFT,
        chain      = GENOCHECK_CHAIN,
    singularity: CROSSMAP_CONTAINER
    shell:
        r"""
        set -euo pipefail
        if [ "{params.needs_lift}" = "True" ]; then
            # CrossMap >=0.6 exposes `CrossMap`; older versions `CrossMap.py`
            if command -v CrossMap >/dev/null 2>&1; then CM=CrossMap;
            else CM=CrossMap.py; fi
            $CM bed {params.chain} {input.bed} {output.bed}
            # CrossMap writes failures to <out>.unmap; guarantee it exists
            touch {output.unmap}
        else
            cp {input.bed} {output.bed}
            : > {output.unmap}
        fi
        """


rule r11c_genocheck_crosswalk:
    """
    Build the hg19 <-> hg38 crosswalk actually used downstream, and restrict to
    transcribed sequence.

      * a site must map to exactly ONE hg38 location, 1 bp wide. Sites with
        multiple mappings are ambiguous and are dropped -- for identity work a
        wrong position is worse than a missing one.
      * the hg38 position must fall in a union-exon interval of the annotation
        (overlapping/alternative exons merged, membership by binary search). A
        site in an intron or intergenic region has no RNA coverage to call from,
        so keeping it would only add missingness.

    Emits the pileup inputs (BED + regions) in hg38, and the crosswalk keyed on
    hg38 chrom:pos that r11i uses to write the result back out in chip coords.
    """
    input:
        lifted = GENOCHECK_OUT + "/sites_hg38_raw.bed",
        unmap  = GENOCHECK_OUT + "/sites_hg38_raw.bed.unmap",
        gtf    = GTF,
    output:
        crosswalk = GENOCHECK_OUT + "/00_genocheck_crosswalk.tsv.gz",
        regions   = GENOCHECK_OUT + "/sites_regions.txt",
        stats     = GENOCHECK_OUT + "/crosswalk.log",
    run:
        import gzip
        import numpy as np
        from collections import Counter

        keep_chroms = set(GENOCHECK_CHROMS)
        counts = {}

        # ---- parse lifted BED; reject non-unique and non-1bp mappings
        name_hits = Counter()
        rows = []
        with open(input.lifted) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                contig38, s, e, name = f[0], int(f[1]), int(f[2]), f[3]
                name_hits[name] += 1
                rows.append((contig38, s, e, name))
        counts["lifted_rows"] = len(rows)

        counts["failed_to_lift"] = sum(1 for _ in open(input.unmap)) \
            if os.path.isfile(input.unmap) else 0

        multi = {n for n, k in name_hits.items() if k > 1}
        counts["dropped_multi_mapping"] = len(multi)

        recs = []
        n_badwidth = n_offchrom = 0
        for contig38, s, e, name in rows:
            if name in multi:
                continue
            if e - s != 1:
                n_badwidth += 1
                continue
            if contig38 not in keep_chroms:
                n_offchrom += 1
                continue
            vid, chrom_orig, bp_chip, a1, a2 = name.split("|")
            recs.append((contig38, e, vid, chrom_orig, int(bp_chip), a1, a2))
        counts["dropped_not_1bp"] = n_badwidth
        counts["dropped_off_target_contig"] = n_offchrom
        counts["unique_mapped_sites"] = len(recs)

        # ---- union exons per contig (hg38; the GTF and the lifted sites agree)
        raw = {}
        opener = gzip.open if str(input.gtf).endswith(".gz") else open
        with opener(input.gtf, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 5)
                if len(f) < 5 or f[2] != "exon" or f[0] not in keep_chroms:
                    continue
                raw.setdefault(f[0], []).append((int(f[3]), int(f[4])))
        starts_by, ends_by = {}, {}
        for contig, ivs in raw.items():
            ivs.sort()
            ms, me = [], []
            cs, ce = ivs[0]
            for s, e in ivs[1:]:
                if s <= ce + 1:
                    ce = max(ce, e)
                else:
                    ms.append(cs); me.append(ce); cs, ce = s, e
            ms.append(cs); me.append(ce)
            starts_by[contig] = np.asarray(ms, dtype=np.int64)
            ends_by[contig] = np.asarray(me, dtype=np.int64)
        counts["exon_intervals_merged"] = int(sum(len(v) for v in starts_by.values()))

        def in_exon(contig, pos):
            s = starts_by.get(contig)
            if s is None:
                return False
            e = ends_by[contig]
            i = int(np.searchsorted(s, pos, side="right")) - 1
            return i >= 0 and pos <= e[i]

        recs = [r for r in recs if in_exon(r[0], r[1])]
        counts["in_union_exon"] = len(recs)

        # final ordering: hg38 chrom-major, position-ascending (mpileup needs it)
        order = {c: i for i, c in enumerate(GENOCHECK_CHROMS)}
        recs.sort(key=lambda r: (order[r[0]], r[1]))

        with gzip.open(output.crosswalk, "wt") as cw, \
                open(output.regions, "w") as rg:
            cw.write("contig38\tpos38\tvariant_id\tchrom_chip\tpos_chip\t"
                     "chip_a1\tchip_a2\n")
            for contig38, pos38, vid, chrom_orig, bp_chip, a1, a2 in recs:
                cw.write(f"{contig38}\t{pos38}\t{vid}\t{chrom_orig}\t"
                         f"{bp_chip}\t{a1}\t{a2}\n")
                rg.write(f"{contig38}:{pos38}-{pos38}\n")

        with open(output.stats, "w") as fh:
            for k, v in counts.items():
                fh.write(f"{k}\t{v}\n")
                print(f"[genocheck][crosswalk] {k}: {v}", file=sys.stderr)

        if counts["lifted_rows"] and \
                counts["failed_to_lift"] > 0.2 * counts["lifted_rows"]:
            print("[genocheck][crosswalk] WARNING: >20% of sites failed to "
                  "lift. Check that the chain file matches the chip's build "
                  "direction (hg19->hg38, not the reverse).", file=sys.stderr)
        if len(recs) < 1000:
            raise ValueError(
                f"only {len(recs)} sites survived liftover + exon restriction; "
                "expected tens of thousands.")


# -----------------------------------------------------------------------------
# r11d / r11e: reference allele at each hg38 site
# -----------------------------------------------------------------------------

rule r11d_genocheck_refbase:
    """
    Reference base at every hg38 site, in batched samtools faidx calls (a
    per-site call would be ~50k process spawns; faidx --region-file is not
    available in older samtools, so regions go in as arguments).
    """
    input:
        regions = GENOCHECK_OUT + "/sites_regions.txt",
        fasta   = GENOME_FA,
    output:
        fa = temp(GENOCHECK_OUT + "/sites_refbase.fa"),
    singularity: SAMTOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        : > {output.fa}
        tr '\n' '\0' < {input.regions} \
          | xargs -0 -n 2000 samtools faidx {input.fasta} >> {output.fa}
        """


rule r11e_genocheck_alleles:
    """
    Crosswalk + reference bases -> the allele file for `bcftools call -C alleles`
    (CHROM POS REF,ALT in hg38 orientation).

    A .bim gives A1/A2 without saying which is reference, and liftover can land
    a site in an inverted region, so orientation is resolved explicitly:
      * REF == A1 -> ALT = A2;  REF == A2 -> ALT = A1
      * neither -> complement both and retry. A successful retry means the site
        is strand-flipped relative to the chip manifest (manifest convention or
        an inverted liftover block). It is recorded in the `flip` column so
        r11i can complement the alleles back to the chip's orientation.
      * still neither -> the site disagrees with the reference; dropped. A large
        count here means the chain or the declared build is wrong.
    """
    input:
        crosswalk = GENOCHECK_OUT + "/00_genocheck_crosswalk.tsv.gz",
        fa        = GENOCHECK_OUT + "/sites_refbase.fa",
    output:
        alleles = GENOCHECK_OUT + "/sites_alleles.tsv",
        bed     = GENOCHECK_OUT + "/sites_hg38.bed",
        table   = GENOCHECK_OUT + "/sites_final.tsv",
        stats   = GENOCHECK_OUT + "/sites_alleles.log",
    run:
        import gzip

        COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

        refbase, key = {}, None
        with open(input.fa) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line[0] == ">":
                    contig, _, rng = line[1:].partition(":")
                    key = (contig, int(rng.split("-")[0]))
                    refbase[key] = ""
                elif key is not None:
                    refbase[key] += line.upper()

        n_in = n_out = n_flip = n_mismatch = n_badref = 0
        out_rows = []
        with gzip.open(input.crosswalk, "rt") as fh:
            next(fh)
            for line in fh:
                contig38, pos38, vid, chrom_chip, pos_chip, a1, a2 = \
                    line.rstrip("\n").split("\t")
                pos38 = int(pos38)
                n_in += 1
                ref = refbase.get((contig38, pos38), "")
                if len(ref) != 1 or ref not in COMP:
                    n_badref += 1
                    continue
                flip = 0
                if ref == a1:
                    alt = a2
                elif ref == a2:
                    alt = a1
                else:
                    c1, c2 = COMP.get(a1), COMP.get(a2)
                    if ref == c1:
                        alt, flip = c2, 1; n_flip += 1
                    elif ref == c2:
                        alt, flip = c1, 1; n_flip += 1
                    else:
                        n_mismatch += 1
                        continue
                out_rows.append((contig38, pos38, ref, alt, vid,
                                 chrom_chip, pos_chip, a1, a2, flip))
                n_out += 1

        with open(output.alleles, "w") as fa_out, \
                open(output.bed, "w") as bed_out, \
                open(output.table, "w") as tab_out:
            tab_out.write("contig38\tpos38\tref38\talt38\tvariant_id\t"
                          "chrom_chip\tpos_chip\tchip_a1\tchip_a2\tflip\n")
            for r in out_rows:
                contig38, pos38, ref, alt = r[0], r[1], r[2], r[3]
                fa_out.write(f"{contig38}\t{pos38}\t{ref},{alt}\n")
                bed_out.write(f"{contig38}\t{pos38 - 1}\t{pos38}\n")
                tab_out.write("\t".join(str(x) for x in r) + "\n")

        stats = {
            "crosswalk_sites_in": n_in,
            "sites_out": n_out,
            "strand_flipped": n_flip,
            "dropped_ref_allele_mismatch": n_mismatch,
            "dropped_non_acgt_reference": n_badref,
        }
        with open(output.stats, "w") as fh:
            for k, v in stats.items():
                fh.write(f"{k}\t{v}\n")
                print(f"[genocheck][alleles] {k}: {v}", file=sys.stderr)

        if n_in and n_mismatch / n_in > 0.05:
            raise ValueError(
                f"{n_mismatch}/{n_in} sites ({100*n_mismatch/n_in:.1f}%) match "
                "the hg38 reference in neither orientation. That is the "
                "signature of a wrong chain file or a mis-declared chip build.")


# -----------------------------------------------------------------------------
# r11f: joint pileup + constrained calling, scattered by hg38 chromosome
# -----------------------------------------------------------------------------

rule r11f_genocheck_call_chrom:
    """
    One joint mpileup across ALL samples per chromosome, so every sample is
    assessed at the same sites and missingness is interpretable.

    mpileup flags:
      -q {MQ}    MAPQ floor; 255 = STAR uniquely-mapped only
      -Q {BQ}    base-quality floor
      -d {DEPTH} per-file depth cap, far above default -- RNA depth at
                 expressed exons is routinely in the thousands
      -B         disable BAQ: it over-penalises reads near splice junctions and
                 systematically loses RNA coverage
      -I         skip indels (unreliable in RNA, absent from the array)
      --ff       exclude unmapped/secondary/QC-fail (+duplicates by default)
      -a AD,DP   per-sample allele depths, for QC
    call flags:
      -m             multiallelic caller
      -C alleles -T  genotype ONLY the given REF,ALT pair; no novel alleles
    """
    input:
        bams = expand("/tmp/data/03_bam_star/{sample}/"
                      "{sample}.Aligned.sortedByCoord.out.patched.md.bam",
                      sample=samples),
        bais = expand("/tmp/data/03_bam_star/{sample}/"
                      "{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
                      sample=samples),
        alleles = GENOCHECK_OUT + "/sites_alleles.tsv",
        bed     = GENOCHECK_OUT + "/sites_hg38.bed",
        fasta   = GENOME_FA,
    output:
        bcf = temp(GENOCHECK_OUT + "/by_chrom/{chrom}.bcf"),
    params:
        mq = GENOCHECK_MIN_MQ, bq = GENOCHECK_MIN_BQ,
        depth = GENOCHECK_MAX_DEPTH, ff = GENOCHECK_FF,
        wd = GENOCHECK_OUT + "/by_chrom",
    threads: 4
    singularity: BCFTOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        P={params.wd}/{wildcards.chrom}

        awk -v c={wildcards.chrom} '$1==c' {input.bed}     > $P.bed
        awk -v c={wildcards.chrom} '$1==c' {input.alleles} > $P.alleles.tsv

        if [ ! -s $P.bed ]; then
            echo "[genocheck][ERR] no sites on {wildcards.chrom}" >&2; exit 1
        fi

        bgzip -f $P.alleles.tsv
        tabix -s1 -b2 -e2 -f $P.alleles.tsv.gz

        # sample names come from the @RG SM tag r03d writes
        printf "%s\n" {input.bams} > $P.bamlist

        bcftools mpileup \
            -b $P.bamlist -f {input.fasta} -R $P.bed \
            -q {params.mq} -Q {params.bq} -d {params.depth} \
            -B -I --ff {params.ff} \
            -a FORMAT/AD,FORMAT/DP -Ou \
        | bcftools call \
            -m -C alleles -T $P.alleles.tsv.gz \
            --threads {threads} -Ob -o {output.bcf}

        bcftools index -f {output.bcf}
        """


# -----------------------------------------------------------------------------
# r11g: confidence filters -> hg38 VCF
# -----------------------------------------------------------------------------

rule r11g_genocheck_filter:
    """
    Genotype level first: any call with DP < MIN_DP or GQ < MIN_GQ becomes
    MISSING rather than the site being dropped. A shallow call is uninformative,
    not wrong, and blanking it keeps the site usable for the samples that do
    have depth -- which matters because RNA depth spans orders of magnitude
    between genes. Then site level: drop sites missing in more than
    MAX_MISSING of samples.

    IDs are set to hg38 CHROM:POS here purely as a join key for r11i, which
    replaces them with the chip's variant IDs.

    NOT applied by default: an allelic-balance filter on hets. Allele-specific
    and mono-allelic expression is pervasive in brain RNA, so a 5:95 het may be
    real biology; blanking those would remove mostly true genotypes. The
    commented line adds it -- test the expression against your bcftools build
    first.
    """
    input:
        bcfs = expand(GENOCHECK_OUT + "/by_chrom/{chrom}.bcf",
                      chrom=GENOCHECK_CHROMS),
    output:
        vcf = GENOCHECK_OUT + "/rnaseq_genotypes.hg38.vcf.gz",
        tbi = GENOCHECK_OUT + "/rnaseq_genotypes.hg38.vcf.gz.tbi",
    params:
        min_dp = GENOCHECK_MIN_DP, min_gq = GENOCHECK_MIN_GQ,
        max_missing = GENOCHECK_MAX_MISSING,
    threads: 4
    singularity: BCFTOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        bcftools concat -a -Ou {input.bcfs} \
        | bcftools +setGT -Ou -- -t q \
              -i 'FMT/DP<{params.min_dp} || FMT/GQ<{params.min_gq}' -n . \
        | bcftools view -i 'F_MISSING<={params.max_missing}' -Ou \
        | bcftools annotate --set-id '%CHROM:%POS' \
              --threads {threads} -Oz -o {output.vcf}

        # stricter alternative -- blank hets with extreme allelic imbalance:
        #   | bcftools +setGT -Ou -- -t q \
        #         -i 'GT="het" && (FMT/AD[:1])/(FMT/DP)<0.10' -n . \

        bcftools index -t -f {output.vcf}
        """


# -----------------------------------------------------------------------------
# r11h / r11i / r11j: PLINK fileset, remapped to chip coordinates
# -----------------------------------------------------------------------------

rule r11h_genocheck_plink_hg38:
    """VCF -> interim PLINK fileset, still in hg38 with CHROM:POS IDs."""
    input:
        vcf = GENOCHECK_OUT + "/rnaseq_genotypes.hg38.vcf.gz",
    output:
        bed = temp(GENOCHECK_OUT + "/_hg38/rna_hg38.bed"),
        bim = temp(GENOCHECK_OUT + "/_hg38/rna_hg38.bim"),
        fam = temp(GENOCHECK_OUT + "/_hg38/rna_hg38.fam"),
    params:
        out_prefix = GENOCHECK_OUT + "/_hg38/rna_hg38",
    threads: 4
    singularity: PLINK2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.out_prefix})
        plink2 --vcf {input.vcf} \
            --double-id --vcf-half-call m --allow-extra-chr \
            --snps-only just-acgt --max-alleles 2 \
            --set-all-var-ids '@:#' \
            --threads {threads} --make-bed --out {params.out_prefix}
        """


rule r11i_genocheck_remap:
    """
    Rewrite the .bim into the chip's coordinate system. The .bed genotype matrix
    and .fam are untouched (variant ORDER is unchanged, only the labels), so
    this is a pure relabelling:

        hg38 chrom:pos  ->  chip chrom, chip bp, chip variant ID

    Alleles are complemented for the sites r11e flagged as strand-flipped, so
    the fileset comes out in the chip's allele orientation and merges cleanly.

    Sample IDs are rewritten to chip IIDs here too if sample_map.tsv exists;
    without it they stay as RNA sample names (kinship still works, the pairs
    are just labelled by RNA name).

    The result is position-UNSORTED, because hg19 and hg38 order variants
    differently. r11j sorts it.
    """
    input:
        bed   = GENOCHECK_OUT + "/_hg38/rna_hg38.bed",
        bim   = GENOCHECK_OUT + "/_hg38/rna_hg38.bim",
        fam   = GENOCHECK_OUT + "/_hg38/rna_hg38.fam",
        table = GENOCHECK_OUT + "/sites_final.tsv",
    output:
        bed = temp(GENOCHECK_OUT + "/_remap/rna_remap.bed"),
        bim = temp(GENOCHECK_OUT + "/_remap/rna_remap.bim"),
        fam = temp(GENOCHECK_OUT + "/_remap/rna_remap.fam"),
        keep = GENOCHECK_OUT + "/_remap/remapped_ids.txt",
        log  = GENOCHECK_OUT + "/remap.log",
    params:
        sample_map = GENOCHECK_SAMPLE_MAP,
    run:
        import shutil

        COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "0": "0"}

        xw = {}
        with open(input.table) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                xw[f[0] + ":" + f[1]] = {
                    "vid": f[4], "chrom": f[5], "bp": f[6], "flip": f[9] == "1"}

        n_in = n_mapped = n_flipped = 0
        kept_ids = []
        with open(input.bim) as fin, open(output.bim, "w") as fout:
            for line in fin:
                chrom, vid, cm, bp, a1, a2 = line.rstrip("\n").split("\t")
                n_in += 1
                rec = xw.get(vid)
                if rec is None:
                    # Not in the crosswalk: keep the row (the .bed matrix must
                    # stay aligned) but mark it so r11j can drop it.
                    fout.write(f"{chrom}\t__DROP__{vid}\t{cm}\t{bp}\t{a1}\t{a2}\n")
                    continue
                if rec["flip"]:
                    a1, a2 = COMP.get(a1, a1), COMP.get(a2, a2)
                    n_flipped += 1
                fout.write(f"{rec['chrom']}\t{rec['vid']}\t{cm}\t"
                           f"{rec['bp']}\t{a1}\t{a2}\n")
                kept_ids.append(rec["vid"])
                n_mapped += 1

        shutil.copy(input.bed, output.bed)

        # fam: optionally rename samples to their chip IIDs
        remap = {}
        if os.path.isfile(params.sample_map) and \
                os.path.getsize(params.sample_map) > 0:
            with open(params.sample_map) as fh:
                for line in fh:
                    f = line.split()
                    if len(f) >= 2:
                        remap[f[0]] = f[1]
        n_renamed = 0
        with open(input.fam) as fin, open(output.fam, "w") as fout:
            for line in fin:
                f = line.split()
                if f[1] in remap:
                    new = remap[f[1]]
                    f[0], f[1] = new, new
                    n_renamed += 1
                fout.write(" ".join(f) + "\n")

        with open(output.keep, "w") as fh:
            fh.write("\n".join(kept_ids) + "\n")

        stats = {"bim_rows_in": n_in, "remapped_to_chip_build": n_mapped,
                 "alleles_complemented": n_flipped,
                 "dropped_not_in_crosswalk": n_in - n_mapped,
                 "samples_renamed_to_chip_iid": n_renamed}
        with open(output.log, "w") as fh:
            for k, v in stats.items():
                fh.write(f"{k}\t{v}\n")
                print(f"[genocheck][remap] {k}: {v}", file=sys.stderr)


rule r11j_genocheck_plink_final:
    """
    Drop the unmapped rows and sort variants into chip-build order.
    This is the deliverable: rnaseq_genotypes.{bed,bim,fam}, in the chip's
    build, with the chip's variant IDs.
    """
    input:
        bed  = GENOCHECK_OUT + "/_remap/rna_remap.bed",
        bim  = GENOCHECK_OUT + "/_remap/rna_remap.bim",
        fam  = GENOCHECK_OUT + "/_remap/rna_remap.fam",
        keep = GENOCHECK_OUT + "/_remap/remapped_ids.txt",
    output:
        bed = GENOCHECK_OUT + "/rnaseq_genotypes.bed",
        bim = GENOCHECK_OUT + "/rnaseq_genotypes.bim",
        fam = GENOCHECK_OUT + "/rnaseq_genotypes.fam",
    params:
        in_prefix  = GENOCHECK_OUT + "/_remap/rna_remap",
        out_prefix = GENOCHECK_OUT + "/rnaseq_genotypes",
    threads: 4
    singularity: PLINK2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        plink2 --bfile {params.in_prefix} \
            --extract {input.keep} \
            --allow-extra-chr --sort-vars \
            --threads {threads} --make-bed --out {params.out_prefix}
        """


# =============================================================================
# KINSHIP AND IDENTITY
# =============================================================================
# Produced unconditionally (needs only the site list):
#   r11k/r11l  RNA x RNA kinship + QC. In a multi-region design this is the
#              donor-structure check: two regions from the same donor must look
#              like duplicates, different donors like unrelated.
#
# Produced only at the IDENTITY tier (.bed/.bim/.fam all present):
#   r11m..r11q merge, cross-dataset kinship, all-pairs IBS, similarity report.
#
# SAMPLE NAMES ARE NOT ASSUMED TO CORRESPOND. Nothing here matches RNA samples
# to array samples by name. Every RNA sample is compared against every array
# sample and ranked by genotype similarity, so the assignment is discovered from
# the data. If sample_map.tsv is supplied it is treated as a HYPOTHESIS: r11q
# reports whether the data-driven best match agrees with it, and flags it when
# it does not. RNA IDs are suffixed '_RNA' before the merge so that identical
# names on both sides cannot collide or be mistaken for each other.
#
# HOW TO READ THESE NUMBERS -- the thresholds are NOT the textbook ones. Allele
# dropout at modest depth turns true hets into apparent homozygotes, which
# biases kinship and concordance DOWNWARD. Expect:
#   * RNA vs RNA, same donor: often near the 0.5 duplicate expectation, because
#     two samples from one donor express the same genes and so drop out at the
#     same sites -- the errors are CORRELATED and partly cancel.
#   * RNA vs array, same person: concordance typically 0.85-0.97 and kinship
#     0.30-0.45 rather than 0.5, because dropout is one-sided. That is a MATCH,
#     not a half-sib.
#   * unrelated pairs sit near 0 kinship in both comparisons; this is the robust
#     part of the signal.
# IBS0 (opposite-homozygote rate) is the sturdiest discriminator: dropout turns
# a het into one of the two homozygotes, which produces IBS1, not IBS0, so IBS0
# stays near zero for a true match even when kinship has been deflated. Judge by
# rank, margin and IBS0 rather than by an absolute cutoff.
# =============================================================================

rule r11k_genocheck_kinship_rna:
    """
    KING-robust kinship among the RNA samples only. --make-king-table also
    reports IBS0 and NSNP, which is what separates a true duplicate (IBS0 ~ 0)
    from a first-degree relative (IBS0 > 0) when dropout has deflated kinship.

    LD pruning is OFF by default: the site set is already small (tens of
    thousands of exonic SNPs) and pruning costs more precision than the LD
    inflation it removes. GENOCHECK_LD_PRUNE=1 enables it.
    """
    input:
        bed = GENOCHECK_OUT + "/rnaseq_genotypes.bed",
    output:
        kin = GENOCHECK_OUT + "/00_kinship_rna.kin0",
    params:
        prefix   = GENOCHECK_OUT + "/rnaseq_genotypes",
        wd       = GENOCHECK_OUT + "/_kinship",
        ld_prune = GENOCHECK_LD_PRUNE,
    threads: 4
    singularity: PLINK2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        T={params.wd}/rna
        if [ "{params.ld_prune}" = "True" ]; then
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --indep-pairwise 200kb 1 0.2 --threads {threads} \
                --out $T.prune >/dev/null
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --extract $T.prune.prune.in \
                --make-king-table --king-table-filter -1 \
                --threads {threads} --out $T >/dev/null
        else
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --make-king-table --king-table-filter -1 \
                --threads {threads} --out $T >/dev/null
        fi
        cp $T.kin0 {output.kin}
        """


rule r11l_genocheck_rna_report:
    """
    Per-sample QC plus the RNA-internal kinship report.

      site_stats    call rate, mean depth, heterozygosity per RNA sample. Read
                    this FIRST -- a sample with a low call rate has too few
                    shared sites for any kinship estimate to mean much, and a
                    low het rate marks the dropout regime where even a correct
                    match looks only partly related.
      kinship_rna_report  every RNA pair, labelled, with the expected
                    same-donor status if donor_map.tsv is present. Unexpected
                    duplicate-level pairs (two "different" donors that are
                    actually one person) and missing ones (two regions of one
                    donor that are not related) both surface here.
    """
    input:
        vcf = GENOCHECK_OUT + "/rnaseq_genotypes.hg38.vcf.gz",
        kin = GENOCHECK_OUT + "/00_kinship_rna.kin0",
    output:
        site_stats = GENOCHECK_OUT + "/00_genocheck_site_stats.tsv",
        report     = GENOCHECK_OUT + "/00_kinship_rna_report.tsv",
    params:
        donor_map = GENOCHECK_DONOR_MAP,
        dup_cut   = GENOCHECK_DUP_KINSHIP,
    run:
        import gzip
        import pandas as pd

        # ---------- per-sample stats, parsed straight from the VCF
        with gzip.open(input.vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    vcf_samples = line.rstrip("\n").split("\t")[9:]
                    break
            n_sites = 0
            called = {s: 0 for s in vcf_samples}
            hets   = {s: 0 for s in vcf_samples}
            dp_sum = {s: 0.0 for s in vcf_samples}
            dp_n   = {s: 0 for s in vcf_samples}
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                n_sites += 1
                fmt = f[8].split(":")
                dp_i = fmt.index("DP") if "DP" in fmt else None
                for s, col in zip(vcf_samples, f[9:]):
                    parts = col.split(":")
                    gt = parts[0]
                    if gt in ("./.", ".", ".|."):
                        continue
                    called[s] += 1
                    if gt in ("0/1", "1/0", "0|1", "1|0"):
                        hets[s] += 1
                    if dp_i is not None and len(parts) > dp_i:
                        v = parts[dp_i]
                        if v not in (".", ""):
                            dp_sum[s] += float(v); dp_n[s] += 1

        stats = pd.DataFrame({
            "sample": vcf_samples,
            "n_sites": n_sites,
            "n_called": [called[s] for s in vcf_samples],
            "call_rate": [called[s] / n_sites if n_sites else 0 for s in vcf_samples],
            "mean_dp": [dp_sum[s] / dp_n[s] if dp_n[s] else 0 for s in vcf_samples],
            "n_het": [hets[s] for s in vcf_samples],
            "het_rate": [hets[s] / called[s] if called[s] else 0 for s in vcf_samples],
        })
        stats.to_csv(output.site_stats, sep="\t", index=False, float_format="%.4f")

        # ---------- RNA-internal kinship
        kin = _genocheck_read_kin0(input.kin)
        kcol = kin.attrs["kinship"]
        kin["label"] = [
            "same_individual" if k >= params.dup_cut else
            "related" if k >= 0.08 else "unrelated" for k in kin[kcol]]

        donors = _genocheck_read_map(params.donor_map)
        if donors:
            kin["same_donor_expected"] = [
                (donors.get(a) is not None and donors.get(a) == donors.get(b))
                for a, b in zip(kin["IID1"], kin["IID2"])]
            kin["agrees_with_donor_map"] = [
                (exp == (lab == "same_individual"))
                for exp, lab in zip(kin["same_donor_expected"], kin["label"])]

        kin.sort_values(kcol, ascending=False).to_csv(
            output.report, sep="\t", index=False, float_format="%.5f")

        print(f"[genocheck][report] {len(vcf_samples)} RNA samples, {n_sites} "
              f"sites, median call rate {stats['call_rate'].median():.3f}",
              file=sys.stderr)
        if donors and "agrees_with_donor_map" in kin:
            bad = int((~kin["agrees_with_donor_map"]).sum())
            print(f"[genocheck][report] {bad} RNA pair(s) disagree with "
                  f"donor_map.tsv -- see {output.report}", file=sys.stderr)


# -----------------------------------------------------------------------------
# IDENTITY tier: requires the complete array fileset
# -----------------------------------------------------------------------------

rule r11m_genocheck_chip_subset:
    """Array fileset restricted to the variants the RNA calls actually cover."""
    input:
        chip_bed = GENOCHECK_CHIP_PREFIX + ".bed",
        bim      = GENOCHECK_OUT + "/rnaseq_genotypes.bim",
    output:
        bed = temp(GENOCHECK_OUT + "/_merge/chip_shared.bed"),
        bim = temp(GENOCHECK_OUT + "/_merge/chip_shared.bim"),
        fam = temp(GENOCHECK_OUT + "/_merge/chip_shared.fam"),
    params:
        chip_prefix = GENOCHECK_CHIP_PREFIX,
        out_prefix  = GENOCHECK_OUT + "/_merge/chip_shared",
        wd          = GENOCHECK_OUT + "/_merge",
    threads: 4
    singularity: PLINK2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        cut -f2 {input.bim} > {params.wd}/shared_ids.txt
        plink2 --bfile {params.chip_prefix} \
            --extract {params.wd}/shared_ids.txt \
            --snps-only just-acgt --max-alleles 2 \
            --allow-extra-chr --sort-vars \
            --threads {threads} --make-bed --out {params.out_prefix}
        """


rule r11n_genocheck_merge:
    """
    Merge RNA + array with the standard plink1.9 two-pass flip procedure:
    attempt the merge; if PLINK reports allele-mismatching variants in .missnp,
    flip the RNA fileset at exactly those variants and retry; if any still
    mismatch, exclude them. Flips surviving r11e/r11i are rare (A/T and C/G were
    excluded, so a flip is resolvable), but a manifest can be internally
    inconsistent, and an unflipped variant reads as a genotype swap at that site.

    RNA sample IDs are suffixed '_RNA' first. This is what makes the comparison
    name-agnostic: the two sides stay distinguishable even when they use the
    same sample names, and no pairing is implied by a shared name.
    """
    input:
        rna_bed  = GENOCHECK_OUT + "/rnaseq_genotypes.bed",
        chip_bed = GENOCHECK_OUT + "/_merge/chip_shared.bed",
    output:
        bed = GENOCHECK_OUT + "/_merge/merged.bed",
        bim = GENOCHECK_OUT + "/_merge/merged.bim",
        fam = GENOCHECK_OUT + "/_merge/merged.fam",
        log = GENOCHECK_OUT + "/merge.log",
    params:
        rna  = GENOCHECK_OUT + "/rnaseq_genotypes",
        chip = GENOCHECK_OUT + "/_merge/chip_shared",
        wd   = GENOCHECK_OUT + "/_merge",
    threads: 4
    singularity: PLINK1_CONTAINER
    shell:
        r"""
        set -euo pipefail
        cd {params.wd}

        awk '{{print $1, $2, $1"_RNA", $2"_RNA"}}' {params.rna}.fam > rna_tag.txt
        plink --bfile {params.rna} --update-ids rna_tag.txt \
              --allow-extra-chr --make-bed --out rna_tagged >/dev/null

        if plink --bfile rna_tagged --bmerge {params.chip} \
                 --allow-extra-chr --make-bed --out merged >merge_pass1.log 2>&1
        then
            echo "merge: clean on first pass" > {output.log}
        else
            if [ -s merged-merge.missnp ]; then
                n=$(wc -l < merged-merge.missnp)
                echo "merge: $n allele-mismatching variants; flipping RNA and retrying" \
                    > {output.log}
                plink --bfile rna_tagged --flip merged-merge.missnp \
                      --allow-extra-chr --make-bed --out rna_flipped >/dev/null
                rm -f merged-merge.missnp
                if plink --bfile rna_flipped --bmerge {params.chip} \
                         --allow-extra-chr --make-bed --out merged \
                         >merge_pass2.log 2>&1
                then
                    echo "merge: clean after flip" >> {output.log}
                else
                    n2=$( [ -s merged-merge.missnp ] && wc -l < merged-merge.missnp || echo 0 )
                    echo "merge: $n2 variants still mismatched after flip; excluding" \
                        >> {output.log}
                    plink --bfile rna_flipped --exclude merged-merge.missnp \
                          --allow-extra-chr --make-bed --out rna_excl >/dev/null
                    plink --bfile rna_excl --bmerge {params.chip} \
                          --allow-extra-chr --make-bed --out merged >/dev/null
                fi
            else
                echo "merge: FAILED with no .missnp -- see merge_pass1.log" \
                    >> {output.log}
                cat merge_pass1.log >&2
                exit 1
            fi
        fi
        awk 'END{{print "variants_in_merged", NR}}' merged.bim >> {output.log}
        awk 'END{{print "samples_in_merged", NR}}' merged.fam >> {output.log}
        """


rule r11o_genocheck_kinship_merged:
    """KING-robust kinship across the merged RNA + array set, all pairs."""
    input:
        bed = GENOCHECK_OUT + "/_merge/merged.bed",
    output:
        kin = GENOCHECK_OUT + "/00_kinship_merged.kin0",
    params:
        prefix   = GENOCHECK_OUT + "/_merge/merged",
        wd       = GENOCHECK_OUT + "/_kinship",
        ld_prune = GENOCHECK_LD_PRUNE,
    threads: 4
    singularity: PLINK2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        T={params.wd}/merged
        if [ "{params.ld_prune}" = "True" ]; then
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --indep-pairwise 200kb 1 0.2 --threads {threads} \
                --out $T.prune >/dev/null
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --extract $T.prune.prune.in \
                --make-king-table --king-table-filter -1 \
                --threads {threads} --out $T >/dev/null
        else
            plink2 --bfile {params.prefix} --allow-extra-chr \
                --make-king-table --king-table-filter -1 \
                --threads {threads} --out $T >/dev/null
        fi
        cp $T.kin0 {output.kin}
        """


rule r11p_genocheck_ibs:
    """
    All-pairs IBS with plink1.9 `--genome full`, which reports IBS0/IBS1/IBS2
    counts. Genotype CONCORDANCE is then IBS2 / (IBS0+IBS1+IBS2) -- the fraction
    of co-called sites where the two samples carry the identical genotype. That
    is the most directly interpretable similarity measure for a sample-swap
    check, and unlike PI_HAT it does not depend on allele-frequency estimates,
    which are unreliable in a set this small.
    """
    input:
        bed = GENOCHECK_OUT + "/_merge/merged.bed",
    output:
        genome = GENOCHECK_OUT + "/_merge/merged.genome",
    params:
        prefix = GENOCHECK_OUT + "/_merge/merged",
    threads: 4
    singularity: PLINK1_CONTAINER
    shell:
        r"""
        set -euo pipefail
        plink --bfile {params.prefix} --genome full --allow-extra-chr \
              --threads {threads} --out {params.prefix} >/dev/null
        test -s {output.genome}
        """


rule r11q_genocheck_identity:
    """
    The similarity report. Joins concordance (from IBS counts) with KING kinship
    and IBS0 rate for every RNA x array pair, then condenses to one row per RNA
    sample: its most similar array sample, the margin to the runner-up, and a
    verdict.

    The MARGIN is the decisive column. A correct assignment sits in a wide gap
    from the next-best candidate; a swap, a mislabelled tube or a contaminated
    sample shows two comparably good partners. Absolute concordance is read
    loosely, because allele dropout depresses it for a correct pairing (see the
    interpretation note above).

    Verdicts: `match` (concordance >= GENOCHECK_MATCH_CONC, IBS0 rate <=
    GENOCHECK_MAX_IBS0, margin >= GENOCHECK_MIN_MARGIN); `ambiguous` (margin
    below threshold -- two candidates fit comparably well); `no_match`
    (nothing in the array set reaches duplicate-level similarity, i.e. the
    donor was never genotyped, or this RNA sample is something else entirely);
    `weak` (passes none of the above cleanly -- inspect by hand, usually a
    low-call-rate sample).

    If sample_map.tsv was supplied, the expected array ID is carried alongside
    and `agrees_with_map` says whether the data agree with it. The map never
    influences the matching itself.
    """
    input:
        genome     = GENOCHECK_OUT + "/_merge/merged.genome",
        kin        = GENOCHECK_OUT + "/00_kinship_merged.kin0",
        site_stats = GENOCHECK_OUT + "/00_genocheck_site_stats.tsv",
    output:
        pairs    = GENOCHECK_OUT + "/00_similarity_pairs.tsv",
        identity = GENOCHECK_OUT + "/00_similarity_report.tsv",
    params:
        sample_map  = GENOCHECK_SAMPLE_MAP,
        match_conc  = GENOCHECK_MATCH_CONC,
        max_ibs0    = GENOCHECK_MAX_IBS0,
        min_margin  = GENOCHECK_MIN_MARGIN,
        dup_cut     = GENOCHECK_DUP_KINSHIP,
    run:
        import pandas as pd

        def origin(iid):
            return "RNA" if str(iid).endswith("_RNA") else "array"

        def strip_tag(iid):
            s = str(iid)
            return s[:-4] if s.endswith("_RNA") else s

        # ---------- IBS -> concordance
        g = pd.read_csv(input.genome, sep=r"\s+")
        need = {"IID1", "IID2", "IBS0", "IBS1", "IBS2"}
        missing = need - set(g.columns)
        if missing:
            raise ValueError(
                f"{input.genome} lacks {sorted(missing)}; `--genome full` is "
                f"required for IBS counts. Columns found: {list(g.columns)}")
        tot = g["IBS0"] + g["IBS1"] + g["IBS2"]
        g["n_compared"] = tot
        g["concordance"] = g["IBS2"] / tot.replace(0, pd.NA)
        g["ibs0_rate_genome"] = g["IBS0"] / tot.replace(0, pd.NA)
        gcols = ["IID1", "IID2", "n_compared", "concordance",
                 "ibs0_rate_genome"]
        if "PI_HAT" in g.columns:
            gcols.append("PI_HAT")
        if "DST" in g.columns:
            gcols.append("DST")
        g = g[gcols]

        # ---------- KING kinship
        k = _genocheck_read_kin0(input.kin)
        kcol = k.attrs["kinship"]
        keep = ["IID1", "IID2", kcol] + \
            (["ibs0_rate"] if "ibs0_rate" in k.columns else [])
        k = k[keep].rename(columns={kcol: "kinship",
                                    "ibs0_rate": "ibs0_rate_king"})

        # pair keys are unordered, so join on a sorted tuple
        for df in (g, k):
            df["_key"] = [tuple(sorted((str(a), str(b))))
                          for a, b in zip(df["IID1"], df["IID2"])]
        m = g.merge(k.drop(columns=["IID1", "IID2"]), on="_key", how="outer")

        # rebuild the pair labels from the key (outer join can leave NaNs)
        m["IID1"] = [k_[0] for k_ in m["_key"]]
        m["IID2"] = [k_[1] for k_ in m["_key"]]
        m = m.drop(columns=["_key"])
        m["origin1"] = m["IID1"].map(origin)
        m["origin2"] = m["IID2"].map(origin)
        m["comparison"] = [
            "RNA-RNA" if a == b == "RNA" else
            "array-array" if a == b == "array" else "RNA-array"
            for a, b in zip(m["origin1"], m["origin2"])]
        m.sort_values("concordance", ascending=False).to_csv(
            output.pairs, sep="\t", index=False, float_format="%.5f")

        # ---------- one row per RNA sample, against ARRAY samples only
        cross = m[m["comparison"] == "RNA-array"].copy()
        cross["rna"] = [a if o1 == "RNA" else b for a, b, o1
                        in zip(cross["IID1"], cross["IID2"], cross["origin1"])]
        cross["arr"] = [b if o1 == "RNA" else a for a, b, o1
                        in zip(cross["IID1"], cross["IID2"], cross["origin1"])]

        smap = _genocheck_read_map(params.sample_map)
        stats = pd.read_csv(input.site_stats, sep="\t").set_index("sample")

        rows = []
        for rid, sub in cross.groupby("rna", sort=True):
            sub = sub.sort_values("concordance", ascending=False)
            top = sub.iloc[0]
            second = sub.iloc[1] if len(sub) > 1 else None
            conc = top["concordance"]
            margin = (conc - second["concordance"]) if second is not None \
                else float("nan")
            ibs0 = top.get("ibs0_rate_king")
            if ibs0 is None or pd.isna(ibs0):
                ibs0 = top.get("ibs0_rate_genome")

            if pd.notna(conc) and conc >= params.match_conc \
                    and (pd.isna(ibs0) or ibs0 <= params.max_ibs0) \
                    and (pd.isna(margin) or margin >= params.min_margin):
                verdict = "match"
            elif pd.notna(margin) and margin < params.min_margin:
                verdict = "ambiguous"
            elif pd.notna(top.get("kinship", float("nan"))) \
                    and top["kinship"] < params.dup_cut:
                verdict = "no_match"
            else:
                verdict = "weak"

            plain = strip_tag(rid)
            expected = smap.get(plain)
            rows.append({
                "rna_sample": plain,
                "best_array_match": top["arr"],
                "concordance": conc,
                "kinship": top.get("kinship", float("nan")),
                "ibs0_rate": ibs0,
                "n_compared": top["n_compared"],
                "runner_up": second["arr"] if second is not None else "NA",
                "runner_up_concordance": second["concordance"]
                if second is not None else float("nan"),
                "margin": margin,
                "verdict": verdict,
                "rna_call_rate": float(stats["call_rate"].get(plain, float("nan"))),
                "expected_array_id": expected if expected else "NA",
                "agrees_with_map": (top["arr"] == expected) if expected else "NA",
            })

        rep = pd.DataFrame(rows).sort_values("rna_sample")
        rep.to_csv(output.identity, sep="\t", index=False, float_format="%.5f")

        counts = rep["verdict"].value_counts().to_dict()
        print(f"[genocheck][identity] verdicts: {counts}", file=sys.stderr)
        if smap:
            bad = rep[rep["agrees_with_map"] == False]
            if len(bad):
                print(f"[genocheck][identity] {len(bad)} sample(s) do NOT match "
                      f"sample_map.tsv: "
                      f"{', '.join(bad['rna_sample'].astype(str))}",
                      file=sys.stderr)
            else:
                print("[genocheck][identity] all samples agree with "
                      "sample_map.tsv", file=sys.stderr)


# =============================================================================
# Reference download instructions (run ONCE before the pipeline)
# =============================================================================
# 1. GENOME FASTA (GTEx-specific: no ALT, no HLA, no decoy contigs)
#    This file is on Google Cloud Storage:
#      gsutil cp gs://gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
#                /tmp/data/00_additional_files/gtex_v11_refs/
#    (gsutil is from the gcloud SDK)
#
#    The same FASTA is served by the Broad's TOPMed public bucket:
#      wget https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
#
# 2. GENCODE v47 GTF
#      cd /tmp/data/00_additional_files/gtex_v11_refs/
#      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
#      gunzip gencode.v47.annotation.gtf.gz
#      mv gencode.v47.annotation.gtf gencode.v47.GRCh38.annotation.gtf
#
# 3. (Optional) Pre-built STAR index from GTEx
#    Building from scratch takes ~30 min on 16 cores. If you want exact
#    parity with GTEx, you can also download their pre-built index:
#      gsutil cp -r gs://gtex-resources/references/star_index_oh75 \
#                /tmp/data/00_additional_files/gtex_v11_refs/
#    BUT note that this index is for sjdbOverhang=75 only, so it only matches
#    if reads are 2x76bp. For other read lengths, build your own.
#
# 4. (Optional, for DROP-OUTRIDER reference panel) GTEx V11 gene counts
#    Used downstream to merge with your samples for OUTRIDER:
#      Download from https://gtexportal.org/home/downloads/adult-gtex/
#      bulk_tissue_expression -- pick the
#      "GTEx_Analysis_v11_RNASeQC_2.4.3_RNAseQC_*_reads.gct.gz" file.
#      Filter to the brain-cortex-related tissues (Brain - Frontal Cortex BA9,
#      Brain - Anterior cingulate cortex BA24, etc.) for OUTRIDER reference.
#
# =============================================================================
# Sample table format:
# =============================================================================
# Tab-separated, columns: name path dataset1 dataset2 extension include_in_analysis
#   name                = unique sample id
#   path                = subdirectory under /tmp/data/ where FASTQs live
#   dataset1            = R1 file basename (without extension)
#   dataset2            = R2 file basename (or "NaN" for SE)
#   extension           = "fq.gz" or "fastq.gz"
#   include_in_analysis = 1 (process) or 0 (skip)
#   gender              = "male"/"female"/"unknown" (for sex inference)
# =============================================================================
