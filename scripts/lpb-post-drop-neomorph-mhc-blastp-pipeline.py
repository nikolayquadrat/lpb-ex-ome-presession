#!/usr/bin/env python3
"""
IZUMO4 neomorph -> neoantigen -> microbial/viral-mimicry pipeline (MOCK / scaffold).
Steps:
  1. input peptide            -- the novel neomorph C-terminus
  2. tile                     -- class I: 8-11-mers; class II: handled by NetMHCIIpan
  3. set alleles              -- SZ07's class I + class II (from arcasHLA)
  4. get strong binders       -- class I (MHCflurry) AND class II (NetMHCIIpan)
  5. BLAST them               -- vs BACTERIAL and VIRAL proteins
  6. find 100% hits           -- full-length exact matches
Install:
    pip install mhcflurry biopython
    mhcflurry-downloads fetch                  # class I models (one-time)
    # class II: NetMHCIIpan is a separate academic download (DTU Health Tech);
    # install it and put `netMHCIIpan` on PATH. (pip alternative: mhcnuggets,
    # less standard.) Class I could also use netMHCpan; MHCflurry is used here
    # only because it is pip-installable.
"""

import io
import os
import re
import sys
import json
import math
import time
import random
import tempfile
import datetime
import subprocess
import argparse
from collections import defaultdict


# ---------------------------------------------------------------------------
# Output plumbing: mirror everything printed to a run log in --out-dir, so every
# mode leaves an on-disk record without touching the ~50 print() call sites.
# ---------------------------------------------------------------------------
class _Tee:
    """Duplicate a stream to a file. Used to capture stdout into the out-dir."""
    def __init__(self, stream, fh):
        self._stream = stream
        self._fh = fh

    def write(self, data):
        self._stream.write(data)
        self._fh.write(data)
        return len(data)

    def flush(self):
        self._stream.flush()
        self._fh.flush()


def _out_path(name):
    """Resolve a filename inside the active output dir (set by main); when no
    out-dir is configured, fall back to the current directory."""
    return os.path.join(OUTPUT_DIR or ".", name)


def _dump_parameters(args, mode):
    """Print every run parameter (CLI args + resolved config constants) so the
    log is a self-contained record. Goes to stdout, which is tee'd to the log."""
    line = "=" * 72
    print(line)
    print("RUN PARAMETERS")
    print(line)
    print(f"  timestamp        : {datetime.datetime.now().isoformat(timespec='seconds')}")
    print(f"  mode             : {mode}")
    print(f"  output_dir       : {OUTPUT_DIR}")
    print("  --- CLI arguments ---")
    for k in sorted(vars(args)):
        print(f"    {k:18}= {getattr(args, k)!r}")
    print("  --- resolved config (effective values this run) ---")
    cfg = [
        ("NEOMORPH_PEPTIDE", NEOMORPH_PEPTIDE),
        ("ARCAS_CLASS_I", ARCAS_CLASS_I),
        ("ARCAS_CLASS_II", ARCAS_CLASS_II),
        ("KMER_LENGTHS_I", KMER_LENGTHS_I),
        ("CLASS_II_WINDOW", CLASS_II_WINDOW),
        ("STRONG_RANK_I", STRONG_RANK_I),
        ("STRONG_RANK_II", STRONG_RANK_II),
        ("BLAST_DB", BLAST_DB),
        ("BLAST_TAXA", BLAST_TAXA),
        ("BLAST_LOCAL_DBS", BLAST_LOCAL_DBS),
        ("BLAST_MATRIX", BLAST_MATRIX),
        ("BLAST_GAPCOSTS", BLAST_GAPCOSTS),
        ("BLAST_EVALUE", BLAST_EVALUE),
        ("BLAST_WORD_SIZE", BLAST_WORD_SIZE),
        ("BLAST_THRESHOLD", BLAST_THRESHOLD),
        ("BLAST_PAUSE_S", BLAST_PAUSE_S),
        ("MAX_MISMATCH", MAX_MISMATCH),
        ("MIN_POSITIVE_FRAC", MIN_POSITIVE_FRAC),
        ("ANCHOR_POSITIONS", ANCHOR_POSITIONS),
        ("REQUIRE_ANCHOR_CONSERVED", REQUIRE_ANCHOR_CONSERVED),
        ("ALLOW_GAPS", ALLOW_GAPS),
        ("PRIORITY_ORGANISMS", PRIORITY_ORGANISMS),
        ("CLADE_TERMS", CLADE_TERMS),
        ("SPECIES_FILE", SPECIES_FILE),
        ("EMAIL", EMAIL),
        ("NETMHCIIPAN_BIN", NETMHCIIPAN_BIN),
    ]
    for k, v in cfg:
        s = repr(v)
        if len(s) > 300:
            s = s[:300] + " ...(truncated)"
        print(f"    {k:24}= {s}")
    print(line)

# ---------------------------------------------------------------------------
# CONFIG  -- edit these
# ---------------------------------------------------------------------------

# (1) neomorph peptide (novel C-terminus).
NEOMORPH_PEPTIDE = "PAKISECRSPAQSRLPRQRDPGAGRSPGFSSCPPPGRSWTKWRQQCTR" # IXUMO4's PAKI + 44aa neomorph sequence

# (3a) TYPED class I alleles. Use his real HLA typing, not these defaults.
#     MHCflurry allele format: "HLA-A*02:01".
ARCAS_CLASS_I = [
    "HLA-A*11:01",
    "HLA-A*31:01",
    "HLA-B*35:01",
    "HLA-B*39:01",
    "HLA-C*04:01",
    "HLA-C*12:03",
]

# (3b) CLASS II alleles from arcasHLA. Give each gene's diploid calls;
#  the converter builds NetMHCIIpan names (DR = beta only; DQ/DP = all a x b pairs).
# DOESN'T WORK ON WINDOWS 
ARCAS_CLASS_II = {
    "DRB1": ["DRB1*01:01"],
    "DQA1": ["DQA1*01:01"],
    "DQB1": ["DQB1*05:01"], # second allele DQB1*05:328 is approximated by DQB1*05:01
    "DPA1": ["DPA1*01:03"],
    "DPB1": ["DPB1*02:01", "DPB1*04:01"],
}

# (2) tiling / windows
KMER_LENGTHS_I = (8, 9, 10, 11)      # class I peptide lengths
CLASS_II_WINDOW = 15                 # class II peptide length NetMHCIIpan slides

# (4) strong-binder thresholds
STRONG_RANK_I = 0.5                  # class I: %Rank (<=0.5 matches online NetMHCpan)
STRONG_RANK_II = 2.0                 # class II: %Rank_EL (<=2 strong, NetMHCIIpan conv.)

# (5/6) BLAST: search BOTH bacteria and viruses; each hit is tagged by taxon.
BLAST_DB = "nr"
# Genomes: no taxon restriction (entrez_query omitted). Remote mode submits
# every binder in ONE batched job -- API-safe. Set entrez to a string here only
# if you want to re-narrow (e.g. "txid5653[ORGN]" for Trypanosomatidae).
# BLAST_TAXA = {"all": None}
BLAST_TAXA = {"all": "bacteria[organism] OR viruses[organism]"}
BLAST_LOCAL_DBS = {"all": "nr"}      # local: name of your combined protein DB
BLAST_MATRIX = "PAM30"               # looser than PAM30 for short peptides (was PAM30)
BLAST_GAPCOSTS = "9 1"               # matches PAM70 (PAM30 wants "9 1")
BLAST_EVALUE = 200000
BLAST_WORD_SIZE = 2
BLAST_THRESHOLD = 16                 # neighborhood word threshold; lower = more sensitive (does not working properly?)
BLAST_PAUSE_S = 12                   # >=10 s: NCBI throttles anonymous heavy use

# --- near-match / conservative-substitution acceptance (loosened from exact-only) ---
# A hit is kept if it spans the whole query with:
#   identities >= qlen - MAX_MISMATCH        (a small budget of non-identical positions)
#   positives  >= ceil(MIN_POSITIVE_FRAC*qlen)  (conservative subs, e.g. K<->R, I<->L, count)
# Optionally require the MHC ANCHOR positions to stay identical/conservative, since a
# substitution there more likely kills binding than one at a TCR-contact position.
# Optional pre-specified pathogen/contamination universe. When --add-species is
# passed, the eukaryotes in this file (+ bacteria/viruses) are OR-COMBINED with
# BLAST_TAXA above (addition, not substitution). Default path = rnaseq-pipe layout.
SPECIES_FILE = r"data\rnaseq-pipe\00_additional_files\contamination\species.txt"
# Higher-taxon (clade) terms covering the contamination list's eukaryotes, for
# --add-clades. Far fewer OR-terms than the ~50-species list, so NCBI's BLAST
# entrez filter accepts it. Preferred for the shuffle-BLAST (combinatorial) null.
CLADE_TERMS = [
    "fungi[organism]",            # Aspergillus, Candida, Cryptococcus, Mucor, Fusarium, Malassezia, ...
    "Apicomplexa[organism]",      # Toxoplasma, Plasmodium, Babesia
    "Kinetoplastida[organism]",   # Trypanosoma, Leishmania
    "Heterolobosea[organism]",    # Naegleria
    "Amoebozoa[organism]",        # Acanthamoeba, Balamuthia
    "Platyhelminthes[organism]",  # Schistosoma, Taenia, Echinococcus
    "Nematoda[organism]",         # Toxocara, Angiostrongylus, Baylisascaris, Gnathostoma
    "Ciliophora[organism]",       # Tetrahymena, Paramecium (decoys)
]
MAX_MISMATCH = 2                     # e.g. a 9-mer -> accept down to 7/9 identity
MIN_POSITIVE_FRAC = 1.0              # every position must be identical OR a conservative sub
ANCHOR_POSITIONS = (2, -1)           # 1-based; -1 = C-terminus. Typical class I anchors (P2, PΩ).
REQUIRE_ANCHOR_CONSERVED = True      # anchors must be identity/positive (set False to disable)
ALLOW_GAPS = False                   # gapped near-matches are non-physiological for a short core
# NCBI contact e-mail for remote BLAST. Intentionally empty here -- pass it at
# run time with --ncbi-email (NCBI wants a real address on programmatic jobs and
# uses it to warn rather than block). Remote BLAST refuses to run without it.
EMAIL = ""
# Output directory for run logs + result files. Set from --out-dir at startup;
# None means "current directory" (see _out_path).
OUTPUT_DIR = None
# organisms to FLAG in the output because of independent evidence for them
# (a pre-specified prior beats anything the broad scan turns up by chance).
PRIORITY_ORGANISMS = ["Homo"]
NETMHCIIPAN_BIN = "netMHCIIpan" # path to the executable
# Human proteome amino-acid frequencies (UniProt/Swiss-Prot, approx %). Used by the
# natural-aa null. Weights need not sum to 100 (random.choices normalizes).
AA_FREQ = {
    'A': 8.25, 'R': 5.53, 'N': 4.06, 'D': 5.46, 'C': 1.38, 'Q': 3.93, 'E': 6.72,
    'G': 7.08, 'H': 2.27, 'I': 5.91, 'L': 9.65, 'K': 5.80, 'M': 2.41, 'F': 3.86,
    'P': 4.74, 'S': 6.65, 'T': 5.36, 'W': 1.10, 'Y': 2.92, 'V': 6.87,
}


# ---------------------------------------------------------------------------
# ALLELE CONVERTERS  (arcasHLA -> tool formats)
# ---------------------------------------------------------------------------
def _two_field_digits(allele):
    """'DRB1*03:01:01' -> '0301' ; 'DQA1*05:01' -> '0501'."""
    after = allele.split("*", 1)[1]
    return "".join(after.split(":")[:2])


def arcas_to_mhcflurry_I(alleles):
    """'A*02:01:01' -> 'HLA-A*02:01' (MHCflurry class I format)."""
    out = []
    for a in alleles:
        a = a.replace("HLA-", "")
        gene, after = a.split("*", 1)
        out.append(f"HLA-{gene}*" + ":".join(after.split(":")[:2]))
    return out


def arcas_to_netmhc2(class_ii):
    """Build NetMHCIIpan allele strings from an arcasHLA class II dict:
        DR  -> 'DRB1_0301'  (beta chain only)
        DQ  -> 'HLA-DQA10501-DQB10201'  (every alpha x beta pairing)
        DP  -> 'HLA-DPA10103-DPB10401'
    """
    out = []
    for gene, alleles in class_ii.items():
        if gene.startswith("DRB"):
            for a in alleles:
                out.append(f"{gene}_{_two_field_digits(a)}")

    def pairs(a_gene, b_gene, prefix):
        for a in class_ii.get(a_gene, []):
            for b in class_ii.get(b_gene, []):
                out.append(f"HLA-{prefix}A1{_two_field_digits(a)}-{prefix}B1{_two_field_digits(b)}")

    pairs("DQA1", "DQB1", "DQ")
    pairs("DPA1", "DPB1", "DP")
    return out


# ---------------------------------------------------------------------------
# progress helper: tqdm if available, else a lightweight inline counter.
# Never a hard dependency (falls back silently if tqdm isn't installed).
# ---------------------------------------------------------------------------
def _progress(iterable, total=None, desc=""):
    try:
        from tqdm import tqdm
        return tqdm(iterable, total=total, desc=desc, unit="it", leave=True)
    except ImportError:
        pass

    def _gen():
        import sys
        n = str(total) if total is not None else "?"
        i = 0
        for item in iterable:
            i += 1
            sys.stderr.write(f"\r  {desc}: {i}/{n}   ")
            sys.stderr.flush()
            yield item
        sys.stderr.write("\n")
    return _gen()


# ---------------------------------------------------------------------------
# (2) TILING (class I)
# ---------------------------------------------------------------------------
def tile(seq, lengths=KMER_LENGTHS_I):
    out = {}
    for k in lengths:
        for i in range(len(seq) - k + 1):
            out.setdefault(seq[i:i + k], i)
    return out                                   # {kmer: offset}


# ---------------------------------------------------------------------------
# (4a) CLASS I  -> strong binders  (MHCflurry)
# ---------------------------------------------------------------------------
def predict_class1(kmers, alleles, strong_rank=STRONG_RANK_I, show_progress=True):
    """Score every k-mer against every allele with MHCflurry; keep rows whose
    percentile rank clears the strong-binder threshold. Uses percentile (not raw
    affinity) to match NetMHCpan's %Rank_EL criterion used online."""
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    peptides = list(kmers)
    rows = []
    alist = (_progress(alleles, total=len(alleles), desc="class I MHC")
             if show_progress else alleles)
    for allele in alist:
        df = predictor.predict(peptides=peptides, alleles=[allele],
                               include_affinity_percentile=True, verbose=0)
        # MHCflurry's EL-like rank column; fall back to affinity_percentile
        rankcol = ("presentation_percentile" if "presentation_percentile" in df.columns
                   else "affinity_percentile")
        for _, r in df.iterrows():
            rank = float(r[rankcol])
            if rank <= strong_rank:
                rows.append({"class": "I", "allele": allele,
                             "peptide": r["peptide"],        # the k-mer
                             "blast_seq": r["peptide"],       # blast the k-mer itself
                             "score": f"rank {rank:.3f}% / {float(r['affinity']):.0f}nM",
                             "offset": kmers[r["peptide"]]})
    return rows


# ---------------------------------------------------------------------------
# (4b) CLASS II -> strong binders  (NetMHCIIpan via subprocess)
# ---------------------------------------------------------------------------
def _run_netmhc2(peptide, allele, window=CLASS_II_WINDOW, binpath=NETMHCIIPAN_BIN):
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        fh.write(f">pep\n{peptide}\n")
        fpath = fh.name
    proc = subprocess.run(
        [binpath, "-a", allele, "-f", fpath, "-length", str(window), "-BA"],
        capture_output=True, text=True, check=True)
    return proc.stdout


def _parse_netmhc2(stdout):
    """Defensive parser: keys on header names, which vary slightly by NetMHCIIpan
    version. Returns list of dicts with peptide, core, rank_el, affinity."""
    def find(header, *cands):
        for c in cands:
            if c in header:
                return header.index(c)
        for i, h in enumerate(header):
            if any(c.lower() in h.lower() for c in cands):
                return i
        return None

    header = None
    idx = {}
    rows = []
    for line in stdout.splitlines():
        s = line.strip()
        if not s or set(s) <= set("-="):
            continue
        if s.startswith("Pos") and "Peptide" in s:
            header = s.split()
            idx = {"pep": find(header, "Peptide"),
                   "core": find(header, "Core"),
                   "rank": find(header, "%Rank_EL", "Rank_EL", "%Rank"),
                   "aff": find(header, "Affinity(nM)", "Affinity")}
            continue
        if header and s[0].isdigit():
            f = s.split()
            try:
                rows.append({
                    "peptide": f[idx["pep"]],
                    "core": f[idx["core"]] if idx["core"] is not None else f[idx["pep"]],
                    "rank_el": float(f[idx["rank"]]),
                    "affinity": float(f[idx["aff"]]) if idx["aff"] is not None else float("nan"),
                })
            except (IndexError, ValueError):
                continue
    return rows


def predict_class2(peptide, alleles, strong_rank=STRONG_RANK_II, show_progress=True):
    rows = []
    alist = (_progress(alleles, total=len(alleles), desc="class II MHC")
             if show_progress else alleles)
    for allele in alist:
        try:
            out = _run_netmhc2(peptide, allele)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"      [warn] NetMHCIIpan failed for {allele}: {e}")
            continue
        for rec in _parse_netmhc2(out):
            if rec["rank_el"] <= strong_rank:
                rows.append({"class": "II", "allele": allele,
                             "peptide": rec["peptide"],       # the 15-mer
                             "blast_seq": rec["core"],         # BLAST the 9-mer core
                             "score": f"rank {rec['rank_el']}%",
                             "offset": None})
    return rows


# ---------------------------------------------------------------------------
# (5) BLAST  (remote NCBI, per taxon)  +  (6) exact-hit filter
# ---------------------------------------------------------------------------
def blast_remote(peptide, entrez, db=BLAST_DB):
    """Single-peptide remote BLAST (kept for ad-hoc use). entrez=None -> no
    taxon filter (all genomes)."""
    from Bio.Blast import NCBIWWW, NCBIXML
    kw = {} if entrez is None else {"entrez_query": entrez}
    handle = NCBIWWW.qblast(
        "blastp",
        db,
        peptide,
        matrix_name=BLAST_MATRIX,
        expect=BLAST_EVALUE,
        word_size=BLAST_WORD_SIZE,
        hitlist_size=250,
        gapcosts=BLAST_GAPCOSTS,
        composition_based_statistics="0",
        threshold=BLAST_THRESHOLD,
        short_query=True, # not sure if it does anything
        **kw)
    return NCBIXML.read(handle)


def parse_organism(hit_def):
    """Extract the source organism from a BLAST hit definition -- nr/RefSeq defs
    carry it in the LAST [square brackets], e.g. '... [Toxoplasma gondii]'."""
    m = re.findall(r"\[([^\[\]]+)\]", hit_def or "")
    return m[-1] if m else "unknown"


# --- optional pre-specified pathogen/contamination universe (--add-species) ---
_EUK_KEYWORDS = ("fung", "yeast", "protozoa", "protist", "helminth",
                 "amoeba", "mould", "mold")


def _parse_species_eukaryotes(path):
    """Read the contamination species list; return the EUKARYOTE organism names
    (fungi/protozoa/helminths/protists). Bacteria & viruses are omitted here
    because they're covered by the broad bacteria[organism]/viruses[organism]
    terms and their short keys don't map to clean names. Classification is by
    section-header keyword, so sentence-like '##' comments don't misfire."""
    names, is_euk = {}, False
    with open(path) as fh:
        for raw in fh:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                if raw.lstrip().startswith("##"):
                    hdr = s.lstrip("#").strip().lower()
                    is_euk = any(k in hdr for k in _EUK_KEYWORDS)
                continue
            if is_euk:
                parts = s.split()
                if len(parts) < 2:
                    continue
                toks = parts[0].split("_")
                nm = (toks[0].capitalize() + " " +
                      (toks[1].lower() if len(toks) > 1 else "")).strip()
                names[nm] = True
    return sorted(names)


def _species_universe(path, include_broad=False):
    """Entrez terms for the listed EUKARYOTES. By default returns ONLY the
    eukaryote clauses (so --add-species is a pure addition to whatever base
    already has bacteria/viruses). Set include_broad=True to also prepend the
    broad bacteria/viruses terms (use when the base is None/unrestricted and you
    want a self-contained universe)."""
    names = _parse_species_eukaryotes(path)
    terms = ([f'"{n}"[organism]' for n in names])
    if include_broad:
        terms = ["bacteria[organism]", "viruses[organism]"] + terms
    return " OR ".join(terms), names


def _species_terms(path):
    """List of per-eukaryote '"Name"[organism]' entrez terms (+ the names)."""
    names = _parse_species_eukaryotes(path)
    return [f'"{n}"[organism]' for n in names], names


def _combine_entrez(base, extra):
    """OR-combine two entrez expressions (either may be None). Preserves base --
    this is what makes --add-species an ADDITION, not a substitution."""
    parts = []
    if base and str(base).strip():
        parts.append(f"({base})")
    if extra and str(extra).strip():
        parts.append(f"({extra})")
    return " OR ".join(parts) if parts else None


def blast_test(peptide, organism, use_local_blast=False):
    """DIAGNOSTIC: BLAST one peptide against a SINGLE organism and print the raw
    outcome. Restricting to one organism shrinks the search to seconds, so this
    isolates 'does BLAST work / is my e-mail+network+params ok' from the slow
    all-genomes search. Prints timing, the exact entrez query, hit count, and the
    first few alignments (or the raw XML head on a parse surprise)."""
    print(f"[test] BLASTing {peptide} vs organism='{organism}'  "
          f"({'local' if use_local_blast else 'remote'})")
    entrez = f'"{organism}"[organism]'
    print(f"[test] entrez_query = {entrez}")
    t0 = time.time()
    try:
        hits = _blast_group([peptide], use_local_blast=use_local_blast,
                            entrez=entrez)
    except SystemExit as ex:
        print(f"[test] BLAST refused: {ex}")
        return
    except Exception as ex:
        print(f"[test] BLAST raised {type(ex).__name__}: {ex}")
        print("[test] -> if this is a timeout, the CODE/CONNECTION is fine but the "
              "public queue is busy; retry, run off-peak, or use --local-blast.")
        return
    dt = time.time() - t0
    total, exact, taxa = hits.get(peptide, (0, 0, set()))
    print(f"[test] returned in {dt:.1f}s: {total} full-length hit(s) "
          f"({exact} exact) across {len(taxa)} taxon label(s)")
    if taxa:
        for org in sorted(taxa):
            print(f"[test]    hit organism: {org}")
    if total == 0:
        print("[test] zero hits is a VALID result (BLAST worked; this peptide simply "
              "has no full-length match in that organism at these settings). The "
              "round-trip succeeded, so the pipeline's BLAST plumbing is working.")
    print(f"[test] SUCCESS: BLAST round-trip completed against '{organism}'.")


def blast_batch_remote(sequences, classmap, entrez=None, db=BLAST_DB):
    """Submit ALL binder sequences in ONE multi-FASTA qblast job (API-safe: one
    job, not N). Records return in submission order; we zip them to the input
    and sanity-check via query length. Returns exact-hit dicts tagged with
    class + organism.
    """
    from Bio.Blast import NCBIWWW, NCBIXML
    if not EMAIL or not str(EMAIL).strip():
        raise SystemExit(
            "[error] remote BLAST needs an NCBI contact e-mail. Pass "
            "--ncbi-email you@example.org (or use --local-blast for a local DB).")
    NCBIWWW.email = EMAIL                      # identify yourself to NCBI
    NCBIWWW.tool = "izumo4_neomorph_pipeline"

    fasta = "".join(f">q{i}\n{s}\n" for i, s in enumerate(sequences))
    kw = {} if entrez is None else {"entrez_query": entrez}
    handle = NCBIWWW.qblast(
        "blastp",
        db,
        fasta,
        matrix_name=BLAST_MATRIX,
        expect=BLAST_EVALUE,
        word_size=BLAST_WORD_SIZE,
        hitlist_size=250,
        gapcosts=BLAST_GAPCOSTS,
        threshold=BLAST_THRESHOLD,
        composition_based_statistics="0",
        short_query=True,
        **kw)
    time.sleep(BLAST_PAUSE_S) # polite settle before parsing

    records = list(NCBIXML.parse(handle))
    if len(records) != len(sequences):
        print(f"      [warn] {len(records)} records for {len(sequences)} queries "
              f"-- order-based mapping may be off")
    out = []
    for seq, rec in zip(sequences, records):
        ql = getattr(rec, "query_letters", len(seq))

        if ql not in (len(seq), 0):
            print(f"      [warn] length mismatch: {seq} (record={ql})")

        # ---------------------------------------------------------
        # DEBUG: what did NCBI BLAST actually return?
        # ---------------------------------------------------------
        raw_hsps = [
            hsp
            for aln in rec.alignments
            for hsp in aln.hsps
        ]

        print(
            f"[BLAST DEBUG] {seq}: "
            f"{len(rec.alignments)} raw alignment(s), "
            f"{len(raw_hsps)} raw HSP(s)"
        )

        # Show a few raw BLAST alignments BEFORE similar_hits()
        for aln in rec.alignments[:10]:
            for hsp in aln.hsps[:3]:
                print(
                    f"    RAW: {aln.accession} "
                    f"align_len={hsp.align_length} "
                    f"id={hsp.identities} "
                    f"pos={getattr(hsp, 'positives', '?')} "
                    f"gaps={getattr(hsp, 'gaps', '?')} "
                    f"E={hsp.expect}"
                )

        # ---------------------------------------------------------
        # Your existing downstream filter
        # ---------------------------------------------------------
        accepted = similar_hits(rec, seq)

        print(
            f"[BLAST DEBUG] {seq}: "
            f"{len(accepted)} hit(s) survive similar_hits()"
        )

        for h in accepted:
            h.update({
                "query": seq,
                "class": classmap.get(seq, "?"),
                "organism": parse_organism(h["subject"]),
            })
            out.append(h)
    return out


def blast_local(peptide, db_path):
    from Bio.Blast import NCBIXML
    proc = subprocess.run(
        ["blastp", "-db", db_path, "-task", "blastp-short",
         "-matrix", BLAST_MATRIX, "-evalue", str(BLAST_EVALUE),
         "-comp_based_stats", "0", "-outfmt", "5", "-max_target_seqs", "50"],
        input=f">q\n{peptide}\n", capture_output=True, text=True, check=True)
    return NCBIXML.read(io.StringIO(proc.stdout))


def exact_hits(record, query):
    """Full-length 100% matches only (identities == len AND alignment spans it)."""
    qlen = len(query)
    hits = []
    for aln in record.alignments:
        for hsp in aln.hsps:
            if hsp.identities == qlen and hsp.align_length == qlen:
                hits.append({"accession": aln.accession,
                             "subject": aln.hit_def,     # full, for organism parse
                             "evalue": hsp.expect})
    return hits


def _anchor_ok(hsp, qlen):
    """Check the MHC anchor positions are identity or conservative (positive) in the
    alignment. Uses hsp.match: ' '=mismatch, '+'=conservative, letter=identity.

    If the match line is missing or not the query length (which happens for some
    BLAST XML, especially remote multi-query jobs), we CANNOT verify the anchors
    from it -- so fall back to the numeric evidence: a hit that already passed the
    identity/positives budget in similar_hits is accepted here rather than
    silently dropped. Set REQUIRE_ANCHOR_CONSERVED=False to skip anchor checks
    entirely."""
    if not REQUIRE_ANCHOR_CONSERVED:
        return True
    match = getattr(hsp, "match", "") or ""
    if len(match) != qlen:
        # Can't map anchors onto positions without a full-length match line.
        # Don't reject on that alone -- the caller already enforced the
        # identity/positives thresholds; accept and let those govern.
        return True
    for a in ANCHOR_POSITIONS:
        idx = (qlen + a) if a < 0 else (a - 1)   # -1 -> last; 1-based -> 0-based
        if not (0 <= idx < qlen):
            continue
        if match[idx] == " ":         # a plain mismatch at an anchor kills it
            return False
    return True


def similar_hits(record, query,
                 max_mismatch=MAX_MISMATCH,
                 min_positive_frac=MIN_POSITIVE_FRAC,
                 allow_gaps=ALLOW_GAPS):
    """Near-matches: full-length, small mismatch budget, conservative substitutions
    allowed (BLAST 'positives'), optionally anchor-preserving. Each hit is tagged
    tier='exact' (identity==qlen) or tier='near'. This is the loosened successor to
    exact_hits -- it catches biochemically similar peptides, not just identical ones."""
    qlen = len(query)
    need_pos = math.ceil(min_positive_frac * qlen)
    hits = []
    for aln in record.alignments:
        for hsp in aln.hsps:
            if hsp.align_length != qlen:
                continue
            if not allow_gaps and getattr(hsp, "gaps", 0):
                continue
            ident = hsp.identities
            pos = getattr(hsp, "positives", ident)   # positives includes identities
            if ident < qlen - max_mismatch:
                continue
            if pos < need_pos:
                continue
            if not _anchor_ok(hsp, qlen):
                continue
            hits.append({"accession": aln.accession,
                         "subject": aln.hit_def,
                         "evalue": hsp.expect,
                         "identity": f"{ident}/{qlen}",
                         "positives": f"{pos}/{qlen}",
                         "tier": "exact" if ident == qlen else "near"})
    return hits


# ---------------------------------------------------------------------------
# ORCHESTRATION
# ---------------------------------------------------------------------------
def run(peptide=NEOMORPH_PEPTIDE, use_local_blast=False):
    print(f"[1] peptide ({len(peptide)} aa): {peptide}")

    alleles_I = arcas_to_mhcflurry_I(ARCAS_CLASS_I)
    alleles_II = arcas_to_netmhc2(ARCAS_CLASS_II)
    print(f"[3] class I alleles : {', '.join(alleles_I)}")
    print(f"    class II alleles: {', '.join(alleles_II)}")

    kmers = tile(peptide)
    print(f"[2] class I tiled -> {len(kmers)} unique k-mers; "
          f"class II -> NetMHCIIpan slides {CLASS_II_WINDOW}-mers internally")

    binders = []
    binders += predict_class1(kmers, alleles_I)
    binders += predict_class2(peptide, alleles_II)
    n1 = sum(b["class"] == "I" for b in binders)
    n2 = sum(b["class"] == "II" for b in binders)
    print(f"[4] strong binders: {n1} class I  +  {n2} class II")
    for b in binders:
        pos = f"@{b['offset']}" if b["offset"] is not None else ""
        print(f"      [{b['class']:<2}] {b['blast_seq']:<12} {b['allele']:<22} "
              f"{b['score']:<10} {pos}")

    # unique sequences to BLAST (class I k-mers + class II cores), remember class
    to_blast = {}
    for b in binders:
        to_blast.setdefault(b["blast_seq"], b["class"])
    seqs = sorted(to_blast)
    mode = "local, per-seq" if use_local_blast else "remote, ONE batched submission"
    print(f"[5] BLASTing {len(seqs)} unique binder sequences vs ALL genomes ({mode}) ...")

    all_exact = []
    if use_local_blast:
        db = BLAST_LOCAL_DBS["all"]
        for seq in seqs:
            try:
                rec = blast_local(seq, db)
            except Exception as e:
                print(f"      [warn] BLAST failed ({seq}): {e}")
                continue
            for h in similar_hits(rec, seq):
                h.update({"query": seq, "class": to_blast[seq],
                          "organism": parse_organism(h["subject"])})
                all_exact.append(h)
    else:
        try:
            all_exact = blast_batch_remote(seqs, to_blast, entrez=BLAST_TAXA["all"])
        except Exception as e:
            print(f"      [warn] batched BLAST failed: {e}")

    # ---- organism-stratified report ----------------------------------------
    # All genomes includes HOST: a SELF (Homo sapiens) hit means the k-mer is NOT
    # unique to the neomorph -> not a mimic (QC flag). Priority organisms (independent
    # prior) escape the look-elsewhere penalty. Hits are now EXACT (100%) or NEAR
    # (conservative substitutions within the mismatch budget, anchors preserved).
    n_exact = sum(h.get("tier") == "exact" for h in all_exact)
    n_near = sum(h.get("tier") == "near" for h in all_exact)
    print(f"[6] full-length hits: {len(all_exact)}  ({n_exact} exact, {n_near} near / conservative)")
    by_org = defaultdict(list)
    for h in all_exact:
        by_org[h["organism"]].append(h)
    for org in sorted(by_org, key=lambda o: -len(by_org[o])):
        hits = by_org[org]
        if org == "Homo sapiens":
            flag = "  <- SELF (also in human proteome; NOT a mimic)"
        elif any(p.lower() in org.lower() for p in PRIORITY_ORGANISMS):
            flag = "  <<< PRIORITY (independent prior -- weigh above chance hits)"
        else:
            flag = ""
        ne = sum(h.get("tier") == "exact" for h in hits)
        print(f"      {org:<30} {len(hits):>3} hit(s) [{ne} exact]{flag}")
        for h in hits:
            tag = h.get("tier", "?")
            idp = f"{h.get('identity','?')} id, {h.get('positives','?')} pos"
            print(f"          [{h['class']:<2}][{tag:<5}] {h['query']:<12} == "
                  f"{h['accession']:<12} ({idp})  {h['subject'][:40]}")

    print("\n>>> Search is now LOOSENED (PAM70, conservative subs, near-matches). This "
          "RAISES the false-positive rate further -- exact and near counts are both "
          "expected background. Run --null on the SAME settings, and for any PRIORITY "
          "or low-%Rank_Pathogen hit, check anchor + TCR-face conservation on SZ07's "
          "actual allele before concluding.")

    # ---- write structured results into the output dir ----------------------
    binders_tsv = _out_path("binders.tsv")
    with open(binders_tsv, "w") as fh:
        fh.write("class\tblast_seq\tallele\tscore\toffset\n")
        for b in binders:
            fh.write(f"{b['class']}\t{b['blast_seq']}\t{b['allele']}\t{b['score']}\t"
                     f"{'' if b['offset'] is None else b['offset']}\n")
    hits_tsv = _out_path("blast_hits.tsv")
    with open(hits_tsv, "w") as fh:
        fh.write("organism\tclass\ttier\tquery\taccession\tidentity\tpositives\t"
                 "subject\tis_self\tis_priority\n")
        for h in all_exact:
            org = h["organism"]
            is_self = "1" if org == "Homo sapiens" else "0"
            is_prio = "1" if any(p.lower() in org.lower() for p in PRIORITY_ORGANISMS) else "0"
            fh.write(f"{org}\t{h['class']}\t{h.get('tier','?')}\t{h['query']}\t"
                     f"{h.get('accession','')}\t{h.get('identity','')}\t"
                     f"{h.get('positives','')}\t{h.get('subject','')}\t{is_self}\t{is_prio}\n")
    print(f"[out] wrote {binders_tsv} ({len(binders)} binders) and "
          f"{hits_tsv} ({len(all_exact)} hits)")
    return binders, all_exact


# ---------------------------------------------------------------------------
# OPTIONAL BUT CRUCIAL: null model (class I binder count vs shuffled background)
# ---------------------------------------------------------------------------
def null_model(peptide=NEOMORPH_PEPTIDE, n_shuffles=200, seed=0):
    rng = random.Random(seed)
    chars = list(peptide)
    alleles_I = arcas_to_mhcflurry_I(ARCAS_CLASS_I)

    real = len(predict_class1(tile(peptide), alleles_I, show_progress=False))
    null = []
    for _ in _progress(range(n_shuffles), total=n_shuffles, desc="binding null"):
        rng.shuffle(chars)
        null.append(len(predict_class1(tile("".join(chars)), alleles_I,
                                       show_progress=False)))
    null.sort()
    import statistics as st
    ge = sum(c >= real for c in null)
    print(f"[null] real class I strong binders: {real}")
    print(f"[null] shuffled mean={st.mean(null):.2f} median={st.median(null)} "
          f"p95={null[int(0.95*len(null))-1]}")
    print(f"[null] p(count >= real) = {ge}/{n_shuffles} = {ge/n_shuffles:.3f}  "
          f"(small = binder-dense; ~0.5 = unremarkable)")
    # (extendable to class II by swapping in predict_class2)


# ---------------------------------------------------------------------------
# BLAST NULL: does a composition-matched shuffle pull the SAME hit cloud?
# ---------------------------------------------------------------------------
# This is the null that answers the MIMICRY question (vs null_model, which only
# tests binder density). It shuffles ONE binder peptide, BLASTs each shuffle vs
# all genomes with the SAME loosened settings, and asks whether the real peptide
# gets more hits -- and from different organisms -- than composition-matched
# scrambles. If the shuffles hit the same GC-rich taxa at the same rate, the
# observed hits are compositional background, not a signal about the sequence.
def _blast_one(seq, use_local_blast=False):
    """BLAST a single peptide vs all genomes; return (total, exact, taxa set)."""
    if use_local_blast:
        rec = blast_local(seq, BLAST_LOCAL_DBS["all"])
        hits = []
        for h in similar_hits(rec, seq):
            h["organism"] = parse_organism(h["subject"])
            hits.append(h)
    else:
        hits = blast_batch_remote([seq], {seq: "I"}, entrez=BLAST_TAXA["all"])
    total = len(hits)
    exact = sum(h.get("tier") == "exact" for h in hits)
    taxa = {h["organism"] for h in hits}
    return total, exact, taxa


def _blast_group(seqs, use_local_blast=False, entrez=None):
    """BLAST a LIST of peptides and return {seq: (total, exact, taxa_set)}.
    Remote path submits the whole list as ONE batched qblast job (not one per
    seq) -- this is what keeps NCBI from throttling. Local loops per seq."""
    out = {s: (0, 0, set()) for s in seqs}
    if use_local_blast:
        for s in seqs:
            t, e, taxa = _blast_one(s, True)
            out[s] = (t, e, taxa)
        return out
    hits = blast_batch_remote(list(seqs), {s: "I" for s in seqs}, entrez=entrez)
    agg = {s: [0, 0, set()] for s in seqs}
    for h in hits:
        q = h["query"]
        agg.setdefault(q, [0, 0, set()])
        agg[q][0] += 1
        if h.get("tier") == "exact":
            agg[q][1] += 1
        agg[q][2].add(h["organism"])
    for s in seqs:
        out[s] = (agg[s][0], agg[s][1], agg[s][2])
    return out


def null_model_blast(peptide, n_shuffles=20, seed=0, use_local_blast=False,
                     batch_size=50):
    """Composition-matched BLAST null for ONE binder peptide (e.g. RQRDPGAGR).
    Shuffle preserves amino-acid composition and only scrambles order, isolating
    'is it THIS sequence' from 'is it the residue composition'.

    The real peptide + all shuffles are BLASTed in ONE batched submission (or a
    few chunks of batch_size), not one request each -- so remote NCBI sees ~1-2
    jobs instead of n_shuffles+1, which avoids the throttling that makes repeated
    single submissions crawl. Use --local-blast for an unthrottled, larger null."""
    import statistics as st
    from collections import Counter
    rng = random.Random(seed)
    chars = list(peptide)

    print(f"[blast-null] testing {peptide}  (n_shuffles={n_shuffles}, "
          f"{'local' if use_local_blast else 'remote, batched'} BLAST)")

    # generate unique composition-matched shuffles up front
    shuffles, seen, attempts = [], set(), 0
    while len(shuffles) < n_shuffles and attempts < n_shuffles * 50:
        rng.shuffle(chars)
        s = "".join(chars)
        attempts += 1
        if s == peptide or s in seen:
            continue
        seen.add(s)
        shuffles.append(s)
    if len(shuffles) < n_shuffles:
        print(f"[blast-null] only {len(shuffles)} distinct shuffles exist for this "
              f"composition; using those")

    # BLAST real + shuffles together, in batches (few requests total)
    all_seqs = [peptide] + shuffles
    results = {}
    batches = [all_seqs[i:i + batch_size] for i in range(0, len(all_seqs), batch_size)]
    for bi, batch in enumerate(_progress(batches, total=len(batches),
                                         desc="blast null batches"), 1):
        try:
            results.update(_blast_group(batch, use_local_blast,
                                        entrez=BLAST_TAXA["all"]))
        except Exception as ex:
            print(f"   [warn] batch {bi}/{len(batches)} failed: {ex}")
    if peptide not in results:
        print("[blast-null] real peptide BLAST failed -- cannot compute null")
        return

    real_total, real_exact, real_taxa = results[peptide]
    print(f"[blast-null] REAL {peptide}: {real_total} hits ({real_exact} exact), "
          f"{len(real_taxa)} distinct taxa")

    null_total, null_exact = [], []
    taxon_counter = Counter()
    for s in shuffles:
        if s not in results:
            continue
        t, e, taxa = results[s]
        null_total.append(t)
        null_exact.append(e)
        taxon_counter.update(taxa)
        print(f"   {s}: {t} hits ({e} exact)")

    if not null_total:
        print("[blast-null] no shuffles completed -- cannot compute null")
        return

    def prank(real, null):
        ge = sum(x >= real for x in null)
        return ge, len(null), ge / len(null)

    gt, nt, pt = prank(real_total, null_total)
    ge, ne, pe = prank(real_exact, null_exact)
    print("\n[blast-null] ---- RESULT ----")
    print(f"[blast-null] total hits: real={real_total}  null mean={st.mean(null_total):.1f}"
          f"  p(null>=real)={gt}/{nt}={pt:.2f}")
    print(f"[blast-null] exact hits: real={real_exact}  null mean={st.mean(null_exact):.1f}"
          f"  p(null>=real)={ge}/{ne}={pe:.2f}")
    print(f"[blast-null] interpretation: p~0.5 => real sits in the middle of the "
          f"shuffled background => hits are COMPOSITIONAL, not sequence-specific.")
    print(f"[blast-null] most common taxa across shuffles (same GC-rich orgs as the "
          f"real run? => compositional attractor):")
    for org, c in taxon_counter.most_common(12):
        print(f"   {c:3d}/{len(null_total)}  {org}")


# ---------------------------------------------------------------------------
# NATURAL-AA NULL: strong-binder count for random peptides drawn from natural
# amino-acid frequencies, matched to the neomorph region length. Unlike the
# shuffle null (which fixes THIS peptide's composition), this varies composition,
# so it answers "is a strong binder here surprising vs a typical peptide at all".
# It is CONSERVATIVE: real frameshift/intron-retention neomorphs are MORE Arg/Pro/
# Gly-skewed than the proteome average (and those residues favour promiscuous class
# I binders), so a proteome-background null gives a LOWER binder mean than a true-
# neomorph null -- i.e. it can only understate how unremarkable the peptide is.
# ---------------------------------------------------------------------------
def _sample_peptide(length, rng):
    aas, wts = zip(*AA_FREQ.items())
    return "".join(rng.choices(aas, weights=wts, k=length))


def _blast_counts(seqs, use_local_blast=False):
    """Return {seq: (total_hits, exact_hits)} for a list of peptides, using the same
    loosened similar_hits filter. Remote path batches all seqs into ONE qblast job."""
    counts = {}
    if not seqs:
        return counts
    if use_local_blast:
        for s in _progress(seqs, total=len(seqs), desc="blast binders"):
            try:
                rec = blast_local(s, BLAST_LOCAL_DBS["all"])
                hs = similar_hits(rec, s)
            except Exception as ex:
                print(f"   [warn] blast failed {s}: {ex}")
                hs = []
            counts[s] = (len(hs), sum(h.get("tier") == "exact" for h in hs))
    else:
        agg = {s: [0, 0] for s in seqs}
        try:
            hits = blast_batch_remote(list(seqs), {s: "I" for s in seqs},
                                      entrez=BLAST_TAXA["all"])
        except Exception as ex:
            print(f"   [warn] batched blast failed: {ex}")
            hits = []
        for h in hits:
            q = h["query"]
            agg.setdefault(q, [0, 0])
            agg[q][0] += 1
            if h.get("tier") == "exact":
                agg[q][1] += 1
        counts = {s: (v[0], v[1]) for s, v in agg.items()}
    return counts


def null_model_natural(peptide, length=None, n_draws=500, seed=0,
                       with_blastp=False, use_local_blast=False, max_blast=50):
    """Rank the real peptide's strong-binder count against random peptides sampled
    from natural aa frequencies at the same length. Uses the SAME predict_class1
    (percentile threshold), so null and real stay on identical criteria.

    with_blastp=True extends the null through the BLAST step: it BLASTs the strong
    binders produced by the random peptides and compares their exact-hit counts to
    the real peptide's strong binders -- i.e. 'do random-composition strong binders
    get BLAST hits like the real one did?'. Heavy (many BLAST queries); random
    binders are de-duped and capped at max_blast."""
    import statistics as st
    rng = random.Random(seed)
    alleles_I = arcas_to_mhcflurry_I(ARCAS_CLASS_I)
    L = length or len(peptide)

    real_rows = predict_class1(tile(peptide), alleles_I, show_progress=False)
    real = len(real_rows)
    null = []
    rand_binder_seqs = []
    for _ in _progress(range(n_draws), total=n_draws, desc="natural-aa null"):
        pep = _sample_peptide(L, rng)
        rows = predict_class1(tile(pep), alleles_I, show_progress=False)
        null.append(len(rows))
        if with_blastp:
            rand_binder_seqs.extend(r["blast_seq"] for r in rows)
    null.sort()
    ge = sum(c >= real for c in null)
    print(f"[nat-null] region length={L}, draws={n_draws}, alleles={len(alleles_I)}")
    print(f"[nat-null] real strong binders: {real}")
    print(f"[nat-null] null mean={st.mean(null):.2f} median={st.median(null)} "
          f"p95={null[int(0.95*len(null))-1]}")
    print(f"[nat-null] p(null >= real) = {ge}/{n_draws} = {ge/n_draws:.3f}")
    print(f"[nat-null] ~0.5 => unremarkable vs a typical peptide. NB this is "
          f"CONSERVATIVE: real neomorphs are more binder-dense, so the true surprise "
          f"is even lower than this p suggests.")

    if not with_blastp:
        return

    # ---- BLAST step: do random-composition strong binders get hits like the real? ----
    real_bind = sorted({r["blast_seq"] for r in real_rows})
    rand_bind = sorted(set(rand_binder_seqs))
    if len(rand_bind) > max_blast:
        rng.shuffle(rand_bind)
        rand_bind = sorted(rand_bind[:max_blast])
    print(f"\n[nat-null+blastp] BLASTing {len(real_bind)} real + {len(rand_bind)} random "
          f"strong binders vs all genomes ({'local' if use_local_blast else 'remote'}) ...")
    real_counts = _blast_counts(real_bind, use_local_blast)
    rand_counts = _blast_counts(rand_bind, use_local_blast)

    n_ex = [v[1] for v in rand_counts.values()]
    if not n_ex:
        print("[nat-null+blastp] no random strong binders to BLAST (raise --n-draws)")
        return
    print(f"[nat-null+blastp] exact hits per RANDOM strong binder: "
          f"mean={st.mean(n_ex):.1f} median={st.median(n_ex)} "
          f"p95={sorted(n_ex)[int(0.95*len(n_ex))-1]}  (n={len(n_ex)})")
    for s, (tot, ex) in sorted(real_counts.items()):
        ge2 = sum(x >= ex for x in n_ex)
        print(f"   real binder {s}: {ex} exact / {tot} total -> "
              f"p(random >= {ex}) = {ge2}/{len(n_ex)} = {ge2/len(n_ex):.2f}")
    print("[nat-null+blastp] real binders' hit counts inside the random distribution "
          "=> BLAST hits are chance for strong binders of this length too.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=("IZUMO4 neomorph neoantigen/mimicry pipeline. Predicts MHC class I "
                     "(MHCflurry) + class II (NetMHCIIpan) binding for the neomorph, "
                     "BLASTs strong binders vs all genomes, and provides three null "
                     "models to test whether binders / BLAST hits exceed chance."),
        epilog=("examples:\n"
                "  full pipeline (bind + BLAST):        python %(prog)s\n"
                "  binding null, shuffle this peptide:  python %(prog)s --null\n"
                "  binding null, natural aa freqs:      python %(prog)s --null-natural-aa\n"
                "     ...also BLAST the random binders: python %(prog)s --null-natural-aa --blastp\n"
                "  BLAST mimicry null for one k-mer:    python %(prog)s --null-with-blastp --peptide RQRDPGAGR\n"
                "  use local BLAST DBs:                 add --local-blast to any BLAST mode\n"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptide", default=NEOMORPH_PEPTIDE,
                    help="peptide to analyse. DEFAULT: the novel C-terminus saved in the "
                         "script (NEOMORPH_PEPTIDE), so binding modes need no --peptide. "
                         "For --null-with-blastp, pass a single binder k-mer (e.g. RQRDPGAGR).")
    ap.add_argument("--local-blast", action="store_true",
                    help="use local BLAST+ against BLAST_LOCAL_DBS instead of remote NCBI")
    ap.add_argument("--null", action="store_true",
                    help="binding null: shuffle THIS peptide (fixes composition; tests order)")
    ap.add_argument("--null-natural-aa", action="store_true",
                    help="binding null: random peptides from natural aa freqs (varies "
                         "composition; the correct 'is a strong binder surprising at all' test). "
                         "Uses the saved neomorph peptide by default.")
    ap.add_argument("--blastp", action="store_true",
                    help="with --null-natural-aa: also BLAST the random-peptide strong binders "
                         "and compare their exact-hit counts to the real peptide's binders")
    ap.add_argument("--null-with-blastp", action="store_true",
                    help="mimicry null: shuffle ONE binder k-mer (--peptide) and BLAST each "
                         "shuffle vs all genomes; do the same GC-rich taxa recur?")
    ap.add_argument("--neomorph-length", type=int, default=None,
                    help="length to sample for --null-natural-aa (default: len of --peptide)")
    ap.add_argument("--n-draws", type=int, default=500,
                    help="random draws for --null-natural-aa (default 500)")
    ap.add_argument("--max-blast", type=int, default=50,
                    help="cap on random strong binders BLASTed under --null-natural-aa --blastp")
    ap.add_argument("--n-shuffles", type=int, default=20,
                    help="shuffles for --null-with-blastp (batched, so remote-safe)")
    ap.add_argument("--batch-size", type=int, default=50,
                    help="peptides per batched qblast submission (--null-with-blastp)")
    ap.add_argument("--add-clades", action="store_true",
                    help="ADD higher-taxon clade terms covering the contamination list "
                         "(fungi + Apicomplexa + Kinetoplastida + helminth/amoeba/ciliate "
                         "phyla) to BLAST_TAXA via OR. NCBI-safe (few terms); preferred for "
                         "the shuffle-BLAST null. Use INSTEAD of --add-species when the full "
                         "species list is rejected as too complex.")
    ap.add_argument("--add-species", action="store_true",
                    help="ADD the per-species contamination universe (bacteria + viruses + "
                         "~48 listed eukaryotes) to BLAST_TAXA via OR. Precise but the long "
                         "OR-list may exceed NCBI's entrez limit; prefer --add-clades if so.")
    ap.add_argument("--species-file", default=SPECIES_FILE,
                    help=f"path to the contamination species list (default: {SPECIES_FILE})")
    ap.add_argument("--out-dir", "-o", default=".",
                    help="directory for the run log and result files (created if "
                         "missing). Every mode writes a timestamped run_<mode>_*.log "
                         "here; the full pipeline also writes binders.tsv + "
                         "blast_hits.tsv. Default: current directory.")
    ap.add_argument("--ncbi-email", default=EMAIL,
                    help="e-mail address sent to NCBI with remote BLAST jobs (NCBI "
                         "uses it to warn rather than block on heavy use). Default: "
                         "the EMAIL set in the script.")
    ap.add_argument("--test-organism", default=None, metavar="ORGANISM",
                    help="DIAGNOSTIC: BLAST a single peptide against ONLY this "
                         "organism (e.g. --test-organism 'Escherichia coli'), print "
                         "the raw result, and exit. Restricting to one organism makes "
                         "the search return in seconds, so this is a quick end-to-end "
                         "check that remote (or --local-blast) BLAST works at all, "
                         "separate from the slow all-genomes default.")
    ap.add_argument("--test-peptide", default="RQRDPGAGR", metavar="SEQ",
                    help="peptide used by --test-organism (default: RQRDPGAGR).")
    a = ap.parse_args()

    # --- output dir + email (module globals every call site reads at call time) ---
    OUTPUT_DIR = os.path.abspath(a.out_dir)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    EMAIL = a.ncbi_email
    SPECIES_FILE = a.species_file

    _mode = ("test_organism" if a.test_organism else
             "null_with_blastp" if a.null_with_blastp else
             "null_natural_aa" if a.null_natural_aa else
             "null" if a.null else "full")
    _stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    _log_path = os.path.join(OUTPUT_DIR, f"run_{_mode}_{_stamp}.log")
    _log_fh = open(_log_path, "w")
    sys.stdout = _Tee(sys.stdout, _log_fh)
    sys.stderr = _Tee(sys.stderr, _log_fh)
    print(f"[out] output dir: {OUTPUT_DIR}")
    print(f"[out] run log:    {_log_path}")

    # optional: OR the pre-specified pathogen universe onto BLAST_TAXA (keeps the
    # body value; every remote call site reads BLAST_TAXA['all'] at call time).
    # --add-clades (few broad taxon terms, NCBI-safe) and/or --add-species (exact
    # ~48-species list). Broad bacteria/viruses are prepended only if the body
    # base is empty, so nothing is duplicated.
    if a.add_clades or a.add_species:
        base = BLAST_TAXA["all"]
        terms, desc = [], []
        if a.add_clades:
            terms += CLADE_TERMS
            desc.append(f"{len(CLADE_TERMS)} clades")
        if a.add_species:
            try:
                sp, snames = _species_terms(a.species_file)
                terms += sp
                desc.append(f"{len(snames)} species")
            except FileNotFoundError:
                print(f"[warn] --add-species: file not found ({a.species_file}); "
                      f"skipping species terms")
        if terms:
            if not (base and str(base).strip()):
                terms = ["bacteria[organism]", "viruses[organism]"] + terms
            BLAST_TAXA["all"] = _combine_entrez(base, " OR ".join(terms))
            eff = BLAST_TAXA["all"] or "ALL (unrestricted)"
            print(f"[universe] added {' + '.join(desc)} to BLAST_TAXA")
            print(f"[universe] effective entrez: {eff[:120]}"
                  f"{'...' if len(eff) > 120 else ''}")

    # full parameter dump (after the universe is resolved, so BLAST_TAXA shows
    # its EFFECTIVE value) -> terminal + log.
    _dump_parameters(a, _mode)

    try:
        if a.test_organism:
            blast_test(a.test_peptide, a.test_organism,
                       use_local_blast=a.local_blast)
        elif a.null_with_blastp:
            # here --peptide should be a SINGLE binder k-mer, e.g. RQRDPGAGR
            null_model_blast(a.peptide, n_shuffles=a.n_shuffles,
                             use_local_blast=a.local_blast, batch_size=a.batch_size)
        elif a.null_natural_aa:
            # --peptide defaults to the saved NEOMORPH region (novel C-term)
            null_model_natural(a.peptide, length=a.neomorph_length, n_draws=a.n_draws,
                               with_blastp=a.blastp, use_local_blast=a.local_blast,
                               max_blast=a.max_blast)
        elif a.null:
            null_model(a.peptide)
        else:
            run(a.peptide, use_local_blast=a.local_blast)
    finally:
        print(f"[out] done; log saved to {_log_path}")
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        _log_fh.close()