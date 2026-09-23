#!/usr/bin/env python3
"""Filter online blastp results for the IZUMO4 neomorph mimicry screen and write a
human-readable tab-delimited table.

Reproduces the pipeline's near-match rule on the EXPORTED results:
  a hit is kept iff, for a full-length (query-spanning) ungapped alignment,
    * no gaps,
    * the alignment spans the whole query (align_len == query_len, q 1..L),
    * <= MAX_MISMATCH non-identical residues (identities >= L - MAX_MISMATCH),
    * EVERY aligned column is identical or a conservative substitution
      (no plain mismatch), and
    * the canonical class I anchor positions P2 and P-Omega (C-terminus) are
      conserved (identical or conservative).

Input: the BLAST XML export (has the aligned strings + the match/midline, which
is what lets us check anchors faithfully). The CSV hit-table is optional and used
only to cross-check counts.

Usage:
    python filter_blast_hits.py RESULT.xml [--csv HITTABLE.csv] \
        [--out kept.tsv] [--max-mismatch 2] [--anchors 2,-1] \
        [--priority Homo] [--peptide RQRDPGAGR]
"""
import argparse, re, sys, xml.etree.ElementTree as ET

def parse_organism(hit_def):
    """nr/RefSeq defs carry the source organism in the LAST [brackets]."""
    m = re.findall(r"\[([^\[\]]+)\]", hit_def or "")
    return m[-1] if m else "unknown"

def anchor_index(a, L):
    """1-based anchor spec -> 0-based column index. Negative counts from the end
    (-1 = last position = P-Omega)."""
    return (L + a) if a < 0 else (a - 1)

def anchors_conserved(midline, qseq, hseq, L, anchors):
    """The class I anchor check. The BLAST MIDLINE encodes each aligned column:
         a LETTER  -> identity           (kept)
         a '+'     -> conservative sub   (kept)
         a SPACE   -> non-conservative mismatch (fails)
    An anchor is 'conserved' iff its column is identity or conservative, i.e. the
    midline character at that position is NOT a space. We also treat a real
    residue-vs-residue identity as conserved when a midline is unavailable."""
    for a in anchors:
        i = anchor_index(a, L)
        if not (0 <= i < L):
            continue
        if midline and len(midline) == L:
            if midline[i] == " ":          # plain mismatch at an anchor -> reject
                return False
        else:                              # no usable midline: fall back to identity
            if i < len(qseq) and i < len(hseq) and qseq[i] != hseq[i]:
                return False
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xml")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--out", default="kept_hits.tsv")
    ap.add_argument("--max-mismatch", type=int, default=2)
    ap.add_argument("--anchors", default="2,-1",
                    help="comma anchor positions, 1-based; -1 = C-terminus")
    ap.add_argument("--priority", default="Homo",
                    help="comma organism substrings to flag PRIORITY")
    ap.add_argument("--peptide", default=None, help="query, for the header note")
    a = ap.parse_args()
    anchors = [int(x) for x in a.anchors.split(",") if x.strip()]
    priority = [p.strip().lower() for p in a.priority.split(",") if p.strip()]

    root = ET.parse(a.xml).getroot()
    def findtext(node, tag):
        el = node.find(f".//{tag}")
        return el.text if el is not None else None
    qlen = int(findtext(root, "BlastOutput_query-len") or 0)
    qdef = findtext(root, "BlastOutput_query-def") or (a.peptide or "query")

    n_hits = n_hsps = 0
    kept = []
    for hit in root.iter("Hit"):
        n_hits += 1
        hit_def = findtext(hit, "Hit_def") or ""
        acc = findtext(hit, "Hit_accession") or ""
        org = parse_organism(hit_def)
        for hsp in hit.iter("Hsp"):
            n_hsps += 1
            ident = int(findtext(hsp, "Hsp_identity") or 0)
            pos   = int(findtext(hsp, "Hsp_positive") or 0)
            gaps  = int(findtext(hsp, "Hsp_gaps") or 0)
            alen  = int(findtext(hsp, "Hsp_align-len") or 0)
            qfrom = int(findtext(hsp, "Hsp_query-from") or 0)
            qto   = int(findtext(hsp, "Hsp_query-to") or 0)
            qseq  = findtext(hsp, "Hsp_qseq") or ""
            hseq  = findtext(hsp, "Hsp_hseq") or ""
            mid   = findtext(hsp, "Hsp_midline") or ""
            ev    = findtext(hsp, "Hsp_evalue") or ""

            # --- the filter ---
            if gaps != 0:                                   continue  # no gaps
            if alen != qlen or qfrom != 1 or qto != qlen:   continue  # full-length span
            if ident < qlen - a.max_mismatch:               continue  # <=MAX_MISMATCH diffs
            if pos != qlen:                                 continue  # every col identity/conservative
            if not anchors_conserved(mid, qseq, hseq, qlen, anchors):
                continue                                              # P2 + P-Omega conserved

            tier = "exact" if ident == qlen else "near"
            is_self = (org == "Homo sapiens")
            is_prio = any(p in org.lower() for p in priority)
            kept.append(dict(organism=org, tier=tier, accession=acc,
                             identity=f"{ident}/{qlen}", positives=f"{pos}/{qlen}",
                             qseq=qseq, hseq=hseq, midline=mid, evalue=ev,
                             hit_def=hit_def,
                             flag=("SELF" if is_self else "PRIORITY" if is_prio else "")))

    # order: exact first, then by organism
    kept.sort(key=lambda h: (h["tier"] != "exact", h["organism"].lower()))

    with open(a.out, "w") as fh:
        fh.write(f"# query\t{qdef}\tlen={qlen}\n")
        fh.write(f"# input_hits\t{n_hits}\tinput_hsps\t{n_hsps}\tkept\t{len(kept)}\n")
        fh.write(f"# rule\tno_gaps; full_length; <= {a.max_mismatch} non-identical; "
                 f"all columns identical/conservative; anchors {a.anchors} conserved\n")
        cols = ["organism","tier","flag","identity","positives","evalue",
                "accession","qseq","hseq","midline","hit_def"]
        fh.write("\t".join(cols) + "\n")
        for h in kept:
            fh.write("\t".join(str(h[c]) for c in cols) + "\n")

    # console summary
    n_exact = sum(h["tier"] == "exact" for h in kept)
    print(f"[in]  {n_hits} hits / {n_hsps} HSPs from {a.xml}")
    print(f"[kept] {len(kept)}  ({n_exact} exact, {len(kept)-n_exact} near)")
    from collections import Counter
    for org, c in Counter(h["organism"] for h in kept).most_common(12):
        print(f"   {c:4d}  {org}")
    print(f"[out] {a.out}")

if __name__ == "__main__":
    main()
