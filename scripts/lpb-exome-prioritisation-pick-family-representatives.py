#!/usr/bin/env python3
"""
pick_representatives.py - Pick one sample per family of related individuals.

Inputs
------
--kin0       KING-format kinship table (.kin0). Whitespace-separated; columns
             include sample IDs and a kinship coefficient. Empty file (no
             relatedness detected) is tolerated.
--fam        PLINK .fam file listing all samples in the cohort. Used as the
             authoritative full sample list, in case some samples don't appear
             in .kin0 (because they weren't related to anyone).
--kinship_threshold   Kinship coefficient above which two samples are
             considered "related" (default 0.0884 = KING's 2nd-degree cutoff).
--out_representatives Output: one sample-ID per line, suitable for
             `bcftools view -S`. Always sorted lexicographically.
--out_log    Plain-text log explaining what was kept, what was dropped, and
             which family each dropped sample belonged to.

Algorithm
---------
1. Read the full sample list from .fam (column 2 = within-family ID).
2. Read related pairs from .kin0 where Kinship >= threshold.
3. Build an undirected graph: nodes = samples, edges = related-pair links.
4. Connected-component analysis identifies families.
5. Within each family of size > 1, pick the first sample (lexicographically).
6. Output the union of (singletons + family representatives).

The greedy/lexicographic choice is deterministic and reproducible. If you
want a more principled choice (e.g., highest-coverage representative),
extend with a sample-quality table.
"""

import argparse
import sys
from collections import defaultdict


def parse_kin0(kin0_path, threshold):
    """
    Parse a KING-format kinship table. Returns set of related-pair tuples
    (a, b) where Kinship >= threshold.
    """
    pairs = set()
    if not kin0_path:
        return pairs

    try:
        with open(kin0_path) as f:
            header = None
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#") or header is None:
                    # Header line -- split and find the kinship column index
                    header = line.lstrip("#").split()
                    continue

                fields = line.split()
                if len(fields) < len(header):
                    continue

                # KING / plink2 .kin0 column names:
                #   FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
                # plink2 sometimes uses lowercase or slight variants.
                # Find the kinship column by matching "kinship" / "Kinship".
                kinship_idx = None
                for i, h in enumerate(header):
                    if h.lower().endswith("kinship"):
                        kinship_idx = i
                        break
                if kinship_idx is None:
                    # Fall back to last column (KING convention)
                    kinship_idx = len(header) - 1

                try:
                    kinship = float(fields[kinship_idx])
                except (ValueError, IndexError):
                    continue
                if kinship < threshold:
                    continue

                # Sample IDs: IID1 (col 1) and IID2 (col 3) in standard KING
                iid1_idx = next((i for i, h in enumerate(header)
                                 if h.upper() in ("IID1", "ID1")), 1)
                iid2_idx = next((i for i, h in enumerate(header)
                                 if h.upper() in ("IID2", "ID2")), 3)
                a, b = fields[iid1_idx], fields[iid2_idx]
                if a == b:
                    continue
                pairs.add((a, b))
    except FileNotFoundError:
        print(f"WARNING: {kin0_path} not found; treating all samples as unrelated",
              file=sys.stderr)
    return pairs


def parse_fam(fam_path):
    """Read sample IDs from PLINK .fam (column 2 = IID)."""
    samples = []
    with open(fam_path) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 2:
                samples.append(fields[1])
    return samples


def find_families(samples, related_pairs):
    """Connected-component analysis. Returns dict family_id -> set of samples."""
    parent = {s: s for s in samples}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]   # path compression
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for a, b in related_pairs:
        if a in parent and b in parent:
            union(a, b)

    families = defaultdict(set)
    for s in samples:
        families[find(s)].add(s)
    return families


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kin0", required=True)
    ap.add_argument("--fam",  required=True)
    ap.add_argument("--kinship_threshold", type=float, default=0.0884)
    ap.add_argument("--out_representatives", required=True)
    ap.add_argument("--out_log", required=True)
    args = ap.parse_args()

    samples = parse_fam(args.fam)
    if not samples:
        sys.exit(f"ERROR: no samples found in {args.fam}")

    related_pairs = parse_kin0(args.kin0, args.kinship_threshold)
    families = find_families(samples, related_pairs)

    # Pick representatives: lexicographically first sample per family.
    representatives = sorted(min(family) for family in families.values())

    # Write outputs
    with open(args.out_representatives, "w") as f:
        for s in representatives:
            f.write(s + "\n")

    with open(args.out_log, "w") as f:
        f.write(f"Cohort size:           {len(samples)}\n")
        f.write(f"Kinship threshold:     {args.kinship_threshold}\n")
        f.write(f"Related-pair count:    {len(related_pairs)}\n")
        f.write(f"Family count:          {len(families)}\n")
        f.write(f"Representative count:  {len(representatives)}\n")
        f.write("\n")
        f.write("Families with > 1 member (kept = representative, dropped = others):\n")
        for fam_id in sorted(families.keys()):
            members = sorted(families[fam_id])
            if len(members) > 1:
                kept = members[0]
                dropped = members[1:]
                f.write(f"  family of {len(members):2d}: keep {kept}; drop {', '.join(dropped)}\n")
        f.write("\n")
        f.write("Singletons (no detected relatives):\n")
        singletons = [m for fam in families.values() for m in fam if len(fam) == 1]
        for s in sorted(singletons):
            f.write(f"  {s}\n")

    print(f"Picked {len(representatives)} representatives "
          f"from {len(samples)} cohort samples "
          f"({len(samples) - len(representatives)} dropped as relatives)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
