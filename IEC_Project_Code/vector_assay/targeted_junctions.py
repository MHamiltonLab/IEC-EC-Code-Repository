#!/usr/bin/env python3
"""Classify LTR soft-clips and summarize the two junction boundaries at supplied target sites."""

import os
import re
from collections import defaultdict

from _assay import OUTPUT_DIR, sam_path, target_sites

U3_5P = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGG"
U5_3P = "AGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGT"
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]


def ident(a, b):
    """best ungapped identity of a aligned anywhere fully inside b (a shorter)."""
    if len(a) > len(b):
        a = a[: len(b)]
    best = 0
    for off in range(0, len(b) - len(a) + 1):
        m = sum(a[j] == b[off + j] for j in range(len(a)))
        best = max(best, m)
    return best / len(a) if a else 0


def classify(clip):
    if len(clip) < 10:
        return ("short", None, 0)
    best = ("none", None, 0)
    for name, ref in (("U3", U3_5P), ("U5", U5_3P)):
        for orient, q in (("+", clip), ("-", rc(clip))):
            qq = q[:30]
            sc = ident(qq, ref)
            if sc > best[2]:
                best = (name, orient, sc)
    return best if best[2] >= 0.85 else ("none", None, best[2])


_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse(cigar, seq):
    p = _CIG.findall(cigar)
    if not p:
        return 0, 0, 0
    lead = int(p[0][0]) if p[0][1] == "S" else 0
    trail = int(p[-1][0]) if p[-1][1] == "S" else 0
    ref = sum(int(n) for n, op in p if op in "MDN=X")
    return lead, trail, ref


WIN = 1500

if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SITES = target_sites()
    TABLE = []

    for lab, raw, chrom, center, gene in SITES:
        f = sam_path(raw)

        u3mol, u5mol = defaultdict(set), defaultdict(set)
        gspan = {"U3": [10**12, 0], "U5": [10**12, 0]}
        for line in open(f):
            if not line.strip() or line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 11 or c[5] == "*" or c[9] == "*":
                continue
            if c[2] != chrom or not (center - WIN <= int(c[3]) <= center + WIN):
                continue
            pos, cigar, seq = int(c[3]), c[5], c[9]
            lead, trail, ref = parse(cigar, seq)
            mb = c[0].split("_molbar")[0]
            if lead >= trail and lead > 0:
                junction, clip = pos, seq[:lead]
            elif trail > 0:
                junction, clip = pos + ref, seq[-trail:]
            else:
                continue
            name, orient, sc = classify(clip)
            if name == "U3":
                u3mol[junction].add(mb)
                gspan["U3"] = [min(gspan["U3"][0], pos), max(gspan["U3"][1], pos + ref)]
            elif name == "U5":
                u5mol[junction].add(mb)
                gspan["U5"] = [min(gspan["U5"][0], pos), max(gspan["U5"][1], pos + ref)]

        def mode(d):
            return max(d.items(), key=lambda kv: len(kv[1])) if d else (None, set())

        j3, m3 = mode(u3mol)
        j5, m5 = mode(u5mol)
        print(f"\n========== {lab}  {chrom}:~{center:,} ==========")
        if j3:
            print(
                f"  5' (U3) junction : chr{chrom[3:]}:{j3:,}   {len(m3)} molecules   "
                f"reads span {gspan['U3'][0]:,}-{gspan['U3'][1]:,} ({gspan['U3'][1] - gspan['U3'][0]} bp genomic)"
            )
        if j5:
            print(
                f"  3' (U5) junction : chr{chrom[3:]}:{j5:,}   {len(m5)} molecules   "
                f"reads span {gspan['U5'][0]:,}-{gspan['U5'][1]:,} ({gspan['U5'][1] - gspan['U5'][0]} bp genomic)"
            )
        if j3 and j5:
            d = abs(j3 - j5)
            verdict = (
                "short separation (<=20 bp)" if d <= 20 else "wide separation (>20 bp)"
            )
            print(f"  >>> 5'-3' separation = {d} bp  ->  {verdict}")
        else:
            d, verdict = "", "incomplete paired-junction evidence"
            print("  >>> Paired junction separation is unavailable")
        TABLE.append(
            {
                "sample": lab,
                "gene": gene,
                "chrom": chrom,
                "U3_5prime_junction": j3 or "",
                "U3_molecules": len(m3) if j3 else 0,
                "U5_3prime_junction": j5 or "",
                "U5_molecules": len(m5) if j5 else 0,
                "junction_separation_bp": d,
                "verdict": verdict,
            }
        )

    cols = [
        "sample",
        "gene",
        "chrom",
        "U3_5prime_junction",
        "U3_molecules",
        "U5_3prime_junction",
        "U5_molecules",
        "junction_separation_bp",
        "verdict",
    ]
    out = os.path.join(OUTPUT_DIR, "targeted_junctions.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in TABLE:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nWrote table: {out}")
