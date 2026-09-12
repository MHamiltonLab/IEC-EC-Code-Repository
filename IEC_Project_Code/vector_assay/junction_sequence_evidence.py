#!/usr/bin/env python3
"""Display sequence-level LTR evidence for a user-supplied SAM file and genomic interval."""

import argparse
import re

U3_5P = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGG"
U5_3P = "AGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGT"
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse(cig):
    p = _CIG.findall(cig)
    if not p:
        return 0, 0, 0
    lead = int(p[0][0]) if p[0][1] == "S" else 0
    trail = int(p[-1][0]) if p[-1][1] == "S" else 0
    ref = sum(int(n) for n, o in p if o in "MDN=X")
    return lead, trail, ref


def best_align(clip, ref):
    """best ungapped placement of clip[:30] in ref -> (identity, offset, aligned_ref)."""
    q = clip[:30]
    best = (0, 0, "")
    if not q:
        return 0, 0, ""
    for o in range(0, len(ref) - len(q) + 1):
        m = sum(q[j] == ref[o + j] for j in range(len(q)))
        if m > best[0]:
            best = (m, o, ref[o : o + len(q)])
    return best[0] / len(q), best[1], best[2]


def classify(clip):
    res = []
    for nm, ref in (("U3", U3_5P), ("U5", U5_3P)):
        for orient, q in (("fwd", clip), ("rev", rc(clip))):
            idy, off, al = best_align(q, ref)
            res.append((idy, nm, orient, off, al, q[:30]))
    res.sort(reverse=True)
    return res[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sam", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="First SAM alignment start to inspect (1-based)",
    )
    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="Last SAM alignment start to inspect (inclusive)",
    )
    args = parser.parse_args()
    if args.start < 1 or args.end < args.start:
        parser.error("Require 1 <= start <= end")
    f = args.sam
    shown = {"U3": 0, "U5": 0}
    print(f"LTR evidence at {args.chrom}:{args.start}-{args.end}\n")
    for line in open(f):
        if not line.strip() or line.startswith("@"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 11 or c[5] == "*" or c[9] == "*":
            continue
        if c[2] != args.chrom or not (args.start <= int(c[3]) <= args.end):
            continue
        pos, cig, seq = int(c[3]), c[5], c[9]
        lead, trail, ref = parse(cig)
        if lead >= trail and lead > 0:
            junction, clip, side = pos, seq[:lead], "LTR-on-left"
        elif trail > 0:
            junction, clip, side = pos + ref, seq[-trail:], "LTR-on-right"
        else:
            continue
        idy, nm, orient, off, al, q = classify(clip)
        if idy < 0.85 or shown[nm] >= 2:
            continue
        shown[nm] += 1
        end = "5' end of provirus" if nm == "U3" else "3' end of provirus"
        print(f"[{nm}]  junction {args.chrom}:{junction:,}   CIGAR={cig}  ({side})")
        print(f"      soft-clip (vector side): {q}")
        print(f"      matches {nm} LTR ({orient}):   {al}")
        print(f"      identity {100 * idy:.0f}%  ->  this read spans the {end}\n")
        if shown["U3"] >= 2 and shown["U5"] >= 2:
            break

    print(f"Displayed U3 reads: {shown['U3']}; U5 reads: {shown['U5']}")
