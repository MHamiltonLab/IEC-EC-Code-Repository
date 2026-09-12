#!/usr/bin/env python3
"""Map SAM soft-clips to a supplied vector and report read-start primer candidates."""

import re
from collections import Counter, defaultdict

import numpy as np
from _assay import VECTOR_FASTA, clip_paths

comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]


def load(p):
    with open(p) as stream:
        lines = list(stream)
    if sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("vector.fasta must contain exactly one reference sequence")
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


K = 12
VEC = VRC = ""
VL = 0
idxF = defaultdict(list)
idxR = defaultdict(list)


def initialize_vector(path=VECTOR_FASTA):
    """Index both orientations of the single supplied vector reference."""
    global VEC, VRC, VL, idxF, idxR
    VEC = load(path)
    VL = len(VEC)
    if VL < K:
        raise ValueError("Vector reference is shorter than the 12-base mapping seed")
    VRC = rc(VEC)
    idxF, idxR = defaultdict(list), defaultdict(list)
    for i in range(VL - K + 1):
        idxF[VEC[i : i + K]].append(i)
        idxR[VRC[i : i + K]].append(i)


_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def clips(path, minlen=20):
    """yield (clip_seq) for soft-clips >= minlen from mapped reads."""
    for line in open(path):
        if not line.strip() or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or f[5] == "*" or f[9] == "*":
            continue
        p = _CIG.findall(f[5])
        seq = f[9]
        if not p:
            continue
        if p[0][1] == "S" and int(p[0][0]) >= minlen:
            yield seq[: int(p[0][0])]
        if p[-1][1] == "S" and int(p[-1][0]) >= minlen:
            yield seq[-int(p[-1][0]) :]


def map_clip(clip):
    """best (strand, vstart, vend, ident) of clip on the vector (+ coords)."""
    best = None
    for strand, (V, idx) in (("+", (VEC, idxF)), ("-", (VRC, idxR))):
        votes = Counter()
        for i in range(len(clip) - K + 1):
            for pos in idx.get(clip[i : i + K], ()):
                votes[pos - i] += 1
        if not votes:
            continue
        off, nv = votes.most_common(1)[0]
        cs, vs = (-off, 0) if off < 0 else (0, off)
        L = min(len(clip) - cs, len(V) - vs)
        if L < K:
            continue
        m = sum(clip[cs + j] == V[vs + j] for j in range(L))
        ident = m / L
        if strand == "+":
            a, b = vs, vs + L
        else:
            a, b = VL - (vs + L), VL - vs
        if best is None or ident * L > best[3] * (best[2] - best[1]):
            best = (strand, a, b, ident)
    return best


def coverage_profile(paths):
    """Pool clip coverage and strand-specific read starts across SAM inputs."""
    cov = np.zeros(VL)
    startF, startR = Counter(), Counter()
    nmap = ntot = 0
    for path in paths:
        for clip in clips(path, minlen=20):
            ntot += 1
            bm = map_clip(clip)
            if bm and bm[3] >= 0.90:
                nmap += 1
                s, a, b, idt = bm
                cov[a:b] += 1
                if s == "+":
                    startF[a] += 1
                else:
                    startR[b - 1] += 1

    return cov, startF, startR, ntot, nmap


if __name__ == "__main__":
    initialize_vector()
    cov, startF, startR, ntot, nmap = coverage_profile(clip_paths())
    print(
        f"clips>=20bp examined: {ntot:,} ; mapped to vector (>=90% id): {nmap:,} "
        f"({100 * nmap / max(ntot, 1):.0f}%)\n"
    )

    print("vector coverage (100-bp bins with >0):")
    for b0 in range(0, VL, 100):
        c = cov[b0 : b0 + 100].sum()
        if c > 0:
            print(f"  {b0:>5d}-{b0 + 100:<5d}: {int(c)}")
    print("\ntop read-start positions (+ strand, vector pos):", startF.most_common(6))
    print("top read-start positions (- strand, vector pos):", startR.most_common(6))
    print(f"\nvector length {VL}; LTR termini ~ 0 (5'R) and ~{VL} (3' end)")

    print(
        "\n=== candidate primer sequences (vector LTR subsequence at read-start peaks) ==="
    )
    print("(+strand read starting at P: the primer's 3' end is just 5' of P)")
    for pos, n in startF.most_common(8):
        if n < 100 or pos < 25:
            continue
        print(
            f"  +  start@{pos:<5d} n={n:<5d}  primer(+) 5'-{VEC[pos - 25 : pos]}-3'   read->{VEC[pos : pos + 12]}..."
        )
    print("(-strand read starting at P: primer on - strand, 3' end at P)")
    for pos, n in startR.most_common(8):
        if n < 100 or pos > VL - 25:
            continue
        print(
            f"  -  start@{pos:<5d} n={n:<5d}  primer(-) 5'-{rc(VEC[pos + 1 : pos + 26])}-3'   read->{rc(VEC[pos - 11 : pos + 1])}..."
        )
