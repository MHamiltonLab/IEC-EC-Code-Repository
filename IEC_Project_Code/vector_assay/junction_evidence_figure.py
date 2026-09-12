#!/usr/bin/env python3
"""Pair a selected-site junction pileup with representative soft-clip sequence matches."""

import argparse
import os
import re
from collections import Counter

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _assay import (
    OUTPUT_DIR,
    figure_style,
    safe_filename,
    sam_path,
    sample_layout,
    save_figure,
    target_sites,
)
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

OUTDIR = OUTPUT_DIR
U3_5P = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGG"
U5_3P = "AGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGT"
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
C_U3, C_U5, C_LTR, TSDc = "#2E78B8", "#C0392B", "#E1B12C", "#F4E5A8"
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse(cig):
    p = _CIG.findall(cig)
    if not p:
        return 0, 0, 0
    lead = int(p[0][0]) if p[0][1] == "S" else 0
    trail = int(p[-1][0]) if p[-1][1] == "S" else 0
    ref = sum(int(n) for n, o in p if o in "MDN=X")
    return lead, trail, ref


def best_align(q, ref):
    q = q[:30]
    best = (0, 0, "")
    if not q:
        return 0, ""
    for o in range(0, len(ref) - len(q) + 1):
        m = sum(q[j] == ref[o + j] for j in range(len(q)))
        if m > best[0]:
            best = (m, o, ref[o : o + len(q)])
    return best[0] / len(q), best[2]


def classify(clip):
    res = []
    for nm, ref in (("U3", U3_5P), ("U5", U5_3P)):
        for q in (clip, rc(clip)):
            idy, seg = best_align(q, ref)
            res.append((idy, nm, q[:30], seg))
    res.sort(key=lambda x: x[0], reverse=True)
    return res[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-index",
        type=int,
        default=1,
        help="1-based row of junction_targets.tsv",
    )
    args = parser.parse_args()
    targets = target_sites()
    if not 1 <= args.target_index <= len(targets):
        parser.error("--target-index must identify a row in junction_targets.tsv")
    LABEL, RAW, CHROM, CENTER, GENE = targets[args.target_index - 1]
    _, PALETTE = sample_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    f = sam_path(RAW)
    reads = {"U3": [], "U5": []}
    examples = {"U3": [], "U5": []}
    for line in open(f):
        if not line.strip() or line.startswith("@"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 11 or c[5] == "*" or c[9] == "*":
            continue
        if c[2] != CHROM or not (CENTER - 300 <= int(c[3]) <= CENTER + 300):
            continue
        pos, cig, seq = int(c[3]), c[5], c[9]
        lead, trail, ref = parse(cig)
        if lead >= trail and lead > 0:
            junction, clip, clo, chi = pos, seq[:lead], pos - lead, pos
        elif trail > 0:
            junction, clip, clo, chi = (
                pos + ref,
                seq[-trail:],
                pos + ref,
                pos + ref + trail,
            )
        else:
            continue
        idy, nm, q, seg = classify(clip)
        if idy < 0.85:
            continue
        reads[nm].append((pos, pos + ref, clo, chi, junction))
        if len(examples[nm]) < 3 and len(clip) >= 25:
            examples[nm].append((q, seg, idy, cig, junction))

    def sub(lst, n=20):
        return (
            lst
            if len(lst) <= n
            else [lst[i] for i in np.linspace(0, len(lst) - 1, n).astype(int)]
        )

    if not reads["U3"] or not reads["U5"]:
        raise ValueError("This paired evidence panel requires both U3 and U5 clips")
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.05, 1.0],
        left=0.04,
        right=0.985,
        top=0.80,
        bottom=0.10,
        wspace=0.12,
    )

    axp = fig.add_subplot(gs[0, 0])
    j3 = Counter(j for *_, j in reads["U3"]).most_common(1)[0][0]
    j5 = Counter(j for *_, j in reads["U5"]).most_common(1)[0][0]
    lo, hi = min(j3, j5) - 95, max(j3, j5) + 95

    def draw(rs, col, y):
        for gs_, ge, clo, chi, jx in rs:
            axp.add_patch(
                Rectangle(
                    (gs_, y), ge - gs_, 0.8, facecolor=col, edgecolor="none", zorder=3
                )
            )
            axp.add_patch(
                Rectangle((clo, y), chi - clo, 0.8, facecolor=C_LTR, zorder=3)
            )
            y += 1
        return y

    y = draw(sub(sorted(reads["U3"])), C_U3, 0)
    y += 1.4
    y = draw(sub(sorted(reads["U5"])), C_U5, y)
    ytop = y
    for j, col in ((j3, C_U3), (j5, C_U5)):
        axp.plot([j, j], [-0.5, ytop + 0.3], color=col, lw=1.6, ls="--", zorder=5)
    a, b = sorted((j3, j5))
    axp.add_patch(
        Rectangle((a, -0.5), b - a, ytop + 0.8, facecolor=TSDc, alpha=0.6, zorder=1)
    )
    axp.text(
        (a + b) / 2,
        ytop + 0.5,
        f"Separation {b - a} bp",
        ha="center",
        fontsize=8.5,
        fontweight="bold",
        color="#7a5c00",
    )
    axp.text(
        lo + 2, -0.35, "5′ (U3) reads", fontsize=8.5, color=C_U3, fontweight="bold"
    )
    axp.text(
        lo + 2,
        (ytop - len(sub(sorted(reads["U5"]))) - 0.6),
        "3′ (U5) reads",
        fontsize=8.5,
        color=C_U5,
        fontweight="bold",
    )
    axp.set_xlim(lo, hi)
    axp.set_ylim(-0.7, ytop + 1.4)
    axp.set_yticks([])
    axp.set_xticks(np.linspace(lo, hi, 4))
    axp.set_xticklabels([f"{int(t):,}" for t in np.linspace(lo, hi, 4)], fontsize=7.5)
    axp.set_xlabel(f"{CHROM} position (bp)", fontsize=8.5)
    axp.set_title(
        f"{LABEL} · {GENE} — read pileup",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
        color=PALETTE[LABEL],
    )
    for sp in ("top", "right", "left"):
        axp.spines[sp].set_visible(False)
    axp.legend(
        handles=[
            Line2D([0], [0], color=C_U3, lw=7, label="host — 5′ (U3)"),
            Line2D([0], [0], color=C_U5, lw=7, label="host — 3′ (U5)"),
            Line2D([0], [0], color=C_LTR, lw=7, label="LTR/vector clip"),
            Line2D([0], [0], color=TSDc, lw=7, label="junction interval"),
        ],
        loc="upper left",
        bbox_to_anchor=(0, 1.0),
        fontsize=7.6,
        frameon=False,
        ncol=2,
    )

    axe = fig.add_subplot(gs[0, 1])
    axe.axis("off")
    axe.set_xlim(0, 1)
    axe.set_ylim(0, 1)
    axe.set_title(
        "Why each read is 5′ or 3′ — the soft-clip sequence",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
    )
    mono = dict(family="DejaVu Sans Mono", fontsize=8.3)
    yy = [0.9]

    def line(txt, color="#222", dy=0.052, bold=False, x=0.01):
        axe.text(
            x,
            yy[0],
            txt,
            color=color,
            **mono,
            va="top",
            fontweight="bold" if bold else "normal",
        )
        yy[0] -= dy

    def block(nm, col, refname, end):
        line(
            f"{'5′' if nm == 'U3' else '3′'} junction reads — clip matches {nm} ({refname})",
            col,
            0.058,
            bold=True,
        )
        for q, seg, idy, cig, junc in examples[nm]:
            mb = "".join("|" if a == c else " " for a, c in zip(q, seg))
            line(f"  read  5'-{q}-3'   ({cig})", col, 0.029)
            line(f"          {mb}        {100 * idy:.0f}% id", "#555", 0.029)
            line(f"  {nm}    5'-{seg}-3'   LTR {end}", "#777", 0.036)
        yy[0] -= 0.010

    block("U3", C_U3, "LTR 5′ end", "U3/5′ terminus")
    block("U5", C_U5, "LTR 3′ end", "U5/3′ terminus")
    line(
        "U3/U5 labels follow the displayed reference-sequence matches.",
        "#222",
        0.05,
        bold=True,
    )
    line(
        "Sequence identity supports the proposed LTR assignment.",
        "#222",
        0.05,
        bold=True,
    )

    fig.suptitle(
        "How 5′ vs 3′ junction reads are identified — pileup + sequence evidence",
        x=0.04,
        ha="left",
        y=0.95,
        fontsize=15,
        fontweight="bold",
        color="#2B2B2B",
    )
    png = os.path.join(
        OUTDIR, f"{safe_filename(RAW)}_{safe_filename(GENE)}_junction_evidence.png"
    )
    save_figure(fig, png)
    plt.close(fig)
    print("Wrote:", png)
