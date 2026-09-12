#!/usr/bin/env python3
"""Display read-level U3/U5 junction pileups for metadata-selected integration sites."""

import math
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


def ident(a, b):
    a = a[: len(b)]
    return max(
        (
            sum(a[j] == b[o + j] for j in range(len(a)))
            for o in range(len(b) - len(a) + 1)
        ),
        default=0,
    ) / max(len(a), 1)


def classify(clip):
    if len(clip) < 10:
        return "none"
    best = ("none", 0)
    for nm, ref in (("U3", U3_5P), ("U5", U5_3P)):
        for q in (clip, rc(clip)):
            s = ident(q[:30], ref)
            if s > best[1]:
                best = (nm, s)
    return best[0] if best[1] >= 0.85 else "none"


_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def parse(cig):
    p = _CIG.findall(cig)
    if not p:
        return 0, 0, 0
    lead = int(p[0][0]) if p[0][1] == "S" else 0
    trail = int(p[-1][0]) if p[-1][1] == "S" else 0
    ref = sum(int(n) for n, o in p if o in "MDN=X")
    return lead, trail, ref


def collect(raw, chrom, center, win=300):
    f = sam_path(raw)
    out = {"U3": [], "U5": []}
    for line in open(f):
        if not line.strip() or line.startswith("@"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 11 or c[5] == "*" or c[9] == "*":
            continue
        if c[2] != chrom or not (center - win <= int(c[3]) <= center + win):
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
        k = classify(clip)
        if k in ("U3", "U5"):
            out[k].append((pos, pos + ref, clo, chi, junction))
    return out


def subsample(lst, n=22):
    if len(lst) <= n:
        return lst
    idx = np.linspace(0, len(lst) - 1, n).astype(int)
    return [lst[i] for i in idx]


if __name__ == "__main__":
    SITES = target_sites()
    _, PALETTE = sample_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 120})
    nrows = math.ceil(len(SITES) / 2)
    figure_height = 2.2 + 4.2 * nrows
    fig, axes = plt.subplots(nrows, 2, figsize=(14, figure_height), squeeze=False)
    fig.subplots_adjust(
        left=0.05, right=0.97, top=1 - 1.8 / figure_height, bottom=0.07, hspace=0.42, wspace=0.13
    )

    for ax, (lab, raw, chrom, center, gene) in zip(axes.ravel(), SITES):
        d = collect(raw, chrom, center)

        def mode(lst):
            cc = Counter(j for *_, j in lst)
            return cc.most_common(1)[0][0] if cc else None

        j3 = mode(d["U3"])
        j5 = mode(d["U5"])
        js = [x for x in (j3, j5) if x]
        if not js:
            ax.set_title(f"{lab} · {gene}: no classified LTR clips")
            ax.axis("off")
            continue
        lo = min(js) - 95
        hi = max(js) + 95

        u5 = subsample(sorted(d["U5"]))
        u3 = subsample(sorted(d["U3"]))
        gap = 1.2

        def draw(reads, col, y):
            for gs, ge, clo, chi, jx in reads:
                ax.add_patch(
                    Rectangle(
                        (gs, y), ge - gs, 0.8, facecolor=col, edgecolor="none", zorder=3
                    )
                )
                ax.add_patch(
                    Rectangle(
                        (clo, y),
                        chi - clo,
                        0.8,
                        facecolor=C_LTR,
                        edgecolor="none",
                        alpha=0.95,
                        zorder=3,
                    )
                )
                y += 1
            return y

        y = draw(u3, C_U3, 0)
        y += gap
        ybreak = y - gap / 2
        y = draw(u5, C_U5, y)
        ytop = y

        for j, c in ((j3, C_U3), (j5, C_U5)):
            if j:
                ax.plot([j, j], [-0.5, ytop + 0.3], color=c, lw=1.6, ls="--", zorder=5)
        if j3 and j5:
            a, b = sorted((j3, j5))
            ax.add_patch(
                Rectangle(
                    (a, -0.5), b - a, ytop + 0.8, facecolor=TSDc, alpha=0.6, zorder=1
                )
            )
            ax.text(
                (a + b) / 2,
                ytop + 1.0,
                f"Separation {b - a} bp",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold",
                color="#7a5c00",
            )
            sub = f"junction separation {b - a} bp"
        else:
            sub = "one LTR class captured"
        ax.text(
            lo + 2,
            ybreak + 0.1,
            "3′ (U5) reads",
            fontsize=8,
            color=C_U5,
            fontweight="bold",
            va="bottom",
        )
        ax.text(
            lo + 2,
            -0.4,
            "5′ (U3) reads",
            fontsize=8,
            color=C_U3,
            fontweight="bold",
            va="bottom",
        )
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.8, ytop + 2.0)
        ax.set_yticks([])
        ax.set_xticks(np.linspace(lo, hi, 5))
        ax.set_xticklabels(
            [f"{int(t):,}" for t in np.linspace(lo, hi, 5)], fontsize=7.5
        )
        ax.set_xlabel(f"{chrom} position (bp)", fontsize=8.5)
        ax.set_title(
            f"{lab} · {gene}  —  {sub}",
            loc="left",
            fontsize=11,
            fontweight="bold",
            color=PALETTE[lab],
            pad=8,
        )
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)

    for unused in axes.ravel()[len(SITES) :]:
        unused.axis("off")

    fig.suptitle(
        "Read-level LTR junction evidence",
        x=0.05,
        ha="left",
        y=1 - 0.1 / figure_height,
        fontsize=16,
        fontweight="bold",
        color="#2B2B2B",
    )
    fig.text(
        0.05,
        1 - 0.65 / figure_height,
        "Each bar = one read · host-mapped portion in the junction colour, "
        "soft-clipped LTR/vector portion in gold",
        ha="left",
        fontsize=9.8,
        color="#666",
    )
    fig.legend(
        handles=[
            Line2D([0], [0], color=C_U3, lw=7, label="host genome — 5′ (U3) reads"),
            Line2D([0], [0], color=C_U5, lw=7, label="host genome — 3′ (U5) reads"),
            Line2D([0], [0], color=C_LTR, lw=7, label="LTR / vector (soft-clip)"),
            Line2D([0], [0], color=TSDc, lw=7, label="junction interval"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 1 - 1.05 / figure_height),
        ncol=4,
        fontsize=8.6,
        frameon=False,
    )
    png = os.path.join(OUTDIR, "junction_pileup.png")
    save_figure(fig, png)
    plt.close(fig)
    print("Wrote:", png)
