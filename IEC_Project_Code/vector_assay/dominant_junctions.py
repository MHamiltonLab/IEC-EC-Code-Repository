#!/usr/bin/env python3
"""Select dominant integration sites and plot LTR-specific junction coverage and read evidence."""

import os
import re
from collections import defaultdict

import matplotlib
import numpy as np
from _assay import (
    GENE_REFERENCE,
    OUTPUT_DIR,
    figure_style,
    load_annotation_engine,
    safe_filename,
    sam_path,
    sample_layout,
    save_figure,
)

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

OUTDIR = OUTPUT_DIR

PSDIR = OUTPUT_DIR / "dominant_junctions"

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


def samfile(stem):
    return sam_path(stem)


def collect(stem, chrom, center, win=320):
    out = {"U3": [], "U5": []}
    for line in open(samfile(stem)):
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
        mb = c[0].split("_molbar")[0]
        if k in ("U3", "U5"):
            out[k].append((pos, pos + ref, clo, chi, junction, mb))
    return out


def modej(lst):
    bj = defaultdict(set)
    for *_, j, mb in lst:
        bj[j].add(mb)
    if not bj:
        return None, 0
    j = max(bj, key=lambda k: len(bj[k]))
    return j, len(bj[j])


def sub(lst, n=26):
    return (
        lst
        if len(lst) <= n
        else [lst[i] for i in np.linspace(0, len(lst) - 1, n).astype(int)]
    )


def generate_junction_reports(sites, palette, output_dir, table_path):
    """Render the common per-site report for dominant or metadata-selected sites."""
    SITES, PAL, PSDIR = sites, palette, output_dir
    PSDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
        }
    )

    TABLE = []

    for lab, stem, chrom, center, gene in SITES:
        d = collect(stem, chrom, center)
        j3, n3 = modej(d["U3"])
        j5, n5 = modej(d["U5"])
        gap = abs(j3 - j5) if (j3 and j5) else ""
        verdict = (
            "short separation (<=20 bp)"
            if (j3 and j5 and gap <= 20)
            else (
                "wide separation (>20 bp)"
                if (j3 and j5)
                else "incomplete paired-junction evidence"
            )
        )
        TABLE.append(
            {
                "sample": lab,
                "gene": gene,
                "chrom": chrom,
                "U3_5prime": j3 or "",
                "U3_mol": n3,
                "U5_3prime": j5 or "",
                "U5_mol": n5,
                "junction_separation_bp": gap,
                "verdict": verdict,
            }
        )

        js = [x for x in (j3, j5) if x]
        if not js:
            continue
        lo = min(js) - 90
        hi = max(js) + 90
        fig, (axc, axp) = plt.subplots(
            2,
            1,
            figsize=(9, 8),
            gridspec_kw={"height_ratios": [1, 1.25], "hspace": 0.40},
        )
        xs = np.arange(lo, hi + 1)

        def cov(lst):
            a = np.zeros(len(xs))
            for gs_, ge, *_ in lst:
                a[max(0, gs_ - lo) : max(0, ge - lo)] += 1
            return a

        c5 = cov(d["U3"])
        c3 = cov(d["U5"])
        mx = max(c5.max(), c3.max(), 1)
        axc.fill_between(xs, 0, c3 / mx, color=C_U5, alpha=0.4, step="mid")
        axc.fill_between(xs, 0, c5 / mx, color=C_U3, alpha=0.4, step="mid")
        for j, col in ((j3, C_U3), (j5, C_U5)):
            if j:
                axc.plot([j, j], [0, 1.08], color=col, lw=2, zorder=5)
        span = hi - lo
        if j3 and j5:
            a, b = sorted((j3, j5))
            mid = (a + b) / 2
            axc.add_patch(
                Rectangle(
                    (a, 0), max(b - a, 0.6), 1.0, facecolor=TSDc, alpha=0.6, zorder=4
                )
            )
            axc.annotate(
                f"5′ (U3)\n{n3:,} mol",
                xy=(j3, 1.08),
                xytext=(lo + span * 0.16, 1.6),
                ha="center",
                fontsize=8.6,
                color=C_U3,
                fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_U3, lw=0.9),
            )
            axc.annotate(
                f"3′ (U5)\n{n5:,} mol",
                xy=(j5, 1.08),
                xytext=(hi - span * 0.16, 1.6),
                ha="center",
                fontsize=8.6,
                color=C_U5,
                fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=C_U5, lw=0.9),
            )
            axc.add_patch(
                Polygon(
                    [(mid - span * 0.05, 2.02), (mid + span * 0.05, 2.02), (mid, 1.72)],
                    closed=True,
                    facecolor="#444",
                    zorder=6,
                )
            )
            axc.text(
                mid,
                2.12,
                "provirus\n(not to scale)",
                ha="center",
                va="bottom",
                fontsize=7.4,
                color="#444",
                fontweight="bold",
            )
            axc.annotate(
                f"Junction separation = {b - a} bp",
                xy=(mid, 0.5),
                xytext=(mid, -0.42),
                ha="center",
                fontsize=9,
                fontweight="bold",
                color="#7a5c00",
                arrowprops=dict(arrowstyle="-", color="#7a5c00", lw=0.9),
            )
        else:
            if j3:
                axc.annotate(
                    f"5′ (U3)\n{n3:,} mol",
                    xy=(j3, 1.08),
                    xytext=(j3, 1.7),
                    ha="center",
                    fontsize=8.6,
                    color=C_U3,
                    fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=C_U3, lw=0.9),
                )
            if j5:
                axc.annotate(
                    f"3′ (U5)\n{n5:,} mol",
                    xy=(j5, 1.08),
                    xytext=(j5, 1.7),
                    ha="center",
                    fontsize=8.6,
                    color=C_U5,
                    fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=C_U5, lw=0.9),
                )
            axc.text(
                (lo + hi) / 2,
                2.0,
                "one LTR class captured",
                ha="center",
                fontsize=8,
                color="#999",
                style="italic",
            )
        axc.set_xlim(lo, hi)
        axc.set_ylim(-0.55, 2.4)
        axc.set_yticks([])
        axc.set_xticks(np.linspace(lo, hi, 5))
        axc.set_xticklabels(
            [f"{int(t):,}" for t in np.linspace(lo, hi, 5)], fontsize=7.5
        )
        axc.set_xlabel(f"{chrom} position (bp)", fontsize=8.5)
        axc.set_title(
            "Coverage view — 5′ (U3) & 3′ (U5) junction separation",
            loc="left",
            fontsize=10.5,
            fontweight="bold",
        )
        for sp in ("top", "right", "left"):
            axc.spines[sp].set_visible(False)

        def draw(rs, col, y):
            for gs_, ge, clo, chi, jx, mb in rs:
                axp.add_patch(
                    Rectangle((gs_, y), ge - gs_, 0.8, facecolor=col, zorder=3)
                )
                axp.add_patch(
                    Rectangle((clo, y), chi - clo, 0.8, facecolor=C_LTR, zorder=3)
                )
                y += 1
            return y

        y = draw(sub(sorted(d["U3"])), C_U3, 0)
        y += 1.3
        ybreak = y - 0.65
        y = draw(sub(sorted(d["U5"])), C_U5, y)
        ytop = max(y, 1)
        for j, col in ((j3, C_U3), (j5, C_U5)):
            if j:
                axp.plot(
                    [j, j], [-0.5, ytop + 0.3], color=col, lw=1.5, ls="--", zorder=5
                )
        if j3 and j5:
            a, b = sorted((j3, j5))
            axp.add_patch(
                Rectangle(
                    (a, -0.5), b - a, ytop + 0.8, facecolor=TSDc, alpha=0.6, zorder=1
                )
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
        if d["U5"]:
            axp.text(
                lo + 2,
                ybreak + 0.3,
                "3′ (U5) reads",
                fontsize=8.5,
                color=C_U5,
                fontweight="bold",
            )
        axp.set_xlim(lo, hi)
        axp.set_ylim(-0.7, ytop + 1.3)
        axp.set_yticks([])
        axp.set_xticks(np.linspace(lo, hi, 5))
        axp.set_xticklabels(
            [f"{int(t):,}" for t in np.linspace(lo, hi, 5)], fontsize=7.5
        )
        axp.set_xlabel(f"{chrom} position (bp)", fontsize=8.5)
        axp.set_title(
            "Read-level pileup (gold = LTR/vector clip)",
            loc="left",
            fontsize=10.5,
            fontweight="bold",
        )
        for sp in ("top", "right", "left"):
            axp.spines[sp].set_visible(False)
        axp.legend(
            handles=[
                Line2D([0], [0], color=C_U3, lw=7, label="host 5′ (U3)"),
                Line2D([0], [0], color=C_U5, lw=7, label="host 3′ (U5)"),
                Line2D([0], [0], color=C_LTR, lw=7, label="LTR/vector"),
                Line2D([0], [0], color=TSDc, lw=7, label="junction interval"),
            ],
            loc="upper center",
            bbox_to_anchor=(0.5, -0.22),
            fontsize=7.6,
            frameon=False,
            ncol=2,
        )
        fig.suptitle(
            f"{lab} · {gene} — integration junction",
            x=0.06,
            ha="left",
            y=0.97,
            fontsize=14.5,
            fontweight="bold",
            color=PAL[lab],
        )
        fig.text(
            0.06, 0.925, f"{chrom}  ·  {verdict}", ha="left", fontsize=9.8, color="#666"
        )
        out = os.path.join(
            PSDIR, f"{safe_filename(lab)}_{safe_filename(gene)}_junction.png"
        )
        save_figure(fig, out)
        plt.close(fig)
        print("Wrote:", os.path.basename(out))

    with open(table_path, "w") as fh:
        cols = [
            "sample",
            "gene",
            "chrom",
            "U3_5prime",
            "U3_mol",
            "U5_3prime",
            "U5_mol",
            "junction_separation_bp",
            "verdict",
        ]
        fh.write("\t".join(cols) + "\n")
        for r in TABLE:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("\nJunction table:")

    for r in TABLE:
        print(
            f"  {r['sample']:5s} {r['gene']:9s} separation={r['junction_separation_bp']!s:>4} {r['verdict']}"
        )


if __name__ == "__main__":
    VA = load_annotation_engine()
    SAMPLES, PAL = sample_layout()
    PSDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    gm = VA.GeneModel(str(GENE_REFERENCE))

    SITES = []

    for stem, lab in SAMPLES.items():
        mol = VA.dedupe_umi(list(VA.parse_sam(samfile(stem))))
        loci = VA.call_sites(mol, window=100)
        if not loci:
            continue
        top = loci[0]
        gene = gm.annotate(top["chrom"], top["pos"])["gene"]
        SITES.append((lab, stem, top["chrom"], top["pos"], gene))

    generate_junction_reports(SITES, PAL, PSDIR, OUTDIR / "dominant_junctions.tsv")
