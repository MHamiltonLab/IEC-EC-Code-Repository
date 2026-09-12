#!/usr/bin/env python3
"""Annotate integration sites by host gene and display gene-level and genome-wide molecule support."""

import os
from collections import Counter, defaultdict

import matplotlib
import numpy as np
from _assay import (
    GENE_REFERENCE,
    OUTPUT_DIR,
    figure_style,
    load_annotation_engine,
    sam_path,
    sample_layout,
    save_figure,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUTDIR = OUTPUT_DIR
INK = "#2B2B2B"
HG19 = {
    "chr1": 249250621,
    "chr2": 243199373,
    "chr3": 198022430,
    "chr4": 191154276,
    "chr5": 180915260,
    "chr6": 171115067,
    "chr7": 159138663,
    "chr8": 146364022,
    "chr9": 141213431,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chrX": 155270560,
    "chrY": 59373566,
}
CHROMS = list(HG19)
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 120,
        "axes.titleweight": "bold",
    }
)


def load_high(stem):
    return VA.dedupe_umi(list(VA.parse_sam(sam_path(stem))))


if __name__ == "__main__":
    VA = load_annotation_engine()
    SAMPLES, PALETTE = sample_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    gm = VA.GeneModel(str(GENE_REFERENCE))
    S = {}
    for stem, lab in SAMPLES.items():
        mol = load_high(stem)
        loci = VA.call_sites(mol, window=100)
        tot = sum(l["support"] for l in loci) or 1
        genes = defaultdict(
            lambda: {
                "support": 0,
                "onc": False,
                "chrom": None,
                "pos": None,
                "region": None,
                "_t": 0,
            }
        )
        for l in loci:
            a = gm.annotate(l["chrom"], l["pos"])
            g = genes[a["gene"]]
            g["support"] += l["support"]
            if l["support"] > g["_t"]:
                g.update(
                    _t=l["support"],
                    chrom=l["chrom"],
                    pos=l["pos"],
                    region=a["region"],
                    onc=a["oncogene"],
                )
        S[lab] = {
            "loci": loci,
            "tot": tot,
            "genes": genes,
            "nmol": len(mol),
            "chrom": Counter(r["chrom"] for r in mol),
        }
    samples = list(SAMPLES.values())

    with open(os.path.join(OUTDIR, "host_gene_support.tsv"), "w") as fh:
        fh.write("sample\tgene\tchrom\tpos\tsupport\tpct\tregion\toncogene\n")
        for s in samples:
            for g, i in sorted(S[s]["genes"].items(), key=lambda kv: -kv[1]["support"]):
                fh.write(
                    f"{s}\t{g}\t{i['chrom']}\t{i['pos']}\t{i['support']}\t{100 * i['support'] / S[s]['tot']:.2f}\t{i['region']}\t{i['onc']}\n"
                )

    print(f"{'sample':6s} {'molecules':>9s} {'loci':>5s}  dominant gene (top4)")
    for s in samples:
        top = sorted(S[s]["genes"].items(), key=lambda kv: -kv[1]["support"])[:4]
        desc = ", ".join(
            f"{g} {100 * i['support'] / S[s]['tot']:.0f}%{'★' if i['onc'] else ''}"
            for g, i in top
        )
        print(f"{s:6s} {S[s]['nmol']:>9d} {len(S[s]['loci']):>5d}  {desc}")

    rows = []
    for s in samples:
        for g, i in sorted(S[s]["genes"].items(), key=lambda kv: -kv[1]["support"])[:5]:
            rows.append((s, g, i))
    rows.sort(key=lambda r: r[2]["support"])
    if not rows:
        raise ValueError(
            "No annotated integration loci were returned by the external engine"
        )
    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(left=0.22, right=0.76, top=0.92, bottom=0.09)
    fm = {
        "CDS exon": "s",
        "UTR exon": "s",
        "exon": "s",
        "intron": "o",
        "intergenic": "D",
    }
    for i, (s, g, info) in enumerate(rows):
        ax.hlines(i, 1, info["support"], color=PALETTE[s], lw=2.3, alpha=0.85)
        ax.scatter(
            info["support"],
            i,
            s=130,
            marker=fm.get(info["region"], "o"),
            color=PALETTE[s],
            edgecolor="white",
            lw=0.8,
            zorder=3,
        )
        ax.text(
            info["support"] * 1.15,
            i,
            f"{info['support']:,} ({100 * info['support'] / S[s]['tot']:.0f}%)",
            va="center",
            fontsize=8.2,
            color=INK,
        )
        if info["onc"]:
            ax.text(info["support"] * 1.04, i + 0.35, "★", fontsize=12, color="#B00000")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(
        [f"{g}  ·  {info['chrom']}:{info['pos'] / 1e6:.2f} Mb" for s, g, info in rows],
        fontsize=8.6,
    )
    for t, (s, g, info) in zip(ax.get_yticklabels(), rows):
        t.set_color("#B00000" if info["onc"] else PALETTE[s])
        t.set_fontweight("bold" if info["onc"] else "normal")
    ax.set_xscale("log")
    ax.set_xlim(1, max(r[2]["support"] for r in rows) * 2.6)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel("Supporting molecules (log)")
    ax.set_title("Top host genes (★ = cancer-gene annotation)", loc="left")
    leg1 = ax.legend(
        handles=[
            Line2D([0], [0], marker="o", ls="", color=PALETTE[s], label=s, ms=9)
            for s in samples
        ],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=8.6,
        title="Sample",
        title_fontsize=9,
        frameon=False,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", ls="", color="#888", label="intron", ms=9),
            Line2D([0], [0], marker="s", ls="", color="#888", label="exon/UTR", ms=9),
            Line2D([0], [0], marker="D", ls="", color="#888", label="intergenic", ms=8),
        ],
        loc="upper left",
        bbox_to_anchor=(1.02, 0.55),
        fontsize=8.6,
        title="Feature",
        title_fontsize=9,
        frameon=False,
    )
    save_figure(fig, f"{OUTDIR}/top_host_genes.png")
    plt.close(fig)

    off = {}
    cum = 0
    for c in CHROMS:
        off[c] = cum
        cum += HG19[c]
    GL = cum
    fig, ax = plt.subplots(figsize=(15, 4.2))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.84, bottom=0.10)
    for i, s in enumerate(samples):
        y = len(samples) - 1 - i
        for c in CHROMS:
            ax.add_patch(
                plt.Rectangle(
                    (off[c], y - 0.16),
                    HG19[c],
                    0.32,
                    facecolor="#F0F0F0",
                    edgecolor="#D5D5D5",
                    lw=0.4,
                )
            )
        xs = [off[l["chrom"]] + l["pos"] for l in S[s]["loci"] if l["chrom"] in off]
        sup = np.array([l["support"] for l in S[s]["loci"] if l["chrom"] in off], float)
        if len(sup):
            ax.scatter(
                xs,
                [y] * len(xs),
                s=12
                + 60 * np.log10(sup) / np.log10(sup.max() if sup.max() > 1 else 10),
                color=PALETTE[s],
                alpha=0.8,
                edgecolor="white",
                lw=0.4,
                zorder=3,
            )
        ax.text(
            -GL * 0.015,
            y,
            s,
            ha="right",
            va="center",
            fontweight="bold",
            color=PALETTE[s],
            fontsize=11,
        )
    for c in CHROMS:
        ax.text(
            off[c] + HG19[c] / 2,
            len(samples) - 0.35,
            c.replace("chr", ""),
            ha="center",
            va="bottom",
            fontsize=7.5,
            color="#888",
        )
    ax.set_xlim(-GL * 0.05, GL)
    ax.set_ylim(-0.6, len(samples) - 0.05)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(
        "Genome-wide lentiviral integration map (marker ∝ log molecule support)",
        loc="left",
        fontsize=13,
        pad=16,
    )
    save_figure(fig, f"{OUTDIR}/genome_integration_map.png")
    plt.close(fig)

    print("\nWrote tables + figures to", OUTDIR)
