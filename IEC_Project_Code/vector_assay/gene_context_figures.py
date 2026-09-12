#!/usr/bin/env python3
"""Aggregate annotated loci by host gene and display transcript context and gene-list flags."""

import math
import os

import matplotlib
import numpy as np
from _assay import OUTPUT_DIR, figure_style, sample_layout, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vectorint_annotate as VA
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

OUTDIR = OUTPUT_DIR
LOCUS_WINDOW = 100
INK = "#2B2B2B"
ONC_RED = "#B00000"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 120,
    }
)

if __name__ == "__main__":
    SAMPLE_LABEL, PALETTE = sample_layout()
    SAMPLE_ORDER = list(SAMPLE_LABEL.values())
    summary_id = os.environ.get("VECTOR_SUMMARY_SAMPLE")
    if summary_id and summary_id not in SAMPLE_LABEL:
        raise ValueError("VECTOR_SUMMARY_SAMPLE must be a sample_id from samples.tsv")
    SUMMARY_LABEL = SAMPLE_LABEL.get(summary_id)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    gm = VA.GeneModel(VA.REFGENE)
    per_sample = {}
    all_rows = []
    for raw, lab in SAMPLE_LABEL.items():
        d = VA.load_sample(raw)
        loci = VA.call_sites(d["HIGH"], window=LOCUS_WINDOW)
        total = sum(l["support"] for l in loci) or 1
        for l in loci:
            l.update(gm.annotate(l["chrom"], l["pos"]))
            l["sample"], l["pct"] = lab, 100 * l["support"] / total
            all_rows.append(l)
        genes = {}
        for l in loci:
            g = genes.setdefault(
                l["gene"],
                {
                    "gene": l["gene"],
                    "support": 0,
                    "oncogene": l["oncogene"],
                    "chrom": l["chrom"],
                    "loci": [],
                },
            )
            g["support"] += l["support"]
            g["loci"].append(l)
        for g in genes.values():
            g["loci"].sort(key=lambda x: -x["support"])
            g["pos"] = g["loci"][0]["pos"]
            g["region"] = g["loci"][0]["region"]
            g["pct"] = 100 * g["support"] / total
            g["n_loci"] = len(g["loci"])
        per_sample[lab] = {
            "genes": sorted(genes.values(), key=lambda x: -x["support"]),
            "loci": loci,
            "total": sum(l["support"] for l in loci),
        }

    samples = [s for s in SAMPLE_ORDER if s in per_sample]

    all_rows.sort(key=lambda r: (SAMPLE_ORDER.index(r["sample"]), -r["support"]))
    with open(os.path.join(OUTDIR, "insertion_sites_annotated.tsv"), "w") as fh:
        cols = [
            "sample",
            "chrom",
            "pos",
            "support",
            "pct",
            "gene",
            "region",
            "dist_tss",
            "strand",
            "oncogene",
        ]
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write(
                "\t".join(
                    f"{r[c]:.2f}" if c == "pct" else str(r.get(c, "")) for c in cols
                )
                + "\n"
            )
    with open(os.path.join(OUTDIR, "insertion_genes_annotated.tsv"), "w") as fh:
        fh.write("sample\tgene\tchrom\tpos\tsupport\tpct\tn_loci\tregion\toncogene\n")
        for s in samples:
            for g in per_sample[s]["genes"]:
                fh.write(
                    f"{s}\t{g['gene']}\t{g['chrom']}\t{g['pos']}\t{g['support']}\t"
                    f"{g['pct']:.2f}\t{g['n_loci']}\t{g['region']}\t{g['oncogene']}\n"
                )

    context_samples = [
        s for s in samples if s != SUMMARY_LABEL and per_sample[s]["genes"]
    ]
    panel_count = len(context_samples) + 1 + int(SUMMARY_LABEL is not None)
    panel_rows = max(1, math.ceil(panel_count / 3))
    fig = plt.figure(figsize=(15.5, 5.5 + 4.5 * panel_rows))
    gs = fig.add_gridspec(
        1 + panel_rows,
        3,
        height_ratios=[1.55] + [1.0] * panel_rows,
        hspace=0.55,
        wspace=0.32,
        left=0.075,
        right=0.975,
        top=0.81,
        bottom=0.055,
    )

    def panel_letter(ax, L):
        ax.text(
            -0.02,
            1.06,
            L,
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="bottom",
            ha="right",
            color=INK,
        )

    axA = fig.add_subplot(gs[0, :])
    rows = []
    for s in samples:
        for g in per_sample[s]["genes"][:5]:
            rows.append({**g, "sample": s})
    rows.sort(key=lambda r: r["support"])
    if not rows:
        raise ValueError(
            "No integration loci were available for the gene-context figure"
        )
    y = np.arange(len(rows))
    feat_marker = {
        "CDS exon": "s",
        "UTR exon": "s",
        "exon": "s",
        "intron": "o",
        "intergenic": "D",
    }
    for yi, r in zip(y, rows):
        col = PALETTE[r["sample"]]
        axA.hlines(yi, 1, r["support"], color=col, lw=2.3, alpha=0.85, zorder=2)
        axA.scatter(
            r["support"],
            yi,
            s=130,
            marker=feat_marker.get(r["region"], "o"),
            color=col,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        extra = f"  ·  {r['n_loci']} loci" if r["n_loci"] > 1 else ""
        axA.text(
            r["support"] * 1.13,
            yi,
            f"{r['support']:,} ({r['pct']:.0f}%){extra}",
            va="center",
            ha="left",
            fontsize=8.2,
            color=INK,
        )
        if r["oncogene"]:
            axA.text(
                r["support"] * 1.04,
                yi + 0.34,
                "★",
                va="center",
                ha="left",
                fontsize=12,
                color=ONC_RED,
            )
    labels = [f"{r['gene']}  ·  {r['chrom']}:{r['pos'] / 1e6:.2f} Mb" for r in rows]
    axA.set_yticks(y)
    axA.set_yticklabels(labels, fontsize=8.6)
    for tick, r in zip(axA.get_yticklabels(), rows):
        tick.set_color(ONC_RED if r["oncogene"] else PALETTE[r["sample"]])
        if r["oncogene"]:
            tick.set_fontweight("bold")
    axA.set_xscale("log")
    axA.set_xlim(1, max(r["support"] for r in rows) * 2.6)
    axA.set_ylim(-0.8, len(rows) - 0.2)
    axA.set_xlabel("Supporting molecules (UMI, log scale)")
    axA.set_title(
        "Top host genes — clonal abundance, genomic feature & cancer-gene flag (★)",
        loc="left",
        pad=10,
    )
    axA.grid(axis="x", color="#ECECEC", lw=0.8, zorder=0)
    samp_handles = [
        Line2D([0], [0], marker="o", ls="", color=PALETTE[s], label=s, ms=9)
        for s in samples
    ]
    feat_handles = [
        Line2D([0], [0], marker="o", ls="", color="#888", label="intron", ms=9),
        Line2D([0], [0], marker="s", ls="", color="#888", label="exon/UTR", ms=9),
        Line2D([0], [0], marker="D", ls="", color="#888", label="intergenic", ms=8),
    ]
    leg1 = axA.legend(
        handles=samp_handles,
        title="Sample",
        fontsize=8.5,
        title_fontsize=9,
        loc="lower right",
        frameon=True,
        framealpha=0.95,
    )
    axA.add_artist(leg1)
    axA.legend(
        handles=feat_handles,
        title="Feature",
        fontsize=8.5,
        title_fontsize=9,
        loc="lower right",
        bbox_to_anchor=(0.81, 0.0),
        frameon=True,
        framealpha=0.95,
    )
    panel_letter(axA, "A")

    def draw_gene_model(
        ax, tx, insertions, color, gene_name, sample, subtitle, chrom, onc
    ):
        txS, txE = tx["txS"], tx["txE"]
        span = txE - txS
        pad = span * 0.06
        cdsS, cdsE = tx["cdsS"], tx["cdsE"]
        has_cds = cdsS < cdsE
        ax.add_line(Line2D([txS, txE], [0, 0], color="#9A9A9A", lw=1.4, zorder=2))
        for xc in np.linspace(txS + span * 0.03, txE - span * 0.03, 14):
            dx = span * 0.012 * (1 if tx["strand"] == "+" else -1)
            ax.plot(
                [xc - dx, xc, xc - dx],
                [0.16, 0, -0.16],
                color="#C2C2C2",
                lw=0.8,
                zorder=2,
            )
        for s, e in zip(tx["exS"], tx["exE"]):
            if has_cds:
                for us, ue in [(s, min(e, cdsS)), (max(s, cdsE), e)]:
                    if ue > us:
                        ax.add_patch(
                            Rectangle(
                                (us, -0.16),
                                ue - us,
                                0.32,
                                facecolor="#6F6F6F",
                                edgecolor="none",
                                zorder=3,
                            )
                        )
                cs, ce = max(s, cdsS), min(e, cdsE)
                if ce > cs:
                    ax.add_patch(
                        Rectangle(
                            (cs, -0.34),
                            ce - cs,
                            0.68,
                            facecolor="#3F3F3F",
                            edgecolor="none",
                            zorder=4,
                        )
                    )
            else:
                ax.add_patch(
                    Rectangle(
                        (s, -0.22),
                        e - s,
                        0.44,
                        facecolor="#6F6F6F",
                        edgecolor="none",
                        zorder=3,
                    )
                )
        smax = max(i["support"] for i in insertions)
        for i in insertions:
            h = 0.7 + 1.3 * (np.log10(i["support"]) / np.log10(max(smax, 10)))
            ax.add_line(
                Line2D(
                    [i["pos"], i["pos"]],
                    [0.34, 0.34 + h],
                    color=color,
                    lw=2.2,
                    zorder=5,
                )
            )
            ax.scatter(
                [i["pos"]],
                [0.34 + h],
                s=150,
                color=color,
                edgecolor="white",
                linewidth=1.0,
                zorder=6,
                marker="v",
            )
            ax.text(
                i["pos"],
                0.34 + h + 0.12,
                f"{i['pct']:.0f}%",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold",
                color=color,
            )
        ax.set_xlim(txS - pad, txE + pad)
        ax.set_ylim(-0.7, 0.34 + 1.3 * 2 + 0.5)
        ax.set_yticks([])
        ticks = np.linspace(txS, txE, 4)
        ax.set_xticks(ticks)
        precision = max(2, int(np.ceil(-np.log10(max(span / 3e6, 1e-9)))))
        ax.set_xticklabels([f"{t / 1e6:.{precision}f}" for t in ticks], fontsize=7.5)
        ax.set_xlabel(f"{chrom} position (Mb)", fontsize=8)
        ax.spines["left"].set_visible(False)
        arrow = "→" if tx["strand"] == "+" else "←"
        title_col = ONC_RED if onc else color
        ax.set_title(
            f"{sample} · {gene_name} {arrow}" + ("  ★" if onc else ""),
            loc="left",
            color=title_col,
            pad=8,
        )
        ax.text(
            0.0,
            0.93,
            subtitle,
            transform=ax.transAxes,
            fontsize=8.6,
            va="top",
            ha="left",
            color="#555",
        )

    next_panel = 0
    for s in context_samples:
        r, c = 1 + next_panel // 3, next_panel % 3
        next_panel += 1
        ax = fig.add_subplot(gs[r, c])
        dom = per_sample[s]["genes"][0]
        tx = gm.transcript_for(dom["chrom"], dom["pos"])
        ins = dom["loci"]
        flag = " · gene-list flag" if dom["oncogene"] else ""
        subtitle = f"{dom['region']} · {dom['support']:,} mol ({dom['pct']:.0f}% of sample){flag}"
        if tx:
            draw_gene_model(
                ax,
                tx,
                ins,
                PALETTE[s],
                dom["gene"],
                s,
                subtitle,
                dom["chrom"],
                dom["oncogene"],
            )
        else:
            ax.axis("off")
            ax.set_title(f"{s} · {dom['gene']}", loc="left", color=PALETTE[s])
            ax.text(
                0,
                0.7,
                f"No transcript overlaps the leading locus.\n{subtitle}",
                transform=ax.transAxes,
            )
        panel_letter(ax, f"B{next_panel}")

    if SUMMARY_LABEL is not None:
        axP = fig.add_subplot(gs[1 + next_panel // 3, next_panel % 3])
        next_panel += 1
        axP.axis("off")
        pl = per_sample[SUMMARY_LABEL]
        n_genes = len({l["gene"] for l in pl["loci"] if l["region"] != "intergenic"})
        onc_genes = sorted({g["gene"] for g in pl["genes"] if g["oncogene"]})
        top = pl["genes"][0] if pl["genes"] else None
        lead = f"{top['pct']:.1f}% ({top['gene']})" if top else "none"
        axP.set_title(
            f"{SUMMARY_LABEL} · support summary",
            loc="left",
            color=PALETTE[SUMMARY_LABEL],
        )
        txt = (
            f"{len(pl['loci']):,} loci · {n_genes:,} overlapping host genes\n\n"
            f"Leading gene: {lead}\n{pl['total']:,} molecules total\n\n"
            f"Genes carrying a reference-list flag: {len(onc_genes)}"
        )
        axP.text(
            0,
            0.9,
            txt,
            transform=axP.transAxes,
            fontsize=9.2,
            va="top",
            linespacing=1.5,
        )

    axK = fig.add_subplot(gs[1 + next_panel // 3, next_panel % 3])
    axK.axis("off")
    axK.set_title("Cancer-gene annotation summary", loc="left", pad=8)
    lines = []
    n_dom_onc = 0
    for s in context_samples:
        dom = per_sample[s]["genes"][0]
        star = " ★" if dom["oncogene"] else ""
        n_dom_onc += dom["oncogene"]
        lines.append(f"{s} → {dom['gene']} ({dom['pct']:.0f}%){star}")
    callout = (
        "Leading annotated gene:\n"
        + "\n".join(lines)
        + f"\n\n{n_dom_onc} of {len(context_samples)} leading genes carry a flag.\n"
        f"★ = membership in the retained {len(VA.CANCER_GENES)}-gene set.\n"
        "Flags include nearest genes for intergenic loci."
    )
    axK.text(
        0,
        0.94,
        callout,
        transform=axK.transAxes,
        fontsize=9.2,
        va="top",
        linespacing=1.45,
    )

    fig.suptitle(
        "Annotated Lentiviral Insertion Sites & Host-Gene Context",
        x=0.075,
        y=0.965,
        ha="left",
        fontsize=21,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.075,
        0.918,
        "Vector integration assay · HIGH_MAPQ molecules · 100-bp adjacent-junction clustering → "
        "host genes · hg19 RefSeq · retained gene-list flags",
        ha="left",
        fontsize=11,
        color="#666666",
    )

    png = os.path.join(OUTDIR, "insertion_sites_annotated.png")
    pdf = os.path.join(OUTDIR, "insertion_sites_annotated.pdf")
    save_figure(fig, png)
    plt.close(fig)
    print(
        "Wrote:",
        png,
        "\n       ",
        pdf,
        "\n       ",
        os.path.join(OUTDIR, "insertion_genes_annotated.tsv"),
    )
