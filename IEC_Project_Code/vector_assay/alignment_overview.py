#!/usr/bin/env python3
"""Summarize HIGH_MAPQ and TRUE_ALL integration alignments and molecule support."""

import os
import re
from collections import Counter, defaultdict

import matplotlib
import numpy as np
from _assay import OUTPUT_DIR, figure_style, sam_path, sample_layout, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = OUTPUT_DIR
JUNCTION_WINDOW = 10
HIGH_MAPQ_CUT = 30

CANON_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_stats(cigar):
    """Return (leading_softclip, trailing_softclip, ref_span, query_len)."""
    parts = _CIG.findall(cigar)
    if not parts:
        return 0, 0, 0, 0
    lead = int(parts[0][0]) if parts[0][1] == "S" else 0
    trail = int(parts[-1][0]) if parts[-1][1] == "S" else 0
    ref_span = sum(int(n) for n, op in parts if op in "MDN=X")
    qlen = sum(int(n) for n, op in parts if op in "MIS=X")
    return lead, trail, ref_span, qlen


def parse_sam(path):
    """Yield dicts of per-read alignment features (mapped reads only)."""
    chrom_len = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    m = dict(
                        tok.split(":", 1)
                        for tok in line.rstrip("\n").split("\t")[1:]
                        if ":" in tok
                    )
                    if "SN" in m and "LN" in m:
                        chrom_len[m["SN"]] = int(m["LN"])
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            flag = int(f[1])
            if flag & 0x4 or f[2] == "*":
                continue
            rname, pos, mapq, cigar = f[2], int(f[3]), int(f[4]), f[5]
            lead, trail, ref_span, qlen = cigar_stats(cigar)
            strand = "-" if (flag & 0x10) else "+"

            if lead >= trail:
                junction = pos
                clip_len = lead
            else:
                junction = pos + ref_span
                clip_len = trail
            molbar = f[0].split("_molbar")[0]
            yield (
                {
                    "chrom": rname,
                    "pos": pos,
                    "mapq": mapq,
                    "strand": strand,
                    "ref_span": ref_span,
                    "clip_len": clip_len,
                    "junction": junction,
                    "molbar": molbar,
                },
                chrom_len,
            )


def dedupe_umi(reads):
    """Keep the highest-MAPQ record for each prefix before ``_molbar``.

    This retains the source assay naming rule; ties retain the first record.
    """
    best = {}
    for r in reads:
        key = r["molbar"]
        if key not in best or r["mapq"] > best[key]["mapq"]:
            best[key] = r
    return list(best.values())


def load_sample(raw):
    """Load HIGH_MAPQ and TRUE_ALL reads (UMI-deduplicated) for one sample."""
    out = {"HIGH": [], "TRUE_ALL": [], "chrom_len": {}}
    for key, _tag in (("HIGH", "_HIGH_MAPQ"), ("TRUE_ALL", "_TRUE_ALL")):
        path = sam_path(raw, key)
        recs = []
        for rec, clen in parse_sam(path):
            recs.append(rec)
            if clen and not out["chrom_len"]:
                out["chrom_len"] = clen
        out[key] = dedupe_umi(recs)
    high_ids = {r["molbar"] for r in out["HIGH"]}
    all_ids = {r["molbar"] for r in out["TRUE_ALL"]}
    if not high_ids <= all_ids:
        raise ValueError(f"{raw}: HIGH_MAPQ molecule IDs must be a subset of TRUE_ALL")
    return out


def call_sites(reads):
    """Cluster adjacent junctions by chromosome, ignoring alignment strand.

    Adjacent junctions at most JUNCTION_WINDOW apart form a chained cluster.
    Support counts unique query-name prefixes; the center is the median.
    """
    by_key = defaultdict(list)
    for r in reads:
        by_key[r["chrom"]].append((r["junction"], r["molbar"]))
    sites = []
    for chrom, items in by_key.items():
        items.sort()
        cluster = [items[0]]
        for j, mb in items[1:]:
            if j - cluster[-1][0] <= JUNCTION_WINDOW:
                cluster.append((j, mb))
            else:
                center = int(np.median([c[0] for c in cluster]))
                sites.append((chrom, ".", center, len({c[1] for c in cluster})))
                cluster = [(j, mb)]
        center = int(np.median([c[0] for c in cluster]))
        sites.append((chrom, ".", center, len({c[1] for c in cluster})))
    return sites


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    SAMPLE_LABEL, PALETTE = sample_layout()
    SAMPLE_ORDER = list(SAMPLE_LABEL.values())
    data = {}
    chrom_len = {}
    for raw, lab in SAMPLE_LABEL.items():
        d = load_sample(raw)
        if not d["TRUE_ALL"] and not d["HIGH"]:
            continue
        d["sites_high"] = call_sites(d["HIGH"])
        d["sites_all"] = call_sites(d["TRUE_ALL"])
        data[lab] = d
        if d["chrom_len"]:
            chrom_len = d["chrom_len"]

    samples = [s for s in SAMPLE_ORDER if s in data]
    if not samples:
        raise ValueError("No mapped integration records were found")

    genome = [(c, chrom_len[c]) for c in CANON_CHROMS if c in chrom_len]
    offset, cum = {}, 0
    for c, L in genome:
        offset[c] = cum
        cum += L
    GENOME_LEN = cum

    summary = []
    for s in samples:
        d = data[s]
        n_high, n_all = len(d["HIGH"]), len(d["TRUE_ALL"])
        sites = d["sites_high"]
        supports = sorted((x[3] for x in sites), reverse=True)
        total = sum(supports) if supports else 0
        top = supports[0] if supports else 0
        chrom_reads = Counter(r["chrom"] for r in d["HIGH"])
        dom_chr, dom_n = chrom_reads.most_common(1)[0] if chrom_reads else ("-", 0)
        summary.append(
            {
                "sample": s,
                "mols_high": n_high,
                "mols_all": n_all,
                "low_mapq": n_all - n_high,
                "sites": len(sites),
                "top_clone_frac": (top / total) if total else 0.0,
                "dom_chr": dom_chr,
                "dom_chr_frac": (dom_n / n_high) if n_high else 0.0,
            }
        )

    with open(os.path.join(OUTDIR, "integration_alignment_summary.tsv"), "w") as fh:
        cols = [
            "sample",
            "mols_high",
            "mols_all",
            "low_mapq",
            "sites",
            "top_clone_frac",
            "dom_chr",
            "dom_chr_frac",
        ]
        fh.write("\t".join(cols) + "\n")
        for row in summary:
            fh.write(
                "\t".join(
                    f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c])
                    for c in cols
                )
                + "\n"
            )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#E6E6E6",
            "grid.linewidth": 0.8,
            "figure.dpi": 120,
        }
    )
    INK = "#2B2B2B"

    fig = plt.figure(figsize=(16, 17))
    gs = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.35, 1.0, 1.0],
        hspace=0.42,
        wspace=0.30,
        left=0.06,
        right=0.975,
        top=0.895,
        bottom=0.055,
    )

    def panel_letter(ax, letter):
        ax.text(
            -0.06,
            1.06,
            letter,
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="bottom",
            ha="right",
            color=INK,
        )

    axA = fig.add_subplot(gs[0, :])
    axA.set_title(
        "Genome-wide lentiviral integration map  (HIGH-MAPQ sites · marker ∝ log molecule support)",
        loc="left",
        pad=12,
    )
    axA.grid(False)

    row_h = 1.0
    for i, s in enumerate(samples):
        y = len(samples) - 1 - i

        for c, L in genome:
            axA.add_patch(
                plt.Rectangle(
                    (offset[c], y - 0.16),
                    L,
                    0.32,
                    facecolor="#F0F0F0",
                    edgecolor="#CFCFCF",
                    lw=0.5,
                    zorder=1,
                )
            )

        sites = data[s]["sites_high"]
        if any(site[0] in offset for site in sites):
            xs, ss = [], []
            for c, strand, j, sup in sites:
                if c in offset:
                    xs.append(offset[c] + j)
                    ss.append(sup)
            ss = np.array(ss, float)
            sizes = 14 + 70 * (
                np.log10(ss) / np.log10(ss.max() if ss.max() > 1 else 10)
            )
            axA.scatter(
                xs,
                np.full(len(xs), y),
                s=sizes,
                color=PALETTE[s],
                alpha=0.78,
                edgecolor="white",
                linewidth=0.4,
                zorder=3,
            )
        axA.text(
            -GENOME_LEN * 0.012,
            y,
            s,
            ha="right",
            va="center",
            fontsize=11,
            fontweight="bold",
            color=PALETTE[s],
        )

    for c, L in genome:
        axA.axvline(offset[c], color="#FFFFFF", lw=0)
        axA.text(
            offset[c] + L / 2,
            len(samples) - 0.35,
            c.replace("chr", ""),
            ha="center",
            va="bottom",
            fontsize=7.5,
            color="#777777",
        )
    for c, L in genome:
        axA.plot(
            [offset[c], offset[c]],
            [-0.5, len(samples) - 0.55],
            color="#DDDDDD",
            lw=0.5,
            zorder=0,
        )
    axA.set_xlim(-GENOME_LEN * 0.05, GENOME_LEN)
    axA.set_ylim(-0.6, len(samples) - 0.1)
    axA.set_yticks([])
    axA.set_xticks([])
    for sp in axA.spines.values():
        sp.set_visible(False)
    panel_letter(axA, "A")

    axB = fig.add_subplot(gs[1, 0])
    x = np.arange(len(samples))
    high = np.array([r["mols_high"] for r in summary], float)
    low = np.array([r["low_mapq"] for r in summary], float)
    axB.bar(
        x, high, color=[PALETTE[s] for s in samples], label="HIGH_MAPQ input", zorder=3
    )
    axB.bar(x, low, bottom=high, color="#BBBBBB", label="TRUE_ALL remainder", zorder=3)
    for xi, (h, l) in enumerate(zip(high, low)):
        axB.text(
            xi,
            h + l,
            f"{int(h + l):,}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=INK,
        )
    axB.set_xticks(x)
    axB.set_xticklabels(samples, rotation=25, ha="right")
    axB.set_ylabel("Integration molecules (UMI)")
    axB.set_title("Molecular support & alignment quality", loc="left")
    axB.legend(fontsize=8, frameon=False, loc="upper right")
    panel_letter(axB, "B")

    axC = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0, 60, 31)
    for s in samples:
        mq = np.array([r["mapq"] for r in data[s]["TRUE_ALL"]])
        if len(mq) == 0:
            continue
        axC.hist(
            mq,
            bins=bins,
            histtype="step",
            linewidth=2.0,
            color=PALETTE[s],
            label=s,
            zorder=3,
        )
    axC.set_yscale("log")
    axC.axvspan(0, HIGH_MAPQ_CUT, color="#F2D5D5", alpha=0.45, zorder=0)
    axC.axvline(HIGH_MAPQ_CUT, color=INK, ls="--", lw=1.2, zorder=4)
    axC.text(
        HIGH_MAPQ_CUT - 1.5,
        axC.get_ylim()[1] * 0.6,
        "← low MAPQ\n   (excluded from HIGH)",
        ha="right",
        va="top",
        fontsize=7.5,
        color="#7A3030",
    )
    axC.set_xlabel("Mapping quality (MAPQ)")
    axC.set_ylabel("Molecule count (log)")
    axC.set_xlim(0, 62)
    axC.set_title("Alignment confidence (TRUE_ALL)", loc="left")
    axC.legend(fontsize=8, frameon=False, loc="upper left")
    panel_letter(axC, "C")

    axD = fig.add_subplot(gs[1, 2])

    chrom_frac = {}
    for s in samples:
        cnt = Counter(r["chrom"] for r in data[s]["HIGH"])
        tot = sum(cnt.values()) or 1
        chrom_frac[s] = {c: cnt[c] / tot for c in cnt}
    top_chroms = sorted(
        {
            c
            for s in samples
            for c, f in chrom_frac[s].items()
            if f >= 0.05 and c in CANON_CHROMS
        },
        key=lambda c: CANON_CHROMS.index(c),
    )
    chr_colors = dict(
        zip(top_chroms, plt.cm.tab20(np.linspace(0, 1, max(len(top_chroms), 1))))
    )
    y = np.arange(len(samples))[::-1]
    for yi, s in zip(y, samples):
        left = 0.0
        for c in top_chroms:
            f = chrom_frac[s].get(c, 0.0)
            if f > 0:
                axD.barh(
                    yi,
                    f,
                    left=left,
                    color=chr_colors[c],
                    edgecolor="white",
                    lw=0.4,
                    zorder=3,
                )
                if f >= 0.08:
                    axD.text(
                        left + f / 2,
                        yi,
                        c.replace("chr", ""),
                        ha="center",
                        va="center",
                        fontsize=7.5,
                        color="white",
                        fontweight="bold",
                    )
                left += f
        other = 1 - left
        if other > 0.001:
            axD.barh(
                yi,
                other,
                left=left,
                color="#DDDDDD",
                edgecolor="white",
                lw=0.4,
                zorder=3,
            )
    axD.set_yticks(y)
    axD.set_yticklabels(samples)
    axD.set_xlim(0, 1)
    axD.set_xlabel("Fraction of integration molecules")
    axD.set_title("Chromosomal integration distribution", loc="left")
    panel_letter(axD, "D")

    axE = fig.add_subplot(gs[2, 0])
    for s in samples:
        sup = sorted((x[3] for x in data[s]["sites_high"]), reverse=True)
        if not sup:
            continue
        axE.plot(
            range(1, len(sup) + 1),
            sup,
            marker="o",
            ms=3.5,
            lw=1.4,
            color=PALETTE[s],
            label=s,
            zorder=3,
        )
    axE.set_xscale("log")
    axE.set_yscale("log")
    axE.set_xlabel("Integration-site rank")
    axE.set_ylabel("Supporting molecules")
    axE.set_title("Integration-site support distribution", loc="left")
    axE.legend(fontsize=8, frameon=False)
    panel_letter(axE, "E")

    axF = fig.add_subplot(gs[2, 1])
    for s in samples:
        mapped = np.array([r["ref_span"] for r in data[s]["HIGH"]])
        clip = np.array([r["clip_len"] for r in data[s]["HIGH"]])
        if len(mapped) == 0:
            continue
        axF.scatter(
            mapped,
            clip,
            s=10,
            color=PALETTE[s],
            alpha=0.35,
            edgecolor="none",
            label=s,
            zorder=3,
        )
    axF.set_xlabel("Genomic mapped length (bp)")
    axF.set_ylabel("Vector soft-clip length (bp)")
    axF.set_title("Chimeric junction architecture", loc="left")
    axF.legend(fontsize=8, frameon=False, markerscale=1.6)
    panel_letter(axF, "F")

    axG = fig.add_subplot(gs[2, 2])
    axG.axis("off")
    axG.set_title("Per-sample summary", loc="left")
    header = ["Sample", "Mols\n(HIGH)", "Sites", "Top clone", "Dom. chr"]
    cell_text = []
    for r in summary:
        cell_text.append(
            [
                r["sample"],
                f"{r['mols_high']:,}",
                str(r["sites"]),
                f"{r['top_clone_frac'] * 100:.0f}%",
                f"{r['dom_chr'].replace('chr', '')} ({r['dom_chr_frac'] * 100:.0f}%)",
            ]
        )
    tbl = axG.table(
        cellText=cell_text, colLabels=header, cellLoc="center", loc="center"
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.7)
    tbl.auto_set_column_width(col=[0])
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if row == 0:
            cell.set_facecolor("#3B3B3B")
            cell.set_text_props(color="white", fontweight="bold")
        else:
            s = samples[row - 1]
            if col == 0:
                cell.set_text_props(color=PALETTE[s], fontweight="bold")
            cell.set_facecolor("#FAFAFA" if row % 2 else "#FFFFFF")
    panel_letter(axG, "G")

    fig.suptitle(
        "Lentiviral Vector Integration & Alignment Landscape",
        x=0.06,
        y=0.975,
        ha="left",
        fontsize=21,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.06,
        0.945,
        "Vector integration assay · HIGH_MAPQ and TRUE_ALL inputs · deduplicated molecule support",
        ha="left",
        fontsize=11.5,
        color="#666666",
    )

    png = os.path.join(OUTDIR, "integration_alignment_landscape.png")
    pdf = os.path.join(OUTDIR, "integration_alignment_landscape.pdf")
    save_figure(fig, png)
    plt.close(fig)

    print("Wrote:")
    print(" ", png)
    print(" ", pdf)
    print(" ", os.path.join(OUTDIR, "integration_alignment_summary.tsv"))
    print()
    print("Summary (UMI-deduplicated molecules):")
    print(
        f"{'sample':10s} {'mols_high':>10s} {'mols_all':>10s} {'sites':>6s} {'top_clone':>10s} {'dom_chr':>10s}"
    )
    for r in summary:
        print(
            f"{r['sample']:10s} {r['mols_high']:>10d} {r['mols_all']:>10d} {r['sites']:>6d} "
            f"{r['top_clone_frac'] * 100:>9.1f}% {r['dom_chr'] + ' ' + format(r['dom_chr_frac'] * 100, '.0f') + '%':>10s}"
        )
