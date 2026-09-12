#!/usr/bin/env python3
"""Inspect candidate primer windows and plot vector coverage and LTR primer geometry."""

import os

import map_vector_softclips as mapper
import matplotlib
import numpy as np
from _assay import OUTPUT_DIR, VECTOR_FASTA, clip_paths, figure_style, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle

if __name__ == "__main__":
    mapper.initialize_vector(VECTOR_FASTA)
    VEC, VL, rc = mapper.VEC, mapper.VL, mapper.rc
    if VL <= 5407:
        raise ValueError(
            "This schematic expects the source vector geometry; review LTR coordinates for your reference"
        )
    cov, startF, startR, _, _ = mapper.coverage_profile(clip_paths())
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    print("=== read-start peaks in 5' LTR (0-181) ===")
    print(
        "+ strand (GSP_U5 reads -> 3'term):",
        sorted(
            [(p, n) for p, n in startF.items() if 60 <= p <= 181], key=lambda x: -x[1]
        )[:8],
    )
    print(
        "- strand (GSP_R  reads -> 5'term):",
        sorted(
            [(p, n) for p, n in startR.items() if 20 <= p <= 120], key=lambda x: -x[1]
        )[:8],
    )

    gsp2_u5_end = 127
    gsp2_r_end = 45

    gsp1_u5 = VEC[gsp2_u5_end - 25 - 22 : gsp2_u5_end - 25]
    gsp1_r = rc(VEC[gsp2_r_end + 1 + 25 : gsp2_r_end + 1 + 25 + 22])
    print(
        f"\nGSP2_U5 (+) read-start {gsp2_u5_end}; GSP1_U5 candidate region just internal (+): 5'-{gsp1_u5}-3'"
    )
    print(
        f"GSP2_R  (-) read-start {gsp2_r_end}; GSP1_R  candidate region just internal (-): 5'-{gsp1_r}-3'"
    )
    secU5 = [(p, n) for p, n in startF.items() if 80 <= p <= 120 and n >= 40]
    secR = [(p, n) for p, n in startR.items() if 60 <= p <= 110 and n >= 40]
    print(
        "secondary + peaks 80-120 (possible GSP1_U5 leakage):",
        sorted(secU5, key=lambda x: -x[1])[:5] or "none",
    )
    print(
        "secondary - peaks 60-110 (possible GSP1_R leakage):",
        sorted(secR, key=lambda x: -x[1])[:5] or "none",
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    C5, C3, CB, PU5, PR = "#7FB3D5", "#F1948A", "#D5DBDB", "#1F6FB2", "#B03A2E"
    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.05, 1.25, 1.15],
        hspace=0.62,
        left=0.06,
        right=0.97,
        top=0.9,
        bottom=0.07,
    )

    axA = fig.add_subplot(gs[0])
    axA.add_patch(Rectangle((0, 0), VL, 0.5, facecolor=CB, edgecolor="none"))
    axA.add_patch(Rectangle((0, 0), 181, 0.5, facecolor=C5, edgecolor="none"))
    axA.add_patch(Rectangle((5407, 0), VL - 5407, 0.5, facecolor=C3, edgecolor="none"))
    axA.text(
        90,
        0.62,
        "5′ LTR",
        ha="center",
        fontsize=8.5,
        color="#1A5276",
        fontweight="bold",
    )
    axA.text(
        5520,
        0.62,
        "3′ LTR",
        ha="center",
        fontsize=8.5,
        color="#922B21",
        fontweight="bold",
    )
    axA.text(
        VL / 2,
        0.25,
        "Vector body",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#555",
    )
    covn = cov / max(cov.max(), 1)
    axA.plot(np.arange(VL), 0.6 + covn * 0.9, color="#34495E", lw=0.7)
    axA.fill_between(np.arange(VL), 0.6, 0.6 + covn * 0.9, color="#34495E", alpha=0.25)
    axA.text(
        VL * 0.5,
        1.62,
        "read LTR-clip coverage",
        ha="center",
        fontsize=8,
        color="#34495E",
    )
    axA.set_xlim(-100, VL + 100)
    axA.set_ylim(0, 1.8)
    axA.set_yticks([])
    axA.set_xticks(range(0, VL + 1, 1000))
    axA.set_xlabel("Vector position (bp; 0-based)", fontsize=9)
    axA.set_title(
        "A. Vector reference and soft-clip coverage",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )

    axB = fig.add_subplot(gs[1])
    axB.set_xlim(-5, 185)
    axB.set_ylim(-1.4, 1.6)
    axB.set_yticks([])
    axB.add_patch(
        Rectangle((0, -0.18), 96, 0.36, facecolor="#AED6F1", edgecolor="#5499C7")
    )
    axB.add_patch(
        Rectangle((96, -0.18), 85, 0.36, facecolor="#A9DFBF", edgecolor="#52BE80")
    )
    axB.text(48, 0, "R", ha="center", va="center", fontweight="bold")
    axB.text(138, 0, "U5", ha="center", va="center", fontweight="bold")

    axB.add_patch(
        FancyArrow(
            71,
            0.55,
            -24,
            0,
            width=0.10,
            head_width=0.26,
            head_length=7,
            length_includes_head=True,
            facecolor=PR,
            edgecolor="none",
        )
    )
    axB.text(59, 0.95, "GSP_R", ha="center", fontsize=9, color=PR, fontweight="bold")
    axB.text(
        59,
        0.78,
        "5′-GGCTTAAGCAGTGGGTTCCCTAGTT-3′",
        ha="center",
        fontsize=6.6,
        color=PR,
        family="DejaVu Sans",
    )

    axB.add_patch(
        FancyArrow(
            102,
            -0.55,
            24,
            0,
            width=0.10,
            head_width=0.26,
            head_length=7,
            length_includes_head=True,
            facecolor=PU5,
            edgecolor="none",
        )
    )
    axB.text(118, -1.0, "GSP_U5", ha="center", fontsize=9, color=PU5, fontweight="bold")
    axB.text(
        118,
        -1.25,
        "5′-GTGTGTGCCCGTCTGTTGTGTGACT-3′",
        ha="center",
        fontsize=6.6,
        color=PU5,
        family="DejaVu Sans",
    )

    axB.add_patch(
        Rectangle(
            (gsp2_r_end + 1 + 25, 0.30),
            22,
            0.12,
            facecolor="none",
            edgecolor=PR,
            ls=":",
            lw=1.0,
        )
    )
    axB.text(
        35, 0.40, "GSP1_R\n(predicted,\ninternal)", ha="center", fontsize=6.6, color=PR
    )
    axB.add_patch(
        Rectangle(
            (80, -0.42), 22, 0.12, facecolor="none", edgecolor=PU5, ls=":", lw=1.0
        )
    )
    axB.text(
        157,
        -0.80,
        "GSP1_U5\n(predicted,\ninternal)",
        ha="center",
        fontsize=6.6,
        color=PU5,
    )
    axB.annotate(
        "← reads out to 5′ host junction",
        xy=(2, 0.55),
        xytext=(2, 1.25),
        fontsize=7.6,
        color=PR,
    )
    axB.annotate(
        "reads out to 3′ host junction →",
        xy=(179, -0.55),
        xytext=(95, 1.25),
        fontsize=7.6,
        color=PU5,
    )
    axB.set_xlabel(
        "position within one LTR (bp); R/U5 identical in both LTRs", fontsize=8.5
    )
    axB.set_title(
        "B. LTR primer map (GSP2 set): two outward-reading primers",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)

    axC = fig.add_subplot(gs[2])
    axC.set_xlim(0, 10)
    axC.set_ylim(-1.1, 1.3)
    axC.axis("off")
    axC.add_patch(
        Rectangle((0, -0.2), 2, 0.4, facecolor="#D7DBDD", edgecolor="#909497")
    )
    axC.text(1, 0, "host genome", ha="center", va="center", fontsize=8)
    axC.add_patch(Rectangle((2, -0.2), 1, 0.4, facecolor=C5, edgecolor="#5499C7"))
    axC.text(2.5, 0, "5′LTR", ha="center", va="center", fontsize=7.5, fontweight="bold")
    axC.add_patch(Rectangle((3, -0.2), 4, 0.4, facecolor=CB, edgecolor="#909497"))
    axC.text(5, 0, "CAR vector body", ha="center", va="center", fontsize=8)
    axC.add_patch(Rectangle((7, -0.2), 1, 0.4, facecolor=C3, edgecolor="#922B21"))
    axC.text(7.5, 0, "3′LTR", ha="center", va="center", fontsize=7.5, fontweight="bold")
    axC.add_patch(
        Rectangle((8, -0.2), 2, 0.4, facecolor="#D7DBDD", edgecolor="#909497")
    )
    axC.text(9, 0, "host genome", ha="center", va="center", fontsize=8)
    axC.add_patch(
        FancyArrow(
            2.5,
            0.45,
            -1.2,
            0,
            width=0.05,
            head_width=0.16,
            head_length=0.25,
            length_includes_head=True,
            facecolor=PR,
            edgecolor="none",
        )
    )
    axC.text(
        2.0,
        0.72,
        "GSP_R → 5′ junction",
        ha="center",
        fontsize=8,
        color=PR,
        fontweight="bold",
    )
    axC.add_patch(
        FancyArrow(
            7.5,
            0.45,
            1.2,
            0,
            width=0.05,
            head_width=0.16,
            head_length=0.25,
            length_includes_head=True,
            facecolor=PU5,
            edgecolor="none",
        )
    )
    axC.text(
        8.0,
        0.72,
        "GSP_U5 → 3′ junction",
        ha="center",
        fontsize=8,
        color=PU5,
        fontweight="bold",
    )
    axC.text(
        5,
        -0.7,
        "Each LTR carries both primers (R/U5 identical); each primer reads outward into flanking host DNA,\n"
        "capturing the 5′ and 3′ integration junctions for sequence-level assessment.",
        ha="center",
        fontsize=8,
        color="#555",
    )
    axC.set_title(
        "C. Assay concept: primers read out of each LTR into the host genome",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )

    fig.suptitle(
        "Vector integration assay: LTR primer map",
        x=0.06,
        ha="left",
        y=0.965,
        fontsize=15,
        fontweight="bold",
        color="#1A1A1A",
    )
    png = os.path.join(OUTPUT_DIR, "vector_primer_map.png")
    save_figure(fig, png)
    plt.close(fig)
    print("\nWrote:", png)
