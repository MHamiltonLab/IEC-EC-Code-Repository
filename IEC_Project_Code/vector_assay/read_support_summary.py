#!/usr/bin/env python3
"""Plot supplied integration fractions and rank/count summaries from per-site read-support tables."""

from pathlib import Path

import matplotlib
import pandas as pd
from _assay import INPUT_DIR, OUTPUT_DIR, figure_style, read_table, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if __name__ == "__main__":
    OUTDIR = OUTPUT_DIR / "read_support"
    OUTDIR.mkdir(parents=True, exist_ok=True)
    figure_style()
    rows = read_table(
        INPUT_DIR / "integration_fractions.tsv",
        ("sample_id", "integration_fraction", "site_counts_tsv"),
    )
    integration_fraction = {
        row["sample_id"]: float(row["integration_fraction"]) for row in rows
    }
    if len(integration_fraction) != len(rows):
        raise ValueError("integration_fractions.tsv: sample identifiers must be unique")
    if any(not 0 <= value <= 1 for value in integration_fraction.values()):
        raise ValueError("Integration fractions must lie between 0 and 1")
    site_files = [
        (row["sample_id"], INPUT_DIR / Path(row["site_counts_tsv"])) for row in rows
    ]
    fig, ax = plt.subplots(figsize=(8, 5))

    samples = list(integration_fraction.keys())
    values = list(integration_fraction.values())

    ax.bar(samples, values)

    ax.set_ylabel("Integration Fraction")
    ax.set_title("Supplied integration fractions")

    plt.xticks(rotation=45)

    plt.tight_layout()

    save_figure(fig, f"{OUTDIR}/integration_fraction_barplot.png")

    plt.close()

    fig, ax = plt.subplots(figsize=(8, 6))

    summary = []

    for sample, f in site_files:
        try:
            df = pd.read_csv(
                f, sep="\t", header=None, names=["chr", "start", "end", "count"]
            )
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=["chr", "start", "end", "count"])

        if df.empty:
            summary.append([sample, 0, 0, 0, 0])
            continue
        if (
            df["count"].isna().any()
            or (df["count"] < 0).any()
            or (df["count"] % 1 != 0).any()
        ):
            raise ValueError(f"{f}: support counts must be nonnegative integers")

        counts = sorted(df["count"], reverse=True)

        ax.plot(
            range(1, len(counts) + 1),
            counts,
            marker="o",
            linewidth=1,
            markersize=3,
            label=sample,
        )

        total = df["count"].sum()
        topclone = df["count"].max()
        unique_sites = len(df)

        clonality = topclone / total if total > 0 else 0

        summary.append([sample, total, topclone, unique_sites, clonality])

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Integration Site Rank")
    ax.set_ylabel("Supporting Reads")

    ax.set_title("Integration-site read support")

    ax.legend(fontsize=7)

    plt.tight_layout()

    save_figure(fig, f"{OUTDIR}/clone_abundance_distribution.png")

    plt.close()

    summary_df = pd.DataFrame(
        summary,
        columns=["sample", "total_reads", "top_clone", "unique_sites", "clonality"],
    )

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.bar(summary_df["sample"], summary_df["clonality"])

    ax.set_ylabel("Largest-site read fraction")
    ax.set_title("Largest-site read fraction by sample")

    plt.xticks(rotation=45)

    plt.tight_layout()

    save_figure(fig, f"{OUTDIR}/clonality_barplot.png")

    plt.close()

    summary_df.to_csv(f"{OUTDIR}/integration_summary.tsv", sep="\t", index=False)

    print("Read-support summaries generated.")
