#!/usr/bin/env python3
"""Plot coverage and molecule-supported junction evidence at metadata-selected sites."""

from _assay import OUTPUT_DIR, sample_layout, target_sites
from dominant_junctions import generate_junction_reports

if __name__ == "__main__":
    _, palette = sample_layout()
    generate_junction_reports(
        target_sites(),
        palette,
        OUTPUT_DIR / "target_junctions",
        OUTPUT_DIR / "target_junction_figures.tsv",
    )
