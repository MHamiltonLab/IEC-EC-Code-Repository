"""Zero-row interfaces for vector-assay metadata. No files are written on import."""

INPUT_SCHEMAS = {
    "samples.tsv": {
        "columns": (
            "sample_id",
            "display_label",
            "high_mapq_sam",
            "true_all_sam",
            "mispriming_sam",
            "color",
        ),
        "rows": [],
    },
    "junction_targets.tsv": {
        "columns": ("sample_id", "chrom", "center", "gene"),
        "rows": [],
    },
}

INPUT_SCHEMAS["integration_fractions.tsv"] = {
    "columns": ("sample_id", "integration_fraction", "site_counts_tsv"),
    "rows": [],
}
