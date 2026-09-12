"""Input interfaces and output settings shared by the vector-assay scripts."""

import csv
import importlib
import os
import re
from functools import lru_cache
from pathlib import Path

INPUT_DIR = Path(
    os.environ.get(
        "VECTOR_INPUT_DIR",
        str(Path(os.environ.get("METHODS_INPUT_DIR", "inputs")) / "vector_assay"),
    )
)
OUTPUT_DIR = Path(os.environ.get("METHODS_OUTPUT_DIR", "results")) / "vector_assay"
VECTOR_FASTA = Path(os.environ.get("VECTOR_FASTA", str(INPUT_DIR / "vector.fasta")))
GENOME_FASTA = Path(os.environ.get("HG19_FASTA", "references/genome.fasta"))
GENE_REFERENCE = Path(os.environ.get("VECTOR_REFGENE", "references/refgene.tsv"))
COLORS = ("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3")


def read_table(path, required):
    """Read a populated TSV with named columns; no study rows are bundled."""
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: populate the empty input interface before running")
    for row in rows:
        if any(not row.get(field, "").strip() for field in required):
            raise ValueError(f"{path}: required fields must be populated")
    return rows


@lru_cache(maxsize=1)
def sample_manifest():
    rows = read_table(INPUT_DIR / "samples.tsv", ("sample_id", "high_mapq_sam"))
    ids = [row["sample_id"] for row in rows]
    labels = [row.get("display_label") or row["sample_id"] for row in rows]
    if len(set(ids)) != len(ids) or len(set(labels)) != len(labels):
        raise ValueError(
            "samples.tsv: sample identifiers and display labels must be unique"
        )
    return {row["sample_id"]: row for row in rows}


def sample_layout():
    labels, palette = {}, {}
    for index, (sample, row) in enumerate(sample_manifest().items()):
        label = row.get("display_label") or sample
        labels[sample] = label
        palette[label] = row.get("color") or COLORS[index % len(COLORS)]
    return labels, palette


def sam_path(sample, kind="HIGH", required=True):
    column = {
        "HIGH": "high_mapq_sam",
        "TRUE_ALL": "true_all_sam",
        "MISPRIMING": "mispriming_sam",
    }[kind]
    value = sample_manifest()[sample].get(column)
    if not value:
        if required:
            raise ValueError(f"samples.tsv: {column} is required for {sample}")
        return None
    path = Path(value)
    if not path.is_absolute():
        path = INPUT_DIR / path
    if not path.is_file():
        raise FileNotFoundError(f"Required SAM input is missing: {path}")
    return path


def clip_paths():
    paths = []
    for sample in sample_manifest():
        paths.append(sam_path(sample))
        extra = sam_path(sample, "MISPRIMING", required=False)
        if extra is not None:
            paths.append(extra)
    return paths


def target_sites():
    labels, _ = sample_layout()
    rows = read_table(
        INPUT_DIR / "junction_targets.tsv", ("sample_id", "chrom", "center", "gene")
    )
    sites = []
    for row in rows:
        center = int(row["center"])
        if center < 1:
            raise ValueError(
                "junction_targets.tsv: center uses positive, 1-based coordinates"
            )
        sites.append(
            (
                labels[row["sample_id"]],
                row["sample_id"],
                row["chrom"],
                center,
                row["gene"],
            )
        )
    return sites


def load_annotation_engine():
    """Load the bundled annotation engine supplied with the methods scripts."""
    return importlib.import_module("vectorint_annotate")


def safe_filename(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip(".") or "sample"


def figure_style():
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(figure, path):
    """Export PDF/SVG with editable text and a PNG preview."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".svg", ".png"):
        figure.savefig(path.with_suffix(suffix), dpi=300, bbox_inches="tight")
