#!/usr/bin/env python3
"""Annotate insertion clusters and summarize read-support diversity."""

import argparse
import csv
import os
import re
from collections import defaultdict
from math import exp, log


def shannon_and_effective(counts):
    s = sum(counts)
    if s <= 0:
        return (0.0, 0.0)
    H = 0.0
    for x in counts:
        if x <= 0:
            continue
        p = x / s
        H -= p * log(p)
    eff = exp(H)
    return (H, eff)


def load_gene_hits(path):
    # bedtools intersect output: A(5) + B(6) => gene name at col 9 (1-based)
    hits = defaultdict(set)  # nameA -> set(gene)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required annotation track is missing: {path}")
    if os.path.getsize(path) == 0:
        return hits
    with open(path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            nameA = parts[3]
            gene = parts[8]
            if gene and gene != ".":
                hits[nameA].add(gene)
    return hits


def load_closest(path):
    closest = {}  # nameA -> (nearest_gene, dist)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required nearest-gene track is missing: {path}")
    if os.path.getsize(path) == 0:
        return closest
    with open(path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            nameA = parts[3]
            gene = parts[8]
            dist = parts[-1]
            if gene == ".":
                closest[nameA] = ("NA", "NA")
            else:
                closest[nameA] = (gene, dist)
    return closest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--exclude-gene-regex",
        default=r"^LINC00486$",
        help="Exclude clusters overlapping or nearest to matching genes; use (?!) to disable.",
    )
    args = ap.parse_args()

    outdir = args.outdir
    tmp = os.path.join(outdir, "tmp")

    clusters_in = os.path.join(outdir, "insertions.clusters.tsv")
    if not os.path.exists(clusters_in):
        raise SystemExit(f"[FATAL] Missing: {clusters_in}")

    hits_genes = load_gene_hits(os.path.join(tmp, "hits.genes.tsv"))
    hits_exons = load_gene_hits(os.path.join(tmp, "hits.exons.tsv"))
    hits_utr = load_gene_hits(os.path.join(tmp, "hits.utr.tsv"))
    closest = load_closest(os.path.join(tmp, "closest.genes.tsv"))

    excl_re = re.compile(args.exclude_gene_regex)

    all_out = os.path.join(outdir, "insertions.clusters.annotated.all.tsv")
    fil_out = os.path.join(outdir, "insertions.clusters.annotated.filtered.tsv")
    per_fil = os.path.join(outdir, "insertions.per_sample.filtered.tsv")

    circos_dir = os.path.join(outdir, "circos")
    os.makedirs(circos_dir, exist_ok=True)

    # Load clusters
    rows = []
    with open(clusters_in, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            nameA = f"{r['sample']}|{r['insertion_id']}"
            r["_nameA"] = nameA

            genes = sorted(hits_genes.get(nameA, set()))
            exons = hits_exons.get(nameA, set())
            utrs = hits_utr.get(nameA, set())

            overlap_genes = ",".join(genes) if genes else ""
            nearest_gene, dist = closest.get(nameA, ("NA", "NA"))

            # Feature class
            if exons:
                fclass = "EXON"
            elif utrs:
                fclass = "UTR"
            elif genes:
                fclass = "INTRON"
            else:
                fclass = "GENOMIC"

            r["feature_class"] = fclass
            r["overlap_genes"] = overlap_genes
            r["nearest_gene"] = nearest_gene
            r["distance_to_nearest_gene_bp"] = dist
            rows.append(r)

    # Write ALL + FILTERED
    out_fields = [
        "sample",
        "insertion_id",
        "chrom",
        "cluster_start_1based",
        "cluster_end_1based",
        "breakpoint_mid_1based",
        "support_total_unique_qname",
        "support_split_unique_qname",
        "support_discordant_unique_qname",
        "vec_contig_mode",
        "vec_pos_mode",
        "feature_class",
        "overlap_genes",
        "nearest_gene",
        "distance_to_nearest_gene_bp",
    ]

    filtered_rows = []
    with open(all_out, "w") as aout, open(fil_out, "w") as fout:
        wa = csv.DictWriter(
            aout, fieldnames=out_fields, delimiter="\t", lineterminator="\n"
        )
        wf = csv.DictWriter(
            fout, fieldnames=out_fields, delimiter="\t", lineterminator="\n"
        )
        wa.writeheader()
        wf.writeheader()

        for r in rows:
            wa.writerow({k: r.get(k, "") for k in out_fields})

            og = r.get("overlap_genes", "")
            ng = r.get("nearest_gene", "")
            drop = False
            if og:
                # any comma-separated gene matches exclude regex
                for g in og.split(","):
                    if excl_re.search(g):
                        drop = True
                        break
            if (not drop) and ng and ng != "NA" and excl_re.search(ng):
                drop = True

            if not drop:
                wf.writerow({k: r.get(k, "") for k in out_fields})
                filtered_rows.append(r)

    # Per-sample filtered summary
    counts_by_sample = defaultdict(list)
    for r in filtered_rows:
        s = r["sample"]
        counts_by_sample[s].append(int(r["support_total_unique_qname"]))

    # Retain samples whose clusters were all excluded.
    all_samples = sorted({r["sample"] for r in rows})
    with open(per_fil, "w") as out:
        out.write(
            "\t".join(
                [
                    "sample",
                    "n_insertions",
                    "total_support",
                    "top_insertion_support",
                    "top_fraction",
                    "shannon_H",
                    "effective_n",
                ]
            )
            + "\n"
        )
        for s in all_samples:
            counts = counts_by_sample.get(s, [])
            total = sum(counts)
            nins = len(counts)
            top = max(counts) if counts else 0
            frac = (top / total) if total > 0 else 0.0
            H, eff = shannon_and_effective(counts)
            out.write(f"{s}\t{nins}\t{total}\t{top}\t{frac:.4f}\t{H:.4f}\t{eff:.3f}\n")

    # Circos tracks from FILTERED
    # links: vecContig vecStart vecEnd genChrom genStart genEnd value
    # hist:  genChrom genStart genEnd value
    by_sample = defaultdict(list)
    for r in filtered_rows:
        by_sample[r["sample"]].append(r)

    all_links_path = os.path.join(circos_dir, "all_samples.links.txt")
    all_hist_path = os.path.join(circos_dir, "all_samples.hist.txt")

    with open(all_links_path, "w") as allL, open(all_hist_path, "w") as allH:
        for s in all_samples:
            s_links = os.path.join(circos_dir, f"{s}.links.txt")
            s_hist = os.path.join(circos_dir, f"{s}.hist.txt")
            rows_s = by_sample.get(s, [])

            # totals for fraction scaling
            total_support = (
                sum(int(r["support_total_unique_qname"]) for r in rows_s)
                if rows_s
                else 0
            )

            with open(s_links, "w") as L, open(s_hist, "w") as Hf:
                for r in rows_s:
                    chrom = r["chrom"]
                    g0 = int(r["cluster_start_1based"]) - 1
                    g1 = int(r["cluster_end_1based"])
                    sup = int(r["support_total_unique_qname"])
                    frac = (sup / total_support) if total_support > 0 else 0.0

                    vec = r["vec_contig_mode"]
                    vpos = r["vec_pos_mode"]
                    try:
                        vpos = int(vpos)
                        v0 = max(0, vpos - 1)
                        v1 = vpos
                    except (TypeError, ValueError):
                        # Omit tracks when the vector position is unavailable.
                        continue

                    # Histogram height is read support; link weight is within-sample support fraction.
                    Hf.write(f"{chrom}\t{g0}\t{g1}\t{sup}\n")
                    allH.write(f"{chrom}\t{g0}\t{g1}\t{sup}\n")

                    # link value = fraction (for thickness mapping)
                    L.write(f"{vec}\t{v0}\t{v1}\t{chrom}\t{g0}\t{g1}\t{frac:.6f}\n")
                    allL.write(f"{vec}\t{v0}\t{v1}\t{chrom}\t{g0}\t{g1}\t{frac:.6f}\n")

    print(f"[OK] Wrote annotated ALL:      {all_out}")
    print(f"[OK] Wrote annotated FILTERED: {fil_out}")
    print(f"[OK] Wrote per-sample FILTERED:{per_fil}")
    print(f"[OK] Wrote circos tracks:      {circos_dir}")


if __name__ == "__main__":
    main()
