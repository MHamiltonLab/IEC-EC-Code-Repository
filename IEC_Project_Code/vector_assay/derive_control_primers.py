#!/usr/bin/env python3
"""Extract gene-specific control primers from 1-based inclusive GTF footprints and an indexed FASTA."""

import os
import sys

from _assay import GENOME_FASTA, INPUT_DIR, OUTPUT_DIR


def resolve_reference(arg=None):
    cands = [arg] if arg else []
    cands += [str(GENOME_FASTA)]
    for c in cands:
        if c and os.path.exists(c):
            return c
    sys.exit(
        "ERROR: reference FASTA not found; pass it as the 2nd argument or set $HG19_FASTA"
    )


def load_fai(fasta):
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        sys.exit(f"ERROR: {fai} missing. Create it with `samtools faidx {fasta}`.")
    idx = {}
    with open(fai) as fh:
        for line in fh:
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            idx[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return idx


def fetch(fasta, idx, contig, start, end):
    """Return uppercase (+)-strand bases for 1-based inclusive [start, end]."""
    if contig not in idx:
        sys.exit(f"ERROR: contig {contig!r} not in reference index")
    length, offset, lb, lw = idx[contig]
    if not 1 <= start <= end <= length:
        raise ValueError(f"Invalid 1-based primer interval: {contig}:{start}-{end}")

    def byte_of(pos):
        p0 = pos - 1
        return offset + (p0 // lb) * lw + (p0 % lb)

    first, last = byte_of(start), byte_of(end)
    with open(fasta, "rb") as fh:
        fh.seek(first)
        raw = fh.read(last - first + 1)
    seq = raw.decode().replace("\n", "").replace("\r", "").upper()
    if len(seq) != end - start + 1:
        sys.exit(
            f"ERROR: extracted {len(seq)} bp, expected {end - start + 1} at {contig}:{start}-{end}"
        )
    return seq


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def parse_attr(attrs, key):

    for field in attrs.split(";"):
        field = field.strip()
        if field.startswith(key + " "):
            return field[len(key) + 1 :].strip().strip('"')
    return ""


def parse_gtf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "primer_bind":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            if strand not in {"+", "-"}:
                raise ValueError("primer_bind features must have + or - strand")
            rows.append(
                {
                    "chrom": chrom,
                    "contig": chrom[3:] if chrom.startswith("chr") else chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "name": parse_attr(attrs, "name"),
                    "gene": parse_attr(attrs, "gene_id"),
                    "function": parse_attr(attrs, "function"),
                }
            )
    return rows


def main():
    gtf = sys.argv[1] if len(sys.argv) > 1 else str(INPUT_DIR / "control_primers.gtf")
    fasta = resolve_reference(sys.argv[2] if len(sys.argv) > 2 else None)
    idx = load_fai(fasta)
    rows = parse_gtf(gtf)
    if not rows:
        sys.exit(f"ERROR: no primer_bind features parsed from {gtf}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = OUTPUT_DIR / "control_primers.tsv"
    header = ["gene", "name", "locus", "strand", "length", "primer_5to3"]
    with open(out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")
        print(f"{'gene':<10}{'locus':<28}{'str':<5}{'len':<5}primer (5'->3')")
        print("-" * 92)
        for r in rows:
            ref = fetch(fasta, idx, r["contig"], r["start"], r["end"])
            primer = ref if r["strand"] == "+" else revcomp(ref)
            locus = f"{r['chrom']}:{r['start']}-{r['end']}"
            out.write(
                "\t".join(
                    [r["gene"], r["name"], locus, r["strand"], str(len(primer)), primer]
                )
                + "\n"
            )
            print(f"{r['gene']:<10}{locus:<28}{r['strand']:<5}{len(primer):<5}{primer}")
    print(f"\nReference: {fasta}")
    print(f"Wrote {out_tsv}")


if __name__ == "__main__":
    main()
