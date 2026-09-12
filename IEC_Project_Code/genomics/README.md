# Genomic analyses

| Script | Function |
| --- | --- |
| [somatic_variant_calling.sh](somatic_variant_calling.sh) | BWA-MEM2 alignment, duplicate marking, base-quality recalibration, capture QC, tumor-only Mutect2 calling, orientation/contamination filtering, and optional Funcotator annotation. |
| [molecular_oncoprint.R](molecular_oncoprint.R) | Annotated gene-mutation oncoprint and Fisher exact analysis of a separately supplied two-by-two count table. |

**Inputs:** Variant calling accepts `SAMPLE R1 R2`; references, known sites, germline allele frequencies, capture intervals, and optional panel-of-normals/Funcotator resources are configured through environment variables. The oncoprint reads `molecular_features.tsv` and `molecular_comparison.tsv`.

**Outputs:** Alignment/QC/variant files under `OUT_ROOT` (default `results/genomics/somatic_variants`); oncoprint figures under `results/genomics/molecular_oncoprint/` with Fisher statistics printed to the console.

