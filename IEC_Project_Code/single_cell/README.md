# Single-cell analyses

| Script | Function |
| --- | --- |
| [car_tissue_repertoire.R](car_tissue_repertoire.R) | CAR/TCR overlays, clonotype sharing, diversity, matched-depth resampling, module scores, and paired pseudobulk expression with Hallmark enrichment. |
| [clonotype_annotation.R](clonotype_annotation.R) | CAR UMI annotation and separate all-clonotyped and CAR-positive repertoire summaries by sample/tissue. |
| [marker_expression.R](marker_expression.R) | Descriptive marker UMAPs and CAR-positive marker-expression distributions. |
| [paired_diversity_tests.R](paired_diversity_tests.R) | Two-sided paired t-tests on patient-level summary values, with statistics, degrees of freedom, mean differences, 95% CIs, and P values. |

**Inputs:** A prepared `single_cell_object.rds`, `single_cell_metadata.tsv`, and, for enrichment, `hallmark_pathways.rds`. Primary QC, integration, and TCR attachment precede these analyses. Configure `CAR_GENE`, `CAR_MIN_UMI`, and clonotype fields near the script headers. The annotation script's `INCLUDE_REGEX` selects sample labels; adapt it to the supplied metadata.

For the comprehensive workflow, set `POSCTRL_SAMPLE` explicitly. Paired subjects are inferred from blood/ileum metadata; `PAIRED_PATIENTS` optionally restricts displays. Pseudobulk DESeq2 uses `~ patient_id + Tissue`. Matched-depth resampling defaults to 1,000 iterations and a 25-cell minimum.

**Outputs:** Tables and vector/raster figures beneath `results/single_cell/<script>/`. The comprehensive script runs the paired-test helper on its per-patient downsampling medians. The helper can also read `paired_diversity_medians.tsv` independently or accept a numeric two-row matrix: rows are tissues, columns are matched patients. Three complete pairs give df = 2. Subsampling iterations are not independent patients, and these three metric tests have no multiplicity correction.

Raw clonotype identifiers may be local to a sample. Pooled tissue summaries require identifiers that represent the same biological clone across samples; otherwise use per-sample summaries. See [input contracts](../docs/input_contracts.md) and [statistical reporting](../docs/statistical_reporting.md).
