# Input contracts

All examples are interfaces, not study data. Input names below are relative to `METHODS_INPUT_DIR` (default `inputs/`) unless stated otherwise. Column names are case-sensitive except where a script explicitly normalizes them. [schemas.R](../inputs/schemas.R) supplies empty objects for the tabular interfaces.

## Clinical and kinetic tables

| File | Grain and principal fields |
| --- | --- |
| `clinical_metadata.tsv` | One participant per row. `Study_ID`, `IEC_enteritis`, baseline covariates, `PD_day`, `Death_day`, `Last_follow_up_day`, and `IEC_Day`. Neurotoxicity analyses additionally use delayed-neurotoxicity flags and onset day. |
| `clinical_outcomes.tsv` | One participant per row. `ID`, `CAR`, `TPN`, `HospitalDays`, `Mortality`, `PostEndoscopyInfections`. `CAR` identifies the infiltrate comparison group in this analysis. Binary fields use `Y`/`N`. |
| `car_measurements.tsv` | Repeated measurements: `Study_ID`, `Day`, `CAR_abs`. CAR abundance is in cells/µL; days are relative to infusion. The scripts recognize additional day/CAR aliases. |
| `lymphocyte_measurements.tsv` | Repeated measurements: `Study_ID`, `Day`, `ALC`, with ALC in thousands of cells/µL. These exact three column names are required. |

`IEC_enteritis` is coded 1 for the case group. Baseline covariates include age, sex, performance status, prior treatment lines, disease characteristics, and laboratory values. Event-day fields are numeric relative times; absence of an event is distinct from a follow-up time. Neurotoxicity indicators use `Yes`/`No` strings. Inspect each script's harmonization block when mapping source fields.

## Tissue imaging

`tissue_cell_counts.tsv` contains one biopsy/measurement row with `Patient_ID`, `Cohort`, `Tissue_Type`, `total_cells`, `CD3`, `CD68`, `CD3_Camelid`, and marker counts. Marker headers such as `CD3-Camelid-GZMB` and `CD3-Camelid-negative-GZMB` distinguish positive and negative compartments.

The supplied interface mixes counts, percentages on a 0–100 scale, and selected fractions on a 0–1 scale. The schema lists the field names; the code's import block explicitly divides relevant percentage columns by 100. In particular, `pct-CD3-Camelid-negative`, `CD103_fraction_tcells`, and the named fraction fields should not be treated as generic 0–100 percentages without checking their definitions.

`patient_labels.tsv` is optional, with `ID` and `Paper_ID` columns. Without it, plotting uses the input identifier directly. The empty fallback contains no mappings. Labels should be unique per input identifier.

## Bulk TCR and single-cell objects

| Input | Interface |
| --- | --- |
| `mixcr/` | TRB clone tables matching the script's filename pattern. The importer selects read-count, CDR3 amino-acid, and V/D/J fields from accepted MiXCR column names. |
| `repertoire_metadata.tsv` | `Sample_ID`, `Paper_ID`, `Day`, `Type`, `Timepoint`, `Pathology`. `Type` defines tissue; pathology values `IEC` and `Control` define the main GI comparison. |
| `repertoire_pairs.tsv` | Optional long-form selection table: `comparison_id`, `Sample`, `title`. Repeated comparison IDs group selected samples into a panel. |
| `single_cell_object.rds` | A prepared Seurat object with RNA counts, a UMAP or other stored reduction, `sample_id` or `orig.ident`, and at least one supported clonotype column. Default clonotype mode is `raw`, using `tcr_raw_clonotype_id`; other modes are configurable. |
| `single_cell_metadata.tsv` | `sample_id`, `subject_id`, `Tissue`. The script harmonizes blood/ileum labels and maps `subject_id` to internal `patient_id`. |
| `hallmark_pathways.rds` | A named list of gene-symbol vectors, or a table with `gs_name`/`gene_symbol` or `pathway`/`gene`. |

Set `POSCTRL_SAMPLE` to a sample represented in both the metadata and object. `PAIRED_PATIENTS` optionally restricts the display subset inferred from paired blood/ileum metadata. Pseudobulk inclusion is determined separately by the script's paired CAR-positive repertoire data. These selections are intentionally explicit.

## Bulk expression and splicing

| File | Interface |
| --- | --- |
| `expression_fpkm.tsv` | Wide table with `GeneID`, `Gene_Name`, `Gene_Biotype`, and one numeric FPKM column per sample. |
| `expression_metadata.tsv` | `Sample_RNAseq`, `Cohort`, `Tissue`, `Timepoint`, `Paper_ID`, `source_id`, and comma-separated gene symbols in `Insertion`. `Cohort` uses `IEC`, `Control`, and optionally `Dup` to identify excluded duplicate samples. |
| `splicing_metadata.tsv` | `sample`, `patient_id`, `group`, `gene`, `chr`, `start`, `end`; optional `paper_sample_id` provides a display label. Coordinate intervals are insertion/target intervals in the configured build. |
| `insertion_coordinates.tsv` | Optional alternate coordinate source for two splicing scripts: `sample`, `gene`, `chrom`, `start`, `end`, or the configured single-position fallback. |
| `*.junction.txt` | Simplified junction tables: chromosome, intron start (0-based), intron end (1-based), read count, and optional annotation. Accepted header aliases are listed in `read_junction_file()`. |
| `bam/` | Coordinate-sorted, indexed BAMs whose sample names match metadata after standard suffix removal and punctuation/case normalization. |
| `GRCh38.gtf.gz` | Compatible gene/transcript/exon annotation containing the gene-name attributes used by the scripts. |
| `hg19ToHg38.over.chain.gz` | Required only when an applicable script is configured for hg19→hg38 liftover. |

Metadata target intervals use 1-based inclusive coordinates when creating `GRanges`. Junction keys retain the convention described above; do not substitute coordinate formats based only on similar field names. The source analysis referenced GRCh38 release 87 annotation; the generic filename does not make alternative releases interchangeable.

## Genomic features and insertion annotation

`molecular_features.tsv` requires `sample_id`, `Grouping`, and `TCL_Genotype` (a delimited gene list). Optional annotation columns are `Product`, `CAR_Detected`, `CD4_CD8`, and `Location`. `Grouping` uses `IEC-EC` and `PTCL` for the retained oncoprint comparison.

`molecular_comparison.tsv` is a separate count table with `Group`, `Outcome`, and integer `n`. Group levels are `G1`/`G2`; outcome levels are `Success`/`Failure`. Populate those categories from a defined comparison; the code contains no fixed observed counts and does not infer the tested endpoint from the oncoprint.

For `annotate_insertions.py --outdir WORK_DIR`, the working directory must already contain:

- `insertions.clusters.tsv`, with `sample`, `insertion_id`, `chrom`, `cluster_start_1based`, `cluster_end_1based`, `breakpoint_mid_1based`, `support_total_unique_qname`, `support_split_unique_qname`, `support_discordant_unique_qname`, `vec_contig_mode`, and `vec_pos_mode`.
- `tmp/hits.genes.tsv`, `tmp/hits.exons.tsv`, and `tmp/hits.utr.tsv`, in BEDTools intersect format comprising a five-column cluster record followed by a six-column annotation record. Field 4 contains `sample|insertion_id`; field 9 contains the gene symbol.
- `tmp/closest.genes.tsv`, in the corresponding nearest-gene format with a final distance column.

Empty annotation tracks represent no hits; missing tracks raise an error. These prepared files are external interfaces. The repository does not supply the unique-read cluster aggregation or BEDTools commands that construct them from extraction output.

The Bash scripts define their own reference/tool variables at the top of each file. The WGS pipeline uses a genome augmented with vector contigs. The exome pipeline requires consistent hg19/b37 FASTA, known-site, germline-frequency, capture-design, and annotation resources.

## Targeted vector integration assay

The targeted assay uses SAM inputs and generic sample/target objects in [vector_assay/input_schemas.py](../vector_assay/input_schemas.py). Read [vector-assay methods and interfaces](vector_assay.md) for file contracts, genomic/reference conventions, and the external annotation-engine requirements. These interfaces are separate from the WGS insertion-cluster tables above.
