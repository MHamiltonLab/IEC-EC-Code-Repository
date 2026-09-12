# Single-cell analyses

[car_tissue_repertoire.R](car_tissue_repertoire.R) compares CAR-positive cells in blood and ileum using an already processed Seurat object. It provides CAR/TCR UMAP overlays, clonotype sharing, repertoire diversity, matched-depth bootstrapping, transcriptional module scores, and patient-paired pseudobulk differential expression with Hallmark enrichment.

**Inputs:** `single_cell_object.rds`, `single_cell_metadata.tsv`, and `hallmark_pathways.rds`. The Seurat object must already contain count data, a dimensional reduction, sample identifiers, and clonotype annotations. The script does not perform primary QC, integration, or TCR attachment.

Set `POSCTRL_SAMPLE` to the positive-control sample identifier. Paired subjects are inferred from blood/ileum metadata; `PAIRED_PATIENTS` optionally supplies a comma-separated subset. `CAR_GENE` defaults to `CILTACELCAR`, with positivity defined by at least one UMI.

**Outputs:** Figures and tables under `results/single_cell/car_tissue_repertoire/`. Pseudobulk DESeq2 uses `~ patient_id + Tissue`. The default repertoire bootstrap uses 1,000 iterations and a 25-cell minimum.

Module-score labels are descriptive transcriptional summaries. 