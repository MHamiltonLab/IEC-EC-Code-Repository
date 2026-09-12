# Bulk RNA expression

[expression_analysis.R](expression_analysis.R) analyzes a gene-by-sample FPKM table. It collapses duplicate gene symbols by the median, applies `log2(FPKM + 1)`, produces expression QC/PCA plots, fits limma models for IEC versus control, and performs Hallmark pathway enrichment. Additional sections compare insertion-associated genes with noninserted samples and within-sample expression reference values.

**Inputs:** `expression_fpkm.tsv` and `expression_metadata.tsv`. The metadata provide sample/cohort labels, insertion-associated genes, and optional tissue-imaging annotations used by downstream comparisons.

**Outputs:** Expression QC, differential-expression and enrichment tables, insertion-associated comparisons, and vector/raster figures under `results/bulk_rna/expression_analysis/`.

The limma analysis operates on supplied FPKM values; it is distinct from the count-based DESeq2 pseudobulk analysis in `single_cell/`. Upstream RNA alignment and quantification are outside this collection. Reference-gene choices and comparison populations remain explicit in the code.

See [input contracts](../docs/input_contracts.md) 