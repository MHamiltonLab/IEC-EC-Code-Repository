# Tissue imaging

COMET-derived cell counts support abundance, phenotype, and tissue-site comparisons.

| Script | Function |
| --- | --- |
| [cd3_normalized_abundance.R](cd3_normalized_abundance.R) | Final COMET abundance/phenotype workflow, denominator audits, compartment comparisons, and editable PDF/EPS figures. |
| [normalization_sensitivity.R](normalization_sensitivity.R) | CD3 versus total-cell normalization, tissue dependence, and repeated-biopsy variability. |
| [model_diagnostics.R](model_diagnostics.R) | Count checks, random-effect comparisons, convergence, and DHARMa residual diagnostics. |

**Inputs:** `tissue_cell_counts.tsv`; optional `patient_labels.tsv` provides display labels.

**Outputs:** Figures and model/summary tables under `results/tissue_imaging/`. The diagnostics script reports to the console and refits its models independently.

The separate sensitivity script explicitly compares CD3 and total-cell normalization. Diagnostics should be aligned with the denominator and model being assessed.
