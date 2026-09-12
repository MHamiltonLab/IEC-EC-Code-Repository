# Clinical analyses

| Script | Function |
| --- | --- |
| [cohort_characteristics.R](cohort_characteristics.R) | Baseline summaries; progression-free and overall survival; nonrelapse mortality and enterocolitis cumulative incidence; reverse Kaplan–Meier follow-up. |
| [enterocolitis_outcomes.R](enterocolitis_outcomes.R) | Nutritional support, hospitalization, mortality, and infection comparisons using Fisher exact tests, Wilcoxon tests, and Poisson regression. |

**Inputs:** `clinical_metadata.tsv` and, for the outcome comparison, `clinical_outcomes.tsv`. These represent participant-level tables. Relative event days are measured from infusion; missing event times and censoring times must remain distinguishable.

**Outputs:** Clinical summary tables, model summaries, survival/cumulative-incidence figures, and outcome plots under `results/clinical/`; some statistics are printed to the console.

The scripts address separate questions and have no required execution order. Group definitions, event coding, and censoring rules are explicit in the code. Outcome-stratified survival plots are descriptive; the code does not implement a time-dependent causal effect analysis.

See [input contracts](../docs/input_contracts.md) for expected fields 
