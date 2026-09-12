# Immune repertoire

[tcr_diversity.R](tcr_diversity.R) reads MiXCR TRB clone tables, collapses sequences by CDR3 amino-acid identity, and computes top-clone fraction, Shannon entropy, inverse Simpson diversity, Pielou evenness, and clonality defined as one minus Pielou evenness.

**Inputs:** Clone tables in `inputs/mixcr/` and `repertoire_metadata.tsv`. The clone-table filename pattern is configured in the script. Optional `repertoire_pairs.tsv` selects samples for paired display panels; each comparison has a generic identifier, sample list, and title.

**Outputs:** Collapsed clone tables, diversity summaries, Wilcoxon comparisons, and repertoire figures under `results/immune_repertoire/tcr_diversity/`.

The script begins after MiXCR processing; raw-read alignment and clone assembly are upstream steps. Sample identifiers must match between metadata and clone filenames. Pair panels use the supplied sample selections and arrange observations by day.

See [input contracts](../docs/input_contracts.md) for field definitions. 
