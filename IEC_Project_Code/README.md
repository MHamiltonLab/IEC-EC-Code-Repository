# Cellular Therapy Analysis Methods

An extended methods code collection for studies of CAR T-cell kinetics, enterocolitis, tissue immune phenotypes, and insertion-associated molecular changes.

The scripts expose the analytical steps, model specifications, and figure-generation logic. They are organized for scientific review and adaptation. Study datasets, patient mappings, reference resources, and prepared analysis objects are supplied separately; this repository does not provide an end-to-end reproduction of manuscript results.

These scripts are not a stand alone data package and are intended to function as an extended methods section detailing analysis of output data.

## Analyses

| Folder | Scope |
| --- | --- |
| [clinical](clinical/) | Cohort characteristics, survival, competing risks, and enterocolitis outcomes |
| [car_kinetics](car_kinetics/) | Early CAR exposure, late persistence, time matching, and neurotoxicity associations |
| [tissue_imaging](tissue_imaging/) | COMET cell abundance, marker phenotypes, normalization sensitivity, and model diagnostics |
| [immune_repertoire](immune_repertoire/) | Bulk TRB clonality and diversity from MiXCR clone tables |
| [single_cell](single_cell/) | CAR-positive single-cell phenotypes, TCR repertoires, and paired pseudobulk expression |
| [bulk_rna](bulk_rna/) | Bulk expression, pathway enrichment, and insertion-associated expression comparisons |
| [splicing](splicing/) | Junction-table, BAM-junction, and transcript/exon-usage analyses |
| [vector_integration](vector_integration/) | Vector–genome read evidence and insertion-cluster annotation |
| [vector_assay](vector_assay/) | Targeted integration-assay alignments, primer maps, host genes, and junction structure |
| [genomics](genomics/) | Somatic variant calling and molecular oncoprints |

## Reading and adapting the code

Each analysis folder has a short README with its scripts, inputs, and principal outputs. [Empty input objects](inputs/schemas.R) describe the table interfaces without supplying observations.

These commands require appropriately populated inputs and the [listed dependencies](docs/dependencies.md). No dependency installation runs automatically.

## Repository scope

The collection contains 32 analysis scripts. Statistical thresholds, model families, and scientific feature names are retained unless this would impact PHI. 


