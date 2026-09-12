# Analysis of Archer vector integration assay from output bam files

Targeted vector-assay methods for alignment quality, primer mapping, host-gene annotation, and LTR junction evidence. The complementary WGS workflow is in [vector_integration](../vector_integration/).

| Script | Function |
| --- | --- |
| [alignment_overview.py](alignment_overview.py) | Compare HIGH_MAPQ and TRUE_ALL molecule support, integration sites, and alignment characteristics. |
| [map_vector_softclips.py](map_vector_softclips.py) | Map soft-clipped sequences to the vector and summarize coverage, read starts, and candidate primers. |
| [primer_discovery.py](primer_discovery.py) | Inspect candidate primer windows and draw the vector/LTR primer schematic. |
| [targeted_junctions.py](targeted_junctions.py) | Classify U3/U5 clips at supplied loci and summarize paired junction boundaries. |
| [dominant_junctions.py](dominant_junctions.py) | Select dominant sites and produce coverage and read-pileup figures. |
| [host_gene_annotation.py](host_gene_annotation.py) | Aggregate integration support by host gene and plot genomic locations. |
| [vectorint_annotate.py](vectorint_annotate.py) | Shared SAM, molecule, locus, and RefGene annotation engine; also prints a locus summary. |
| [gene_context_figures.py](gene_context_figures.py) | Plot gene-level support, transcript structure, and retained cancer-gene flags. |
| [target_junction_figures.py](target_junction_figures.py) | Produce standalone coverage/pileup reports for metadata-selected sites. |
| [junction_pileups.py](junction_pileups.py) | Compare selected sites in a panel using read-count junction modes. |
| [junction_evidence_figure.py](junction_evidence_figure.py) | Pair a selected-site pileup with representative LTR sequence matches. |
| [read_support_summary.py](read_support_summary.py) | Plot supplied integration fractions and per-site read-support summaries. |
| [derive_control_primers.py](derive_control_primers.py) | Extract strand-aware control-primer sequences from GTF footprints and indexed FASTA. |
| [junction_sequence_evidence.py](junction_sequence_evidence.py) | Show representative LTR matches within a supplied genomic interval. |

**Inputs:** External SAM files, sample/target metadata, vector FASTA, and compatible genomic references. [input_schemas.py](input_schemas.py) contains empty metadata objects. [primer_reference.tsv](primer_reference.tsv) retains only the supplied reagent sequences and reference positions; no sample observations are included.

**Outputs:** Tables and PNG/PDF/SVG figures under `results/vector_assay/`; mapping and sequence-evidence scripts also report to the console. PDF and SVG exports preserve editable text.

**Dependencies:** Python 3, NumPy, Matplotlib, and pandas for the read-support summary. The supplied **`vectorint_annotate.py`** module is included in this folder.

Run scripts from the repository root. `VECTOR_INPUT_DIR` overrides `inputs/vector_assay/`; `METHODS_OUTPUT_DIR` controls the output root. See [vector-assay methods and interfaces](../docs/vector_assay.md) for reference assumptions, thresholds, and input contracts.
