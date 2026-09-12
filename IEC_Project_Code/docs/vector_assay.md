# Vector-assay methods and interfaces

This folder preserves the analysis stages in the supplied targeted integration-assay scripts. It complements the WGS evidence-extraction workflow. The assays use different inputs and support definitions; their outputs should not be concatenated as equivalent measurements.

## Input objects

The default input directory is `inputs/vector_assay/`. Set `VECTOR_INPUT_DIR` to override it, or use `METHODS_INPUT_DIR` to change the parent input directory. [input_schemas.py](../vector_assay/input_schemas.py) defines empty column/row objects without embedding sample data.

| Input | Fields or resource |
| --- | --- |
| `samples.tsv` | Required `sample_id`, `high_mapq_sam`; optional `display_label`, `true_all_sam`, `mispriming_sam`, `color`. The alignment overview also requires `true_all_sam`. Each row identifies one analysis sample; row order determines display order. |
| `junction_targets.tsv` | `sample_id`, `chrom`, `center`, `gene`. Target sample IDs must occur in `samples.tsv`; `center` is a positive 1-based genomic coordinate. |
| `vector.fasta` | One vector reference sequence. `VECTOR_FASTA` overrides the filename. |
| `control_primers.gtf` | `primer_bind` features with 1-based inclusive coordinates, `+`/`-` strand, and `name`/`gene_id` attributes. |
| Genomic FASTA and `.fai` | Default `references/genome.fasta`; override with `HG19_FASTA` or the control-primer script's second positional argument. |
| Gene annotation | Default `references/refgene.tsv`; override with `VECTOR_REFGENE`. UCSC extended refGene/genePred table with a leading bin column and at least 16 fields; plain text or gzip. |
| `integration_fractions.tsv` | `sample_id`, `integration_fraction`, `site_counts_tsv`. Fractions are supplied values on a 0–1 scale; counts filenames resolve relative to the vector input directory. |
| Per-site read counts | Headerless TSV with four columns: chromosome, start, end, count. The read-support summary uses counts only and does not reinterpret interval coordinates. |

SAM filenames in the manifest resolve relative to the vector input directory unless explicitly absolute. Explicit filenames replace the source's participant-specific prefixes and batch-specific glob patterns. Optional mispriming SAM files are pooled with HIGH_MAPQ files for primer mapping; supply the intended assay outputs to avoid counting duplicated input content twice.

The two-row [primer reference](../vector_assay/primer_reference.tsv) contains GSP_U5 and GSP_R sequences, reported strand, and vector positions. Both supplied sequences are 25 bases long. The malformed free-text footer was removed. The source identifies both as GSP2 primers; this file is a reagent/reference description, not a patient table or independent validation of primer identity. Position values retain the source notation and must be checked against the supplied vector reference before adaptation.

## Alignment overview

The standalone SAM parser skips unmapped alignments and extracts CIGAR reference span, leading/trailing soft-clip length, MAPQ, and strand. It assigns the junction to the larger-clip side, favoring the leading side on ties: `POS` for leading clips, and `POS + reference_span` for trailing clips. These are the source's SAM-based boundary conventions, not an implicit conversion to BED coordinates.

Deduplication keeps the highest-MAPQ record for each query-name prefix before `_molbar`; ties keep the first record. This rule depends on the assay's read naming. It is not a general UMI error-correction algorithm, and records without that token retain their complete query name as the key.

Site calling groups by chromosome, ignoring strand, and greedily links adjacent junction coordinates separated by at most 10 bp. The cluster center is the integer median; support counts distinct query-name prefixes. Chaining can produce a cluster wider than 10 bp. The dominant-site fraction divides its support by summed HIGH_MAPQ site support.

HIGH_MAPQ and TRUE_ALL are supplied file categories. The code does not recreate their upstream filtering. The histogram marks MAPQ 30; the composition panel reports HIGH_MAPQ molecules and the TRUE_ALL count remainder. A new input check requires HIGH_MAPQ molecule identifiers to be a subset of TRUE_ALL before reporting that remainder. These counts and fractions describe molecular support, not independently established cell fractions.

## Vector clips and primer candidates

Soft clips of at least 20 bases are mapped in both orientations by 12-mer offset voting and ungapped sequence comparison. The mapper selects the highest identity-times-overlap score and accepts identity of at least 0.90. It pools per-base coverage and orientation-specific read starts across the configured samples. Vector mapping intervals use 0-based, half-open coordinates; reverse read starts use the interval's last included base.

Candidate reporting inspects the eight largest read-start peaks per orientation, requires at least 100 clips, and retrieves adjacent 25-base sequence windows. The schematic preserves the source reference geometry: the first LTR ends at 181, the second begins at 5407, the R/U5 boundary is 96, and selected GSP2 read-start positions are 127 and 45. The proposed 22-base GSP1 windows and secondary-peak thresholds are heuristic candidate definitions retained from the source. They are not recovered or experimentally confirmed GSP1 reagent sequences. Review these constants when using a different vector reference.

The R-side candidate rectangle now uses the same interval as the candidate sequence calculation. All plotting scripts export PDF/SVG with text objects and PNG previews; platform-specific font loading has been removed.

## Junction evidence

The targeted-site script selects SAM alignment starts within ±1,500 bp of each input center. The dominant-site script requests 100-bp site clustering from the bundled engine, takes its first returned locus, and inspects alignment starts within ±320 bp. The supplied engine sorts loci by decreasing support. `target_junction_figures.py` uses metadata-selected centers with the same ±320-bp collection, molecule-count modes, and shared report function.

U3/U5 classification compares up to the first 30 bases of each clip and its reverse complement against the retained LTR-terminal sequences. Targeted and dominant-site summaries require clips of at least 10 bases and identity of at least 0.85. For each class, the modal boundary maximizes the number of distinct query-name prefixes. Paired separation is the absolute difference between these two modal boundaries. The original 20-bp threshold is retained as a short/wide separation label.

The display text now reports junction separation without asserting that every short interval proves a target-site duplication or excludes a host deletion. That interpretation requires sequence and orientation assessment. The source figure previously printed a no-deletion statement regardless of the measured gap; the measurement algorithm is unchanged. Missing one or both LTR classes is labeled incomplete evidence, without assuming a specific cause.

Coverage and per-site report panels use alignment records; their modal-boundary support uses distinct prefixes. At most 26 evenly spaced records per LTR class are shown. The separate `junction_pileups.py` display uses a ±300-bp window, modes counted from alignment records, and up to 22 records per class. These modes can differ from molecule-supported modes when records are duplicated.

`junction_evidence_figure.py` also uses a ±300-bp window and alignment-record modes, displaying up to 20 records per class and three sequence examples per class. Examples require clips of at least 25 bases; pileup classification has no separate 10-base minimum. `--target-index` selects a 1-based row of the target metadata. This paired display requires both LTR classes.

The console sequence-evidence script displays up to two matches per class at identity ≥0.85, using a user-supplied interval; it also retains the source behavior without the 10-base minimum. Fixed participant examples and prewritten conclusions have been removed from both sequence displays. Exact-match rows illustrate the classification rule without asserting that every accepted partial match is uniquely assignable.

## Host-gene annotation engine

The supplied **[vectorint_annotate.py](../vector_assay/vectorint_annotate.py)** module is included and used by the host-gene and dominant-site analyses. It implements query-name-prefix deduplication and the same chromosome-wise chained clustering described above, with configurable windows and support-sorted output. The separate WGS annotator retains its own input schema.

| Required interface | Expected behavior from the calling code |
| --- | --- |
| `parse_sam(path)` | Iterable of alignment dictionaries, including `chrom` and fields consumed by deduplication. |
| `dedupe_umi(records)` | Deduplicated alignment dictionaries. |
| `call_sites(records, window=100)` | Locus dictionaries with `chrom`, `pos`, and `support`; dominant-site selection assumes descending support. |
| `GeneModel(reference)` | Gene-model object constructed from an extended refGene/genePred table, including a leading bin column. |
| `GeneModel.annotate(chrom, pos)` | Dictionary containing `gene`, `region`, and `oncogene`. |

The source reference filename had a GTF-like extension, but the parser actually expects extended refGene/genePred columns: transcript name, chromosome, strand, transcript/CDS bounds, comma-delimited exon starts/ends, and gene symbol at the positions coded in `GeneModel`. A nine-column GTF will not populate this model; an empty model now raises an explicit input error.

Feature selection ranks CDS exon above UTR/exon above intron, with a 0.5 coding-transcript bonus. The separate representative-transcript function prefers coding transcripts, then the longest span, so the transcript drawn need not be the one selected by feature ranking. The retained gene-list flag is membership of the uppercased gene symbol in the source-curated set. The source describes this set as COSMIC-CGC-based, but supplies no release or derivation record; it is not presented as a current or complete census.

For an overlapping transcript, `dist_tss` is strand-aware distance from its TSS. For an intergenic locus, the same field instead holds signed genomic distance to the nearest transcript boundary, with no maximum-distance cutoff. The flag can therefore refer to a nearest gene without the insertion overlapping that gene. Gene-list membership or proximity alone is not an inference of insertional oncogenesis.

GenePred transcript/exon/CDS intervals use 0-based half-open endpoints. The source directly passes its SAM-derived junction coordinate to these interval comparisons; this convention has been preserved. Boundary-adjacent annotations require explicit coordinate review before adaptation, rather than assuming that renaming the resource resolves a possible one-base offset. The genome display retains hg19 chromosome lengths; alignments and annotation resources must use a compatible build.

Gene-level support sums all annotated loci assigned to a gene. Displayed location, feature, and cancer-gene flag come from that gene's most-supported locus. Percentages divide gene support by total called-locus support within the sample. `gene_context_figures.py` additionally writes locus/gene tables and draws transcript-context panels for samples selected from metadata. Set `VECTOR_SUMMARY_SAMPLE` to a sample ID to display its descriptive support summary instead of a transcript panel. The source's fixed sample classification and prewritten interpretation have been replaced by measured counts and fractions.

## Supplied fractions and read-support summaries

`read_support_summary.py` preserves a separate workflow based on headerless per-site count tables. It ranks site counts, sums support, reports the number of rows/sites, and calculates `clonality = largest site count / total count`, using zero when total count is zero. These are read-support summaries, not the SAM engine's deduplicated molecule counts. Empty count files produce zero summary values.

The integration-fraction bar plot reads an external fraction for each sample. The source supplied observed values without the numerator/denominator calculation; no new calculation has been inferred. Define those quantities when populating `integration_fractions.tsv`, and do not substitute the largest-site fraction or HIGH_MAPQ/TRUE_ALL ratio merely because each is a fraction. The original observed values have been removed.

## Control-primer extraction

The control-primer script reads indexed FASTA intervals without invoking SAMtools. Coordinates are 1-based inclusive, and negative-strand features are reverse-complemented. The source strips `chr` prefixes to match unprefixed reference contigs; supply a compatible FASTA index. Interval bounds and strand are now checked before extraction. The returned sequence is the gene-specific footprint; universal adapters are not encoded in these GTF features.
