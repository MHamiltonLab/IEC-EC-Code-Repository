# Vector integration by WGS

| Script | Function |
| --- | --- |
| [extract_vector_junctions.sh](extract_vector_junctions.sh) | Optional fastp trimming, BWA alignment, CAR-related read selection, duplicate marking, and extraction/binning of vector–genome evidence. |
| [annotate_insertions.py](annotate_insertions.py) | Annotate supplied insertion clusters, apply a gene-exclusion rule, summarize read-support diversity, and write Circos tracks. |

**Inputs:** The extraction script accepts `SAMPLE R1 R2`, a combined genome/vector reference (`REF`), and a GATK 4 jar (`GATK_JAR`). Tools resolve from `PATH` or explicit environment settings. Defaults retain MAPQ ≥20, 10-bp bins, and vector contigs named `ciltacel`, `axicel`, or `CAR_VECTOR`.

**Outputs:** Extraction writes sample-level alignment, QC, evidence, and bin tables beneath `OUT_ROOT` (default `results/vector_integration`). Annotation reads/writes the folder supplied by `--outdir`.

The scripts have different input schemas: annotation requires precomputed unique-read insertion clusters and BEDTools annotation tracks. The supplied extraction bins are evidence-row counts, and cannot be passed directly to annotation. The intermediate aggregation/annotation preparation is outside this collection.

See [input contracts](../docs/input_contracts.md) for the exact annotation interface 
