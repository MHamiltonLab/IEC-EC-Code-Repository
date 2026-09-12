# Splicing and transcript usage

| Script | Function |
| --- | --- |
| [junction_table_analysis.R](junction_table_analysis.R) | Junction-table import, gene-level splice summaries, insertion-proximal comparisons, and focused transcript annotation. |
| [bam_junction_analysis.R](bam_junction_analysis.R) | Direct extraction and comparison of splice junctions from indexed BAMs. |
| [transcript_and_exon_usage.R](transcript_and_exon_usage.R) | Transcript-aware junction annotation plus exon-bin and transcript-exon coverage usage. |

**Inputs:** `splicing_metadata.tsv`, a compatible `GRCh38.gtf.gz`, and either `*.junction.txt` reports or indexed BAMs in `inputs/bam/`. Optional insertion-coordinate and hg19→hg38 chain files support the first two scripts. The transcript/exon script expects hg38 coordinates directly.

**Outputs:** Sample-matching audits, gene/junction summaries, reference comparisons, transcript annotations, and coverage figures under `results/splicing/`.

Each script builds its own objects from its inputs. The table-based script expects simplified junction reports rather than raw STAR `SJ.out.tab`. Sample names must match after removing standard file suffixes and normalizing punctuation/case; no study-specific prefix translation is applied.

Gene targets and reference-sample rules are visible near the beginning of each script.
