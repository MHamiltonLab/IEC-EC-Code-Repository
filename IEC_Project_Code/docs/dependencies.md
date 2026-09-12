# Dependencies

Includes both required packages and packages used only by optional branches. Scripts declare their imports explicitly. 

## R

The code uses contemporary R syntax, including the native pipe in some files. The single-cell workflow expects Seurat-compatible count layers and metadata; verify compatibility with the prepared object. Base R and selected standard namespaces are omitted from the table below.

| Analysis folder | Referenced packages |
| --- | --- |
| `bulk_rna` | `dplyr`, `fgsea`, `ggrepel`, `limma`, `msigdbr`, `pheatmap`, `purrr`, `RColorBrewer`, `readr`, `showtext`, `sysfonts`, `tibble`, `tidyr`, `tidyverse` |
| `car_kinetics` | `broom`, `Cairo`, `car`, `dplyr`, `emmeans`, `forcats`, `ggplot2`, `gridExtra`, `lme4`, `lmerTest`, `mgcv`, `purrr`, `readr`, `scales`, `stringr`, `tableone`, `tibble`, `tidyr` |
| `clinical` | `broom`, `cmprsk`, `dplyr`, `forcats`, `ggplot2`, `gridExtra`, `mgcv`, `patchwork`, `purrr`, `readr`, `scales`, `stringr`, `survival`, `survminer`, `tableone`, `tibble`, `tidyr` |
| `genomics` | `circlize`, `ComplexHeatmap`, `dplyr`, `readr`, `stringr`, `wesanderson` |
| `immune_repertoire` | `dplyr`, `ggplot2`, `RColorBrewer`, `readr`, `showtext`, `stringr`, `sysfonts`, `systemfonts`, `tibble`, `tidyverse` |
| `single_cell` | `DESeq2`, `dplyr`, `fgsea`, `ggplot2`, `ggprism`, `ggrepel`, `Matrix`, `patchwork`, `scales`, `Seurat`, `SeuratObject`, `tibble`, `tidyr` |
| `splicing` | `data.table`, `dplyr`, `GenomeInfoDb`, `GenomicAlignments`, `GenomicRanges`, `ggplot2`, `ggrepel`, `IRanges`, `pheatmap`, `purrr`, `readr`, `Rsamtools`, `rtracklayer`, `S4Vectors`, `showtext`, `stringr`, `tibble`, `tidyr` |
| `tissue_imaging` | `broom.mixed`, `DHARMa`, `dplyr`, `emmeans`, `forcats`, `ggplot2`, `glmmTMB`, `lme4`, `performance`, `RColorBrewer`, `readr`, `scales`, `showtext`, `stringr`, `tibble`, `tidyr`, `tidyverse` |

`tidyverse` imports its standard component packages. Some plotting branches use Cairo, system fonts, or optional text-rendering helpers. `PLOT_FONT_FILE`, where used, accepts an explicitly supplied font file; source code contains no platform-specific font locations. The bulk RNA script obtains Hallmark sets through `msigdbr`, while the single-cell script reads a prepared Hallmark object.

## Command-line tools

| Workflow | Tools and resources |
| --- | --- |
| Vector evidence extraction | Bash, BWA, SAMtools, fastp when trimming is enabled, Java, a GATK 4 jar, and standard Unix text tools. Indexed combined genome/vector reference. |
| Insertion annotation | Python 3 standard library. Precomputed insertion clusters and BEDTools-format overlap/nearest-gene tracks. |
| Somatic variant calling | Bash, BWA-MEM2, SAMtools, BCFtools, GATK 4, and standard Unix text tools. FASTA, known-site, germline-frequency, capture-design, and optional annotation resources. |

## Targeted vector-assay Python

The `vector_assay` scripts use Python 3, NumPy, and Matplotlib; `read_support_summary.py` additionally uses pandas. SAM parsing and indexed-FASTA extraction use the standard library where present in the supplied code. The supplied `vectorint_annotate.py` engine is included and shared by the annotation-dependent analyses. It requires an external extended refGene/genePred table; compatible reference resources and their conventions are described in [vector-assay methods](vector_assay.md). 
