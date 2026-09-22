#!/usr/bin/env Rscript
# Single-cell CAR and TCR analysis
# Compare CAR-positive blood and ileum cells, repertoire diversity, module scores, and pseudobulk expression.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(patchwork)
  library(ggrepel)
  library(grid)
})

HAS_GGPRISM <- requireNamespace("ggprism", quietly = TRUE)
HAS_DESEQ2  <- requireNamespace("DESeq2",  quietly = TRUE)
HAS_FGSEA   <- requireNamespace("fgsea",   quietly = TRUE)

# Inputs and analysis parameters
BASE_DIR    <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
QC_DIR      <- Sys.getenv("QC_DIR",   unset=BASE_DIR)

INPUT_RDS   <- Sys.getenv("INPUT_RDS",   unset=file.path(QC_DIR, "single_cell_object.rds"))
SAMPLE_META <- Sys.getenv("SAMPLE_META", unset=file.path(QC_DIR, "single_cell_metadata.tsv"))
OUTDIR      <- Sys.getenv("OUTDIR",      unset=file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "single_cell", "car_tissue_repertoire"))

CAR_GENE    <- Sys.getenv("CAR_GENE", unset="CILTACELCAR")
CAR_MIN_UMI <- as.integer(Sys.getenv("CAR_MIN_UMI", unset="1"))

CLONE_ID_MODE <- Sys.getenv("CLONE_ID_MODE", unset="raw")

BOOT_B      <- as.integer(Sys.getenv("BOOT_B", unset="1000"))
SEED        <- as.integer(Sys.getenv("SEED", unset="1"))

PT_SIZE_ALL   <- as.numeric(Sys.getenv("PT_SIZE_ALL",  unset="0.25"))
PT_SIZE_CAR   <- as.numeric(Sys.getenv("PT_SIZE_CAR",  unset="0.35"))
PT_ALPHA_ALL  <- as.numeric(Sys.getenv("PT_ALPHA_ALL", unset="0.70"))
PT_ALPHA_CAR  <- as.numeric(Sys.getenv("PT_ALPHA_CAR", unset="0.95"))

OUT_EXT     <- Sys.getenv("OUT_EXT", unset="eps")
TOP_HM_N    <- as.integer(Sys.getenv("TOP_HM_N", unset="40"))

HALLMARK_RDS <- Sys.getenv("HALLMARK_RDS", unset=file.path(QC_DIR, "hallmark_pathways.rds"))

CELLTYPE_COL <- Sys.getenv("CELLTYPE_COL", unset="")

POSCTRL_SAMPLE <- Sys.getenv("POSCTRL_SAMPLE", unset = "")
if (!nzchar(POSCTRL_SAMPLE)) {
  stop("Set POSCTRL_SAMPLE to the positive-control sample_id in single_cell_metadata.tsv.")
}

TOP_TCR_N <- as.integer(Sys.getenv("TOP_TCR_N", unset="20"))

HIST_BINS    <- as.integer(Sys.getenv("HIST_BINS",    unset="100"))
HIST_X_TRANS <- tolower(Sys.getenv("HIST_X_TRANS",    unset="log1p"))
ECDF_ALPHA   <- as.numeric(Sys.getenv("ECDF_ALPHA",   unset="0.9"))

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"plots"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"tables"),showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"logs"),  showWarnings=FALSE, recursive=TRUE)

cat("=== START SINGLE-CELL ANALYSIS ===\n")
cat("QC_DIR       :", QC_DIR, "\n")
cat("INPUT_RDS    :", INPUT_RDS, "\n")
cat("SAMPLE_META  :", SAMPLE_META, "\n")
cat("OUTDIR       :", OUTDIR, "\n")
cat("HALLMARK_RDS :", HALLMARK_RDS, "\n")
cat("CAR_GENE     :", CAR_GENE, "  CAR_MIN_UMI:", CAR_MIN_UMI, "\n")
cat("CLONE_ID_MODE:", CLONE_ID_MODE, "\n")
cat("BOOT_B       :", BOOT_B, " SEED:", SEED, "\n")
cat("OUT_EXT      :", OUT_EXT, " (PDF & PNG always written too)\n")
cat("CELLTYPE_COL :", ifelse(nchar(CELLTYPE_COL)>0, CELLTYPE_COL, "<auto>"), "\n")
cat("POSCTRL_SAMPLE:", POSCTRL_SAMPLE, "\n")
cat("TOP_TCR_N    :", TOP_TCR_N, "\n")
cat("HIST_BINS    :", HIST_BINS, " HIST_X_TRANS:", HIST_X_TRANS, " ECDF_ALPHA:", ECDF_ALPHA, "\n")
cat("====================================\n\n")

stopifnot(file.exists(INPUT_RDS))
stopifnot(file.exists(SAMPLE_META))
stopifnot(file.exists(HALLMARK_RDS))

safe_family <- function(fam_pref = "Arial") {
  ps_ok <- tryCatch({
    fams <- names(grDevices::psFonts())
    fam_pref %in% fams
  }, error = function(e) FALSE)
  if (ps_ok) return(fam_pref)
  "Helvetica"
}

FONT_FAMILY <- safe_family(Sys.getenv("PLOT_FONT_FAMILY", unset = "Arial"))

theme_prism_safe <- function(base_size=18, family=NULL) {
  fam <- safe_family(if (is.null(family)) FONT_FAMILY else family)
  if (HAS_GGPRISM) ggprism::theme_prism(base_size = base_size, base_family = fam)
  else theme_classic(base_size = base_size, base_family = fam)
}

theme_umap <- function(base_size = 20, legend_pos = "right", arrows = FALSE, family = NULL) {
  fam <- safe_family(if (is.null(family)) FONT_FAMILY else family)
  th <- theme_classic(base_size = base_size, base_family = fam) +
    theme(
      axis.ticks   = element_blank(),
      axis.text    = element_blank(),
      axis.title   = element_text(hjust = 0),
      legend.position = legend_pos,
      plot.title   = element_text(face = "bold"),
      panel.spacing = unit(0.6, "lines")
    )
  if (arrows) {
    th <- th +
      theme(
        axis.line.x = element_line(linewidth = 0.6, arrow = arrow(length = unit(0.18, "cm"), ends = "last")),
        axis.line.y = element_line(linewidth = 0.6, arrow = arrow(length = unit(0.18, "cm"), ends = "last"))
      )
  } else {
    th <- th + theme(axis.line.x = element_blank(), axis.line.y = element_blank())
  }
  th
}

save_plot_eps <- function(p, filename_base, width, height) {
  out <- filename_base
  if (!grepl("\\.eps$", out, ignore.case=TRUE)) out <- paste0(out, ".eps")
  if (capabilities("cairo")) {
    ggsave(out, plot=p, width=width, height=height, bg="white",
           device=function(...) grDevices::cairo_ps(..., onefile=FALSE, fallback_resolution=600))
  } else {
    ggsave(out, plot=p, width=width, height=height, bg="white",
           device=function(...) grDevices::postscript(..., onefile=FALSE, paper="special", horizontal=FALSE))
  }
  invisible(out)
}

save_plot_all <- function(p, filename_base, width, height, png_dpi=600) {

  if (capabilities("cairo")) {
    ggsave(paste0(filename_base, ".png"), plot=p, width=width, height=height,
           dpi=png_dpi, bg="white",
           device=function(...) grDevices::png(..., type="cairo"))
  } else {
    ggsave(paste0(filename_base, ".png"), plot=p, width=width, height=height,
           dpi=png_dpi, bg="white")
  }

  if (capabilities("cairo")) {
    ggsave(paste0(filename_base, ".pdf"), plot=p, width=width, height=height, bg="white",
           device=grDevices::cairo_pdf)
  } else {
    ggsave(paste0(filename_base, ".pdf"), plot=p, width=width, height=height, bg="white",
           device="pdf")
  }

  save_plot_eps(p, filename_base, width, height)
}

save_plot_any <- function(p, filename_base, width, height, ext="eps", png_dpi=600) {
  save_plot_all(p, filename_base, width=width, height=height, png_dpi=png_dpi)
  invisible(NULL)
}

scale_clonefreq <- function(...) {
  scale_color_gradientn(
    colours = c("#1f77b4", "#00bfc4", "#f6d743"),
    trans = scales::pseudo_log_trans(sigma = 0.02),
    ...
  )
}

cloneType_levels <- c(
  "Hyperexpanded (100 < X <= 500)",
  "Large (20 < X <= 100)",
  "Medium (5 < X <= 20)",
  "Small (1 < X <= 5)",
  "Single (X == 1)"
)
cloneType_cols <- c(
  "Hyperexpanded (100 < X <= 500)" = "#F0E442",
  "Large (20 < X <= 100)"          = "#E69F00",
  "Medium (5 < X <= 20)"           = "#D55E00",
  "Small (1 < X <= 5)"             = "#0072B2",
  "Single (X == 1)"                = "#332288"
)
tissue_cols <- c("Blood"="#E56AA6", "Ileum"="#2AA1B1", "PosCtrl"="#8B5CF6")

axis_guides <- guides(x = ggplot2::guide_axis(cap = TRUE),
                      y = ggplot2::guide_axis(cap = TRUE))

subset_cells_safe <- function(obj, cells, label="(subset)") {
  if (is.numeric(cells)) stop("subset_cells_safe received NUMERIC cells (indices). Convert to barcodes first.")
  cells <- as.character(cells)
  cells2 <- intersect(cells, Cells(obj))
  if (length(cells2) == 0) {
    cat("\n[ERROR] subset_cells_safe", label, "\n")
    cat("  Provided cells:", length(cells), "\n")
    cat("  Matched cells :", length(cells2), "\n")
    cat("  Example provided:", paste(head(cells, 3), collapse=", "), "\n")
    cat("  Example Cells(obj):", paste(head(Cells(obj), 3), collapse=", "), "\n")
    stop("No cells matched Cells(obj).")
  }
  subset(obj, cells=cells2)
}

get_counts_safe <- function(obj, assay="RNA") {
  m <- tryCatch(GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m)) return(m)
  ass <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(ass), error=function(e) character(0))
  if (length(layers) == 0) stop("No counts slot or layers found for assay=", assay)
  ly <- layers[grepl("count", layers, ignore.case=TRUE)]
  if (length(ly) == 0) ly <- layers[1]
  mat <- tryCatch(SeuratObject::LayerData(ass, layer=ly[1]), error=function(e) NULL)
  if (is.null(mat)) stop("Failed to extract counts from assay=", assay, " layers.")
  if (!inherits(mat,"dgCMatrix")) mat <- as(mat,"dgCMatrix")
  mat
}

gene_umi_from_assay <- function(obj, gene, assay="RNA") {
  cells_all <- Cells(obj)

  m <- tryCatch(GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m) && gene %in% rownames(m)) {
    v <- as.numeric(m[gene, , drop=FALSE]); names(v) <- colnames(m)
    out <- setNames(numeric(length(cells_all)), cells_all); out[names(v)] <- v; return(out)
  }

  ass <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(ass), error=function(e) character(0))
  if (length(layers)==0) return(setNames(numeric(length(cells_all)), cells_all))
  pref <- c("counts","counts_raw","raw_counts")
  layers_use <- intersect(pref, layers)
  if (length(layers_use)==0) layers_use <- layers[grepl("count", layers, ignore.case=TRUE)]
  out <- setNames(numeric(length(cells_all)), cells_all)
  for (ly in layers_use) {
    mat <- tryCatch(SeuratObject::LayerData(ass, layer=ly), error=function(e) NULL)
    if (is.null(mat)) next
    if (!inherits(mat,"dgCMatrix")) mat <- as(mat,"dgCMatrix")
    if (!(gene %in% rownames(mat))) next
    v <- as.numeric(mat[gene, , drop=FALSE]); names(v) <- colnames(mat)
    common <- intersect(names(v), names(out))
    if (length(common)>0) out[common] <- out[common] + v[common]
  }
  out
}

calc_diversity <- function(freqs) {
  freqs <- freqs[is.finite(freqs) & freqs > 0]
  if (length(freqs)==0) return(list(n_clones=0L, shannon=NA, inv_simpson=NA, clonality=NA))
  n_clones <- length(freqs)
  shannon <- -sum(freqs * log(freqs))
  inv_simpson <- 1 / sum(freqs^2)
  clonality <- if (n_clones>1) 1 - (shannon / log(n_clones)) else NA
  list(n_clones=as.integer(n_clones), shannon=shannon, inv_simpson=inv_simpson, clonality=clonality)
}

cloneType_bin <- function(clone_size) {
  dplyr::case_when(
    is.na(clone_size) ~ NA_character_,
    clone_size > 100  ~ "Hyperexpanded (100 < X <= 500)",
    clone_size > 20   ~ "Large (20 < X <= 100)",
    clone_size > 5    ~ "Medium (5 < X <= 20)",
    clone_size > 1    ~ "Small (1 < X <= 5)",
    clone_size == 1   ~ "Single (X == 1)",
    TRUE ~ NA_character_
  )
}

paired_stats <- function(x_blood, x_ileum) {
  ok <- is.finite(x_blood) & is.finite(x_ileum)
  x_blood <- x_blood[ok]; x_ileum <- x_ileum[ok]
  n <- length(x_blood)
  if (n < 2) {
    return(tibble(
      n=n,
      wilcox_exact_2s=NA_real_,
      wilcox_asym_2s=NA_real_,
      sign_exact_2s=NA_real_,
      median_delta=NA_real_,
      mean_delta=NA_real_,
      r_effect=NA_real_
    ))
  }
  d <- x_ileum - x_blood
  w_exact <- suppressWarnings(wilcox.test(x_ileum, x_blood, paired=TRUE, alternative="two.sided", exact=TRUE)$p.value)
  w_asym  <- suppressWarnings(wilcox.test(x_ileum, x_blood, paired=TRUE, alternative="two.sided", exact=FALSE)$p.value)
  npos <- sum(d > 0); ntie <- sum(d == 0); nuse <- n - ntie
  sign_p <- if (nuse >= 1) {
    2 * min(pbinom(npos, size=nuse, prob=0.5), 1 - pbinom(npos-1, size=nuse, prob=0.5))
  } else NA_real_
  r_eff <- NA_real_
  z <- tryCatch({
    wt <- suppressWarnings(wilcox.test(x_ileum, x_blood, paired=TRUE, exact=FALSE))
    as.numeric(qnorm(wt$p.value/2, lower.tail=FALSE)) * sign(mean(d))
  }, error=function(e) NA_real_)
  if (is.finite(z)) r_eff <- z / sqrt(n)
  tibble(
    n=n,
    wilcox_exact_2s=w_exact,
    wilcox_asym_2s=w_asym,
    sign_exact_2s=sign_p,
    median_delta=median(d, na.rm=TRUE),
    mean_delta=mean(d, na.rm=TRUE),
    r_effect=r_eff
  )
}

load_hallmark_pathways <- function(rds_path) {
  x <- readRDS(rds_path)
  if (is.list(x) && length(x)>0 && all(vapply(x, is.character, logical(1)))) return(x)
  if (is.data.frame(x) && all(c("gs_name","gene_symbol") %in% colnames(x))) return(split(x$gene_symbol, x$gs_name))
  if (is.data.frame(x) && all(c("pathway","gene") %in% colnames(x))) return(split(x$gene, x$pathway))
  stop("Unrecognized hallmark RDS format in: ", rds_path)
}

detect_celltype_col <- function(md) {
  if (nchar(CELLTYPE_COL) > 0 && CELLTYPE_COL %in% colnames(md)) return(CELLTYPE_COL)
  candidates <- c(
    "predicted.celltype.l2","predicted.celltype.l1","predicted.celltype.l3",
    "predicted.celltype","Azimuth.celltype","azimuth_celltype","celltype","CellType"
  )
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) > 0) return(hit[1])
  fuzzy <- grep("predicted\\.celltype|azimuth|celltype", colnames(md), ignore.case=TRUE, value=TRUE)
  if (length(fuzzy) > 0) return(fuzzy[1])
  NA_character_
}

consensus_val <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing=TRUE))[1]
}

facet_BP <- function(sample_id_vec, tissue_vec, posctrl_sample) {
  out <- ifelse(sample_id_vec == posctrl_sample, "PosCtrl", as.character(tissue_vec))
  out <- ifelse(out %in% c("Blood", "PosCtrl"), out, NA_character_)
  factor(out, levels = c("Blood","PosCtrl"))
}

cat("Loading object...\n")
# Prepared single-cell object
obj <- readRDS(INPUT_RDS)

assay_names <- tryCatch(Assays(obj), error=function(e) character(0))
cat("Cells:", ncol(obj), " Assays:", paste(assay_names, collapse=","), "\n")

if (!"sample_id" %in% colnames(obj@meta.data)) {
  if ("orig.ident" %in% colnames(obj@meta.data)) obj$sample_id <- as.character(obj$orig.ident)
  else stop("No sample_id or orig.ident in meta.data")
}

sample_meta <- read.delim(SAMPLE_META, stringsAsFactors=FALSE) %>%
  rename(patient_id = subject_id) %>%
  mutate(sample_id = as.character(sample_id),
         Tissue_raw = as.character(Tissue)) %>%
  mutate(Tissue = dplyr::case_when(
    Tissue_raw %in% c("Blood","blood","PB","PBMC") ~ "Blood",
    Tissue_raw %in% c("Ileum","ileum","TI","Terminal Ileum","Terminal_ileum") ~ "Ileum",
    TRUE ~ Tissue_raw
  )) %>%
  select(sample_id, patient_id, Tissue) %>%
  distinct(sample_id, .keep_all=TRUE)

cat("\n[DEBUG] Tissue unique in sample_metadata:\n")
print(sort(unique(sample_meta$Tissue)))
cat("\n[DEBUG] sample_metadata Tissue counts:\n")
print(table(sample_meta$Tissue))

keep_samples <- unique(sample_meta$sample_id)
obj <- subset(obj, subset = sample_id %in% keep_samples)

md <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  mutate(sample_id = as.character(sample_id)) %>%
  select(-any_of(c("patient_id", "Tissue"))) %>%
  left_join(sample_meta, by="sample_id") %>%
  tibble::column_to_rownames("cell_barcode")
obj@meta.data <- as.data.frame(md)

stopifnot(all(c("patient_id","Tissue","sample_id") %in% colnames(obj@meta.data)))
obj$patient_id <- as.character(obj@meta.data$patient_id)
obj$Tissue     <- factor(as.character(obj@meta.data$Tissue), levels = c("Blood","Ileum"))

# Paired subjects and positive control
paired_patients <- sample_meta %>%
  group_by(patient_id) %>%
  filter(all(c("Blood", "Ileum") %in% Tissue)) %>%
  pull(patient_id) %>% unique()
selected_patients <- Sys.getenv("PAIRED_PATIENTS", unset = "")
if (nzchar(selected_patients)) {
  paired_patients <- intersect(strsplit(selected_patients, ",", fixed = TRUE)[[1]], paired_patients)
}
pc_patient <- unique(sample_meta$patient_id[sample_meta$sample_id == POSCTRL_SAMPLE])
if (length(pc_patient) != 1L || is.na(pc_patient)) {
  stop("POSCTRL_SAMPLE must resolve to one subject in single_cell_metadata.tsv.")
}

obj$patient_id <- factor(obj$patient_id,
                         levels = unique(c(paired_patients, sort(unique(obj$patient_id)))))

cat("\n[DEBUG] Tissue mapped counts (per cell after join):\n")
print(table(obj$Tissue, useNA="ifany"))
cat("After filtering: cells=", ncol(obj), "\n")

red <- if ("umap" %in% names(obj@reductions)) "umap" else {
  u <- grep("umap", names(obj@reductions), ignore.case=TRUE, value=TRUE)
  if (length(u)>0) u[1] else names(obj@reductions)[1]
}
cat("Using reduction:", red, "\n")

assay_counts <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cat("Using assay for CAR UMI:", assay_counts, "\n")
car_v <- gene_umi_from_assay(obj, gene=CAR_GENE, assay=assay_counts)
obj$CILTACELCAR_umi <- as.numeric(car_v)
obj$CAR_pos <- obj$CILTACELCAR_umi >= CAR_MIN_UMI
cat("CAR+ cells:", sum(obj$CAR_pos, na.rm=TRUE), "/", ncol(obj), "\n")

md0 <- obj@meta.data
pick_clone_col <- function(mode) {
  mode <- tolower(mode)
  if (mode == "ctstrict" && "CTstrict" %in% colnames(md0)) return("CTstrict")
  if (mode == "raw"     && "tcr_raw_clonotype_id" %in% colnames(md0)) return("tcr_raw_clonotype_id")
  if (mode == "cdr3"    && "tcr_cdr3_aa" %in% colnames(md0)) return("tcr_cdr3_aa")
  if (mode == "auto") {
    if ("CTstrict" %in% colnames(md0)) return("CTstrict")
    if ("tcr_raw_clonotype_id" %in% colnames(md0)) return("tcr_raw_clonotype_id")
    if ("tcr_cdr3_aa" %in% colnames(md0)) return("tcr_cdr3_aa")
  }
  NA_character_
}
clone_col <- pick_clone_col(CLONE_ID_MODE)
if (is.na(clone_col)) {
  stop("No clonotype column found for CLONE_ID_MODE=", CLONE_ID_MODE,
       " (need one of: CTstrict, tcr_raw_clonotype_id, tcr_cdr3_aa).")
}
cat("Using clonotype column:", clone_col, "\n")

obj$CTstrict  <- as.character(obj@meta.data[[clone_col]])
obj$has_clone <- !is.na(obj$CTstrict) & obj$CTstrict != "" & obj$CTstrict != "NA"
cat("CAR+ clonotyped cells:", sum(obj$CAR_pos & obj$has_clone, na.rm=TRUE), "\n")

cells_car <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  filter(CAR_pos, has_clone, Tissue %in% c("Blood","Ileum"), !is.na(patient_id)) %>%
  mutate(sample_id = as.character(sample_id))

if (nrow(cells_car) == 0) stop("No CAR+ clonotyped cells after filtering (Blood/Ileum).")

# CAR-positive clonotype abundance
per_clone <- cells_car %>%
  group_by(sample_id, patient_id, Tissue, CTstrict) %>%
  summarise(clone_size = n(), .groups="drop")

denom <- per_clone %>%
  group_by(sample_id) %>%
  summarise(n_car_clonotyped = sum(clone_size), .groups="drop")

per_clone <- per_clone %>%
  left_join(denom, by="sample_id") %>%
  mutate(
    clone_freq = clone_size / n_car_clonotyped,
    cloneType = factor(cloneType_bin(clone_size), levels=cloneType_levels)
  )

write.table(per_clone, file.path(OUTDIR,"tables","per_clone_CARpos.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

obj$clone_size_car <- NA_real_
obj$clone_freq_car <- NA_real_
obj$cloneType_car  <- NA_character_

cells_car2 <- cells_car %>%
  left_join(per_clone %>% select(sample_id, CTstrict, clone_size, clone_freq, cloneType),
            by=c("sample_id","CTstrict"))

obj$clone_size_car[cells_car2$cell_barcode] <- cells_car2$clone_size
obj$clone_freq_car[cells_car2$cell_barcode] <- cells_car2$clone_freq
obj$cloneType_car[cells_car2$cell_barcode]  <- as.character(cells_car2$cloneType)
obj$cloneType_car <- factor(obj$cloneType_car, levels=cloneType_levels)

seq_cols_pref <- c(
  "tcr_cdr3_aa","tcr_cdr3_nt",
  "tcr_trb_cdr3_aa","tcr_trb_cdr3_nt","tcr_trb_v_gene","tcr_trb_j_gene",
  "tcr_tra_cdr3_aa","tcr_tra_cdr3_nt","tcr_tra_v_gene","tcr_tra_j_gene"
)
seq_cols <- intersect(seq_cols_pref, colnames(obj@meta.data))

ct_seq_annot <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  filter(has_clone) %>%
  group_by(CTstrict) %>%
  summarise(
    across(all_of(seq_cols), consensus_val),
    .groups="drop"
  )

top_car_sample <- per_clone %>%
  arrange(sample_id, desc(clone_size)) %>%
  group_by(sample_id) %>%
  mutate(rank_in_sample = row_number()) %>%
  slice_head(n = TOP_TCR_N) %>%
  ungroup() %>%
  left_join(ct_seq_annot, by="CTstrict") %>%
  select(sample_id, patient_id, Tissue, CTstrict, clone_size, clone_freq, rank_in_sample, all_of(seq_cols))

write.table(top_car_sample,
            file.path(OUTDIR,"tables","tcr_top_sequences_CARpos_by_sample.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

blood_patient <- cells_car %>%
  filter(Tissue=="Blood") %>%
  group_by(patient_id, CTstrict) %>%
  summarise(blood_clone_size = n(), .groups="drop")

blood_denom <- blood_patient %>%
  group_by(patient_id) %>%
  summarise(blood_total = sum(blood_clone_size), .groups="drop")

blood_patient <- blood_patient %>%
  left_join(blood_denom, by="patient_id") %>%
  mutate(blood_clone_freq = blood_clone_size / blood_total) %>%
  group_by(patient_id) %>%
  arrange(desc(blood_clone_size), .by_group=TRUE) %>%
  mutate(blood_rank = row_number()) %>%
  ungroup()

ileum_top_vs_blood <- per_clone %>%
  filter(Tissue=="Ileum") %>%
  arrange(sample_id, desc(clone_size)) %>%
  group_by(sample_id) %>%
  mutate(rank_in_ileum = row_number()) %>%
  slice_head(n = TOP_TCR_N) %>%
  ungroup() %>%
  left_join(ct_seq_annot, by="CTstrict") %>%
  left_join(blood_patient, by=c("patient_id","CTstrict")) %>%
  mutate(
    present_in_blood = !is.na(blood_clone_size),
    blood_clone_freq = ifelse(is.na(blood_clone_freq), 0, blood_clone_freq)
  ) %>%
  select(sample_id, patient_id, CTstrict, clone_size, clone_freq, rank_in_ileum,
         present_in_blood, blood_clone_size, blood_clone_freq, blood_rank,
         Tissue, all_of(seq_cols))

write.table(ileum_top_vs_blood,
            file.path(OUTDIR,"tables","ileum_top_vs_blood_match_TOPN.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

match_summary <- ileum_top_vs_blood %>%
  group_by(patient_id, sample_id) %>%
  summarise(
    n_top_ileum = n(),
    n_matched_in_blood = sum(present_in_blood, na.rm=TRUE),
    frac_matched_in_blood = ifelse(n_top_ileum>0, n_matched_in_blood/n_top_ileum, NA_real_),
    median_blood_rank_matched = suppressWarnings(median(blood_rank[present_in_blood], na.rm=TRUE)),
    .groups="drop"
  )
write.table(match_summary,
            file.path(OUTDIR,"tables","ileum_top_vs_blood_match_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

p1 <- DimPlot(obj, reduction=red, group.by="Tissue",
              pt.size=PT_SIZE_ALL, raster=FALSE, alpha=PT_ALPHA_ALL) +
  scale_color_manual(values=tissue_cols, na.value="grey80") +
  labs(title="All cells: Blood vs Ileum", x="UMAP1", y="UMAP2") +
  theme_umap(base_size=24, legend_pos="right", arrows=TRUE) +
  theme(legend.title = element_blank()) +
  axis_guides
save_plot_all(p1, file.path(OUTDIR,"plots","UMAP_all_cells_Tissue"), width=12, height=10)

emb_all <- as.data.frame(Embeddings(obj, reduction = red))
colnames(emb_all)[1:2] <- c("UMAP_1","UMAP_2")
emb_all$cell_barcode <- rownames(emb_all)

emb_all <- emb_all %>%
  left_join(obj@meta.data %>% tibble::rownames_to_column("cell_barcode"), by = "cell_barcode") %>%
  mutate(
    facet_BP = facet_BP(sample_id, Tissue, POSCTRL_SAMPLE),
    is_cilta = as.logical(CAR_pos)
  ) %>%
  filter(!is.na(facet_BP))

p_cilta_BP_onTop <- ggplot() +
  geom_point(
    data = emb_all %>% filter(!is_cilta),
    aes(x=UMAP_1, y=UMAP_2),
    color = "grey85", alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL
  ) +
  geom_point(
    data = emb_all %>% filter(is_cilta),
    aes(x=UMAP_1, y=UMAP_2),
    color = "#6A1B9A", alpha = PT_ALPHA_CAR, size = PT_SIZE_CAR
  ) +
  facet_wrap(~ facet_BP, nrow = 1, scales = "fixed") +
  labs(
    title = "CILTACELCAR+ (RNA) — Blood vs Positive control (positives on top)",
    x="UMAP1", y="UMAP2"
  ) +
  theme_umap(base_size=22, legend_pos="none", arrows=FALSE) +
  theme(strip.text = element_text(size=18, face="bold")) +
  axis_guides
save_plot_all(p_cilta_BP_onTop, file.path(OUTDIR,"plots","UMAP_CILTApos_faceted_Blood_PosCtrl_onTop"),
              width = 14, height = 6)

p_cilta_BP_normal <- ggplot(emb_all, aes(x=UMAP_1, y=UMAP_2, color = is_cilta)) +
  geom_point(alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL) +
  scale_color_manual(values = c(`FALSE`="grey80", `TRUE`="#6A1B9A"), name="CILTA+") +
  facet_wrap(~ facet_BP, nrow = 1, scales = "fixed") +
  labs(
    title = "CILTACELCAR+ (RNA) — Blood vs Positive control (normal order)",
    x="UMAP1", y="UMAP2"
  ) +
  theme_umap(base_size=22, legend_pos="right", arrows=FALSE) +
  theme(strip.text = element_text(size=18, face="bold")) +
  axis_guides
save_plot_all(p_cilta_BP_normal, file.path(OUTDIR,"plots","UMAP_CILTApos_faceted_Blood_PosCtrl_normal"),
              width = 14, height = 6)

obj_car <- subset(obj, subset = CAR_pos & has_clone & Tissue %in% c("Blood","Ileum"))

p2 <- DimPlot(obj_car, reduction=red, group.by="cloneType_car",
              pt.size=PT_SIZE_CAR, raster=FALSE, alpha=PT_ALPHA_CAR) +
  scale_color_manual(values=cloneType_cols, na.value="grey80", name="Clone category") +
  labs(title="CAR+ clonotyped: clone-size category", x="UMAP1", y="UMAP2") +
  theme_umap(base_size=22, arrows=TRUE) +
  axis_guides
save_plot_all(p2, file.path(OUTDIR,"plots","UMAP_CARpos_cloneType"), width=12, height=10)

obj$Tissue3 <- factor(ifelse(obj$sample_id == POSCTRL_SAMPLE, "PosCtrl", as.character(obj$Tissue)),
                      levels=c("Blood","Ileum","PosCtrl"))
obj_car3 <- subset(obj, subset = CAR_pos & has_clone & Tissue3 %in% c("Blood","Ileum","PosCtrl"))

um <- as.data.frame(Embeddings(obj_car3, reduction=red))
colnames(um)[1:2] <- c("UMAP_1","UMAP_2")
um$cell_barcode <- rownames(um)
um <- um %>%
  left_join(obj_car3@meta.data %>% tibble::rownames_to_column("cell_barcode"), by="cell_barcode")

p4 <- ggplot(um, aes(x=UMAP_1, y=UMAP_2, color=clone_freq_car)) +
  geom_point(size=PT_SIZE_CAR, alpha=PT_ALPHA_CAR) +
  scale_clonefreq(name="Clone freq\n(within sample)", na.value="grey80") +
  facet_wrap(~Tissue3, nrow=1) +
  labs(title="CAR+ clonotyped: clone frequency (Blood vs Ileum vs PosCtrl)", x="UMAP1", y="UMAP2") +
  theme_umap(base_size=20, arrows=FALSE) +
  theme(strip.text = element_text(size=18, face="bold")) +
  axis_guides
save_plot_all(p4, file.path(OUTDIR,"plots","UMAP_CARpos_clone_frequency_Tissue3_faceted"), width=18, height=8)

emb_all_s <- as.data.frame(Embeddings(obj, reduction = red))
colnames(emb_all_s)[1:2] <- c("UMAP_1", "UMAP_2")
emb_all_s$cell_barcode <- rownames(emb_all_s)

emb_all_s <- emb_all_s %>%
  dplyr::left_join(
    obj@meta.data %>% tibble::rownames_to_column("cell_barcode"),
    by = "cell_barcode"
  ) %>%
  dplyr::mutate(
    sample_id = as.factor(as.character(sample_id)),
    is_cilta  = as.logical(CAR_pos)
  )

p_cilta_bySample_onTop <- ggplot() +
  geom_point(
    data = dplyr::filter(emb_all_s, !is_cilta),
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey85", alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL
  ) +
  geom_point(
    data = dplyr::filter(emb_all_s,  is_cilta),
    aes(x = UMAP_1, y = UMAP_2),
    color = "#6A1B9A", alpha = PT_ALPHA_CAR, size = PT_SIZE_CAR
  ) +
  facet_wrap(~ sample_id, scales = "free", ncol = 4) +
  labs(title = "CILTACELCAR+ (RNA) — by sample (positives on top)", x = "UMAP1", y = "UMAP2") +
  theme_umap(base_size = 18, legend_pos = "none", arrows=FALSE) +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  axis_guides
save_plot_all(p_cilta_bySample_onTop, file.path(OUTDIR, "plots", "UMAP_CILTApos_bySample_onTop"),
              width = 16, height = 12)

p_cilta_bySample_normal <- ggplot(emb_all_s, aes(x = UMAP_1, y = UMAP_2, color = is_cilta)) +
  geom_point(alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL) +
  scale_color_manual(values = c(`FALSE` = "grey80", `TRUE` = "#6A1B9A"), name="CILTA+") +
  facet_wrap(~ sample_id, scales = "free", ncol = 4) +
  labs(title = "CILTACELCAR+ (RNA) — by sample (normal order)", x = "UMAP1", y = "UMAP2") +
  theme_umap(base_size = 18, legend_pos = "right", arrows=FALSE) +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  axis_guides
save_plot_all(p_cilta_bySample_normal, file.path(OUTDIR, "plots", "UMAP_CILTApos_bySample_normal"),
              width = 16, height = 12)

obj_car_bySample <- subset(obj, subset = CAR_pos & has_clone & Tissue %in% c("Blood","Ileum"))

um_s <- as.data.frame(Embeddings(obj_car_bySample, reduction = red))
colnames(um_s)[1:2] <- c("UMAP_1", "UMAP_2")
um_s$cell_barcode <- rownames(um_s)
um_s <- um_s %>%
  dplyr::left_join(
    obj_car_bySample@meta.data %>% tibble::rownames_to_column("cell_barcode"),
    by = "cell_barcode"
  ) %>%
  dplyr::mutate(sample_id = as.factor(as.character(sample_id)))

p_cf_bySample <- ggplot(um_s, aes(x = UMAP_1, y = UMAP_2, color = clone_freq_car)) +
  geom_point(size = PT_SIZE_CAR, alpha = PT_ALPHA_CAR) +
  scale_clonefreq(name = "Clone freq\n(within sample)", na.value = "grey80") +
  facet_wrap(~ sample_id, scales = "free", ncol = 4) +
  labs(title = "CAR+ clonotyped — clone frequency (by sample)", x = "UMAP1", y = "UMAP2") +
  theme_umap(base_size = 18, arrows = FALSE) +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  axis_guides
save_plot_all(p_cf_bySample, file.path(OUTDIR, "plots", "UMAP_CARpos_clone_frequency_bySample"),
              width = 16, height = 12)

emb_all <- as.data.frame(Embeddings(obj, reduction = red))
colnames(emb_all)[1:2] <- c("UMAP_1","UMAP_2")
emb_all$cell_barcode <- rownames(emb_all)

emb_all <- emb_all %>%
  dplyr::left_join(
    obj@meta.data %>% tibble::rownames_to_column("cell_barcode"),
    by = "cell_barcode"
  ) %>%
  dplyr::mutate(
    patient_id = factor(as.character(patient_id), levels = unique(c(paired_patients, setdiff(as.character(patient_id), paired_patients)))),
    Tissue     = factor(as.character(Tissue),     levels = c("Blood","Ileum")),
    is_cilta   = as.logical(CAR_pos)
  )

emb_pairs <- emb_all %>% dplyr::filter(patient_id %in% paired_patients, Tissue %in% c("Blood","Ileum"))

p_cilta_pairs_onTop <- ggplot() +
  geom_point(
    data = dplyr::filter(emb_pairs, !is_cilta),
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey85", alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL
  ) +
  geom_point(
    data = dplyr::filter(emb_pairs,  is_cilta),
    aes(x = UMAP_1, y = UMAP_2),
    color = "#6A1B9A", alpha = PT_ALPHA_CAR, size = PT_SIZE_CAR
  ) +
  facet_grid(rows = vars(Tissue), cols = vars(patient_id), scales = "free") +
  labs(title = "CILTACELCAR+ (RNA) — paired patients (Blood top, Ileum bottom)",
       x = "UMAP1", y = "UMAP2") +
  theme_umap(base_size = 18, legend_pos = "none", arrows = FALSE) +
  theme(strip.text = element_text(size = 10, face = "bold"))
save_plot_all(p_cilta_pairs_onTop, file.path(OUTDIR,"plots","UMAP_CILTApos_pairs_grid_onTop"),
              width = 14, height = 8)

p_cilta_pairs_normal <- ggplot(emb_pairs, aes(x = UMAP_1, y = UMAP_2, color = is_cilta)) +
  geom_point(alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL) +
  scale_color_manual(values = c(`FALSE`="grey80", `TRUE`="#6A1B9A"), name = "CILTA+") +
  facet_grid(rows = vars(Tissue), cols = vars(patient_id), scales = "free") +
  labs(title = "CILTACELCAR+ (RNA) — paired patients (normal order)",
       x = "UMAP1", y = "UMAP2") +
  theme_umap(base_size = 18, legend_pos = "right", arrows = FALSE) +
  theme(strip.text = element_text(size = 10, face = "bold"))
save_plot_all(p_cilta_pairs_normal, file.path(OUTDIR,"plots","UMAP_CILTApos_pairs_grid_normal"),
              width = 14, height = 8)

emb_pc <- emb_all %>% dplyr::filter(as.character(patient_id) == pc_patient)
if (nrow(emb_pc) > 0) {
  p_cilta_pc_onTop <- ggplot() +
    geom_point(
      data = dplyr::filter(emb_pc, !is_cilta),
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey85", alpha = PT_ALPHA_ALL, size = PT_SIZE_ALL
    ) +
    geom_point(
      data = dplyr::filter(emb_pc,  is_cilta),
      aes(x = UMAP_1, y = UMAP_2),
      color = "#6A1B9A", alpha = PT_ALPHA_CAR, size = PT_SIZE_CAR
    ) +
    coord_equal() +
    labs(title = paste0("CILTACELCAR+ (RNA) — Positive control (", pc_patient, ")"),
         x = "UMAP1", y = "UMAP2") +
    theme_umap(base_size = 18, legend_pos = "none", arrows = TRUE)
  save_plot_all(p_cilta_pc_onTop, file.path(OUTDIR,"plots","UMAP_CILTApos_PosCtrl_only"),
                width = 8, height = 6)
}

obj_car_pairs <- subset(obj, subset = CAR_pos & has_clone & patient_id %in% paired_patients & Tissue %in% c("Blood","Ileum"))
if (ncol(obj_car_pairs) > 0) {
  um_pairs <- as.data.frame(Embeddings(obj_car_pairs, reduction = red))
  colnames(um_pairs)[1:2] <- c("UMAP_1","UMAP_2")
  um_pairs$cell_barcode <- rownames(um_pairs)
  um_pairs <- um_pairs %>%
    dplyr::left_join(obj_car_pairs@meta.data %>% tibble::rownames_to_column("cell_barcode"), by="cell_barcode") %>%
    dplyr::mutate(
      patient_id = factor(as.character(patient_id), levels = paired_patients),
      Tissue     = factor(as.character(Tissue), levels = c("Blood","Ileum"))
    )

  p_cf_pairs <- ggplot(um_pairs, aes(x = UMAP_1, y = UMAP_2, color = clone_freq_car)) +
    geom_point(size = PT_SIZE_CAR, alpha = PT_ALPHA_CAR) +
    scale_clonefreq(name = "Clone freq\n(within sample)", na.value = "grey80") +
    facet_grid(rows = vars(Tissue), cols = vars(patient_id), scales = "free") +
    labs(title = "CAR+ clonotyped — clone frequency (paired; Blood top)",
         x = "UMAP1", y = "UMAP2") +
    theme_umap(base_size = 18, arrows = FALSE) +
    theme(strip.text = element_text(size = 10, face = "bold"))
  save_plot_all(p_cf_pairs, file.path(OUTDIR,"plots","UMAP_CARpos_clone_frequency_pairs_grid"),
                width = 14, height = 8)
}

obj_car_pc <- subset(obj, subset = CAR_pos & has_clone & patient_id == pc_patient)
if (ncol(obj_car_pc) > 0) {
  um_pc <- as.data.frame(Embeddings(obj_car_pc, reduction = red))
  colnames(um_pc)[1:2] <- c("UMAP_1","UMAP_2")
  um_pc$cell_barcode <- rownames(um_pc)
  um_pc <- um_pc %>%
    dplyr::left_join(obj_car_pc@meta.data %>% tibble::rownames_to_column("cell_barcode"), by="cell_barcode")

  p_cf_pc <- ggplot(um_pc, aes(x = UMAP_1, y = UMAP_2, color = clone_freq_car)) +
    geom_point(size = PT_SIZE_CAR, alpha = PT_ALPHA_CAR) +
    scale_clonefreq(name = "Clone freq\n(PC) ", na.value = "grey80") +
    coord_equal() +
    labs(title = paste0("CAR+ clonotyped — clone frequency (PC ", pc_patient, ")"),
         x = "UMAP1", y = "UMAP2") +
    theme_umap(base_size = 18, arrows = TRUE)
  save_plot_all(p_cf_pc, file.path(OUTDIR,"plots","UMAP_CARpos_clone_frequency_PosCtrl_only"),
                width = 8, height = 6)
}

summarize_umi <- function(val, grp, thr) {
  tibble(value = as.numeric(val), group = as.character(grp)) %>%
    group_by(group) %>%
    summarise(
      n              = dplyr::n(),
      n_nonzero      = sum(value > 0, na.rm=TRUE),
      frac_nonzero   = ifelse(n > 0, n_nonzero / n, NA_real_),
      n_pos_thresh   = sum(value >= thr, na.rm=TRUE),
      frac_pos_thresh= ifelse(n > 0, n_pos_thresh / n, NA_real_),
      mean           = mean(value, na.rm=TRUE),
      sd             = sd(value, na.rm=TRUE),
      median         = median(value, na.rm=TRUE),
      mad            = mad(value, na.rm=TRUE),
      p95            = as.numeric(quantile(value, 0.95, na.rm=TRUE, names=FALSE, type=7)),
      p99            = as.numeric(quantile(value, 0.99, na.rm=TRUE, names=FALSE, type=7)),
      p999           = as.numeric(quantile(value, 0.999, na.rm=TRUE, names=FALSE, type=7)),
      max            = suppressWarnings(max(value, na.rm=TRUE)),
      .groups="drop"
    )
}

md_diag <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  mutate(
    Tissue3 = factor(ifelse(sample_id == POSCTRL_SAMPLE, "PosCtrl", as.character(Tissue)),
                     levels=c("Blood","Ileum","PosCtrl")),
    CILTA_UMI = as.numeric(CILTACELCAR_umi),
    CAR_pos   = as.logical(CAR_pos)
  )

diag_tissue <- summarize_umi(md_diag$CILTA_UMI, md_diag$Tissue3, thr=CAR_MIN_UMI) %>%
  rename(Tissue3 = group)
write.table(diag_tissue, file.path(OUTDIR,"tables","ciltaRNA_diagnostic_by_Tissue3.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

diag_sample <- summarize_umi(md_diag$CILTA_UMI, md_diag$sample_id, thr=CAR_MIN_UMI) %>%
  rename(sample_id = group)
write.table(diag_sample, file.path(OUTDIR,"tables","ciltaRNA_diagnostic_by_sample.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

p_hist_tissue <- ggplot(md_diag, aes(x=CILTA_UMI)) +
  geom_histogram(bins=HIST_BINS, fill="#4F6A92", alpha=0.9, color="white", linewidth=0.1) +
  geom_vline(xintercept=CAR_MIN_UMI, linetype="dashed", color="#E31A1C", linewidth=0.8) +
  facet_wrap(~Tissue3, scales="free_y", nrow=1) +
  labs(title=paste0("CILTACELCAR RNA UMI by Tissue3 (threshold=", CAR_MIN_UMI, ")"),
       x="UMI", y="Cells") +
  theme_bw(base_size=12) +
  theme(strip.text = element_text(size=12, face="bold"))
if (HIST_X_TRANS == "log1p") {
  p_hist_tissue <- p_hist_tissue + scale_x_continuous(trans=scales::log1p_trans())
}
save_plot_any(p_hist_tissue, file.path(OUTDIR,"plots","HIST_CILTA_RNA_by_Tissue3"), width=14, height=5, ext=OUT_EXT)

p_hist_sample <- ggplot(md_diag, aes(x=CILTA_UMI)) +
  geom_histogram(bins=HIST_BINS, fill="#4F6A92", alpha=0.9, color="white", linewidth=0.1) +
  geom_vline(xintercept=CAR_MIN_UMI, linetype="dashed", color="#E31A1C", linewidth=0.8) +
  facet_wrap(~sample_id, scales="free_y") +
  labs(title=paste0("CILTACELCAR RNA UMI by sample (threshold=", CAR_MIN_UMI, ")"),
       x="UMI", y="Cells") +
  theme_bw(base_size=11) +
  theme(strip.text = element_text(size=9, face="bold"))
if (HIST_X_TRANS == "log1p") {
  p_hist_sample <- p_hist_sample + scale_x_continuous(trans=scales::log1p_trans())
}
save_plot_any(p_hist_sample, file.path(OUTDIR,"plots","HIST_CILTA_RNA_by_sample"), width=16, height=10, ext=OUT_EXT)

p_ecdf <- ggplot(md_diag, aes(x=CILTA_UMI, color=Tissue3)) +
  stat_ecdf(geom="step", alpha=ECDF_ALPHA) +
  geom_vline(xintercept=CAR_MIN_UMI, linetype="dashed", color="#E31A1C", linewidth=0.8) +
  labs(title="ECDF of CILTACELCAR RNA UMI by Tissue3", x="UMI", y="ECDF") +
  theme_bw(base_size=12) +
  theme(legend.position="right")
if (HIST_X_TRANS == "log1p") {
  p_ecdf <- p_ecdf + scale_x_continuous(trans=scales::log1p_trans())
}
save_plot_any(p_ecdf, file.path(OUTDIR,"plots","ECDF_CILTA_RNA_by_Tissue3"), width=10, height=7, ext=OUT_EXT)

car_counts_tissue <- md_diag %>% group_by(Tissue3) %>% summarise(n_pos=sum(CAR_pos, na.rm=TRUE), .groups="drop")
p_bar_tissue <- ggplot(car_counts_tissue, aes(x=Tissue3, y=n_pos, fill=Tissue3)) +
  geom_col(width=0.65) +
  geom_text(aes(label=n_pos), vjust=-0.2, size=4) +
  scale_fill_manual(values=tissue_cols) +
  labs(title=paste0("# CAR+ cells (", CAR_GENE, " ≥ ", CAR_MIN_UMI, ") by Tissue3"), x=NULL, y="# CAR+") +
  theme_bw(base_size=12) + theme(legend.position="none")
save_plot_any(p_bar_tissue, file.path(OUTDIR,"plots","BAR_CARpos_by_Tissue3"), width=8, height=6, ext=OUT_EXT)

car_counts_sample <- md_diag %>% group_by(sample_id) %>% summarise(n_pos=sum(CAR_pos, na.rm=TRUE), .groups="drop")
p_bar_sample <- ggplot(car_counts_sample, aes(x=sample_id, y=n_pos)) +
  geom_col(fill="#6A1B9A", width=0.65) +
  geom_text(aes(label=n_pos), vjust=-0.2, size=3.2) +
  labs(title=paste0("# CAR+ cells (", CAR_GENE, " ≥ ", CAR_MIN_UMI, ") by sample"), x="sample_id", y="# CAR+") +
  theme_bw(base_size=10) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
save_plot_any(p_bar_sample, file.path(OUTDIR,"plots","BAR_CARpos_by_sample"), width=14, height=7, ext=OUT_EXT)

cat("\n[DIAG] Wrote CILTACELCAR RNA diagnostics tables & plots.\n")

ct_col <- detect_celltype_col(obj_car@meta.data)
if (!is.na(ct_col) && ct_col %in% colnames(obj_car@meta.data)) {
  cat("Detected celltype column:", ct_col, "\n")
  obj_car$celltype_detected <- as.character(obj_car@meta.data[[ct_col]])

  ct_comp <- obj_car@meta.data %>%
    tibble::rownames_to_column("cell_barcode") %>%
    filter(!is.na(patient_id), Tissue %in% c("Blood","Ileum")) %>%
    group_by(patient_id, Tissue, celltype_detected) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(patient_id, Tissue) %>%
    mutate(frac = n/sum(n)) %>%
    ungroup()

  write.table(ct_comp, file.path(OUTDIR,"tables","celltype_fractions.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  p_ct <- ggplot(ct_comp, aes(x=Tissue, y=frac)) +
    geom_boxplot(outlier.shape=NA, alpha=0.25, width=0.55) +
    geom_point(aes(color=patient_id), size=3.0, position=position_jitter(width=0.06), alpha=0.95) +
    geom_line(aes(group=patient_id), color="grey55", linewidth=0.5, alpha=0.8) +
    facet_wrap(~celltype_detected, scales="free_y") +
    scale_y_continuous(labels=percent_format(accuracy=1)) +
    labs(title=paste0("Celltype fractions (", ct_col, "): Blood vs Ileum"), x="", y="Fraction") +
    theme_prism_safe(base_size=14) +
    theme(legend.position="none")
  save_plot_eps(p_ct, file.path(OUTDIR,"plots","celltype_fractions_faceted"), width=18, height=11)

  ct_bp <- obj_car@meta.data %>%
    tibble::rownames_to_column("cell_barcode") %>%
    filter(!is.na(patient_id)) %>%
    mutate(facet_BP = facet_BP(sample_id, Tissue, POSCTRL_SAMPLE)) %>%
    filter(!is.na(facet_BP)) %>%
    group_by(patient_id, facet_BP, celltype_detected) %>%
    summarise(n = n(), .groups="drop") %>%
    group_by(patient_id, facet_BP) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

  p_ct_bp <- ggplot(ct_bp, aes(x=facet_BP, y=frac)) +
    geom_boxplot(outlier.shape=NA, alpha=0.25, width=0.55) +
    geom_point(aes(color=patient_id),
               size=3.0, position=position_jitter(width=0.06), alpha=0.95) +
    geom_line(aes(group=patient_id), color="grey55", linewidth=0.5, alpha=0.8) +
    facet_wrap(~celltype_detected, scales="free_y") +
    scale_y_continuous(labels=percent_format(accuracy=1)) +
    labs(title=paste0("Celltype fractions (", ct_col, "): Blood vs Positive control"),
         x="", y="Fraction", color="Patient") +
    theme_prism_safe(base_size=14) +
    theme(legend.position="none",
          strip.text = element_text(size=12, face="bold"))
  save_plot_any(p_ct_bp, file.path(OUTDIR,"plots","celltype_fractions_faceted_Blood_PosCtrl"),
                width=18, height=11, ext=OUT_EXT)

} else {
  cat("No celltype column detected; skipping celltype fraction stats.\n")
}

div_by_sample <- per_clone %>%
  group_by(sample_id, patient_id, Tissue) %>%
  summarise(
    n_cells_clonotyped = sum(clone_size),
    n_clones = dplyr::n(),
    shannon = calc_diversity(clone_freq)$shannon,
    inv_simpson = calc_diversity(clone_freq)$inv_simpson,
    clonality = calc_diversity(clone_freq)$clonality,
    .groups="drop"
  )
write.table(div_by_sample, file.path(OUTDIR,"tables","diversity_by_sample_CARpos.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

div_patient <- div_by_sample %>%
  group_by(patient_id, Tissue) %>%
  summarise(
    shannon = weighted.mean(shannon, w=n_cells_clonotyped, na.rm=TRUE),
    inv_simpson = weighted.mean(inv_simpson, w=n_cells_clonotyped, na.rm=TRUE),
    clonality = weighted.mean(clonality, w=n_cells_clonotyped, na.rm=TRUE),
    n_clones = sum(n_clones, na.rm=TRUE),
    n_cells_clonotyped = sum(n_cells_clonotyped, na.rm=TRUE),
    .groups="drop"
  )

paired <- div_patient %>%
  pivot_wider(names_from=Tissue, values_from=c(shannon, inv_simpson, clonality, n_clones, n_cells_clonotyped)) %>%
  filter(!is.na(shannon_Blood) & !is.na(shannon_Ileum))

pvals <- bind_rows(
  paired_stats(paired$shannon_Blood, paired$shannon_Ileum) %>% mutate(metric="shannon"),
  paired_stats(paired$inv_simpson_Blood, paired$inv_simpson_Ileum) %>% mutate(metric="inv_simpson"),
  paired_stats(paired$clonality_Blood, paired$clonality_Ileum) %>% mutate(metric="clonality"),
  paired_stats(paired$n_clones_Blood, paired$n_clones_Ileum) %>% mutate(metric="n_clones")
) %>% select(metric, everything())

write.table(pvals, file.path(OUTDIR,"tables","paired_tests_metrics.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

div_long <- div_patient %>%
  filter(patient_id %in% paired$patient_id) %>%
  pivot_longer(cols=c(shannon, inv_simpson, clonality, n_clones),
               names_to="metric", values_to="value")

plot_metric <- function(metric_name, ylab_txt) {
  df <- div_long %>% filter(metric == metric_name)
  ggplot(df, aes(x=Tissue, y=value)) +
    geom_boxplot(outlier.shape=NA, alpha=0.25, width=0.55) +
    geom_point(aes(color=patient_id), size=4, position=position_jitter(width=0.06), alpha=0.95) +
    geom_line(aes(group=patient_id), color="grey55", linewidth=0.7, alpha=0.8) +
    labs(title=paste0("CAR+ clonotyped: ", metric_name, " (paired patients)"),
         x="", y=ylab_txt) +
    theme_prism_safe(base_size=18) +
    theme(legend.position="none")
}
pA <- plot_metric("shannon", "Shannon")
pB <- plot_metric("inv_simpson", "Inverse Simpson")
pC <- plot_metric("clonality", "Clonality (1 - Pielou)")
pD <- plot_metric("n_clones", "Number of clones")
save_plot_eps((pA | pB) / (pC | pD),
              file.path(OUTDIR,"plots","diversity_boxplots_paired"),
              width=16, height=11)

set.seed(SEED)
# Matched-depth repertoire bootstrap
MIN_BOOT_CELLS <- as.integer(Sys.getenv("MIN_BOOT_CELLS", unset="25"))

boot_counts <- cells_car2 %>%
  count(patient_id, Tissue, name="n_cells") %>%
  tidyr::pivot_wider(names_from=Tissue, values_from=n_cells, values_fill=0) %>%
  mutate(
    n_target = pmin(Blood, Ileum),
    eligible_min_cells = (Blood >= MIN_BOOT_CELLS) & (Ileum >= MIN_BOOT_CELLS)
  ) %>%
  arrange(desc(eligible_min_cells), patient_id)

write.table(boot_counts,
            file.path(OUTDIR,"tables","bootstrap_patient_cell_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

paired_pids <- unique(paired$patient_id)
boot_counts <- boot_counts %>%
  mutate(paired_in_metrics = patient_id %in% paired_pids)

boot_excluded <- boot_counts %>%
  filter(paired_in_metrics, !eligible_min_cells)

if (nrow(boot_excluded) > 0) {
  write.table(boot_excluded,
              file.path(OUTDIR,"tables","bootstrap_patient_excluded.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}

boot_eligible_pids <- boot_counts %>%
  filter(paired_in_metrics, eligible_min_cells) %>%
  pull(patient_id)

cat("\n[BOOT] Paired patients (from diversity metrics):", length(paired_pids), "\n")
cat("[BOOT] MIN_BOOT_CELLS:", MIN_BOOT_CELLS, "\n")
cat("[BOOT] Eligible for bootstrap:", length(boot_eligible_pids), " -> ",
    paste(boot_eligible_pids, collapse=", "), "\n")
if (nrow(boot_excluded) > 0) {
  cat("[BOOT] Excluded (insufficient CAR+ clonotyped cells in Blood/Ileum):\n")
  print(boot_excluded)
}

boot_patient <- function(pid) {
  df <- cells_car2 %>% filter(patient_id==pid)
  bl <- df %>% filter(Tissue=="Blood")
  il <- df %>% filter(Tissue=="Ileum")
  if (nrow(bl) < MIN_BOOT_CELLS || nrow(il) < MIN_BOOT_CELLS) return(NULL)
  n_target <- min(nrow(bl), nrow(il))
  one_side <- function(dat, tissue_label) {
    replicate(BOOT_B, {
      idx <- sample.int(nrow(dat), size=n_target, replace=FALSE)
      s <- dat$CTstrict[idx]
      freqs <- as.numeric(table(s)); freqs <- freqs / sum(freqs)
      div <- calc_diversity(freqs)
      c(shannon=div$shannon, inv_simpson=div$inv_simpson, clonality=div$clonality, n_clones=div$n_clones)
    }) |>
      t() |>
      as.data.frame() |>
      mutate(Tissue=tissue_label, patient_id=pid, n_target=n_target,
             n_blood=nrow(bl), n_ileum=nrow(il))
  }
  dplyr::bind_rows(one_side(bl, "Blood_downsampled"), one_side(il, "Ileum_downsampled"))
}

boot_df <- dplyr::bind_rows(lapply(boot_eligible_pids, boot_patient))

if (nrow(boot_df) > 0) {
  write.table(boot_df, file.path(OUTDIR,"tables","bootstrap_downsampled_metrics.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  boot_med <- boot_df %>%
    group_by(patient_id, Tissue) %>%
    summarise(
      shannon = median(shannon, na.rm=TRUE),
      inv_simpson = median(inv_simpson, na.rm=TRUE),
      clonality = median(clonality, na.rm=TRUE),
      n_clones = median(n_clones, na.rm=TRUE),
      n_target = first(n_target),
      n_blood = first(n_blood),
      n_ileum = first(n_ileum),
      .groups="drop"
    ) %>%
    mutate(
      Tissue_plot = dplyr::recode(Tissue,
                                  "Blood_downsampled"="Blood",
                                  "Ileum_downsampled"="Ileum",
                                  .default=Tissue),
      Tissue_plot = factor(Tissue_plot, levels=c("Blood","Ileum"))
    )
  write.table(boot_med, file.path(OUTDIR,"tables","bootstrap_downsampled_medians.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  source(file.path("single_cell", "paired_diversity_tests.R"), local = TRUE)
  boot_tests <- paired_diversity_tests(boot_med)
  write.table(boot_tests, file.path(OUTDIR, "tables", "downsampled_paired_t_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(boot_tests, digits = 8)

  boot_long_med <- boot_med %>%
    select(patient_id, Tissue_plot, shannon, inv_simpson, clonality) %>%
    pivot_longer(cols=c(shannon, inv_simpson, clonality),
                 names_to="metric", values_to="value") %>%
    mutate(
      metric = factor(metric, levels=c("shannon","inv_simpson","clonality"),
                      labels=c("Shannon","Inverse Simpson","Clonality (1 - Pielou)"))
    )

  plot_boot_metric <- function(metric_label, ylab_txt = NULL) {
    if (is.null(ylab_txt)) {
      ylab_txt <- switch(metric_label,
                         "Shannon"                  = "Shannon diversity",
                         "Inverse Simpson"          = "Inverse Simpson",
                         "Clonality (1 - Pielou)"   = "Clonality (1 - Pielou)",
                         metric_label)
    }
    df <- boot_long_med %>% dplyr::filter(metric == metric_label)
    metric_key <- switch(metric_label, "Shannon" = "shannon",
                         "Inverse Simpson" = "inv_simpson",
                         "Clonality (1 - Pielou)" = "clonality")
    test_row <- boot_tests[boot_tests$metric == metric_key, , drop = FALSE]
    test_label <- if (nrow(test_row) && test_row$status == "ok") {
      paste0("Paired t-test: P = ", formatC(test_row$p_value, format = "g", digits = 3),
             "; n = ", test_row$n_pairs, " pairs")
    } else "Paired t-test: not estimable"
    ggplot(df, aes(x = Tissue_plot, y = value)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.55) +
      geom_point(aes(color = patient_id),
                 size = 3.2, alpha = 0.95,
                 position = position_jitter(width = 0.05, height = 0)) +
      geom_line(aes(group = patient_id),
                color = "grey55", linewidth = 0.6, alpha = 0.85) +
      labs(x = "", y = ylab_txt, subtitle = test_label) +
      theme_prism_safe(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        strip.text  = element_text(size = 11, face = "bold")
      )
  }

  p_bs1 <- plot_boot_metric("Shannon")
  p_bs2 <- plot_boot_metric("Inverse Simpson")
  p_bs3 <- plot_boot_metric("Clonality (1 - Pielou)")

  p_boot_med <- (p_bs1 | p_bs2 | p_bs3) +
    plot_annotation(
      title = paste0("Matched-depth subsampling (median across ", BOOT_B, " iterations)"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  save_plot_eps(p_boot_med, file.path(OUTDIR,"plots","bootstrap_downsampled_median_paired_boxplots"),
                width = 18, height = 6)
  ggsave(file.path(OUTDIR,"plots","bootstrap_downsampled_median_paired_boxplots.png"),
         plot = p_boot_med, width = 18, height = 6, dpi = 600, bg = "white")

} else {
  cat("\n[BOOT] No bootstrap output produced (no eligible patients).\n")
}

DefaultAssay(obj_car) <- if ("SCT" %in% Assays(obj_car)) "SCT" else DefaultAssay(obj_car)

GS_naive  <- list(c("CCR7","IL7R","TCF7","LEF1","LTB"))
GS_cyto   <- list(c("NKG7","GNLY","PRF1","GZMB","CTSW"))
GS_exh    <- list(c("PDCD1","LAG3","HAVCR2","TIGIT","TOX"))
GS_prolif <- list(c("MKI67","TOP2A","TYMS","HMGB2","STMN1"))
GS_trm    <- list(c("ITGAE","CXCR6","ZNF683","RGS1","ITGA1"))
GS_treg   <- list(c("FOXP3","IL2RA","CTLA4","IKZF2"))

obj_car <- AddModuleScore(obj_car, features=GS_naive,  name="Naive",  search=TRUE)
obj_car <- AddModuleScore(obj_car, features=GS_cyto,   name="Cytotox",search=TRUE)
obj_car <- AddModuleScore(obj_car, features=GS_exh,    name="Exhaust",search=TRUE)
obj_car <- AddModuleScore(obj_car, features=GS_prolif, name="Prolif", search=TRUE)
obj_car <- AddModuleScore(obj_car, features=GS_trm,    name="TRM",    search=TRUE)
obj_car <- AddModuleScore(obj_car, features=GS_treg,   name="Treg",   search=TRUE)

ms <- obj_car@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  select(cell_barcode, patient_id, Tissue, Naive1, Cytotox1, Exhaust1, Prolif1, TRM1, Treg1) %>%
  mutate(state = c("Naive","Cytotox","Exhaust","Prolif","TRM","Treg")[max.col(select(., Naive1, Cytotox1, Exhaust1, Prolif1, TRM1, Treg1),
                                                                              ties.method="first")])
obj_car$state <- ms$state[match(rownames(obj_car@meta.data), ms$cell_barcode)]

state_comp <- ms %>%
  group_by(patient_id, Tissue, state) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(patient_id, Tissue) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

write.table(state_comp, file.path(OUTDIR,"tables","CARpos_state_fractions.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

state_wide <- state_comp %>%
  select(patient_id, Tissue, state, frac) %>%
  pivot_wider(names_from=Tissue, values_from=frac) %>%
  filter(!is.na(Blood) & !is.na(Ileum))

state_stats <- state_wide %>%
  group_by(state) %>%
  summarise(
    n = n(),
    wilcox_exact_2s = suppressWarnings(wilcox.test(Ileum, Blood, paired=TRUE, exact=TRUE)$p.value),
    wilcox_asym_2s  = suppressWarnings(wilcox.test(Ileum, Blood, paired=TRUE, exact=FALSE)$p.value),
    median_delta = median(Ileum - Blood),
    .groups="drop"
  ) %>%
  mutate(p_adj = p.adjust(wilcox_exact_2s, method="BH"))

write.table(state_stats, file.path(OUTDIR,"tables","CARpos_state_fraction_paired_stats.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

p_state <- ggplot(state_comp, aes(x=Tissue, y=frac)) +
  geom_boxplot(outlier.shape=NA, alpha=0.25, width=0.55) +
  geom_point(aes(color=patient_id), size=3.5, position=position_jitter(width=0.06), alpha=0.95) +
  geom_line(aes(group=patient_id), color="grey55", linewidth=0.6, alpha=0.8) +
  facet_wrap(~state, scales="free_y", ncol=3) +
  scale_y_continuous(labels=percent_format(accuracy=1)) +
  labs(title="CAR+ clonotyped state (dominant module score): Blood vs Ileum", x="", y="Fraction") +
  theme_prism_safe(base_size=16) +
  theme(legend.position="none")
save_plot_eps(p_state, file.path(OUTDIR,"plots","CARpos_state_fractions_faceted"), width=16, height=10)

# Patient-paired pseudobulk expression
pb_cells <- obj_car@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  filter(patient_id %in% paired$patient_id)

obj_pb <- subset_cells_safe(obj_car, pb_cells$cell_barcode, label="pseudobulk subset")

assay_pb <- if ("RNA" %in% Assays(obj_pb)) "RNA" else DefaultAssay(obj_pb)
counts <- get_counts_safe(obj_pb, assay=assay_pb)
meta   <- obj_pb@meta.data %>% tibble::rownames_to_column("cell_barcode")

meta$grp <- paste(meta$patient_id, meta$Tissue, sep="__")
cat("Pseudobulking counts...\n")

grp_levels <- unique(meta$grp)
pb_mat <- sapply(grp_levels, function(g) {
  cells <- meta$cell_barcode[meta$grp == g]
  Matrix::rowSums(counts[, cells, drop=FALSE])
})
pb_mat <- as.matrix(pb_mat)
colnames(pb_mat) <- grp_levels

coldata <- tibble(sample = grp_levels) %>%
  separate(sample, into=c("patient_id","Tissue"), sep="__", remove=FALSE) %>%
  mutate(Tissue = factor(Tissue, levels=c("Blood","Ileum"))) %>%
  as.data.frame()

pid_levels <- unique(coldata$patient_id)
pid_levels <- pid_levels[order(as.integer(gsub("\\D+","", pid_levels)), pid_levels)]
coldata$patient_id <- factor(coldata$patient_id, levels=pid_levels)

ord <- order(coldata$Tissue, coldata$patient_id)
coldata <- coldata[ord, , drop=FALSE]
rownames(coldata) <- coldata$sample
pb_mat <- pb_mat[, coldata$sample, drop=FALSE]

keepg <- rowSums(pb_mat) >= 10
pb_mat <- pb_mat[keepg, , drop=FALSE]
cat("Genes kept for DE:", nrow(pb_mat), "\n")

if (!HAS_DESEQ2) stop("DESeq2 not available in this R environment.")
cat("Running DESeq2 with design ~ patient_id + Tissue ...\n")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(pb_mat),
                                      colData = coldata,
                                      design = ~ patient_id + Tissue)
dds <- DESeq2::DESeq(dds, quiet=TRUE)
res <- DESeq2::results(dds, contrast=c("Tissue","Ileum","Blood"))

de_res <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(p_adj = padj,
         avg_log2FC = log2FoldChange)

write.table(de_res, file.path(OUTDIR,"tables","DE_pseudobulk_Ileum_vs_Blood_CARpos.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

if (HAS_FGSEA) {
  cat("Running OFFLINE GSEA (Hallmark) from cached RDS...\n")
  pathways <- load_hallmark_pathways(HALLMARK_RDS)
  ranks <- de_res$stat
  names(ranks) <- de_res$gene
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing=TRUE)

  fg <- NULL
  fg <- tryCatch({
    if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
      fgsea::fgseaMultilevel(pathways=pathways, stats=ranks)
    } else {
      fgsea::fgsea(pathways=pathways, stats=ranks, nperm=10000)
    }
  }, error=function(e) {
    stop("fgsea failed: ", conditionMessage(e))
  })

  fg <- as.data.frame(fg) %>% arrange(padj)
  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(x, collapse=";"), FUN.VALUE=character(1))
  }

  write.table(fg, file.path(OUTDIR,"tables","GSEA_Hallmark_fgsea.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  fg2 <- fg %>% filter(!is.na(padj))
  top_pos <- fg2 %>% arrange(desc(NES)) %>% slice_head(n=15)
  top_neg <- fg2 %>% arrange(NES) %>% slice_head(n=15)
  top_fg <- bind_rows(top_pos, top_neg) %>%
    mutate(pathway = factor(pathway, levels=pathway[order(NES)]))

  p_gsea <- ggplot(top_fg, aes(x=pathway, y=NES)) +
    geom_col(alpha=0.85) +
    coord_flip() +
    labs(title="GSEA (Hallmark) on DESeq2 stat: Ileum vs Blood (CAR+)",
         x="", y="NES") +
    theme_prism_safe(base_size=14)
  save_plot_eps(p_gsea, file.path(OUTDIR,"plots","GSEA_Hallmark_top_pos_neg"), width=12, height=9)
} else {
  cat("GSEA skipped: fgsea not installed.\n")
}

cat("\n=== SessionInfo ===\n")
print(sessionInfo())
cat("\n=== DONE ===\n")
