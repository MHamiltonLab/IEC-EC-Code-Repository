#!/usr/bin/env Rscript
# Single-cell marker expression
# Display T-cell and myeloid markers on UMAPs and compare CAR-positive marker distributions descriptively.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(patchwork)
  library(grid)
})

HAS_GGH4X      <- requireNamespace("ggh4x", quietly = TRUE)
HAS_GGPRISM    <- requireNamespace("ggprism", quietly = TRUE)
HAS_WES        <- requireNamespace("wesanderson", quietly = TRUE)
HAS_BREWER     <- requireNamespace("RColorBrewer", quietly = TRUE)

BASE_DIR    <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
QC_DIR      <- Sys.getenv("QC_DIR",   unset=BASE_DIR)

INPUT_RDS   <- Sys.getenv("INPUT_RDS", unset=file.path(QC_DIR, "single_cell_object.rds"))
SAMPLE_META <- Sys.getenv("SAMPLE_META", unset=file.path(QC_DIR, "single_cell_metadata.tsv"))
OUTDIR      <- Sys.getenv("OUTDIR", unset=file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "single_cell", "marker_expression"))

CAR_GENE    <- Sys.getenv("CAR_GENE", unset="CILTACELCAR")
CAR_MIN_UMI <- as.integer(Sys.getenv("CAR_MIN_UMI", unset="1"))

PT_SIZE_ALL   <- as.numeric(Sys.getenv("PT_SIZE_ALL",  unset="0.25"))
PT_ALPHA_ALL  <- as.numeric(Sys.getenv("PT_ALPHA_ALL", unset="0.70"))
PT_SIZE_CAR   <- as.numeric(Sys.getenv("PT_SIZE_CAR",  unset="0.35"))
PT_ALPHA_CAR  <- as.numeric(Sys.getenv("PT_ALPHA_CAR", unset="0.95"))

OUT_EXT     <- Sys.getenv("OUT_EXT", unset="eps")

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"plots"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"tables"),showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUTDIR,"logs"),  showWarnings=FALSE, recursive=TRUE)

cat("=== START CAR marker validation ===\n")
cat("QC_DIR     :", QC_DIR, "\n")
cat("INPUT_RDS  :", INPUT_RDS, "\n")
cat("SAMPLE_META:", SAMPLE_META, "\n")
cat("OUTDIR     :", OUTDIR, "\n")
cat("CAR_GENE   :", CAR_GENE, " CAR_MIN_UMI:", CAR_MIN_UMI, "\n")
cat("OUT_EXT    :", OUT_EXT, "\n")
cat("===================================\n\n")

stopifnot(file.exists(INPUT_RDS))
stopifnot(file.exists(SAMPLE_META))

theme_prism_safe <- function(base_size=18, family="Arial") {
  if (HAS_GGPRISM) ggprism::theme_prism(base_size = base_size, base_family = family)
  else theme_classic(base_size = base_size, base_family = family)
}

guide_axis_trunc_safe <- function() {
  if (!HAS_GGH4X) return(NULL)
  ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(3, "cm")
  )
}

theme_umap_arrow <- function(base_size=22, legend_pos="right", family="Arial") {
  axg <- guide_axis_trunc_safe()
  p <- theme_classic(base_size=base_size, base_family=family) +
    theme(
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      axis.title = element_text(hjust = 0),
      axis.line.x = element_line(
        linewidth = 1.1,
        arrow = arrow(length = unit(0.32, "cm"), ends = "last")
      ),
      axis.line.y = element_line(
        linewidth = 1.1,
        arrow = arrow(length = unit(0.32, "cm"), ends = "last")
      ),
      legend.position = legend_pos,
      plot.title = element_text(face="bold")
    )
  if (!is.null(axg)) p <- p + guides(x = axg, y = axg)
  p
}

save_plot_eps <- function(p, filename, width, height) {
  out <- filename
  if (!grepl("\\.eps$", out, ignore.case=TRUE)) out <- paste0(out, ".eps")
  if (capabilities("cairo")) {
    ggsave(out, plot=p, width=width, height=height,
           device=function(...) grDevices::cairo_ps(..., onefile=FALSE, fallback_resolution=600),
           dpi=600)
  } else {
    ggsave(out, plot=p, width=width, height=height,
           device=function(...) grDevices::postscript(..., onefile=FALSE, paper="special", horizontal=FALSE),
           dpi=600)
  }
  invisible(out)
}

wes_or_fallback <- function(n, name="Zissou1") {
  if (HAS_WES) {
    return(wesanderson::wes_palette(name, n=n, type="discrete"))
  }
  if (HAS_BREWER) {
    n2 <- max(3, min(8, n))
    return(RColorBrewer::brewer.pal(n2, "Set2")[seq_len(n)])
  }
  grDevices::rainbow(n)
}

tissue_cols <- c("Blood"="#E56AA6", "Ileum"="#2AA1B1")

car_cols <- c("FALSE"="grey82", "TRUE"=wes_or_fallback(2, "Darjeeling1")[2])

get_counts <- function(obj, assay="RNA") {
  m <- tryCatch(GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m)) return(m)

  ass <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(ass), error=function(e) character(0))
  if (length(layers)==0) stop("No counts slot or layers found for assay=", assay)

  out <- NULL
  use <- layers[grepl("count", layers, ignore.case=TRUE)]
  if (length(use)==0) use <- layers
  for (ly in use) {
    mat <- tryCatch(SeuratObject::LayerData(ass, layer=ly), error=function(e) NULL)
    if (is.null(mat)) next
    if (!inherits(mat,"dgCMatrix")) mat <- as(mat,"dgCMatrix")
    out <- if (is.null(out)) mat else out + mat
  }
  out
}

gene_umi <- function(obj, gene, assay="RNA") {
  m <- get_counts(obj, assay=assay)
  if (!gene %in% rownames(m)) return(setNames(rep(0, ncol(m)), colnames(m)))
  v <- Matrix::colSums(m[gene,,drop=FALSE])
  setNames(as.numeric(v), names(v))
}

choose_assay_for_features <- function(obj, genes, preferred=c("SCT","RNA")) {
  preferred <- preferred[preferred %in% Assays(obj)]
  if (length(preferred)==0) return(DefaultAssay(obj))
  hits <- sapply(preferred, function(a) sum(genes %in% rownames(obj[[a]])))
  preferred[which.max(hits)]
}

keep_present <- function(obj, genes) {
  all_feats <- unique(unlist(lapply(Assays(obj), function(a) rownames(obj[[a]]))))
  genes[genes %in% all_feats]
}

cat("Loading object...\n")
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

keep_samples <- unique(sample_meta$sample_id)
obj <- subset(obj, subset = sample_id %in% keep_samples)

md <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  mutate(sample_id = as.character(sample_id)) %>%
  select(-any_of(c("patient_id", "Tissue"))) %>%
  left_join(sample_meta, by="sample_id") %>%
  tibble::column_to_rownames("cell_barcode")
obj@meta.data <- as.data.frame(md)

cat("\n[DIAGNOSTICS] Tissue mapped counts:\n")
print(table(obj$Tissue, useNA="ifany"))

red <- if ("umap" %in% names(obj@reductions)) "umap" else {
  u <- grep("umap", names(obj@reductions), ignore.case=TRUE, value=TRUE)
  if (length(u)>0) u[1] else names(obj@reductions)[1]
}
cat("Using reduction:", red, "\n")

assay_counts <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cat("Using assay for CAR UMI:", assay_counts, "\n")
obj$CAR_umi <- gene_umi(obj, CAR_GENE, assay=assay_counts)
obj$CAR_pos <- obj$CAR_umi >= CAR_MIN_UMI
cat("CAR+ cells:", sum(obj$CAR_pos, na.rm=TRUE), "/", ncol(obj), "\n")

t_markers <- c("CD3D","CD3E","TRAC","TRBC1","TRBC2","CD247","LCK","LTB","IL7R")
m_markers <- c("LYZ","LST1","S100A8","S100A9","FCGR3A","MS4A7","LGALS3")

t_markers_ok <- keep_present(obj, t_markers)
m_markers_ok <- keep_present(obj, m_markers)

cat("\nT markers present:", paste(t_markers_ok, collapse=", "), "\n")
cat("Myeloid markers present:", paste(m_markers_ok, collapse=", "), "\n")
if (length(t_markers_ok) < 3) cat("[WARN] Few T markers found; check gene symbols / assay.\n")
if (length(m_markers_ok) < 3) cat("[WARN] Few myeloid markers found; check gene symbols / assay.\n")

assay_expr <- choose_assay_for_features(obj, unique(c(t_markers_ok, m_markers_ok, CAR_GENE)))
DefaultAssay(obj) <- assay_expr
cat("Using assay for expression plots:", assay_expr, "\n")

if (length(t_markers_ok) >= 3) obj <- AddModuleScore(obj, features=list(t_markers_ok), name="TcellScore", search=TRUE)
if (length(m_markers_ok) >= 3) obj <- AddModuleScore(obj, features=list(m_markers_ok), name="MyeloidScore", search=TRUE)

p_tissue <- DimPlot(obj, reduction=red, group.by="Tissue",
                    pt.size=PT_SIZE_ALL, raster=FALSE, alpha=PT_ALPHA_ALL) +
  scale_color_manual(values=tissue_cols, na.value="grey85") +
  labs(title="All cells: Blood vs Ileum", x="UMAP1", y="UMAP2") +
  theme_umap_arrow(base_size=28, family="Arial") +
  theme(legend.title=element_blank())
save_plot_eps(p_tissue, file.path(OUTDIR,"plots","UMAP_all_cells_Tissue"), width=12, height=10)

p_car <- DimPlot(obj, reduction=red, group.by="CAR_pos",
                 pt.size=PT_SIZE_ALL, raster=FALSE, alpha=PT_ALPHA_ALL) +
  scale_color_manual(values=car_cols, na.value="grey85") +
  labs(title=paste0("All cells: CAR+ calls (", CAR_GENE, " UMI ≥ ", CAR_MIN_UMI, ")"),
       x="UMAP1", y="UMAP2") +
  theme_umap_arrow(base_size=26, family="Arial") +
  theme(legend.title=element_blank())
save_plot_eps(p_car, file.path(OUTDIR,"plots","UMAP_all_cells_CARpos"), width=12, height=10)

grad_t <- c("grey95", tissue_cols[["Ileum"]])
grad_m <- c("grey95", "#8C7AA9")

plot_feature_panel <- function(features, title, cols, ncol=3) {
  feats <- features[features %in% rownames(obj[[assay_expr]])]
  if (length(feats)==0) return(NULL)
  FeaturePlot(
    obj, features=feats, reduction=red,
    cols=cols, order=TRUE, pt.size=PT_SIZE_ALL,
    min.cutoff="q05", max.cutoff="q95", ncol=ncol, raster=FALSE
  ) &
    theme_umap_arrow(base_size=18, family="Arial") &
    theme(plot.title = element_text(face="bold", size=18)) +
    plot_annotation(title = title)
}

p_car_gene <- NULL
if (CAR_GENE %in% rownames(obj[[assay_expr]])) {
  p_car_gene <- FeaturePlot(
    obj, features=CAR_GENE, reduction=red,
    cols=c("grey95", tissue_cols[["Blood"]]),
    order=TRUE, pt.size=PT_SIZE_ALL,
    min.cutoff="q05", max.cutoff="q95", raster=FALSE
  ) +
    labs(title=paste0("CAR transgene signal: ", CAR_GENE), x="UMAP1", y="UMAP2") +
    theme_umap_arrow(base_size=26, family="Arial")
  save_plot_eps(p_car_gene, file.path(OUTDIR,"plots","UMAP_Feature_CARgene_allcells"), width=12, height=10)
}

p_tcells <- plot_feature_panel(t_markers_ok[1:min(6,length(t_markers_ok))],
                               "T-cell markers (all cells)", cols=grad_t, ncol=3)
if (!is.null(p_tcells)) {
  save_plot_eps(p_tcells, file.path(OUTDIR,"plots","UMAP_Feature_Tcell_markers_allcells"),
                width=16, height=10)
}

p_myel <- plot_feature_panel(m_markers_ok[1:min(6,length(m_markers_ok))],
                             "Myeloid markers (all cells)", cols=grad_m, ncol=3)
if (!is.null(p_myel)) {
  save_plot_eps(p_myel, file.path(OUTDIR,"plots","UMAP_Feature_Myeloid_markers_allcells"),
                width=16, height=10)
}

score_feats <- c()
if ("TcellScore1" %in% colnames(obj@meta.data)) score_feats <- c(score_feats, "TcellScore1")
if ("MyeloidScore1" %in% colnames(obj@meta.data)) score_feats <- c(score_feats, "MyeloidScore1")

if (length(score_feats) > 0) {
  p_scores <- FeaturePlot(
    obj, features=score_feats, reduction=red,
    cols=c("grey95", "#2AA1B1"), order=TRUE,
    pt.size=PT_SIZE_ALL, raster=FALSE, ncol=length(score_feats)
  ) &
    theme_umap_arrow(base_size=20, family="Arial") +
    plot_annotation(title="Module scores (all cells): T-cell vs Myeloid")
  save_plot_eps(p_scores, file.path(OUTDIR,"plots","UMAP_ModuleScores_Tcell_vs_Myeloid"),
                width=14, height=6)
}

obj_car <- subset(obj, subset = CAR_pos & Tissue %in% c("Blood","Ileum"))
cat("\nCAR+ cells (Blood/Ileum):", ncol(obj_car), "\n")
if (ncol(obj_car) == 0) stop("No CAR+ cells after filtering to Blood/Ileum.")

DefaultAssay(obj_car) <- assay_expr

fill_cols <- tissue_cols

obj_car$Tissue <- factor(obj_car$Tissue, levels=c("Blood","Ileum"))

mk_vln <- function(features, title, fname, ncol=3) {
  feats <- features[features %in% rownames(obj_car[[assay_expr]])]
  if (length(feats)==0) return(NULL)

  p <- VlnPlot(
    obj_car, features=feats, group.by="Tissue",
    pt.size=0, combine=TRUE, ncol=ncol
  ) &
    theme_prism_safe(base_size=16, family="Arial") &
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title   = element_text(face="bold"),
      axis.text.x  = element_text(angle=15, hjust=1)
    )

  p <- p & scale_fill_manual(values=fill_cols, drop=FALSE)

  p <- p + plot_annotation(title = title)
  save_plot_eps(p, file.path(OUTDIR,"plots", fname), width=16, height=10)
  p
}

mk_vln(t_markers_ok[1:min(9,length(t_markers_ok))],
       "CAR+ cells: T-cell markers by Tissue (Blood vs Ileum)",
       "Vln_CARpos_Tcell_markers_by_Tissue", ncol=3)

mk_vln(m_markers_ok[1:min(9,length(m_markers_ok))],
       "CAR+ cells: Myeloid markers (negative control) by Tissue",
       "Vln_CARpos_Myeloid_markers_by_Tissue", ncol=3)

score_feats_car <- c()
if ("TcellScore1" %in% colnames(obj_car@meta.data)) score_feats_car <- c(score_feats_car, "TcellScore1")
if ("MyeloidScore1" %in% colnames(obj_car@meta.data)) score_feats_car <- c(score_feats_car, "MyeloidScore1")
if (length(score_feats_car) > 0) {
  mk_vln(score_feats_car,
         "CAR+ cells: module scores by Tissue (T-cell vs Myeloid)",
         "Vln_CARpos_ModuleScores_by_Tissue", ncol=2)
}

cat("\n=== SessionInfo ===\n")
print(sessionInfo())
cat("\n=== DONE marker validation ===\n")
