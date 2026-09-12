#!/usr/bin/env Rscript
# Bulk RNA expression analysis
# Analyze log-transformed FPKM with limma, pathway enrichment, and insertion-associated expression comparisons.

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
meta_path <- file.path(data_dir, "expression_metadata.tsv")
fpkm_path <- file.path(data_dir, "expression_fpkm.tsv")

out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "bulk_rna", "expression_analysis")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(123)

has_repel <- requireNamespace("ggrepel", quietly = TRUE)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

foggy_sf <- c("#4C6A87", "#7B99B6", "#9CB7CE", "#C0CEDD",
              "#9AA6B2", "#72808E", "#B8A9B4", "#C9BCC6", "#A6B8BE")

get_muted_palette <- function(n){
  base <- foggy_sf
  if(n <= length(base)) return(base[seq_len(n)])
  return(colorRampPalette(base)(n))
}
scale_fill_muted  <- function(n = NULL){ scale_fill_manual(values = get_muted_palette(ifelse(is.null(n), 9, n))) }
scale_color_muted <- function(n = NULL){ scale_color_manual(values = get_muted_palette(ifelse(is.null(n), 9, n))) }

base_font     <- "Arial"
fallback_font <- "Helvetica"

has_showtext <- requireNamespace("showtext", quietly = TRUE) &&
  requireNamespace("sysfonts", quietly = TRUE)
if (has_showtext) {
  library(showtext); library(sysfonts)
  arial_candidates <- Sys.getenv("PLOT_FONT_FILE", unset = "")
  arial_path <- arial_candidates[file.exists(arial_candidates)][1]
  if (!is.na(arial_path)) {
    sysfonts::font_add("Arial", arial_path)
    showtext::showtext_auto(enable = TRUE)
    message("Using Arial via showtext from: ", arial_path)
  } else {
    message("Arial TTF not found; will fallback to Helvetica/sans if needed.")
  }
} else {
  message("showtext/sysfonts not installed; proceeding without them.")
}

test_family <- tryCatch({
  grDevices::pdf(NULL, family = base_font); grDevices::dev.off(); base_font
}, error = function(e) fallback_font)
base_font <- test_family

theme_set(
  theme_classic(base_family = base_font) +
    theme(
      plot.title  = element_text(face = "bold"),
      legend.position = "right",
      axis.text   = element_text(color = "gray20"),
      axis.title  = element_text(color = "gray20"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text  = element_text(face = "bold")
    )
)

pdf_device_simple <- function(filename, width = 6, height = 4, family = base_font, ...) {
  fam <- if (tolower(family) %in% c("arial")) base_font else family
  grDevices::pdf(file = filename, width = width, height = height,
                 family = fam, useDingbats = FALSE, ...)
}
eps_device_simple <- function(filename, width = 6, height = 4, family = "Helvetica", ...) {
  grDevices::postscript(file = filename, width = width, height = height,
                        onefile = FALSE, horizontal = FALSE, paper = "special",
                        family = family, ...)
}

save_plot_all <- function(plot, file_base, width = 6, height = 4) {

  pdf_file <- file.path(plot_dir, paste0(file_base, ".pdf"))
  pdf_device_simple(filename = pdf_file, width = width, height = height, family = base_font)
  if (exists("showtext_begin", where = asNamespace("showtext"), inherits = FALSE)) {
    try(showtext::showtext_begin(), silent = TRUE)
  }
  print(plot)
  if (exists("showtext_end", where = asNamespace("showtext"), inherits = FALSE)) {
    try(showtext::showtext_end(), silent = TRUE)
  }
  grDevices::dev.off()

  eps_file <- file.path(plot_dir, paste0(file_base, ".eps"))
  eps_device_simple(filename = eps_file, width = width, height = height, family = "Helvetica")
  print(plot)
  grDevices::dev.off()
}

save_pheatmap_all <- function(pheat, file_base, width = 7, height = 6) {
  pdf_file <- file.path(plot_dir, paste0(file_base, ".pdf"))
  pdf_device_simple(filename = pdf_file, width = width, height = height, family = base_font)
  if (exists("showtext_begin", where = asNamespace("showtext"), inherits = FALSE)) {
    try(showtext::showtext_begin(), silent = TRUE)
  }
  grid::grid.newpage(); print(pheat)
  if (exists("showtext_end", where = asNamespace("showtext"), inherits = FALSE)) {
    try(showtext::showtext_end(), silent = TRUE)
  }
  grDevices::dev.off()

  eps_file <- file.path(plot_dir, paste0(file_base, ".eps"))
  eps_device_simple(filename = eps_file, width = width, height = height, family = "Helvetica")
  grid::grid.newpage(); print(pheat)
  grDevices::dev.off()
}

message("Loading metadata...")
# Expression and sample inputs
meta <- readr::read_tsv(meta_path, show_col_types = FALSE)

message("Loading FPKM matrix...")
fpkm <- readr::read_tsv(fpkm_path, show_col_types = FALSE)

sample_cols <- intersect(meta$Sample_RNAseq, colnames(fpkm))
if(length(sample_cols) == 0){
  stop("No overlapping sample columns found between metadata and FPKM.")
}
meta <- meta %>% filter(Sample_RNAseq %in% sample_cols) %>%
  arrange(match(Sample_RNAseq, sample_cols))
fpkm <- fpkm %>% select(GeneID, Gene_Name, Gene_Biotype, all_of(sample_cols))

for (nm in c("COMET_perc_CD3", "COMET_perc_total")) {
  if (nm %in% colnames(meta)) meta[[nm]] <- suppressWarnings(as.numeric(meta[[nm]]))
}

meta_nd <- meta %>% filter(is.na(Cohort) | Cohort != "Dup")
samples_nd <- meta_nd$Sample_RNAseq

# Gene collapse and transformation
expr_long <- fpkm %>%
  select(Gene_Name, all_of(sample_cols)) %>%
  group_by(Gene_Name) %>%
  summarise(across(everything(),
                   ~ if(all(is.na(.x))) NA_real_ else median(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  filter(!is.na(Gene_Name) & Gene_Name != "")

expr_mat <- expr_long %>% column_to_rownames("Gene_Name") %>% as.matrix()

log_expr <- log2(expr_mat + 1)

expr_df <- as.data.frame(log_expr) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample_RNAseq", values_to = "logFPKM") %>%
  left_join(meta_nd, by = "Sample_RNAseq")

p_violin <- ggplot(expr_df, aes(x = Sample_RNAseq, y = logFPKM, fill = Cohort)) +
  geom_violin(trim = TRUE, color = "gray30", linewidth = 0.3) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 5, color = "black") +
  coord_flip() +
  scale_fill_muted() +
  labs(title = "Expression distributions (log2(FPKM+1))", x = "Sample", y = "log2(FPKM+1)")
save_plot_all(p_violin, "qc_violin_logFPKM", width = 7.5, height = 8)

expr_nd <- log_expr[, samples_nd, drop = FALSE]
pca <- prcomp(t(expr_nd), center = TRUE, scale. = FALSE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Sample_RNAseq") %>%
  left_join(meta_nd, by = "Sample_RNAseq")
var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Cohort, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_muted() +
  labs(title = "PCA of samples (log2(FPKM+1))",
       x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"))
save_plot_all(p_pca, "qc_pca_all", width = 6.5, height = 5.2)

dist_m <- as.matrix(dist(t(expr_nd)))
annot_df <- meta_nd %>% select(Sample_RNAseq, Cohort, Tissue, Timepoint) %>% column_to_rownames("Sample_RNAseq")
pheat <- pheatmap(
  dist_m, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
  annotation_row = annot_df, annotation_col = annot_df,
  main = "Sample distance (Euclidean on log2(FPKM+1)) [No Duplicates]",
  silent = TRUE
)
save_pheatmap_all(pheat, "qc_sample_distance", width = 8, height = 7)

meta_de <- meta_nd %>%
  filter(Cohort %in% c("IEC", "Control")) %>%
  mutate(Cohort = factor(Cohort, levels = c("Control", "IEC")))
samples_de <- meta_de$Sample_RNAseq
expr_de <- log_expr[, samples_de, drop = FALSE]

# limma differential expression
design <- model.matrix(~ Cohort, data = meta_de)
colnames(design) <- c("Intercept", "IEC_vs_Control")

fit <- lmFit(expr_de, design)
fit <- eBayes(fit, robust = TRUE)

de_tbl <- topTable(fit, coef = "IEC_vs_Control", number = Inf, sort.by = "P") %>%
  rownames_to_column("Gene") %>% as_tibble()
readr::write_tsv(de_tbl, file.path(table_dir, "DE_IEC_vs_Control.tsv"))

de_counts <- de_tbl %>%
  filter(adj.P.Val < 0.05) %>%
  summarise(Up_in_IEC = sum(logFC > 0), Down_in_IEC = sum(logFC < 0), .groups = "drop")
readr::write_tsv(de_counts, file.path(table_dir, "DE_IEC_vs_Control_counts.tsv"))

de_tbl <- de_tbl %>% mutate(
  Signif = case_when(adj.P.Val < 0.05 & abs(logFC) >= 1 ~ "FDR<0.05 & |logFC|>=1",
                     adj.P.Val < 0.05 ~ "FDR<0.05",
                     TRUE ~ "NS"),
  negLogP = -log10(P.Value)
)
label_up   <- de_tbl %>% filter(adj.P.Val < 0.05, logFC > 0)  %>% arrange(adj.P.Val) %>% slice_head(n = 10) %>% pull(Gene)
label_down <- de_tbl %>% filter(adj.P.Val < 0.05, logFC < 0)  %>% arrange(adj.P.Val) %>% slice_head(n = 10) %>% pull(Gene)
label_genes <- unique(c(label_up, label_down))

p_volcano <- ggplot(de_tbl, aes(x = logFC, y = negLogP, color = Signif)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c(
    "FDR<0.05 & |logFC|>=1" = "#66c2a5",
    "FDR<0.05"              = "#8da0cb",
    "NS"                    = "#bdbdbd"
  )) +
  labs(title = "IEC vs Control (limma on log2(FPKM+1))", x = "log2 fold-change (IEC vs Control)", y = "-log10(P)") +
  theme(legend.title = element_blank())

if(length(label_genes) > 0){
  lab_df <- de_tbl %>% filter(Gene %in% label_genes)
  if (has_repel) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data = lab_df, aes(label = Gene),
        size = 3, max.overlaps = Inf, box.padding = 0.25, point.padding = 0.2, seed = 123,
        segment.color = "grey60", segment.size = 0.3
      )
  } else {
    p_volcano <- p_volcano +
      geom_text(data = lab_df, aes(label = Gene), size = 2.7, vjust = -0.3)
  }
}
save_plot_all(p_volcano, "DE_volcano_IEC_vs_Control", width = 6.8, height = 5.4)

de_sig <- de_tbl %>% filter(adj.P.Val < 0.05)
n_up <- min(25, sum(de_sig$logFC > 0)); n_dn <- min(25, sum(de_sig$logFC < 0))
if(n_up + n_dn >= 2){
  top_up <- de_sig %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% slice_head(n = n_up)
  top_dn <- de_sig %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% slice_head(n = n_dn)
  de_top <- bind_rows(top_up, top_dn)
} else {
  de_top <- de_tbl %>% arrange(adj.P.Val) %>% slice_head(n = min(50, nrow(.)))
}
hm_genes <- unique(de_top$Gene)
hm_mat <- expr_de[hm_genes, , drop = FALSE]
ann_col <- meta_de %>% select(Sample_RNAseq, Cohort, Tissue, Timepoint) %>% column_to_rownames("Sample_RNAseq")
ann_row <- de_top %>% select(Gene, logFC) %>%
  mutate(Direction = ifelse(logFC > 0, "Up_in_IEC", "Down_in_IEC")) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames("Gene") %>% select(Direction)
pheat_top <- pheatmap(
  hm_mat, scale = "row",
  clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
  annotation_col = ann_col, annotation_row = ann_row,
  show_rownames = TRUE, main = "Top DE genes (balanced up/down)",
  silent = TRUE
)
save_pheatmap_all(pheat_top, "DE_TopBalanced_heatmap", width = 7.8, height = 9)

stopifnot("t" %in% names(de_tbl), "Gene" %in% names(de_tbl))
ranks <- de_tbl %>%
  dplyr::filter(is.finite(t)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(Gene, stat = t) %>%
  tibble::deframe()

if (exists("pathways")) rm(pathways)

msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct()

hallmark_pathways <- split(x = msig_h$gene_symbol, f = msig_h$gs_name)

hallmark_pathways <- lapply(hallmark_pathways, unique)

cat("Unique Hallmark set names: ", length(unique(msig_h$gs_name)), "\n")
cat("Length of hallmark_pathways: ", length(hallmark_pathways), "\n")
print(head(names(hallmark_pathways)))

names(ranks) <- toupper(trimws(names(ranks)))
hallmark_pathways <- lapply(hallmark_pathways, function(g) unique(toupper(trimws(g))))

overlap_counts <- vapply(hallmark_pathways, function(g) sum(g %in% names(ranks)), integer(1))
cat("[fgsea] pathways with >= 15 overlapping genes: ", sum(overlap_counts >= 15), "\n")

minSize <- 15; maxSize <- 500
hallmark_ok <- hallmark_pathways[overlap_counts >= minSize & overlap_counts <= maxSize]
cat("[fgsea] pathways passing size filter: ", length(hallmark_ok), "\n")

if (length(hallmark_ok) == 0) {
  for (alt_min in c(12, 10, 8)) {
    hall_ok2 <- hallmark_pathways[overlap_counts >= alt_min & overlap_counts <= maxSize]
    cat("[fgsea] trying minSize =", alt_min, " -> N sets:", length(hall_ok2), "\n")
    if (length(hall_ok2) > 0) { hallmark_ok <- hall_ok2; minSize <- alt_min; break }
  }
  if (length(hallmark_ok) == 0) stop("No Hallmark pathways have enough overlap; check gene IDs.")
}

ranks_sorted <- sort(ranks, decreasing = TRUE)

fgsea_res <- fgsea::fgseaMultilevel(
  pathways = hallmark_ok,
  stats    = ranks_sorted,
  minSize  = minSize,
  maxSize  = maxSize,
  nproc    = 0
)

if (nrow(fgsea_res) == 0) {
  message("[fgsea] Multilevel returned 0 rows; trying fgseaSimple (permutations).")
  fgsea_res <- fgsea::fgseaSimple(
    pathways = hallmark_ok,
    stats    = ranks_sorted,
    nperm    = 10000,
    minSize  = minSize,
    maxSize  = maxSize,
    nproc    = 0
  )
}

fgsea_res <- fgsea_res %>% dplyr::arrange(padj, dplyr::desc(NES))
readr::write_tsv(fgsea_res, file.path(table_dir, "GSEA_Hallmark_full.tsv"))
print(fgsea_res %>% dplyr::select(pathway, size, NES, pval, padj) %>% head(10))

fgsea_res <- fgsea_res %>%
  dplyr::arrange(padj, dplyr::desc(NES))

readr::write_tsv(fgsea_res, file.path(table_dir, "GSEA_Hallmark_full.tsv"))

print(fgsea_res %>% dplyr::select(pathway, NES, pval, padj) %>% head(10))

gsea_fdr <- 0.05
topN     <- 15

to_plot_sig <- fgsea_res %>%
  filter(padj < gsea_fdr) %>%
  arrange(padj) %>%
  slice_head(n = topN)

p_gsea <- ggplot(to_plot_sig, aes(y = reorder(pathway, NES), x = NES)) +
  geom_col(fill = "#4C6A87") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  labs(title = "GSEA (Hallmark): IEC vs Control",
       x = "Normalized Enrichment Score", y = "Pathway") +
  coord_cartesian(clip = "off")

save_plot_all(p_gsea, "GSEA_Hallmark_bar", width = 7.2, height = 5.0)

readr::write_tsv(to_plot_sig, file.path(table_dir, "GSEA_Hallmark_bar_shown.tsv"))

 to_plot_all <- fgsea_res %>%
   arrange(padj) %>%
   slice_head(n = topN) %>%
   mutate(Sig = if_else(padj < gsea_fdr, "FDR < 0.05", "NS"))

 p_gsea2 <- ggplot(to_plot_all, aes(y = reorder(pathway, NES), x = NES, fill = Sig)) +
   geom_col() +
   scale_fill_manual(values = c("FDR < 0.05" = "#4C6A87", "NS" = "#C0CEDD")) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
   labs(title = "GSEA (Hallmark): IEC vs Control",
        x = "Normalized Enrichment Score", y = "Pathway", fill = NULL)
 save_plot_all(p_gsea2, "GSEA_Hallmark_bar_2color", width = 7.2, height = 5.0)
 readr::write_tsv(to_plot_all, file.path(table_dir, "GSEA_Hallmark_bar_shown_2color.tsv"))

row_lc <- tolower(rownames(log_expr))
cilta_idx <- which(trimws(row_lc) == "ciltacel")
if(length(cilta_idx) == 1){
  cilta_row <- rownames(log_expr)[cilta_idx]

  car_ic <- expr_df %>%
    filter(Gene == cilta_row, Cohort %in% c("IEC", "Control")) %>%
    mutate(Cohort = factor(Cohort, levels = c("Control", "IEC")))
  w_p <- tryCatch({ wilcox.test(logFPKM ~ Cohort, data = car_ic, exact = FALSE)$p.value }, error = function(e) NA_real_)
  p_lab <- if(is.na(w_p)) "Wilcoxon p = NA" else paste0("Wilcoxon p = ", formatC(w_p, format = "e", digits = 2))

  p_cilta_ic <- ggplot(car_ic, aes(x = Cohort, y = logFPKM, fill = Cohort)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.85) +
    scale_fill_manual(values = c("Control" = "#9CB7CE", "IEC" = "#4C6A87")) +
    labs(title = "Cilta-cel expression: IEC vs Control", x = "Cohort", y = "log2(FPKM+1)") +
    annotate("text", x = 1.5, y = max(car_ic$logFPKM, na.rm = TRUE) * 1.05, label = p_lab, size = 3.6, vjust = 0)
  save_plot_all(p_cilta_ic, "CAR_ciltacel_IEC_vs_Control", width = 5.2, height = 4.6)

  tibble(Test = "Wilcoxon rank-sum", Gene = cilta_row, Group1 = "Control", Group2 = "IEC", p_value = w_p) %>%
    write_tsv(file.path(table_dir, "CAR_ciltacel_Wilcoxon.tsv"))
} else {
  warning("Gene 'ciltacel' not found in expression matrix; skipping IHC correlations.")
}

meta_ins <- meta_nd %>%
  mutate(Insertion = ifelse(is.na(Insertion), NA_character_, Insertion)) %>%
  mutate(Insertion_genes = strsplit(ifelse(is.na(Insertion), "", Insertion), ",")) %>%
  mutate(Insertion_genes = lapply(Insertion_genes, function(x) {
    x <- trimws(x); x <- x[nchar(x) > 0]; unique(x)
  }))

all_inserted_genes <- unique(unlist(meta_ins$Insertion_genes))
all_inserted_genes <- all_inserted_genes[all_inserted_genes %in% rownames(log_expr)]

expr_df_ins <- expr_df %>%
  mutate(InsertionGene = ifelse(Gene %in% all_inserted_genes, "Inserted", "Other")) %>%
  mutate(CohortIC = case_when(Cohort %in% c("IEC", "Control") ~ Cohort, TRUE ~ NA_character_)) %>%
  filter(!is.na(CohortIC))

wilcox_grouped <- expr_df_ins %>%
  group_by(CohortIC) %>%
  group_modify(~{
    x <- .x %>% filter(InsertionGene == "Inserted") %>% pull(logFPKM)
    y <- .x %>% filter(InsertionGene == "Other")    %>% pull(logFPKM)
    p <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    tibble(p_value = p)
  }) %>% ungroup() %>%
  mutate(p_label = paste0("Wilcoxon p = ", formatC(p_value, format = "e", digits = 2)))
readr::write_tsv(wilcox_grouped, file.path(table_dir, "Insertion_vs_Other_Wilcoxon_byCohort.tsv"))

p_ins_box <- ggplot(expr_df_ins, aes(x = InsertionGene, y = logFPKM, fill = InsertionGene)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.95) +
  scale_fill_manual(values = c("Inserted" = "#4C6A87", "Other" = "#C0CEDD")) +
  facet_wrap(~ CohortIC, scales = "free_y") +
  labs(
    title = "Expression of genes with insertion vs all other genes (within cohort)",
    subtitle = "Inserted = union of genes listed in any sample’s Insertion field; Other = all remaining expressed genes.\nEach facet shows the distribution within that cohort.",
    x = "",
    y = "log2(FPKM+1)"
  ) +
  theme(legend.position = "none") +
  geom_text(
    data = wilcox_grouped,
    aes(x = 1.5, y = Inf, label = paste0("Wilcoxon p = ", formatC(p_value, format = "e", digits = 2))),
    vjust = 1.2, inherit.aes = FALSE, size = 3.6
  )

save_plot_all(p_ins_box, "Insertion_vs_Other_boxplots", width = 7.2, height = 4.8)

p_ins_box <- p_ins_box + geom_text(data = wilcox_grouped, aes(x = 1.5, y = Inf, label = p_label),
                                   vjust = 1.2, inherit.aes = FALSE, size = 3.6)
save_plot_all(p_ins_box, "Insertion_vs_Other_boxplots", width = 7.2, height = 4.8)

ins_map <- meta_ins %>%
  select(Sample_RNAseq, Cohort, Insertion_genes) %>%
  filter(lengths(Insertion_genes) > 0) %>%
  unnest(Insertion_genes) %>%
  rename(Gene = Insertion_genes) %>%
  filter(Gene %in% rownames(log_expr))
ic_samples <- meta_nd %>% filter(Cohort %in% c("IEC", "Control")) %>% pull(Sample_RNAseq)

pergene_stats <- list()
for(g in unique(ins_map$Gene)){
  inserted_samples <- ins_map %>% filter(Gene == g) %>% pull(Sample_RNAseq) %>% unique()
  inserted_ic <- intersect(inserted_samples, ic_samples)
  noninsert_ic <- setdiff(ic_samples, inserted_samples)
  x <- as.numeric(log_expr[g, inserted_ic, drop = TRUE])
  y <- as.numeric(log_expr[g, noninsert_ic, drop = TRUE])
  if(length(x) >= 1 & length(y) >= 2){
    p <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    median_diff <- median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    pergene_stats[[g]] <- tibble(Gene = g, n_inserted = length(x), n_noninsert = length(y),
                                 median_inserted = median(x, na.rm = TRUE),
                                 median_noninsert = median(y, na.rm = TRUE),
                                 median_diff = median_diff, p_value = p)
  } else {
    pergene_stats[[g]] <- tibble(Gene = g, n_inserted = length(x), n_noninsert = length(y),
                                 median_inserted = NA_real_, median_noninsert = NA_real_,
                                 median_diff = NA_real_, p_value = NA_real_)
  }
}
pergene_df <- bind_rows(pergene_stats) %>% mutate(p_adj = p.adjust(p_value, method = "BH")) %>% arrange(p_adj)
readr::write_tsv(pergene_df, file.path(table_dir, "Insertion_PerGene_Wilcoxon.tsv"))

build_gene_df <- function(g){
  inserted_samples <- ins_map %>% filter(Gene == g) %>% pull(Sample_RNAseq) %>% unique()
  group <- ifelse(colnames(log_expr) %in% inserted_samples, "Inserted", "Non-insert")
  tibble(Gene = g, Sample_RNAseq = colnames(log_expr), Group = group,
         logFPKM = as.numeric(log_expr[g, colnames(log_expr)])) %>%
    filter(Sample_RNAseq %in% ic_samples) %>%
    left_join(meta_nd %>% select(Sample_RNAseq, Cohort), by = "Sample_RNAseq")
}
plot_df <- purrr::map_dfr(unique(ins_map$Gene), build_gene_df)
if(nrow(plot_df) > 0){
  keep_genes <- pergene_df %>% filter(n_inserted >= 1, n_noninsert >= 2) %>% pull(Gene)
  plot_df2 <- plot_df %>% filter(Gene %in% keep_genes) %>%
    mutate(Gene = factor(Gene, levels = keep_genes))
  p_labels <- pergene_df %>% filter(Gene %in% keep_genes) %>%
    transmute(Gene, p_label = paste0("BH-FDR = ", ifelse(is.na(p_adj), "NA", formatC(p_adj, format = "e", digits = 2))))
  ymax_by_gene <- plot_df2 %>% group_by(Gene) %>% summarise(ymax = max(logFPKM, na.rm = TRUE), .groups = "drop")
  p_lab_df <- left_join(p_labels, ymax_by_gene, by = "Gene")

  p_pergene <- ggplot(plot_df2, aes(x = Group, y = logFPKM, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("Inserted" = "#4C6A87", "Non-insert" = "#C0CEDD")) +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(title = "Inserted vs Non-inserted (per gene) in IEC + Control", x = "", y = "log2(FPKM+1)") +
    theme(legend.position = "none") +
    geom_text(data = p_lab_df, aes(x = 1.5, y = ymax * 1.05, label = p_label),
              inherit.aes = FALSE, size = 3.2, vjust = 0)
  save_plot_all(p_pergene, "Insertion_PerGene_boxplots", width = 8.5, height = 7.5)
}

expr_ic_wide <- log_expr[, ic_samples, drop = FALSE]
sample_medians_ic <- apply(expr_ic_wide, 2, median, na.rm = TRUE)

pergene_vsother_stats <- list(); pergene_vsother_plot_list <- list()
for (g in all_inserted_genes) {
  if (!g %in% rownames(expr_ic_wide)) next
  x_vec <- as.numeric(expr_ic_wide[g, ])
  y_vec <- as.numeric(sample_medians_ic)
  if (sum(is.finite(x_vec) & is.finite(y_vec)) >= 3) {
    p <- tryCatch(wilcox.test(x_vec, y_vec, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
    med_diff <- median(x_vec - y_vec, na.rm = TRUE)
  } else { p <- NA_real_; med_diff <- NA_real_ }
  pergene_vsother_stats[[g]] <- tibble(Gene = g, n_samples = length(ic_samples),
                                       median_gene = median(x_vec, na.rm = TRUE),
                                       median_other = median(y_vec, na.rm = TRUE),
                                       median_diff = med_diff, p_value = p)
  pergene_vsother_plot_list[[g]] <- tibble(
    Sample_RNAseq = rep(ic_samples, times = 2),
    Group = factor(rep(c("InsertedGene", "OtherGenes"), each = length(ic_samples)),
                   levels = c("InsertedGene", "OtherGenes")),
    logFPKM = c(x_vec, y_vec),
    GeneFacet = g
  )
}
if (length(pergene_vsother_stats) > 0) {
  pergene_vsother_df <- bind_rows(pergene_vsother_stats) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>% arrange(p_adj)
  readr::write_tsv(pergene_vsother_df, file.path(table_dir, "Insertion_PerGene_vsOther_PairedWilcoxon.tsv"))
  plot_df_vsother <- bind_rows(pergene_vsother_plot_list)
  ord <- pergene_vsother_df$Gene
  plot_df_vsother$GeneFacet <- factor(plot_df_vsother$GeneFacet, levels = ord)
  p_labels_vsother <- pergene_vsother_df %>%
    transmute(GeneFacet = factor(Gene, levels = ord),
              p_label = paste0("BH-FDR = ", ifelse(is.na(p_adj), "NA", formatC(p_adj, format = "e", digits = 2))))
  ymax_vsother <- plot_df_vsother %>% group_by(GeneFacet) %>% summarise(ymax = max(logFPKM, na.rm = TRUE), .groups = "drop")
  p_lab_df_vsother <- left_join(p_labels_vsother, ymax_vsother, by = "GeneFacet")

  p_vsother <- ggplot(plot_df_vsother, aes(x = Group, y = logFPKM, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.65) +
    scale_fill_manual(values = c("InsertedGene" = "#4C6A87", "OtherGenes" = "#C0CEDD")) +
    facet_wrap(~ GeneFacet, scales = "free_y") +
    labs(title = "Per-gene: Inserted gene vs Other genes (IEC+Control; paired by sample)", x = "", y = "log2(FPKM+1)") +
    theme(legend.position = "none") +
    geom_text(data = p_lab_df_vsother, aes(x = 1.5, y = ymax * 1.05, label = p_label),
              inherit.aes = FALSE, size = 3.0, vjust = 0)
  save_plot_all(p_vsother, "Insertion_PerGene_vsOther_boxplots", width = 9.5, height = 8.5)
}

comp_build <- function(g){
  inserted_samples_for_g <- ins_map %>% filter(Gene == g) %>% pull(Sample_RNAseq) %>% unique()
  comparators <- setdiff(ic_samples, inserted_samples_for_g)
  if(length(comparators) == 0) return(NULL)
  df_nonins <- tibble(Gene = g, Sample_RNAseq = comparators,
                      logFPKM = as.numeric(log_expr[g, comparators]),
                      Set = "Non-insert")
  df_ins <- tibble(Gene = g, Sample_RNAseq = inserted_samples_for_g,
                   logFPKM = as.numeric(log_expr[g, inserted_samples_for_g]),
                   Set = "Inserted")
  bind_rows(df_nonins, df_ins)
}
comp_plot_df <- map_dfr(unique(ins_map$Gene), comp_build) %>%
  left_join(meta_nd %>% select(Sample_RNAseq, Cohort), by = "Sample_RNAseq")

if(nrow(comp_plot_df) > 0){
  p_comp <- ggplot(comp_plot_df, aes(x = Gene, y = logFPKM, color = Set)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.8, size = 1.8) +
    scale_color_manual(values = c("Inserted" = "#4C6A87", "Non-insert" = "#C0CEDD")) +
    labs(title = "Inserted sample vs comparator distribution (IEC+Control)", x = "Gene", y = "log2(FPKM+1)")
  save_plot_all(p_comp, "Insertion_vs_Comparator_jitter", width = 7.4, height = 4.8)
}

z_list <- list()
for(i in seq_len(nrow(ins_map))){
  smp <- ins_map$Sample_RNAseq[i]
  coh <- ins_map$Cohort[i]
  g   <- ins_map$Gene[i]

  x <- as.numeric(log_expr[g, smp])
  inserted_samples_for_g <- ins_map %>% filter(Gene == g) %>% pull(Sample_RNAseq) %>% unique()
  comparators <- setdiff(ic_samples, inserted_samples_for_g)

  if(length(comparators) >= 2){
    vec <- as.numeric(log_expr[g, comparators])
    mu  <- mean(vec, na.rm = TRUE)
    sdv <- stats::sd(vec, na.rm = TRUE)
    z   <- ifelse(is.finite(sdv) && sdv > 0, (x - mu)/sdv, NA_real_)
    z_list[[length(z_list)+1]] <- tibble(Sample_RNAseq = smp, Cohort = coh, Gene = g,
                                         Expression_logFPKM = x, Mean_nonInsert = mu, SD_nonInsert = sdv, Zscore = z)
  } else {
    z_list[[length(z_list)+1]] <- tibble(Sample_RNAseq = smp, Cohort = coh, Gene = g,
                                         Expression_logFPKM = x, Mean_nonInsert = NA_real_, SD_nonInsert = NA_real_, Zscore = NA_real_)
  }
}
z_table <- bind_rows(z_list) %>%
  left_join(meta_nd %>% select(Sample_RNAseq, source_id, Paper_ID, Tissue, Timepoint), by = "Sample_RNAseq") %>%
  arrange(Gene, desc(Zscore))
readr::write_tsv(z_table, file.path(table_dir, "InsertionImpact_Zscores.tsv"))

if(nrow(z_table) > 0){
  p_z <- ggplot(z_table, aes(x = reorder(paste0(Gene, " (", Sample_RNAseq, ")"), Zscore), y = Zscore, fill = Cohort)) +
    geom_col() +
    coord_flip() +
    scale_fill_muted() +
    labs(title = "Insertion impact (Z-scores) vs IEC+Control non-inserted peers",
         x = "Gene (Sample)", y = "Z-score")
  save_plot_all(p_z, "InsertionImpact_Zscores_bar", width = 7.8, height = 6.8)
}

if(length(cilta_idx) == 1){
  meta_corr <- meta_nd %>%
    filter(Cohort %in% c("IEC", "Control")) %>%
    mutate(Paper_ID = ifelse(is.na(Paper_ID), "NA", Paper_ID))
  cilta_vec <- as.numeric(log_expr[rownames(log_expr)[cilta_idx], meta_corr$Sample_RNAseq])
  df_corr_base <- meta_corr %>%
    transmute(Sample_RNAseq, Paper_ID, Cohort, ciltacel_logFPKM = cilta_vec,
              COMET_perc_CD3 = .data[["COMET_perc_CD3"]],
              COMET_perc_total = .data[["COMET_perc_total"]])

  compute_and_plot_corr <- function(df, yvar, file_stub){
    d <- df %>% select(Sample_RNAseq, Paper_ID, Cohort, ciltacel_logFPKM, !!sym(yvar)) %>%
      rename(Y = !!sym(yvar)) %>%
      filter(is.finite(ciltacel_logFPKM), is.finite(Y))
    if(nrow(d) < 3){
      warning("Not enough paired points for correlation: ", yvar)
      return(invisible(NULL))
    }
    pear  <- cor.test(d$ciltacel_logFPKM, d$Y, method = "pearson")
    spear <- cor.test(d$ciltacel_logFPKM, d$Y, method = "spearman", exact = FALSE)
    stats_tbl <- tibble(Metric = yvar, N = nrow(d),
                        Pearson_r = unname(pear$estimate), Pearson_p = pear$p.value,
                        Spearman_rho = unname(spear$estimate), Spearman_p = spear$p.value)
    readr::write_tsv(stats_tbl, file.path(table_dir, paste0("COR_stats_", file_stub, ".tsv")))

    n_pat <- d %>% distinct(Paper_ID) %>% nrow()
    pal_pat <- get_muted_palette(n_pat)
    ord_ids <- d %>% distinct(Paper_ID) %>% arrange(Paper_ID) %>% pull()
    names(pal_pat) <- ord_ids

    lab_txt <- paste0(
      "Pearson r = ", formatC(stats_tbl$Pearson_r, digits = 2, format = "f"),
      " (p = ", formatC(stats_tbl$Pearson_p, format = "e", digits = 2), ")\n",
      "Spearman \u03C1 = ", formatC(stats_tbl$Spearman_rho, digits = 2, format = "f"),
      " (p = ", formatC(stats_tbl$Spearman_p, format = "e", digits = 2), ")"
    )

    p <- ggplot(d, aes(x = ciltacel_logFPKM, y = Y, color = Paper_ID, shape = Cohort)) +
      geom_point(size = 2.6, alpha = 0.9) +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "grey20", linewidth = 0.7) +
      scale_color_manual(values = pal_pat) +
      scale_shape_manual(values = c("Control" = 16, "IEC" = 17)) +
      labs(title = paste0("ciltacel vs ", yvar, " (IEC & Control; no duplicates)"),
           x = "ciltacel log2(FPKM+1)", y = yvar, color = "Patient ID") +
      annotate("text", x = Inf, y = Inf, label = lab_txt, hjust = 1.02, vjust = 1.2, size = 3.4)

    save_plot_all(p, paste0("COR_ciltacel_", file_stub), width = 6.8, height = 5.4)
  }

  if("COMET_perc_CD3" %in% names(df_corr_base))   compute_and_plot_corr(df_corr_base, "COMET_perc_CD3",   "vs_COMET_perc_CD3")
  if("COMET_perc_total" %in% names(df_corr_base)) compute_and_plot_corr(df_corr_base, "COMET_perc_total", "vs_COMET_perc_total")
}

if (!exists("meta_ins")) {
  stopifnot(exists("meta_nd"))
  meta_ins <- meta_nd %>%
    dplyr::mutate(Insertion = ifelse(is.na(Insertion), "", Insertion)) %>%
    dplyr::mutate(Insertion_genes = strsplit(Insertion, ",")) %>%
    dplyr::mutate(Insertion_genes = lapply(Insertion_genes, function(x) {
      x <- trimws(x); x <- x[nchar(x) > 0]; unique(x)
    }))
}

if (!exists("all_inserted_genes")) {
  stopifnot(exists("log_expr"))
  all_inserted_genes <- unique(unlist(meta_ins$Insertion_genes))

  all_inserted_genes <- intersect(all_inserted_genes, rownames(log_expr))
}

if (!exists("expr_df")) {
  expr_df <- as.data.frame(log_expr) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Sample_RNAseq", values_to = "logFPKM") %>%
    dplyr::left_join(meta_nd, by = "Sample_RNAseq")
}

expr_df_ins <- expr_df %>%
  dplyr::mutate(InsertionGene = ifelse(Gene %in% all_inserted_genes, "Inserted", "Other")) %>%
  dplyr::mutate(CohortIC = dplyr::case_when(Cohort %in% c("IEC", "Control") ~ Cohort, TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(CohortIC))

ic_samples <- meta_nd %>%
  dplyr::filter(Cohort %in% c("IEC","Control")) %>%
  dplyr::pull(Sample_RNAseq) %>%
  unique()

med_by_sample <- expr_df_ins %>%
  dplyr::filter(Sample_RNAseq %in% ic_samples) %>%
  dplyr::group_by(CohortIC, Sample_RNAseq, InsertionGene) %>%
  dplyr::summarise(median_logFPKM = median(logFPKM, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = InsertionGene, values_from = median_logFPKM)

readr::write_tsv(med_by_sample, file.path(table_dir, "Insertion_MediansPerSample.tsv"))

paired_by_cohort <- med_by_sample %>%
  dplyr::group_by(CohortIC) %>%
  dplyr::group_modify(~{
    d <- .x %>% dplyr::filter(is.finite(Inserted), is.finite(Other))
    n <- nrow(d)
    if (n >= 2) {
      wt  <- wilcox.test(d$Inserted, d$Other, paired = TRUE, exact = FALSE,
                         conf.int = TRUE, conf.level = 0.95)
      wtg <- wilcox.test(d$Inserted, d$Other, paired = TRUE, exact = FALSE,
                         alternative = "greater")
      tibble::tibble(
        n_pairs     = n,
        W           = unname(wt$statistic),
        p_two_sided = wt$p.value,
        p_greater   = wtg$p.value,
        HL_est      = unname(wt$estimate),
        CI_lower    = wt$conf.int[1],
        CI_upper    = wt$conf.int[2],
        median_diff = stats::median(d$Inserted - d$Other, na.rm = TRUE)
      )
    } else {
      tibble::tibble(
        n_pairs = n, W = NA_real_, p_two_sided = NA_real_, p_greater = NA_real_,
        HL_est = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, median_diff = NA_real_
      )
    }
  }) %>%
  dplyr::ungroup()

readr::write_tsv(paired_by_cohort, file.path(table_dir, "Insertion_PairedWilcoxon_ByCohort.tsv"))

d_all <- med_by_sample %>% dplyr::filter(is.finite(Inserted), is.finite(Other))
if (nrow(d_all) >= 2) {
  wt_all <- wilcox.test(d_all$Inserted, d_all$Other, paired = TRUE, exact = FALSE,
                        conf.int = TRUE, conf.level = 0.95)
  wt_all_g <- wilcox.test(d_all$Inserted, d_all$Other, paired = TRUE, exact = FALSE,
                          alternative = "greater")
  overall_out <- tibble::tibble(
    n_pairs     = nrow(d_all),
    W           = unname(wt_all$statistic),
    p_two_sided = wt_all$p.value,
    p_greater   = wt_all_g$p.value,
    HL_est      = unname(wt_all$estimate),
    CI_lower    = wt_all$conf.int[1],
    CI_upper    = wt_all$conf.int[2],
    median_diff = stats::median(d_all$Inserted - d_all$Other, na.rm = TRUE)
  )
} else {
  overall_out <- tibble::tibble(
    n_pairs = nrow(d_all), W = NA_real_, p_two_sided = NA_real_, p_greater = NA_real_,
    HL_est = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, median_diff = NA_real_
  )
}
readr::write_tsv(overall_out, file.path(table_dir, "Insertion_PairedWilcoxon_Overall.tsv"))

plot_ms <- med_by_sample %>%
  tidyr::pivot_longer(c(Inserted, Other), names_to = "Group", values_to = "Median") %>%
  dplyr::mutate(Group = factor(Group, levels = c("Other","Inserted")))

lab_df <- paired_by_cohort %>%
  dplyr::mutate(
    label = paste0(
      "n = ", n_pairs,
      "\nWilcoxon W = ", ifelse(is.na(W), "NA", formatC(W, digits = 3, format = "f")),
      "\np (two-sided) = ", ifelse(is.na(p_two_sided), "NA", formatC(p_two_sided, format = "e", digits = 2)),
      "\nHL = ", ifelse(is.na(HL_est), "NA", sprintf("%.2f", HL_est)),
      " [", ifelse(is.na(CI_lower), "NA", sprintf("%.2f", CI_lower)),
      ", ", ifelse(is.na(CI_upper), "NA", sprintf("%.2f", CI_upper)), "]"
    )
  )

p_paired_medians <- ggplot(plot_ms, aes(x = Group, y = Median, group = Sample_RNAseq)) +
  geom_line(alpha = 0.35, color = "grey50") +
  geom_point(aes(color = Group), size = 2.6, alpha = 0.95) +
  scale_color_manual(values = c("Inserted" = "#4C6A87", "Other" = "#C0CEDD")) +
  facet_wrap(~ CohortIC, scales = "free_y") +
  labs(title = "Per-sample medians: Inserted genes vs Other genes",
       x = NULL, y = "Median log2(FPKM+1)", color = NULL) +
  theme(legend.position = "bottom") +

  geom_text(
    data = lab_df,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.2, size = 3.2
  )

save_plot_all(p_paired_medians, "Insertion_Medians_Paired_byCohort", width = 7.6, height = 5.4)
