#!/usr/bin/env Rscript
# CAR and clonotype annotation
# Annotate CAR UMIs and summarize all-clonotyped versus CAR-positive repertoires by sample and tissue.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

msg  <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
warn <- function(...) cat("[WARN] ", sprintf(...), "\n", sep = "")
die  <- function(...) { cat("[ERROR] ", sprintf(...), "\n", sep = ""); quit(status = 2) }

QC_DIR <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "single_cell", "clonotype_annotation"))
INPUT_RDS <- Sys.getenv("INPUT_RDS", unset = file.path(QC_DIR, "single_cell_object.rds"))
CAR_GENE      <- Sys.getenv("CAR_GENE", unset = "CILTACELCAR")
CAR_MIN_UMI   <- as.numeric(Sys.getenv("CAR_MIN_UMI", unset = "1"))
CLONE_ID_MODE <- Sys.getenv("CLONE_ID_MODE", unset = "raw")
INCLUDE_REGEX <- Sys.getenv("INCLUDE_REGEX", unset = "(_IL_|_BL_)")

if (OUT_DIR == "") die("OUT_DIR is empty")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

msg("QC_DIR        = %s", QC_DIR)
msg("OUT_DIR       = %s", OUT_DIR)
msg("INPUT_RDS     = %s", INPUT_RDS)
msg("CAR_GENE      = %s (CAR_MIN_UMI=%s)", CAR_GENE, CAR_MIN_UMI)
msg("CLONE_ID_MODE = %s", CLONE_ID_MODE)

if (!file.exists(INPUT_RDS)) die("Missing INPUT_RDS: %s", INPUT_RDS)

get_layers <- function(obj, assay = "RNA") {
  if (!assay %in% names(obj@assays)) return(character(0))
  ass <- obj[[assay]]
  tryCatch(SeuratObject::Layers(ass), error = function(e) character(0))
}

get_gene_umi_from_layers <- function(obj, assay = "RNA", gene, prefer_prefix = "counts") {
  if (!assay %in% names(obj@assays)) {
    die("Assay '%s' not present. Available: %s", assay, paste(names(obj@assays), collapse = ", "))
  }
  all_cells <- colnames(obj)
  out <- setNames(rep(0, length(all_cells)), all_cells)

  layers <- get_layers(obj, assay = assay)
  if (length(layers) == 0) {

    mat <- tryCatch({
      Seurat::GetAssayData(obj, assay = assay, layer = "counts")
    }, error = function(e) {
      die("Failed GetAssayData(assay=%s, layer=counts): %s", assay, conditionMessage(e))
    })
    if (!gene %in% rownames(mat)) return(out)
    v <- Matrix::colSums(mat[gene, , drop = FALSE])
    out[names(v)] <- as.numeric(v)
    return(out)
  }

  layers_use <- layers[grepl(prefer_prefix, layers, fixed = TRUE)]
  if (length(layers_use) == 0) layers_use <- layers

  ass <- obj[[assay]]
  for (ly in layers_use) {
    mat <- tryCatch({
      SeuratObject::LayerData(ass, layer = ly)
    }, error = function(e) {
      warn("LayerData failed for layer=%s: %s", ly, conditionMessage(e))
      NULL
    })
    if (is.null(mat)) next
    if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
    if (!gene %in% rownames(mat)) next

    cols <- colnames(mat)
    v <- Matrix::colSums(mat[gene, , drop = FALSE])
    out[cols] <- out[cols] + as.numeric(v)
  }
  out
}

calc_diversity <- function(freqs) {
  freqs <- freqs[is.finite(freqs) & freqs > 0]
  if (length(freqs) == 0) {
    return(list(n_clones=0L, shannon=NA_real_, simpson=NA_real_, inv_simpson=NA_real_, pielou=NA_real_, clonality=NA_real_))
  }
  n_clones <- length(freqs)
  shannon <- -sum(freqs * log(freqs))
  simpson <- 1 - sum(freqs^2)
  inv_simpson <- 1 / sum(freqs^2)
  pielou <- if (n_clones > 1) shannon / log(n_clones) else NA_real_
  clonality <- if (n_clones > 1) 1 - (shannon / log(n_clones)) else NA_real_
  list(n_clones=as.integer(n_clones), shannon=shannon, simpson=simpson, inv_simpson=inv_simpson, pielou=pielou, clonality=clonality)
}

obj <- readRDS(INPUT_RDS)
msg("Loaded object: %d cells", ncol(obj))

md <- obj@meta.data
md$cell_barcode <- rownames(md)

if (!"sample_id" %in% colnames(md)) {
  if ("orig.ident" %in% colnames(md)) {
    md$sample_id <- as.character(md$orig.ident)
    warn("sample_id missing; using orig.ident")
  } else {
    die("No sample_id or orig.ident found in meta.data")
  }
}

if (nchar(INCLUDE_REGEX) > 0) {
  keep <- grepl(INCLUDE_REGEX, md$sample_id)
  if (any(!keep)) {
    warn("Dropping %d cells whose sample_id does not match INCLUDE_REGEX=%s", sum(!keep), INCLUDE_REGEX)
    md <- md[keep, , drop = FALSE]
  }
}

md$tissue <- NA_character_
md$tissue[grepl("_IL_", md$sample_id)] <- "IL"
md$tissue[grepl("_BL_", md$sample_id)] <- "BL"

msg("Extracting CAR gene UMIs for %s from RNA layers...", CAR_GENE)
car_umi <- get_gene_umi_from_layers(obj, assay = "RNA", gene = CAR_GENE, prefer_prefix = "counts")
car_umi <- car_umi[md$cell_barcode]

md$CAR_gene <- CAR_GENE
md$CAR_umi  <- as.numeric(car_umi)
md$CAR_pos  <- md$CAR_umi >= CAR_MIN_UMI

msg("CAR+ cells: %d / %d (>= %s UMI)", sum(md$CAR_pos, na.rm = TRUE), nrow(md), CAR_MIN_UMI)

needed <- c("tcr_raw_clonotype_id", "tcr_cdr3_aa", "tcr_clone_n")
present <- intersect(needed, colnames(md))
msg("TCR columns present: %s", paste(present, collapse = ", "))

if (!"tcr_raw_clonotype_id" %in% colnames(md) && !"tcr_cdr3_aa" %in% colnames(md)) {
  die("No tcr_raw_clonotype_id or tcr_cdr3_aa found. This RDS does not include attached TCR annotations.")
}

if (tolower(CLONE_ID_MODE) == "cdr3") {
  md$CTstrict <- as.character(md$tcr_cdr3_aa)
  md$clone_id_source <- "tcr_cdr3_aa"
} else {
  md$CTstrict <- as.character(md$tcr_raw_clonotype_id)
  md$clone_id_source <- "tcr_raw_clonotype_id"
}

md$has_clone <- !is.na(md$CTstrict) & md$CTstrict != "" & md$CTstrict != "NA"
msg("Cells with clonotype assigned (CTstrict): %d / %d", sum(md$has_clone), nrow(md))

compute_clone_tables <- function(df, subset_label) {
  df <- df %>% mutate(subset = subset_label)

  denom <- df %>%
    filter(has_clone) %>%
    group_by(sample_id) %>%
    summarise(n_clonotyped = dplyr::n(), .groups = "drop")

  per_clone <- df %>%
    filter(has_clone) %>%
    group_by(sample_id, tissue, CTstrict) %>%
    summarise(clone_size = dplyr::n(), .groups = "drop") %>%
    left_join(denom, by = "sample_id") %>%
    mutate(clone_freq = ifelse(n_clonotyped > 0, clone_size / n_clonotyped, NA_real_),
           subset = subset_label)

  per_cell <- df %>%
    left_join(
      per_clone %>% select(sample_id, CTstrict, clone_size, clone_freq, subset),
      by = c("sample_id", "CTstrict", "subset")
    )

  list(per_clone = per_clone, per_cell = per_cell)
}

res_all <- compute_clone_tables(md, subset_label = "ALL")
res_car <- compute_clone_tables(md %>% filter(CAR_pos), subset_label = sprintf("CARpos_%s_umi_ge_%s", CAR_GENE, CAR_MIN_UMI))

div_by_sample <- function(per_clone_df) {
  per_clone_df %>%
    group_by(sample_id, tissue, subset) %>%
    summarise(
      n_cells_clonotyped = sum(clone_size),
      n_clones = dplyr::n(),
      shannon = calc_diversity(clone_freq)$shannon,
      simpson = calc_diversity(clone_freq)$simpson,
      inv_simpson = calc_diversity(clone_freq)$inv_simpson,
      pielou = calc_diversity(clone_freq)$pielou,
      clonality = calc_diversity(clone_freq)$clonality,
      .groups = "drop"
    )
}

div_sample <- bind_rows(
  div_by_sample(res_all$per_clone),
  div_by_sample(res_car$per_clone)
)

div_by_tissue <- function(df_cells, subset_label) {
  df_cells <- df_cells %>% filter(has_clone)
  if (nrow(df_cells) == 0) return(tibble())

  per_clone <- df_cells %>%
    group_by(tissue, CTstrict) %>%
    summarise(clone_size = dplyr::n(), .groups = "drop") %>%
    group_by(tissue) %>%
    mutate(clone_freq = clone_size / sum(clone_size)) %>%
    ungroup() %>%
    mutate(subset = subset_label)

  per_clone %>%
    group_by(tissue, subset) %>%
    summarise(
      n_cells_clonotyped = sum(clone_size),
      n_clones = dplyr::n(),
      shannon = calc_diversity(clone_freq)$shannon,
      simpson = calc_diversity(clone_freq)$simpson,
      inv_simpson = calc_diversity(clone_freq)$inv_simpson,
      pielou = calc_diversity(clone_freq)$pielou,
      clonality = calc_diversity(clone_freq)$clonality,
      .groups = "drop"
    )
}

div_tissue <- bind_rows(
  div_by_tissue(md, "ALL"),
  div_by_tissue(md %>% filter(CAR_pos), sprintf("CARpos_%s_umi_ge_%s", CAR_GENE, CAR_MIN_UMI))
)

topN <- 25
top_clones_car <- res_car$per_clone %>%
  arrange(sample_id, desc(clone_size)) %>%
  group_by(sample_id) %>%
  slice_head(n = topN) %>%
  ungroup()

gz_tsv <- function(df, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(div_sample, file.path(OUT_DIR, "car_clonality_summary_by_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(div_tissue, file.path(OUT_DIR, "car_clonality_summary_by_tissue.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(top_clones_car, file.path(OUT_DIR, "car_clonality_top_clones_by_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

per_cell_out <- bind_rows(
  res_all$per_cell %>% mutate(context = "ALL"),
  res_car$per_cell %>% mutate(context = "CARpos")
) %>%
  select(cell_barcode, sample_id, tissue, CAR_gene, CAR_umi, CAR_pos,
         clone_id_source, CTstrict, has_clone, clone_size, clone_freq, context, subset)

gz_tsv(per_cell_out, file.path(OUT_DIR, "car_clonality_per_cell.tsv.gz"))

obj2 <- obj

obj2 <- subset(obj2, cells = md$cell_barcode)
obj2@meta.data <- md
saveRDS(obj2, file.path(OUT_DIR, "car_clonality_annotated_object.rds"))

sink(file.path(OUT_DIR, "sessionInfo.car_clonality.txt"))
print(sessionInfo())
sink()

msg("Wrote outputs to %s", OUT_DIR)

pdf(file.path(OUT_DIR, "car_clonality_diversity_plots.pdf"), width = 11, height = 8.5)

p1 <- div_sample %>%
  mutate(is_car = grepl("^CARpos_", subset)) %>%
  ggplot(aes(x = sample_id, y = clonality, shape = tissue)) +
  geom_point(size = 2) +
  facet_wrap(~ is_car, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Clonality by sample", y = "Clonality (1 - normalized Shannon)", x = "Sample")
print(p1)

p2 <- div_sample %>%
  mutate(is_car = grepl("^CARpos_", subset)) %>%
  ggplot(aes(x = sample_id, y = shannon, shape = tissue)) +
  geom_point(size = 2) +
  facet_wrap(~ is_car, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Shannon diversity by sample", y = "Shannon", x = "Sample")
print(p2)

p3 <- div_tissue %>%
  mutate(is_car = grepl("^CARpos_", subset)) %>%
  ggplot(aes(x = tissue, y = clonality)) +
  geom_point(size = 3) +
  facet_wrap(~ is_car, scales = "free_y") +
  theme_bw() +
  labs(title = "Clonality by tissue (pooled)", y = "Clonality", x = "Tissue")
print(p3)

dev.off()

make_rank_df <- function(per_clone_df) {
  per_clone_df %>%
    group_by(sample_id, tissue) %>%
    arrange(desc(clone_freq), .by_group = TRUE) %>%
    mutate(rank = row_number(), cum_freq = cumsum(clone_freq)) %>%
    ungroup()
}
rank_car <- make_rank_df(res_car$per_clone)

pdf(file.path(OUT_DIR, "car_clonality_clone_rank_plots.pdf"), width = 11, height = 8.5)

p_rank <- rank_car %>%
  ggplot(aes(x = rank, y = cum_freq, group = sample_id)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tissue) +
  theme_bw() +
  labs(title = sprintf("CAR+ clone rank curves (%s >= %s UMI)", CAR_GENE, CAR_MIN_UMI),
       x = "Clone rank", y = "Cumulative frequency")
print(p_rank)

dev.off()

msg("Plots written. DONE.")
