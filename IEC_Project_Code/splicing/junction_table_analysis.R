#!/usr/bin/env Rscript
# Junction-table splice analysis
# Evaluate gene-level splicing, insertion-proximal junctions, and transcript-aware annotations.

args <- commandArgs(trailingOnly = TRUE)
default_base_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
base_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_base_dir
base_dir <- normalizePath(path.expand(base_dir), mustWork = TRUE)

meta_path <- file.path(base_dir, "splicing_metadata.tsv")
insertion_path <- file.path(base_dir, "insertion_coordinates.tsv")
gtf_path <- file.path(base_dir, "GRCh38.gtf.gz")
chain_path <- file.path(base_dir, "hg19ToHg38.over.chain.gz")

# Coordinate conventions
coord_source <- "metadata"
coord_build  <- "hg38"

junction_pattern <- "\\.junction\\.txt$"

base_font_size    <- 14
plot_title_size   <- 16
axis_title_size   <- 14
axis_text_size    <- 12
legend_title_size <- 13
legend_text_size  <- 12
strip_text_size   <- 13
heatmap_text_size <- 12

out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "splicing", "junction_table_analysis")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
log_dir   <- file.path(out_dir, "logs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Junction inclusion and comparison parameters
local_pad_bp            <- 1000
min_junction_reads      <- 1
min_gene_total_reads    <- 5
min_controls_per_gene   <- 2
label_top_junc_n        <- 12
use_gene_body_window    <- TRUE
absent_ref_detect_min_n    <- 2
absent_ref_detect_min_rate <- 0.25

prefer_metadata_controls <- TRUE

meta_sample_col   <- "sample"
meta_patient_col  <- "patient_id"
meta_group_col    <- "group"
meta_gene_col     <- "gene"
meta_chr_col      <- "chr"
meta_start_col    <- "start"
meta_end_col      <- "end"

ins_sample_col    <- "sample"
ins_gene_col      <- "gene"
ins_chr_col       <- "chrom"
ins_start_col     <- "start"
ins_end_col       <- "end"
ins_pos_col       <- "position"

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(pheatmap)
})

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}
base_font <- "Helvetica"

theme_set(
  theme_classic(base_size = base_font_size, base_family = base_font) +
    theme(
      plot.title  = element_text(face = "bold", size = plot_title_size),
      legend.position = "right",
      axis.text   = element_text(color = "gray20", size = axis_text_size),
      axis.title  = element_text(color = "gray20", size = axis_title_size),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text  = element_text(face = "bold", size = strip_text_size)
    )
)

foggy_sf <- c("#4C6A87", "#7B99B6", "#9CB7CE", "#C0CEDD",
              "#9AA6B2", "#72808E", "#B8A9B4", "#C9BCC6", "#A6B8BE")
get_muted_palette <- function(n){
  if (n <= length(foggy_sf)) return(foggy_sf[seq_len(n)])
  colorRampPalette(foggy_sf)(n)
}

pdf_device_live_text <- function(filename, width = 6, height = 4, ...) {
  grDevices::pdf(
    file = filename, width = width, height = height,
    family = base_font, useDingbats = FALSE,
    version = "1.4", colormodel = "srgb", ...
  )
}
eps_device_live_text <- function(filename, width = 6, height = 4, ...) {
  grDevices::postscript(
    file = filename, width = width, height = height,
    onefile = FALSE, horizontal = FALSE, paper = "special",
    family = base_font, colormodel = "srgb", ...
  )
}
render_to_device <- function(open_device, draw) {
  open_device()
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(NULL)
}
save_plot_all <- function(plot, file_base, width = 6, height = 4) {
  pdf_file <- file.path(plot_dir, paste0(file_base, ".pdf"))
  render_to_device(
    function() pdf_device_live_text(pdf_file, width = width, height = height),
    function() print(plot)
  )

  eps_file <- file.path(plot_dir, paste0(file_base, ".eps"))
  render_to_device(
    function() eps_device_live_text(eps_file, width = width, height = height),
    function() print(plot)
  )
}
save_pheatmap_all <- function(pheat, file_base, width = 7, height = 6) {
  pdf_file <- file.path(plot_dir, paste0(file_base, ".pdf"))
  render_to_device(
    function() pdf_device_live_text(pdf_file, width = width, height = height),
    function() {
      grid::grid.newpage()
      print(pheat)
    }
  )

  eps_file <- file.path(plot_dir, paste0(file_base, ".eps"))
  render_to_device(
    function() eps_device_live_text(eps_file, width = width, height = height),
    function() {
      grid::grid.newpage()
      print(pheat)
    }
  )
}

msg <- function(...) cat(paste0(..., "\n"))

read_delim_flex <- function(path) {
  data.table::fread(path, sep = "\t", header = TRUE, data.table = TRUE)
}

norm_names <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
norm_chr   <- function(x) gsub("^chr", "", as.character(x), ignore.case = TRUE)
safe_num   <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Generic sample matching
canonical_sample_core <- function(x) {
  s <- basename(as.character(x))
  s <- sub("(?i)(\\.aligned\\.sortedbycoord\\.out)?\\.bam$", "", s, perl = TRUE)
  s <- sub("(?i)\\.junction\\.txt$", "", s, perl = TRUE)
  toupper(gsub("[^A-Za-z0-9]+", "_", s))
}

empirical_p_2sided <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  if (!is.finite(x) || length(ref) < 2) return(NA_real_)
  med <- median(ref, na.rm = TRUE)
  more_extreme <- sum(abs(ref - med) >= abs(x - med), na.rm = TRUE)
  (more_extreme + 1) / (length(ref) + 1)
}
make_z <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  if (!is.finite(x) || length(ref) < 2) return(NA_real_)
  mu <- mean(ref, na.rm = TRUE); sdv <- stats::sd(ref, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(NA_real_)
  (x - mu)/sdv
}

msg("Loading metadata...")
meta <- as_tibble(read_delim_flex(meta_path))
names(meta) <- norm_names(names(meta))

needed_meta <- c(meta_sample_col, meta_group_col, meta_gene_col)
missing_meta <- base::setdiff(needed_meta, names(meta))
if (length(missing_meta) > 0) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))

meta <- meta %>%
  mutate(
    sample_raw = .data[[meta_sample_col]],
    sample_core = canonical_sample_core(.data[[meta_sample_col]]),
    group = as.character(.data[[meta_group_col]]),
    patient_id = if (meta_patient_col %in% names(meta)) as.character(.data[[meta_patient_col]]) else NA_character_,
    gene_raw = as.character(.data[[meta_gene_col]]),
    chrom_in = if (meta_chr_col %in% names(meta)) norm_chr(.data[[meta_chr_col]]) else NA_character_,
    start_in = if (meta_start_col %in% names(meta)) safe_num(.data[[meta_start_col]]) else NA_real_,
    end_in   = if (meta_end_col %in% names(meta)) safe_num(.data[[meta_end_col]]) else NA_real_
  ) %>%
  mutate(gene_raw = ifelse(gene_raw %in% c("NA", "", "Na", "na"), NA_character_, gene_raw))

meta_expanded <- meta %>%
  mutate(gene = strsplit(ifelse(is.na(gene_raw), "", gene_raw), ",")) %>%
  unnest(gene, keep_empty = TRUE) %>%
  mutate(gene = trimws(gene), gene = na_if(gene, ""))

write_tsv(meta_expanded, file.path(table_dir, "metadata_expanded_gene_rows.tsv"))

ins_tbl <- tibble()
if (coord_source == "insertions") {
  msg("Loading insertion table...")
  ins_tbl <- as_tibble(read_delim_flex(insertion_path))
  names(ins_tbl) <- norm_names(names(ins_tbl))

  if (!(ins_sample_col %in% names(ins_tbl) && ins_gene_col %in% names(ins_tbl) && ins_chr_col %in% names(ins_tbl))) {
    stop("Insertion file is missing required configured columns.")
  }

  ins_tbl <- ins_tbl %>%
    mutate(
      sample_raw = .data[[ins_sample_col]],
      sample_core = canonical_sample_core(.data[[ins_sample_col]]),
      gene_raw = as.character(.data[[ins_gene_col]]),
      chrom_in = norm_chr(.data[[ins_chr_col]])
    )

  if (ins_start_col %in% names(ins_tbl) && ins_end_col %in% names(ins_tbl)) {
    ins_tbl <- ins_tbl %>% mutate(start_in = safe_num(.data[[ins_start_col]]), end_in = safe_num(.data[[ins_end_col]]))
  } else if (ins_pos_col %in% names(ins_tbl)) {
    ins_tbl <- ins_tbl %>% mutate(start_in = safe_num(.data[[ins_pos_col]]), end_in = safe_num(.data[[ins_pos_col]]))
  } else {
    stop("Insertion file must have either start/end columns or a single position column.")
  }

  ins_tbl <- ins_tbl %>%
    mutate(gene = strsplit(ifelse(is.na(gene_raw), "", gene_raw), ",")) %>%
    unnest(gene, keep_empty = TRUE) %>%
    mutate(gene = trimws(gene), gene = na_if(gene, "")) %>%
    dplyr::select(sample_raw, sample_core, gene, chrom_in, start_in, end_in)

  write_tsv(ins_tbl, file.path(table_dir, "insertions_expanded_gene_rows.tsv"))
}

msg("Loading junction files...")
junction_files <- list.files(path = base_dir, pattern = junction_pattern, full.names = TRUE)
if (length(junction_files) == 0) stop("No junction files found with pattern: ", junction_pattern)

jfile_tbl <- tibble(
  file = junction_files,
  file_base = basename(junction_files),
  sample_core = canonical_sample_core(basename(junction_files))
)

match_audit <- jfile_tbl %>%
  left_join(meta_expanded %>% distinct(sample_raw, sample_core), by = "sample_core") %>%
  mutate(match_status = ifelse(is.na(sample_raw), "unmatched", "matched"))
write_tsv(match_audit, file.path(table_dir, "sample_matching_audit_junction_files.tsv"))

msg("\n=== Junction file matching audit ===")
print(dplyr::count(match_audit, match_status))
unmatched_files <- match_audit %>% dplyr::filter(match_status == "unmatched")
if (nrow(unmatched_files) > 0) {
  msg("Warning: unmatched junction files detected. These files will still be read, but group/coordinate annotation may be missing.")
  print(unmatched_files)
}

coord_tbl <- tibble()
if (coord_source == "metadata") {
  coord_tbl <- meta_expanded %>%
    transmute(sample_raw, sample_core, group, patient_id, gene, chrom_in, start_in, end_in)
} else {
  coord_tbl <- meta_expanded %>%
    dplyr::select(sample_raw, sample_core, group, patient_id) %>% distinct() %>%
    left_join(ins_tbl, by = c("sample_raw", "sample_core"))
}

if (coord_build == "hg19") {
  if (!file.exists(chain_path)) stop("coord_build='hg19' but chain file not found: ", chain_path)
  msg("Running liftover hg19 -> hg38...")
  chain <- import.chain(chain_path)
  gr_hg19 <- GRanges(
    seqnames = paste0("chr", norm_chr(coord_tbl$chrom_in)),
    ranges = IRanges(coord_tbl$start_in, coord_tbl$end_in),
    sample_raw = coord_tbl$sample_raw, sample_core = coord_tbl$sample_core,
    group = coord_tbl$group, patient_id = coord_tbl$patient_id, gene = coord_tbl$gene
  )
  lo <- liftOver(gr_hg19, chain)
  lo_n <- lengths(lo)
  lo_first <- unlist(lo[lo_n >= 1], use.names = FALSE)
  coord_tbl <- tibble(
    sample_raw = mcols(lo_first)$sample_raw,
    sample_core = mcols(lo_first)$sample_core,
    group = mcols(lo_first)$group,
    patient_id = mcols(lo_first)$patient_id,
    gene = mcols(lo_first)$gene,
    chrom = norm_chr(as.character(seqnames(lo_first))),
    start = start(lo_first),
    end = end(lo_first),
    liftover_n = as.integer(lo_n[lo_n >= 1])
  )
} else {
  coord_tbl <- coord_tbl %>% transmute(
    sample_raw, sample_core, group, patient_id, gene,
    chrom = norm_chr(chrom_in), start = start_in, end = end_in,
    liftover_n = NA_integer_
  )
}

coord_tbl <- coord_tbl %>%
  mutate(has_interval = !is.na(chrom) & is.finite(start) & is.finite(end)) %>%
  mutate(start2 = ifelse(has_interval, pmin(start, end), start),
         end2   = ifelse(has_interval, pmax(start, end), end),
         start = start2, end = end2) %>%
  dplyr::select(-start2, -end2)
write_tsv(coord_tbl, file.path(table_dir, "target_coordinates_used.tsv"))

msg("\n=== Coordinate availability summary ===")
print(dplyr::count(coord_tbl, has_interval))

# Junction-table import
read_junction_file <- function(file) {
  dt <- fread(file, sep = "\t", header = TRUE, data.table = TRUE)
  names(dt) <- norm_names(names(dt))
  pick <- function(cands) { hit <- base::intersect(cands, names(dt)); if (length(hit)) hit[1] else NA_character_ }
  chrom_col <- pick(c("chrom", "chr", "chromosome"))
  start_col <- pick(c("intron_st_0_based", "intron_st_0_based_", "intron_st_0_b", "intron_st", "start"))
  end_col   <- pick(c("intron_end_1_based", "intron_end_1_based_", "intron_end_1_b", "intron_end", "end"))
  count_col <- pick(c("read_count", "count", "reads", "nreads"))
  anno_col  <- pick(c("annotation", "annot"))
  req <- c(chrom_col, start_col, end_col, count_col)
  if (any(is.na(req))) stop("Could not parse junction columns in: ", basename(file))

  tibble(
    file = basename(file),
    sample_core = canonical_sample_core(basename(file)),
    chrom = norm_chr(dt[[chrom_col]]),
    start = safe_num(dt[[start_col]]),
    end   = safe_num(dt[[end_col]]),
    count = safe_num(dt[[count_col]]),
    annotation = if (!is.na(anno_col)) trimws(as.character(dt[[anno_col]])) else NA_character_
  ) %>%
    dplyr::filter(is.finite(start), is.finite(end), is.finite(count), count >= min_junction_reads) %>%
    mutate(junction_id = paste0(chrom, ":", start, ":", end))
}

junctions <- map_dfr(junction_files, read_junction_file)
write_tsv(junctions, file.path(table_dir, "all_junction_rows_raw.tsv"))

junctions <- junctions %>%
  left_join(meta_expanded %>% distinct(sample_core, sample_raw, group, patient_id), by = "sample_core") %>%
  mutate(sample = coalesce(sample_raw, sample_core))

msg("Importing GTF and mapping junctions to genes...")
gtf <- import(gtf_path)
genes <- gtf[gtf$type == "gene"]
if (!"gene_name" %in% names(mcols(genes))) stop("GTF lacks gene_name.")
genes <- genes[!is.na(mcols(genes)$gene_name)]
seqlevelsStyle(genes) <- "UCSC"

genes_df <- tibble(
  gene = mcols(genes)$gene_name,
  chrom = norm_chr(as.character(seqnames(genes))),
  gene_start = start(genes),
  gene_end   = end(genes)
) %>% distinct(gene, chrom, gene_start, gene_end)
write_tsv(genes_df, file.path(table_dir, "genes_from_gtf.tsv"))

gj <- GRanges(
  seqnames = paste0("chr", junctions$chrom),
  ranges = IRanges(start = junctions$start + 1, end = junctions$end),
  junction_id = junctions$junction_id
)
seqlevelsStyle(gj) <- "UCSC"
hits <- findOverlaps(gj, genes, ignore.strand = TRUE)
map_df <- tibble(
  junction_id = mcols(gj)$junction_id[queryHits(hits)],
  gene = mcols(genes)$gene_name[subjectHits(hits)]
) %>% distinct(junction_id, gene)

junctions_g <- junctions %>% inner_join(map_df, by = "junction_id") %>%
  left_join(genes_df, by = c("gene", "chrom"))
write_tsv(junctions_g, file.path(table_dir, "junction_rows_gene_mapped.tsv"))

msg("Mapped junction rows: ", nrow(junctions_g), "; unique genes: ", n_distinct(junctions_g$gene))

analysis_rows <- junctions_g %>%
  left_join(coord_tbl %>% dplyr::select(sample_core, gene, group, patient_id, chrom_target = chrom, start_target = start, end_target = end, has_interval),
            by = c("sample_core", "gene")) %>%
  mutate(
    is_target_sample_gene = !is.na(group.y) | !is.na(group.x),
    group = coalesce(group.y, group.x),
    patient_id = coalesce(patient_id.y, patient_id.x),
    chrom_target = ifelse(is.na(chrom_target), chrom, chrom_target),
    has_interval = ifelse(is.na(has_interval), FALSE, has_interval),
    junction_len = end - start + 1,
    spans_target_interval = ifelse(has_interval & chrom == chrom_target, start <= end_target & end >= start_target, FALSE),
    donor_near_target    = ifelse(has_interval & chrom == chrom_target, abs(start - start_target) <= local_pad_bp | abs(start - end_target) <= local_pad_bp, FALSE),
    acceptor_near_target = ifelse(has_interval & chrom == chrom_target, abs(end - start_target) <= local_pad_bp | abs(end - end_target) <= local_pad_bp, FALSE),
    near_target = spans_target_interval | donor_near_target | acceptor_near_target,
    annotation_simple = case_when(
      is.na(annotation) ~ "unknown",
      str_detect(tolower(annotation), "annotated") ~ "annotated",
      str_detect(tolower(annotation), "novel") ~ "novel",
      TRUE ~ tolower(annotation)
    )
  )

if (use_gene_body_window) {
  analysis_rows <- analysis_rows %>%
    mutate(in_gene_body = (start <= gene_end & end >= gene_start)) %>%
    dplyr::filter(in_gene_body)
}

usage_rows <- analysis_rows %>%
  group_by(sample, sample_core, group, patient_id, gene) %>%
  mutate(gene_total_reads = sum(count, na.rm = TRUE), usage = count / pmax(gene_total_reads, 1)) %>%
  ungroup()

sample_gene <- usage_rows %>%
  group_by(sample, sample_core, group, patient_id, gene) %>%
  summarise(
    gene_total_reads = sum(count, na.rm = TRUE),
    n_junctions = n_distinct(junction_id),
    annotated_reads = sum(count[annotation_simple == "annotated"], na.rm = TRUE),
    novel_reads     = sum(count[annotation_simple == "novel"], na.rm = TRUE),
    novel_fraction  = novel_reads / pmax(gene_total_reads, 1),
    local_reads     = sum(count[near_target], na.rm = TRUE),
    local_fraction  = local_reads / pmax(gene_total_reads, 1),
    local_novel_reads = sum(count[near_target & annotation_simple == "novel"], na.rm = TRUE),
    local_novel_fraction = local_novel_reads / pmax(local_reads, 1),
    entropy = -sum(usage * log(usage + 1e-6), na.rm = TRUE),
    has_interval = any(has_interval),
    start_target = dplyr::first(
      start_target[has_interval & !is.na(start_target)],
      default = NA_real_
    ),
    end_target = dplyr::first(
      end_target[has_interval & !is.na(end_target)],
      default = NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(gene_total_reads >= min_gene_total_reads)

write_tsv(sample_gene, file.path(table_dir, "sample_gene_metrics.tsv"))

msg("\n=== sample_gene summary ===")
print(sample_gene %>% summarise(
  n_rows = n(),
  n_genes = n_distinct(gene),
  n_samples = n_distinct(sample),
  median_gene_reads = median(gene_total_reads, na.rm = TRUE),
  median_novel_fraction = median(novel_fraction, na.rm = TRUE)
))

case_tbl <- coord_tbl %>%
  dplyr::filter(!is.na(gene), !toupper(group) %in% c("CONTROL")) %>%
  distinct(sample_core, sample_raw, gene, group)

control_meta <- meta_expanded %>%
  dplyr::filter(toupper(group) %in% c("CONTROL")) %>%
  distinct(sample_core, sample_raw, group)

single_case <- list()
per_junction <- list()
groupwise <- list()

for (g in unique(case_tbl$gene)) {
  sg <- sample_gene %>% dplyr::filter(gene == g)
  if (nrow(sg) == 0) next

  case_samples <- base::intersect(case_tbl$sample_core[case_tbl$gene == g], sg$sample_core)
  if (length(case_samples) == 0) next

  if (prefer_metadata_controls) {
    ctrl_samples <- base::intersect(control_meta$sample_core, sg$sample_core)
  } else {
    ctrl_samples <- base::setdiff(sg$sample_core, case_samples)
  }
  if (length(ctrl_samples) < min_controls_per_gene) next

  for (s in case_samples) {
    one <- sg %>% dplyr::filter(sample_core == s)
    ref <- sg %>% dplyr::filter(sample_core %in% ctrl_samples)
    if (nrow(one) == 0 || nrow(ref) < min_controls_per_gene) next
    single_case[[length(single_case) + 1]] <- tibble(
      gene = g,
      sample = one$sample[1],
      group = one$group[1],
      n_controls = nrow(ref),
      gene_total_reads = one$gene_total_reads[1],
      entropy = one$entropy[1],
      entropy_z = make_z(one$entropy[1], ref$entropy),
      entropy_emp_p = empirical_p_2sided(one$entropy[1], ref$entropy),
      novel_fraction = one$novel_fraction[1],
      novel_fraction_z = make_z(one$novel_fraction[1], ref$novel_fraction),
      novel_fraction_emp_p = empirical_p_2sided(one$novel_fraction[1], ref$novel_fraction),
      local_fraction = one$local_fraction[1],
      local_fraction_z = make_z(one$local_fraction[1], ref$local_fraction),
      local_fraction_emp_p = empirical_p_2sided(one$local_fraction[1], ref$local_fraction),
      local_novel_fraction = one$local_novel_fraction[1],
      local_novel_fraction_z = make_z(one$local_novel_fraction[1], ref$local_novel_fraction),
      local_novel_fraction_emp_p = empirical_p_2sided(one$local_novel_fraction[1], ref$local_novel_fraction)
    )
  }

  if (length(case_samples) >= 2) {
    cases <- sg %>% dplyr::filter(sample_core %in% case_samples)
    ctrls <- sg %>% dplyr::filter(sample_core %in% ctrl_samples)
    gtest <- function(x, y) tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    groupwise[[length(groupwise) + 1]] <- tibble(
      gene = g,
      n_cases = nrow(cases),
      n_controls = nrow(ctrls),
      entropy_p = gtest(cases$entropy, ctrls$entropy),
      novel_fraction_p = gtest(cases$novel_fraction, ctrls$novel_fraction),
      local_fraction_p = gtest(cases$local_fraction, ctrls$local_fraction),
      local_novel_fraction_p = gtest(cases$local_novel_fraction, ctrls$local_novel_fraction)
    )
  }

  ug_observed <- usage_rows %>% dplyr::filter(gene == g)
  top_j <- ug_observed %>% group_by(junction_id) %>% summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>% slice_head(n = max(25, label_top_junc_n)) %>% pull(junction_id)

  comparison_samples <- sg %>%
    dplyr::filter(sample_core %in% c(case_samples, ctrl_samples)) %>%
    dplyr::select(sample_core, sample) %>%
    distinct(sample_core, .keep_all = TRUE)
  junction_meta <- ug_observed %>%
    dplyr::filter(junction_id %in% top_j) %>%
    group_by(junction_id) %>%
    summarise(
      annotation_simple = dplyr::first(annotation_simple),
      chrom = dplyr::first(chrom),
      junction_start = dplyr::first(start),
      junction_end = dplyr::first(end),
      .groups = "drop"
    )
  ug <- tidyr::crossing(
    sample_core = comparison_samples$sample_core,
    junction_id = top_j
  ) %>%
    left_join(comparison_samples, by = "sample_core") %>%
    left_join(junction_meta, by = "junction_id") %>%
    left_join(
      ug_observed %>%
        group_by(sample_core, junction_id) %>%
        summarise(
          count = dplyr::first(count),
          usage = dplyr::first(usage),
          .groups = "drop"
        ),
      by = c("sample_core", "junction_id")
    ) %>%
    mutate(
      count = replace_na(count, 0),
      usage = replace_na(usage, 0)
    )

  for (jid in unique(ug$junction_id)) {
    dd <- ug %>% dplyr::filter(junction_id == jid)
    ref_usage <- dd %>% dplyr::filter(sample_core %in% ctrl_samples) %>% pull(usage)
    ref_count <- dd %>% dplyr::filter(sample_core %in% ctrl_samples) %>% pull(count)
    ref_logcount <- log1p(ref_count)
    for (s in case_samples) {
      obs <- dd %>% dplyr::filter(sample_core == s)
      if (nrow(obs) == 0) next
      target_s <- coord_tbl %>%
        dplyr::filter(sample_core == s, gene == g, has_interval) %>%
        slice_head(n = 1)
      has_target <- nrow(target_s) == 1
      same_chrom <- has_target && obs$chrom[1] == target_s$chrom[1]
      spans_target <- same_chrom &&
        obs$junction_start[1] <= target_s$end[1] &&
        obs$junction_end[1] >= target_s$start[1]
      near_target <- spans_target || (
        same_chrom && (
          abs(obs$junction_start[1] - target_s$start[1]) <= local_pad_bp ||
          abs(obs$junction_start[1] - target_s$end[1]) <= local_pad_bp ||
          abs(obs$junction_end[1] - target_s$start[1]) <= local_pad_bp ||
          abs(obs$junction_end[1] - target_s$end[1]) <= local_pad_bp
        )
      )
      per_junction[[length(per_junction) + 1]] <- tibble(
        gene = g,
        sample = obs$sample[1],
        junction_id = jid,
        annotation = obs$annotation_simple[1],
        annotation_simple = obs$annotation_simple[1],
        chrom = obs$chrom[1],
        start = obs$junction_start[1],
        end = obs$junction_end[1],
        junction_label = paste0(
          obs$chrom[1], ":", obs$junction_start[1], ":", obs$junction_end[1]
        ),
        inserted_usage = obs$usage[1],
        inserted_count = obs$count[1],
        case_usage = obs$usage[1],
        case_count = obs$count[1],
        case_logcount = log1p(obs$count[1]),
        case_detect = obs$count[1] > 0,
        n_controls_total = length(ref_usage),
        n_controls_with_junction = sum(ref_usage > 0, na.rm = TRUE),
        ref_n = length(ref_usage),
        ref_detect_n = sum(ref_usage > 0, na.rm = TRUE),
        ref_detect_rate = mean(ref_usage > 0, na.rm = TRUE),
        ref_mean_count = mean(ref_count, na.rm = TRUE),
        ref_median_count = median(ref_count, na.rm = TRUE),
        control_mean_usage = mean(ref_usage, na.rm = TRUE),
        control_median_usage = median(ref_usage, na.rm = TRUE),
        ref_mean_usage = mean(ref_usage, na.rm = TRUE),
        ref_median_usage = median(ref_usage, na.rm = TRUE),
        usage_z = make_z(obs$usage[1], ref_usage),
        usage_emp_p = empirical_p_2sided(obs$usage[1], ref_usage),
        logcount_z = make_z(log1p(obs$count[1]), ref_logcount),
        logcount_emp_p = empirical_p_2sided(log1p(obs$count[1]), ref_logcount),
        absent_in_case_supported_elsewhere =
          obs$count[1] == 0 &&
          sum(ref_usage > 0, na.rm = TRUE) >= absent_ref_detect_min_n &&
          mean(ref_usage > 0, na.rm = TRUE) >= absent_ref_detect_min_rate,
        absent_loss_score = ifelse(
          obs$count[1] == 0 &&
            sum(ref_usage > 0, na.rm = TRUE) >= absent_ref_detect_min_n &&
            mean(ref_usage > 0, na.rm = TRUE) >= absent_ref_detect_min_rate,
          mean(ref_usage > 0, na.rm = TRUE) * mean(ref_usage, na.rm = TRUE),
          0
        ),
        expressed_only_in_case =
          obs$count[1] > 0 && sum(ref_usage > 0, na.rm = TRUE) == 0,
        gained_novel_score = ifelse(
          obs$count[1] > 0 &&
            sum(ref_usage > 0, na.rm = TRUE) == 0 &&
            obs$annotation_simple[1] == "novel",
          obs$usage[1],
          0
        ),
        near_target = near_target,
        spans_target_interval = spans_target
      )
    }
  }
}

single_case_df <- bind_rows(single_case)
if (nrow(single_case_df) > 0) {
  single_case_df <- single_case_df %>% mutate(
    entropy_fdr = p.adjust(entropy_emp_p, method = "BH"),
    novel_fraction_fdr = p.adjust(novel_fraction_emp_p, method = "BH"),
    local_fraction_fdr = p.adjust(local_fraction_emp_p, method = "BH"),
    local_novel_fraction_fdr = p.adjust(local_novel_fraction_emp_p, method = "BH")
  )
}
write_tsv(single_case_df, file.path(table_dir, "single_case_disruption_results.tsv"))

groupwise_df <- bind_rows(groupwise)
if (nrow(groupwise_df) > 0) {
  groupwise_df <- groupwise_df %>% mutate(
    entropy_fdr = p.adjust(entropy_p, method = "BH"),
    novel_fraction_fdr = p.adjust(novel_fraction_p, method = "BH"),
    local_fraction_fdr = p.adjust(local_fraction_p, method = "BH"),
    local_novel_fraction_fdr = p.adjust(local_novel_fraction_p, method = "BH")
  )
}
write_tsv(groupwise_df, file.path(table_dir, "groupwise_disruption_results.tsv"))

per_junction_df <- bind_rows(per_junction)
if (nrow(per_junction_df) > 0) per_junction_df <- per_junction_df %>% mutate(usage_fdr = p.adjust(usage_emp_p, method = "BH"))
write_tsv(per_junction_df, file.path(table_dir, "per_junction_case_control_results.tsv"))

msg("\n=== Top single-case results ===")
if (nrow(single_case_df) > 0) {
  print(single_case_df %>% arrange(local_fraction_fdr, entropy_fdr) %>% slice_head(n = 20))
} else {
  msg("No single-case results were produced. Check sample matching and control availability.")
}

ann_df <- usage_rows %>%
  group_by(sample, annotation_simple) %>%
  summarise(reads = sum(count), .groups = "drop") %>%
  group_by(sample) %>% mutate(frac = reads / sum(reads)) %>% ungroup()

p_ann <- ggplot(ann_df, aes(x = reorder(sample, frac, FUN = sum), y = frac, fill = annotation_simple)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("annotated" = "#4C6A87", "novel" = "#9CB7CE", "unknown" = "#C0CEDD")) +
  labs(title = "Junction annotation composition by sample", x = "Sample", y = "Fraction of junction reads", fill = NULL)
save_plot_all(p_ann, "qc_annotation_fraction_by_sample", width = 7.8, height = 6.8)

sample_gene_plot <- sample_gene %>% mutate(GroupClass = ifelse(toupper(group) == "CONTROL", "Control", "Inserted/Case"))
for (metric in c("novel_fraction", "entropy", "local_fraction", "local_novel_fraction")) {
  plot_df <- sample_gene_plot %>% dplyr::filter(is.finite(.data[[metric]]))

  if (nrow(plot_df) == 0) next

  p <- ggplot(plot_df, aes(x = GroupClass, y = .data[[metric]], fill = GroupClass)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, size = 1.7, alpha = 0.75) +
    scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
    labs(title = paste0("Sample-gene ", metric, " by group"), x = NULL, y = metric) +
    theme(legend.position = "none")

  save_plot_all(p, paste0("global_", metric, "_by_group"), width = 5.8, height = 4.8)
}

plot_genes_sc <- character(0)
plot_genes_gp <- character(0)

if (exists("single_case_df") &&
    nrow(single_case_df) > 0 &&
    all(c("gene", "local_fraction_fdr", "entropy_fdr") %in% colnames(single_case_df))) {
  plot_genes_sc <- single_case_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
}

if (exists("groupwise_df") &&
    nrow(groupwise_df) > 0 &&
    all(c("gene", "local_fraction_fdr", "entropy_fdr") %in% colnames(groupwise_df))) {
  plot_genes_gp <- groupwise_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
}

plot_genes <- unique(c(plot_genes_sc, plot_genes_gp, case_tbl$gene))
plot_genes <- plot_genes[!is.na(plot_genes)]

for (g in plot_genes) {
  sg <- sample_gene %>% dplyr::filter(gene == g)
  if (nrow(sg) == 0) next
  sg2 <- sg %>% mutate(GroupClass = ifelse(toupper(group) == "CONTROL", "Control", "Inserted/Case"))

  sg_long <- sg2 %>%
    dplyr::select(sample, GroupClass, entropy, novel_fraction, local_fraction, local_novel_fraction) %>%
    pivot_longer(cols = c(entropy, novel_fraction, local_fraction, local_novel_fraction), names_to = "metric", values_to = "value")

  p_metrics <- ggplot(sg_long, aes(x = GroupClass, y = value, fill = GroupClass)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.75) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
    labs(title = paste0(g, ": splice disruption metrics"), x = NULL, y = NULL) +
    theme(legend.position = "none")
  save_plot_all(p_metrics, paste0("gene_metrics_", g), width = 8.2, height = 5.6)

  ug <- usage_rows %>% dplyr::filter(gene == g)
  if (nrow(ug) == 0) next
  top_j <- ug %>% group_by(junction_id) %>% summarise(total = sum(count), .groups = "drop") %>% arrange(desc(total)) %>% slice_head(n = label_top_junc_n) %>% pull(junction_id)
  hm_df <- ug %>% dplyr::filter(junction_id %in% top_j) %>% dplyr::select(junction_id, sample, usage) %>% pivot_wider(names_from = sample, values_from = usage, values_fill = 0)
  if (nrow(hm_df) >= 2 && ncol(hm_df) >= 3) {
    hm_mat <- as.data.frame(hm_df)
    rownames(hm_mat) <- hm_mat$junction_id
    hm_mat$junction_id <- NULL
    hm_mat <- as.matrix(hm_mat)
    ann_col <- sg2 %>% distinct(sample, GroupClass) %>% as.data.frame()
    rownames(ann_col) <- ann_col$sample
    ann_col$sample <- NULL
    ann_col <- ann_col[colnames(hm_mat), , drop = FALSE]
    pheat <- pheatmap(hm_mat, scale = "row", annotation_col = ann_col,
                      clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                      main = paste0(g, ": top junction usage (row-scaled)"),
                      fontsize = heatmap_text_size,
                      fontsize_row = heatmap_text_size,
                      fontsize_col = heatmap_text_size,
                      angle_col = 45,
                      silent = TRUE)
    save_pheatmap_all(pheat, paste0("gene_heatmap_", g), width = 8.2, height = 7.2)
  }

  pj <- per_junction_df %>% dplyr::filter(gene == g)
  if (nrow(pj) > 0) {
    keep_j <- pj %>% group_by(junction_id) %>% summarise(max_abs_z = max(abs(usage_z), na.rm = TRUE), .groups = "drop") %>% arrange(desc(max_abs_z)) %>% slice_head(n = label_top_junc_n) %>% pull(junction_id)
    pjb <- pj %>% dplyr::filter(junction_id %in% keep_j) %>% mutate(label = ifelse(near_target | spans_target_interval, paste0(junction_id, " *"), junction_id))
    p_bar <- ggplot(pjb, aes(x = reorder(label, inserted_usage), y = inserted_usage, fill = sample)) +
      geom_col(position = position_dodge(width = 0.8)) + coord_flip() +
      scale_fill_manual(values = get_muted_palette(max(3, n_distinct(pjb$sample)))) +
      labs(title = paste0(g, ": inserted-sample junction usage (* near/spans target interval)"), x = "Junction", y = "Usage", fill = "Sample")
    save_plot_all(p_bar, paste0("gene_junction_bar_", g), width = 9.0, height = 6.2)
  }
}

priority_tbl <- if (nrow(single_case_df) > 0) {
  single_case_df %>%
    mutate(score = pmax(abs(local_fraction_z), 0, na.rm = TRUE) + pmax(abs(entropy_z), 0, na.rm = TRUE) + pmax(abs(novel_fraction_z), 0, na.rm = TRUE)) %>%
    arrange(local_fraction_fdr, entropy_fdr, desc(score))
} else {
  tibble()
}
write_tsv(priority_tbl, file.path(table_dir, "priority_single_case_ranking.tsv"))

msg("\n=== Priority single-case ranking ===")
if (nrow(priority_tbl) > 0) print(priority_tbl %>% slice_head(n = 25))

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
msg("\nDone. Outputs written to: ", normalizePath(out_dir))

has_repel <- requireNamespace("ggrepel", quietly = TRUE)

display_lut <- meta_expanded %>%
  mutate(
    display_label = dplyr::coalesce(
      if ("paper_sample_id" %in% colnames(meta_expanded)) as.character(paper_sample_id) else NA_character_,
      if ("patient_id" %in% colnames(meta_expanded)) as.character(patient_id) else NA_character_,
      as.character(sample_raw)
    )
  ) %>%
  distinct(sample_core, sample_raw, display_label)

sample_gene <- sample_gene %>%
  left_join(display_lut %>% dplyr::select(sample_core, display_label), by = "sample_core") %>%
  mutate(display_label = dplyr::coalesce(display_label, sample))

usage_rows <- usage_rows %>%
  left_join(display_lut %>% dplyr::select(sample_core, display_label), by = "sample_core") %>%
  mutate(display_label = dplyr::coalesce(display_label, sample))

per_junction_df <- per_junction_df %>%
  left_join(display_lut %>% dplyr::select(sample_raw, display_label), by = c("sample" = "sample_raw")) %>%
  mutate(display_label = dplyr::coalesce(display_label, sample))

if (exists("single_case_df") && nrow(single_case_df) > 0) {
  single_case_df <- single_case_df %>%
    left_join(display_lut %>% dplyr::select(sample_raw, display_label), by = c("sample" = "sample_raw")) %>%
    mutate(display_label = dplyr::coalesce(display_label, sample))
}

case_control_class <- function(x) {
  ifelse(toupper(x) == "CONTROL", "Control", "Inserted/Case")
}

plot_genes_sc <- character(0)
plot_genes_gp <- character(0)

if (exists("single_case_df") &&
    nrow(single_case_df) > 0 &&
    all(c("gene", "local_fraction_fdr", "entropy_fdr") %in% colnames(single_case_df))) {
  plot_genes_sc <- single_case_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
}

if (exists("groupwise_df") &&
    nrow(groupwise_df) > 0 &&
    all(c("gene", "local_fraction_fdr", "entropy_fdr") %in% colnames(groupwise_df))) {
  plot_genes_gp <- groupwise_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
}

plot_genes <- unique(c(plot_genes_sc, plot_genes_gp, case_tbl$gene))
plot_genes <- plot_genes[!is.na(plot_genes)]

if (length(plot_genes) == 0 && exists("single_case_df") && nrow(single_case_df) > 0) {
  plot_genes <- unique(single_case_df$gene)
}

plot_gene_audit <- tibble::tibble(gene = unique(case_tbl$gene)) %>%
  dplyr::filter(!is.na(gene), nzchar(gene)) %>%
  dplyr::mutate(
    present_in_sample_gene = gene %in% unique(sample_gene$gene),
    present_in_usage_rows = gene %in% unique(usage_rows$gene),
    present_in_per_junction_results = gene %in% unique(per_junction_df$gene),
    scheduled_for_plotting = gene %in% plot_genes
  )
readr::write_tsv(
  plot_gene_audit,
  file.path(table_dir, "metadata_case_gene_plot_audit.tsv")
)
msg("\n=== Metadata case-gene plot audit ===")
print(plot_gene_audit)

sample_gene_plot <- sample_gene %>%
  mutate(GroupClass = case_control_class(group))

for (metric in c("novel_fraction", "entropy", "local_fraction", "local_novel_fraction")) {
  plot_df <- sample_gene_plot %>% dplyr::filter(is.finite(.data[[metric]]))
  if (nrow(plot_df) == 0) next

  p <- ggplot(plot_df, aes(x = GroupClass, y = .data[[metric]], fill = GroupClass)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.15, size = 1.7, alpha = 0.75) +
    scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
    labs(title = paste0("Sample-gene ", metric, " by group"), x = NULL, y = metric) +
    theme(legend.position = "none")

  save_plot_all(p, paste0("global_", metric, "_by_group"), width = 5.8, height = 4.8)
}

for (g in plot_genes) {

  sg <- sample_gene %>%
    dplyr::filter(gene == g) %>%
    mutate(
      GroupClass = case_control_class(group),
      Display = dplyr::coalesce(display_label, sample)
    )

  if (nrow(sg) == 0) next

  sg_long <- sg %>%
    dplyr::select(sample, Display, GroupClass, entropy, novel_fraction, local_fraction, local_novel_fraction) %>%
    pivot_longer(
      cols = c(entropy, novel_fraction, local_fraction, local_novel_fraction),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(is.finite(value))

  if (nrow(sg_long) > 0) {
    p_metrics <- ggplot(sg_long, aes(x = GroupClass, y = value, fill = GroupClass)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.95) +
      geom_jitter(width = 0.12, size = 1.4, alpha = 0.45, color = "grey35") +
      geom_point(
        data = sg_long %>% dplyr::filter(GroupClass == "Inserted/Case"),
        color = "#4C6A87", size = 2.8, alpha = 0.95
      ) +
      facet_wrap(~ metric, scales = "free_y") +
      scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
      labs(
        title = paste0(g, ": splice disruption metrics (case labeled)"),
        x = NULL, y = NULL
      ) +
      theme(legend.position = "none")

    if (has_repel) {
      p_metrics <- p_metrics +
        ggrepel::geom_text_repel(
          data = sg_long %>% dplyr::filter(GroupClass == "Inserted/Case"),
          aes(label = Display),
          size = 4.2,
          seed = 123,
          max.overlaps = Inf,
          box.padding = 0.2,
          point.padding = 0.15,
          segment.color = "grey60",
          segment.size = 0.3
        )
    } else {
      p_metrics <- p_metrics +
        geom_text(
          data = sg_long %>% dplyr::filter(GroupClass == "Inserted/Case"),
          aes(label = Display),
          size = 4.0,
          vjust = -0.4
        )
    }

    save_plot_all(p_metrics, paste0("gene_metrics_labeled_", g), width = 9.0, height = 5.8)
  }

  ug_observed <- usage_rows %>%
    dplyr::filter(gene == g)

  if (nrow(ug_observed) == 0) next

  top_j <- ug_observed %>%
    group_by(junction_id) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = label_top_junc_n) %>%
    pull(junction_id)

  sample_labels <- sg %>%
    dplyr::select(sample, GroupClass, Display) %>%
    distinct(sample, .keep_all = TRUE)

  ug <- tidyr::crossing(
    sample = sample_labels$sample,
    junction_id = top_j
  ) %>%
    left_join(sample_labels, by = "sample") %>%
    left_join(
      ug_observed %>%
        group_by(sample, junction_id) %>%
        summarise(
          count = dplyr::first(count),
          usage = dplyr::first(usage),
          .groups = "drop"
        ),
      by = c("sample", "junction_id")
    ) %>%
    mutate(
      count = replace_na(count, 0),
      usage = replace_na(usage, 0)
    )

  hm_df <- ug %>%
    dplyr::select(junction_id, sample, Display, GroupClass, usage) %>%
    distinct(junction_id, sample, .keep_all = TRUE) %>%
    mutate(HeatmapSample = ifelse(GroupClass == "Inserted/Case", paste0("[CASE] ", Display), Display)) %>%
    dplyr::select(junction_id, HeatmapSample, usage) %>%
    pivot_wider(names_from = HeatmapSample, values_from = usage, values_fill = 0)

  if (nrow(hm_df) >= 2 && ncol(hm_df) >= 3) {
    hm_mat <- as.data.frame(hm_df)
    rownames(hm_mat) <- hm_mat$junction_id
    hm_mat$junction_id <- NULL
    hm_mat <- as.matrix(hm_mat)

    sample_order <- ug %>%
      distinct(Display, GroupClass) %>%
      mutate(HeatmapSample = ifelse(GroupClass == "Inserted/Case", paste0("[CASE] ", Display), Display)) %>%
      arrange(desc(GroupClass == "Inserted/Case"), HeatmapSample) %>%
      pull(HeatmapSample)

    sample_order <- base::intersect(sample_order, colnames(hm_mat))
    hm_mat <- hm_mat[, sample_order, drop = FALSE]

    ann_col <- ug %>%
      distinct(Display, GroupClass) %>%
      mutate(HeatmapSample = ifelse(GroupClass == "Inserted/Case", paste0("[CASE] ", Display), Display)) %>%
      dplyr::select(HeatmapSample, GroupClass) %>%
      as.data.frame()

    rownames(ann_col) <- ann_col$HeatmapSample
    ann_col$HeatmapSample <- NULL
    ann_col <- ann_col[colnames(hm_mat), , drop = FALSE]

    pheat <- pheatmap(
      hm_mat,
      scale = "row",
      annotation_col = ann_col,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      main = paste0(g, ": top junction usage (case columns marked with [CASE])"),
      fontsize = heatmap_text_size,
      fontsize_row = heatmap_text_size,
      fontsize_col = heatmap_text_size,
      angle_col = 45,
      silent = TRUE
    )

    save_pheatmap_all(pheat, paste0("gene_heatmap_labeled_", g), width = 9.5, height = 7.4)
  }

  top_box_j <- ug %>%
    group_by(junction_id) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = label_top_junc_n) %>%
    pull(junction_id)

  box_df <- ug %>%
    dplyr::filter(junction_id %in% top_box_j) %>%
    mutate(JunctionLabel = junction_id)

  near_map <- per_junction_df %>%
    dplyr::filter(gene == g) %>%
    group_by(junction_id) %>%
    summarise(any_near = any(near_target | spans_target_interval), .groups = "drop")

  box_df <- box_df %>%
    left_join(near_map, by = "junction_id") %>%
    mutate(JunctionLabel = ifelse(!is.na(any_near) & any_near, paste0(junction_id, " *"), junction_id))

  p_box <- ggplot(box_df, aes(x = reorder(JunctionLabel, usage, FUN = median), y = usage)) +
    geom_boxplot(
      data = box_df %>% dplyr::filter(GroupClass == "Control"),
      fill = "#C0CEDD",
      outlier.shape = NA,
      width = 0.7
    ) +
    geom_jitter(
      data = box_df %>% dplyr::filter(GroupClass == "Control"),
      width = 0.12,
      size = 1.2,
      alpha = 0.45,
      color = "grey45"
    ) +
    geom_point(
      data = box_df %>% dplyr::filter(GroupClass == "Inserted/Case"),
      color = "#4C6A87",
      size = 2.8
    ) +
    coord_flip() +
    labs(
      title = paste0(g, ": case vs controls junction-usage comparison (* near/spans target interval)"),
      x = "Junction",
      y = "Within-gene junction usage"
    )

  if (has_repel) {
    p_box <- p_box +
      ggrepel::geom_text_repel(
        data = box_df %>% dplyr::filter(GroupClass == "Inserted/Case"),
        aes(label = Display),
        size = 4.2,
        seed = 123,
        max.overlaps = Inf,
        box.padding = 0.15,
        point.padding = 0.15,
        segment.color = "grey60",
        segment.size = 0.3,
        direction = "y"
      )
  } else {
    p_box <- p_box +
      geom_text(
        data = box_df %>% dplyr::filter(GroupClass == "Inserted/Case"),
        aes(label = Display),
        hjust = -0.1,
        size = 4.0
      )
  }

  save_plot_all(p_box, paste0("gene_junction_case_vs_controls_", g), width = 10.2, height = 6.8)

  pj <- per_junction_df %>% dplyr::filter(gene == g)
  if (nrow(pj) > 0) {
    keep_j <- pj %>%
      group_by(junction_id) %>%
      summarise(max_abs_z = max(abs(usage_z), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_abs_z)) %>%
      slice_head(n = label_top_junc_n) %>%
      pull(junction_id)

    pjb <- pj %>%
      dplyr::filter(junction_id %in% keep_j) %>%
      mutate(
        label = ifelse(near_target | spans_target_interval, paste0(junction_id, " *"), junction_id),
        Display = dplyr::coalesce(display_label, sample)
      )

    p_bar <- ggplot(pjb, aes(x = reorder(label, inserted_usage), y = inserted_usage, fill = Display)) +
      geom_col(position = position_dodge(width = 0.8)) +
      coord_flip() +
      scale_fill_manual(values = get_muted_palette(max(3, n_distinct(pjb$Display)))) +
      labs(
        title = paste0(g, ": case sample junction usage (* near/spans target interval)"),
        x = "Junction",
        y = "Usage",
        fill = "Patient / sample"
      )

    save_plot_all(p_bar, paste0("gene_junction_bar_labeled_", g), width = 10.0, height = 6.5)
  }
}

cat("\nPlotting complete.\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

msg("\n============================================================")
msg("Transcript-aware annotation for MAFG/SIRT7 + entropy plot")
msg("============================================================")

needed_objs_tx <- c("per_junction_df", "sample_gene", "case_tbl", "table_dir", "plot_dir")
missing_objs_tx <- needed_objs_tx[!vapply(needed_objs_tx, exists, logical(1))]
if (length(missing_objs_tx) > 0) {
  stop("Missing required objects for transcript-aware analysis: ",
       paste(missing_objs_tx, collapse = ", "))
}

if (!exists("gtf")) {
  if (!exists("gtf_path")) stop("Neither gtf nor gtf_path exists.")
  msg("GTF object not found in memory; importing: ", gtf_path)
  gtf <- rtracklayer::import(gtf_path)
}

# Focused transcript annotations
target_tx_genes <- c("MAFG", "SIRT7")
boundary_tol_bp <- 2

safe_mcol <- function(gr, colname, default = NA_character_) {
  if (colname %in% names(S4Vectors::mcols(gr))) {
    as.character(S4Vectors::mcols(gr)[[colname]])
  } else {
    rep(default, length(gr))
  }
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

classify_junction_to_transcript <- function(j_start, j_end, exon_tbl, tol = 2) {

  left_hits <- exon_tbl %>%
    dplyr::filter(abs(exon_end - j_start) <= tol)

  right_hits <- exon_tbl %>%
    dplyr::filter(abs(exon_start - j_end) <= tol)

  if (nrow(left_hits) == 0 && nrow(right_hits) == 0) {
    return(tibble(
      event_type = "unmapped_to_transcript_exons",
      left_exon_number = NA_character_,
      right_exon_number = NA_character_,
      left_exon_rank_genomic = NA_integer_,
      right_exon_rank_genomic = NA_integer_,
      skipped_exon_numbers = NA_character_,
      skipped_exon_rank_genomic = NA_character_,
      n_skipped_exons = 0L
    ))
  }

  if (nrow(left_hits) > 0 && nrow(right_hits) > 0) {

    out <- list()

    for (i in seq_len(nrow(left_hits))) {
      for (j in seq_len(nrow(right_hits))) {

        l <- left_hits[i, ]
        r <- right_hits[j, ]

        if (r$exon_rank_genomic == l$exon_rank_genomic + 1) {
          event_type <- "canonical_adjacent_junction"
          skipped <- exon_tbl[0, ]
        } else if (r$exon_rank_genomic > l$exon_rank_genomic + 1) {
          event_type <- "exon_skipping_candidate"
          skipped <- exon_tbl %>%
            dplyr::filter(
              exon_rank_genomic > l$exon_rank_genomic,
              exon_rank_genomic < r$exon_rank_genomic
            )
        } else {
          event_type <- "noncanonical_reverse_or_overlapping_boundary"
          skipped <- exon_tbl[0, ]
        }

        out[[length(out) + 1]] <- tibble(
          event_type = event_type,
          left_exon_number = as.character(l$exon_number),
          right_exon_number = as.character(r$exon_number),
          left_exon_rank_genomic = as.integer(l$exon_rank_genomic),
          right_exon_rank_genomic = as.integer(r$exon_rank_genomic),
          skipped_exon_numbers = collapse_unique(as.character(skipped$exon_number)),
          skipped_exon_rank_genomic = collapse_unique(as.character(skipped$exon_rank_genomic)),
          n_skipped_exons = nrow(skipped)
        )
      }
    }

    return(dplyr::bind_rows(out))
  }

  if (nrow(left_hits) > 0 && nrow(right_hits) == 0) {
    return(tibble(
      event_type = "alternative_acceptor_candidate",
      left_exon_number = collapse_unique(as.character(left_hits$exon_number)),
      right_exon_number = NA_character_,
      left_exon_rank_genomic = suppressWarnings(min(left_hits$exon_rank_genomic, na.rm = TRUE)),
      right_exon_rank_genomic = NA_integer_,
      skipped_exon_numbers = NA_character_,
      skipped_exon_rank_genomic = NA_character_,
      n_skipped_exons = 0L
    ))
  }

  if (nrow(left_hits) == 0 && nrow(right_hits) > 0) {
    return(tibble(
      event_type = "alternative_donor_candidate",
      left_exon_number = NA_character_,
      right_exon_number = collapse_unique(as.character(right_hits$exon_number)),
      left_exon_rank_genomic = NA_integer_,
      right_exon_rank_genomic = suppressWarnings(min(right_hits$exon_rank_genomic, na.rm = TRUE)),
      skipped_exon_numbers = NA_character_,
      skipped_exon_rank_genomic = NA_character_,
      n_skipped_exons = 0L
    ))
  }
}

calc_cds_impact <- function(tx_id, event_row, exon_tbl, cds_tbl) {

  tx_cds <- cds_tbl %>%
    dplyr::filter(transcript_id == tx_id)

  if (nrow(tx_cds) == 0) {
    return(tibble(
      affected_exon_numbers = NA_character_,
      cds_bases_affected = NA_integer_,
      cds_bases_affected_mod3 = NA_integer_,
      cds_overlap_class = "no_CDS_annotation_for_transcript",
      predicted_functional_impact = "noncoding_or_unannotated_CDS"
    ))
  }

  event_type <- event_row$event_type[1]

  affected_exons <- exon_tbl[0, ]

  if (event_type == "exon_skipping_candidate" &&
      !is.na(event_row$skipped_exon_rank_genomic[1])) {

    ranks <- as.integer(strsplit(event_row$skipped_exon_rank_genomic[1], ";")[[1]])
    affected_exons <- exon_tbl %>%
      dplyr::filter(exon_rank_genomic %in% ranks)

  } else if (event_type %in% c("canonical_adjacent_junction",
                               "alternative_acceptor_candidate",
                               "alternative_donor_candidate",
                               "splice_boundary_change_CDS_impact_uncertain")) {

    ranks <- c(event_row$left_exon_rank_genomic[1], event_row$right_exon_rank_genomic[1])
    ranks <- ranks[is.finite(ranks)]
    affected_exons <- exon_tbl %>%
      dplyr::filter(exon_rank_genomic %in% ranks)
  }

  if (nrow(affected_exons) == 0) {
    return(tibble(
      affected_exon_numbers = NA_character_,
      cds_bases_affected = 0L,
      cds_bases_affected_mod3 = NA_integer_,
      cds_overlap_class = "no_mapped_affected_exon",
      predicted_functional_impact = "unmapped_event_functional_impact_uncertain"
    ))
  }

  cds_bases <- 0L

  for (i in seq_len(nrow(affected_exons))) {
    ex_start <- affected_exons$exon_start[i]
    ex_end <- affected_exons$exon_end[i]

    overlaps <- tx_cds %>%
      dplyr::mutate(
        ov_start = pmax(cds_start, ex_start),
        ov_end = pmin(cds_end, ex_end),
        ov_width = pmax(0, ov_end - ov_start + 1)
      )

    cds_bases <- cds_bases + sum(overlaps$ov_width, na.rm = TRUE)
  }

  cds_mod <- ifelse(cds_bases > 0, cds_bases %% 3L, NA_integer_)

  cds_overlap_class <- dplyr::case_when(
    cds_bases == 0 ~ "UTR_or_nonCDS_event",
    cds_bases > 0 ~ "CDS_overlapping_event",
    TRUE ~ "unknown"
  )

  predicted_functional_impact <- dplyr::case_when(
    event_type == "exon_skipping_candidate" & cds_bases > 0 & cds_mod == 0 ~
      "predicted_in_frame_CDS_skip",
    event_type == "exon_skipping_candidate" & cds_bases > 0 & cds_mod != 0 ~
      "predicted_frameshift_CDS_skip_possible_NMD",
    event_type == "exon_skipping_candidate" & cds_bases == 0 ~
      "predicted_UTR_or_nonCDS_skip",
    event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate",
                      "canonical_adjacent_junction") & cds_bases > 0 ~
      "splice_boundary_change_CDS_impact_uncertain",
    event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate",
                      "canonical_adjacent_junction") & cds_bases == 0 ~
      "splice_boundary_change_UTR_or_nonCDS",
    TRUE ~ "functional_impact_uncertain"
  )

  tibble(
    affected_exon_numbers = collapse_unique(as.character(affected_exons$exon_number)),
    cds_bases_affected = as.integer(cds_bases),
    cds_bases_affected_mod3 = as.integer(cds_mod),
    cds_overlap_class = cds_overlap_class,
    predicted_functional_impact = predicted_functional_impact
  )
}

msg("\nBuilding transcript-aware exon and CDS models...")

tx_records <- gtf[gtf$type == "transcript"]
exon_records <- gtf[gtf$type == "exon"]
cds_records <- gtf[gtf$type == "CDS"]

GenomeInfoDb::seqlevelsStyle(exon_records) <- "UCSC"
if (length(cds_records) > 0) GenomeInfoDb::seqlevelsStyle(cds_records) <- "UCSC"
if (length(tx_records) > 0) GenomeInfoDb::seqlevelsStyle(tx_records) <- "UCSC"

exon_model <- tibble(
  gene = safe_mcol(exon_records, "gene_name"),
  gene_id = safe_mcol(exon_records, "gene_id"),
  transcript_id = safe_mcol(exon_records, "transcript_id"),
  transcript_name = safe_mcol(exon_records, "transcript_name"),
  transcript_biotype = dplyr::coalesce(
    safe_mcol(exon_records, "transcript_biotype"),
    safe_mcol(exon_records, "transcript_type")
  ),
  exon_number_raw = safe_mcol(exon_records, "exon_number"),
  chrom = norm_chr(as.character(GenomicRanges::seqnames(exon_records))),
  strand = as.character(GenomicRanges::strand(exon_records)),
  exon_start = GenomicRanges::start(exon_records),
  exon_end = GenomicRanges::end(exon_records)
) %>%
  dplyr::filter(gene %in% target_tx_genes, !is.na(transcript_id), transcript_id != "") %>%
  dplyr::mutate(
    exon_number = dplyr::if_else(
      is.na(exon_number_raw) | exon_number_raw == "",
      NA_character_,
      exon_number_raw
    )
  ) %>%
  dplyr::group_by(gene, transcript_id, chrom, strand, exon_start, exon_end) %>%
  dplyr::summarise(
    gene_id = dplyr::first(gene_id),
    transcript_name = dplyr::first(transcript_name),
    transcript_biotype = dplyr::first(transcript_biotype),
    exon_number = dplyr::first(exon_number),
    .groups = "drop"
  ) %>%
  dplyr::group_by(gene, transcript_id) %>%
  dplyr::arrange(exon_start, exon_end, .by_group = TRUE) %>%
  dplyr::mutate(
    exon_rank_genomic = dplyr::row_number(),
    exon_number = dplyr::coalesce(exon_number, as.character(exon_rank_genomic))
  ) %>%
  dplyr::ungroup()

cds_model <- tibble(
  gene = safe_mcol(cds_records, "gene_name"),
  transcript_id = safe_mcol(cds_records, "transcript_id"),
  chrom = norm_chr(as.character(GenomicRanges::seqnames(cds_records))),
  strand = as.character(GenomicRanges::strand(cds_records)),
  cds_start = GenomicRanges::start(cds_records),
  cds_end = GenomicRanges::end(cds_records)
) %>%
  dplyr::filter(gene %in% target_tx_genes, !is.na(transcript_id), transcript_id != "")

readr::write_tsv(exon_model, file.path(table_dir, "transcript_exon_model_MAFG_SIRT7.tsv"))
readr::write_tsv(cds_model, file.path(table_dir, "transcript_cds_model_MAFG_SIRT7.tsv"))

altered_junctions <- per_junction_df %>%
  dplyr::filter(gene %in% target_tx_genes) %>%
  dplyr::filter(
    absent_in_case_supported_elsewhere |
      expressed_only_in_case |
      abs(usage_z) >= 1.5 |
      abs(logcount_z) >= 1.5
  ) %>%
  dplyr::mutate(
    alteration_class = dplyr::case_when(
      absent_in_case_supported_elsewhere ~ "lost_in_case_supported_elsewhere",
      expressed_only_in_case ~ "gained_case_only",
      usage_z <= -1.5 | logcount_z <= -1.5 ~ "reduced_in_case",
      usage_z >= 1.5  | logcount_z >= 1.5  ~ "increased_in_case",
      TRUE ~ "altered"
    )
  ) %>%
  dplyr::select(
    gene, sample, display_label, junction_id, junction_label,
    chrom, start, end, annotation_simple,
    alteration_class,
    case_count, case_usage,
    ref_detect_n, ref_detect_rate,
    ref_mean_count, ref_mean_usage,
    usage_z, logcount_z,
    absent_in_case_supported_elsewhere,
    absent_loss_score,
    expressed_only_in_case,
    gained_novel_score,
    near_target, spans_target_interval
  ) %>%
  dplyr::distinct()

readr::write_tsv(altered_junctions, file.path(table_dir, "altered_junctions_for_transcript_annotation_MAFG_SIRT7.tsv"))

msg("\nAltered junctions selected for transcript-aware annotation:")
print(altered_junctions %>% dplyr::count(gene, alteration_class))

tx_annotation_rows <- list()

for (i in seq_len(nrow(altered_junctions))) {

  jrow <- altered_junctions[i, ]
  g <- jrow$gene
  j_start <- jrow$start
  j_end <- jrow$end

  tx_ids <- exon_model %>%
    dplyr::filter(gene == g) %>%
    dplyr::pull(transcript_id) %>%
    unique()

  for (tx in tx_ids) {
    ex_tx <- exon_model %>%
      dplyr::filter(gene == g, transcript_id == tx) %>%
      dplyr::arrange(exon_rank_genomic)

    if (nrow(ex_tx) == 0) next

    class_tbl <- classify_junction_to_transcript(j_start, j_end, ex_tx, tol = boundary_tol_bp)

    for (k in seq_len(nrow(class_tbl))) {

      ctbl <- class_tbl[k, ]
      impact_tbl <- calc_cds_impact(tx, ctbl, ex_tx, cds_model)

      tx_annotation_rows[[length(tx_annotation_rows) + 1]] <- dplyr::bind_cols(
        jrow,
        tibble(
          transcript_id = tx,
          transcript_name = dplyr::first(ex_tx$transcript_name),
          transcript_biotype = dplyr::first(ex_tx$transcript_biotype),
          strand = dplyr::first(ex_tx$strand),
          n_exons_transcript = nrow(ex_tx)
        ),
        ctbl,
        impact_tbl
      )
    }
  }
}

tx_junction_annotation <- dplyr::bind_rows(tx_annotation_rows)
if (nrow(tx_junction_annotation) == 0) {

  tx_junction_annotation <- altered_junctions[0, ] %>%
    dplyr::mutate(
      transcript_id = character(),
      transcript_name = character(),
      transcript_biotype = character(),
      strand = character(),
      n_exons_transcript = integer(),
      event_type = character(),
      left_exon_number = character(),
      right_exon_number = character(),
      left_exon_rank_genomic = integer(),
      right_exon_rank_genomic = integer(),
      skipped_exon_numbers = character(),
      skipped_exon_rank_genomic = character(),
      n_skipped_exons = integer(),
      affected_exon_numbers = character(),
      cds_bases_affected = integer(),
      cds_bases_affected_mod3 = integer(),
      cds_overlap_class = character(),
      predicted_functional_impact = character()
    )
  msg("No altered MAFG/SIRT7 junctions met the transcript-aware thresholds.")
}

tx_junction_annotation <- tx_junction_annotation %>%
  dplyr::mutate(
    transcript_event_score = dplyr::case_when(
      event_type == "exon_skipping_candidate" ~ 100,
      event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate") ~ 50,
      event_type == "canonical_adjacent_junction" ~ 25,
      TRUE ~ 0
    ) +
      dplyr::coalesce(abs(usage_z), 0) +
      dplyr::coalesce(abs(logcount_z), 0) +
      dplyr::coalesce(absent_loss_score, 0) +
      dplyr::coalesce(gained_novel_score, 0)
  ) %>%
  dplyr::arrange(gene, dplyr::desc(transcript_event_score))

readr::write_tsv(
  tx_junction_annotation,
  file.path(table_dir, "transcript_aware_junction_annotation_MAFG_SIRT7.tsv")
)

tx_event_summary <- tx_junction_annotation %>%
  dplyr::group_by(
    gene, transcript_id, transcript_name, transcript_biotype,
    event_type, predicted_functional_impact
  ) %>%
  dplyr::summarise(
    n_junctions = dplyr::n_distinct(junction_id),
    n_lost = sum(alteration_class == "lost_in_case_supported_elsewhere", na.rm = TRUE),
    n_gained = sum(alteration_class == "gained_case_only", na.rm = TRUE),
    max_abs_usage_z = suppressWarnings(max(abs(usage_z), na.rm = TRUE)),
    max_abs_logcount_z = suppressWarnings(max(abs(logcount_z), na.rm = TRUE)),
    max_absent_loss_score = suppressWarnings(max(absent_loss_score, na.rm = TRUE)),
    max_gained_novel_score = suppressWarnings(max(gained_novel_score, na.rm = TRUE)),
    affected_exon_numbers = collapse_unique(affected_exon_numbers),
    skipped_exon_numbers = collapse_unique(skipped_exon_numbers),
    max_cds_bases_affected = suppressWarnings(max(cds_bases_affected, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    max_abs_usage_z = ifelse(is.infinite(max_abs_usage_z), NA_real_, max_abs_usage_z),
    max_abs_logcount_z = ifelse(is.infinite(max_abs_logcount_z), NA_real_, max_abs_logcount_z),
    max_absent_loss_score = ifelse(is.infinite(max_absent_loss_score), NA_real_, max_absent_loss_score),
    max_gained_novel_score = ifelse(is.infinite(max_gained_novel_score), NA_real_, max_gained_novel_score),
    max_cds_bases_affected = ifelse(is.infinite(max_cds_bases_affected), NA_real_, max_cds_bases_affected)
  ) %>%
  dplyr::arrange(
    gene,
    dplyr::desc(n_lost + n_gained),
    dplyr::desc(max_absent_loss_score),
    dplyr::desc(max_gained_novel_score),
    dplyr::desc(max_abs_usage_z)
  )

readr::write_tsv(
  tx_event_summary,
  file.path(table_dir, "transcript_event_summary_MAFG_SIRT7.tsv")
)

msg("\nTop transcript-aware event summary:")
print(tx_event_summary %>% dplyr::slice_head(n = 40))

top_tx_events <- tx_junction_annotation %>%
  dplyr::filter(
    event_type != "unmapped_to_transcript_exons"
  ) %>%
  dplyr::group_by(gene, junction_id, alteration_class) %>%
  dplyr::slice_max(order_by = transcript_event_score, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(gene, dplyr::desc(transcript_event_score)) %>%
  dplyr::select(
    gene, sample, display_label, alteration_class,
    junction_label, annotation_simple,
    case_count, case_usage,
    ref_detect_n, ref_detect_rate,
    usage_z, logcount_z,
    absent_loss_score, gained_novel_score,
    transcript_id, transcript_name, transcript_biotype,
    event_type,
    left_exon_number, right_exon_number,
    skipped_exon_numbers, n_skipped_exons,
    affected_exon_numbers,
    cds_bases_affected, cds_bases_affected_mod3,
    cds_overlap_class,
    predicted_functional_impact,
    transcript_event_score
  )

readr::write_tsv(
  top_tx_events,
  file.path(table_dir, "top_transcript_aware_events_MAFG_SIRT7.tsv")
)

msg("\nTop transcript-aware junction events:")
print(top_tx_events %>% dplyr::slice_head(n = 30))

msg("\nGenerating entropy summary figure across genes...")

case_lookup_entropy <- case_tbl %>%
  dplyr::select(gene, sample_core) %>%
  dplyr::distinct() %>%
  dplyr::mutate(is_case_gene = TRUE)

entropy_plot_df <- sample_gene %>%
  dplyr::left_join(case_lookup_entropy, by = c("gene", "sample_core")) %>%
  dplyr::mutate(
    is_case_gene = ifelse(is.na(is_case_gene), FALSE, is_case_gene),
    Role = ifelse(is_case_gene, "Inserted case", "All other samples")
  ) %>%
  dplyr::filter(gene %in% unique(case_tbl$gene)) %>%
  dplyr::filter(is.finite(entropy))

gene_order_entropy <- entropy_plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    case_entropy = dplyr::first(entropy[Role == "Inserted case"]),
    ref_median_entropy = median(entropy[Role == "All other samples"], na.rm = TRUE),
    entropy_delta = case_entropy - ref_median_entropy,
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(abs(entropy_delta))) %>%
  dplyr::pull(gene)

entropy_plot_df$gene <- factor(entropy_plot_df$gene, levels = gene_order_entropy)

p_entropy_all <- ggplot(entropy_plot_df, aes(x = gene, y = entropy)) +
  geom_point(
    data = entropy_plot_df %>% dplyr::filter(Role == "All other samples"),
    color = "black",
    size = 2,
    position = position_jitter(width = 0.15, height = 0)
  ) +
  geom_point(
    data = entropy_plot_df %>% dplyr::filter(Role == "Inserted case"),
    color = "red",
    size = 3
  ) +
  coord_flip() +
  labs(
    title = "Splice-junction entropy by target gene",
    subtitle = "Black = all other samples; red = vector-inserted case sample",
    x = "Target gene",
    y = "Junction usage entropy"
  )

save_plot_all(
  p_entropy_all,
  "entropy_by_gene_case_vs_allother",
  width = 7.5,
  height = 5.8
)

entropy_summary_table <- entropy_plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    case_sample = dplyr::first(sample[Role == "Inserted case"]),
    case_label = dplyr::first(display_label[Role == "Inserted case"]),
    case_entropy = dplyr::first(entropy[Role == "Inserted case"]),
    ref_median_entropy = median(entropy[Role == "All other samples"], na.rm = TRUE),
    ref_mean_entropy = mean(entropy[Role == "All other samples"], na.rm = TRUE),
    entropy_delta = case_entropy - ref_median_entropy,
    n_reference = sum(Role == "All other samples"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(abs(entropy_delta)))

readr::write_tsv(
  entropy_summary_table,
  file.path(table_dir, "entropy_by_gene_case_vs_allother_summary.tsv")
)

msg("\nEntropy summary:")
print(entropy_summary_table)

msg("\nTranscript annotation complete.")
msg("New transcript-aware outputs:")
msg("  - ", file.path(table_dir, "transcript_exon_model_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "transcript_cds_model_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "altered_junctions_for_transcript_annotation_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "transcript_aware_junction_annotation_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "transcript_event_summary_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "top_transcript_aware_events_MAFG_SIRT7.tsv"))
msg("  - ", file.path(table_dir, "entropy_by_gene_case_vs_allother_summary.tsv"))
msg("New entropy figure:")
msg("  - ", file.path(plot_dir, "entropy_by_gene_case_vs_allother.pdf"))
msg("  - ", file.path(plot_dir, "entropy_by_gene_case_vs_allother.eps"))

writeLines(
  c(
    "Figure export settings",
    "PDF device: grDevices::pdf",
    "EPS device: grDevices::postscript",
    paste0("Font family: ", base_font),
    paste0("Base font size: ", base_font_size, " pt"),
    "showtext: disabled",
    "Vector figures retain selectable text"
  ),
  con = file.path(log_dir, "figure_export_settings.txt")
)
