#!/usr/bin/env Rscript
# BAM splice-junction analysis
# Extract and compare splice-junction evidence at insertion-associated genes.

args <- commandArgs(trailingOnly = TRUE)
default_base_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
base_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_base_dir
base_dir <- normalizePath(path.expand(base_dir), mustWork = TRUE)

bam_dir        <- file.path(base_dir, "bam")
meta_path <- file.path(base_dir, "splicing_metadata.tsv")
insertion_path <- file.path(base_dir, "insertion_coordinates.tsv")
gtf_path <- file.path(base_dir, "GRCh38.gtf.gz")
chain_path <- file.path(base_dir, "hg19ToHg38.over.chain.gz")

# Coordinate conventions
coord_source <- "metadata"
coord_build  <- "hg38"

only_matched_bams <- TRUE

# BAM extraction and comparison parameters
fetch_pad_bp <- 500
local_pad_bp <- 1000
min_gene_total_reads <- 5
min_controls_per_gene <- 2
label_top_junc_n <- 12

base_font_size    <- 14
plot_title_size   <- 16
axis_title_size   <- 14
axis_text_size    <- 12
legend_title_size <- 13
legend_text_size  <- 12
strip_text_size   <- 13
heatmap_text_size <- 12

out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "splicing", "bam_junction_analysis")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
log_dir   <- file.path(out_dir, "logs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

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
  library(Rsamtools)
  library(GenomicAlignments)
  library(pheatmap)
})

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}
base_font <- "Helvetica"

theme_set(
  theme_classic(base_size = base_font_size, base_family = base_font) + theme(
    plot.title  = element_text(face = "bold", size = plot_title_size),
    legend.position = "right",
    axis.text  = element_text(color = "gray20", size = axis_text_size),
    axis.title = element_text(color = "gray20", size = axis_title_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = strip_text_size)
  )
)

foggy_sf <- c("#4C6A87", "#7B99B6", "#9CB7CE", "#C0CEDD",
              "#9AA6B2", "#72808E", "#B8A9B4", "#C9BCC6", "#A6B8BE")
get_muted_palette <- function(n){ if (n <= length(foggy_sf)) foggy_sf[seq_len(n)] else colorRampPalette(foggy_sf)(n) }

pdf_device_live_text <- function(filename, width = 6, height = 4, ...) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    family = base_font,
    useDingbats = FALSE,
    version = "1.4",
    colormodel = "srgb",
    ...
  )
}
eps_device_live_text <- function(filename, width = 6, height = 4, ...) {
  grDevices::postscript(
    file = filename,
    width = width,
    height = height,
    onefile = FALSE,
    horizontal = FALSE,
    paper = "special",
    family = base_font,
    colormodel = "srgb",
    ...
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
read_delim_flex <- function(path) data.table::fread(path, sep = "\t", header = TRUE, data.table = TRUE)
norm_names <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
norm_chr <- function(x) gsub("^chr", "", as.character(x), ignore.case = TRUE)
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
# Generic sample matching
canonical_sample_core <- function(x) {
  s <- basename(as.character(x))
  s <- sub("(?i)(\\.aligned\\.sortedbycoord\\.out)?\\.bam$", "", s, perl = TRUE)
  s <- sub("(?i)\\.junction\\.txt$", "", s, perl = TRUE)
  toupper(gsub("[^A-Za-z0-9]+", "_", s))
}

get_bam_seqlevels <- function(bam_path) {
  hdr <- Rsamtools::scanBamHeader(bam_path)[[1]]$targets
  names(hdr)
}

match_seqname_to_bam <- function(seqname, bam_seqlevels) {
  s <- as.character(seqname)
  core <- gsub("^chr", "", s, ignore.case = TRUE)
  candidates <- unique(c(
    s,
    core,
    paste0("chr", core),
    if (core %in% c("M", "MT")) c("MT", "M", "chrM", "chrMT") else character(0)
  ))
  hit <- candidates[candidates %in% bam_seqlevels]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

convert_gr_to_bam_seqstyle <- function(gr, bam_seqlevels) {
  if (length(gr) == 0) return(gr)
  new_seq <- vapply(
    as.character(GenomicRanges::seqnames(gr)),
    match_seqname_to_bam,
    character(1),
    bam_seqlevels = bam_seqlevels
  )
  keep <- !is.na(new_seq)
  if (!any(keep)) return(gr[0])

  out <- GenomicRanges::GRanges(
    seqnames = new_seq[keep],
    ranges = GenomicRanges::ranges(gr)[keep],
    strand = GenomicRanges::strand(gr)[keep]
  )
  S4Vectors::mcols(out) <- S4Vectors::mcols(gr)[keep, , drop = FALSE]
  names(out) <- names(gr)[keep]
  out
}

empirical_p_2sided <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  if (!is.finite(x) || length(ref) < 2) return(NA_real_)
  med <- median(ref, na.rm = TRUE)
  more_extreme <- sum(abs(ref - med) >= abs(x - med), na.rm = TRUE)
  (more_extreme + 1)/(length(ref) + 1)
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
  ins_tbl <- as_tibble(read_delim_flex(insertion_path))
  names(ins_tbl) <- norm_names(names(ins_tbl))
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
    stop("Insertion file must have start/end or position columns.")
  }
  ins_tbl <- ins_tbl %>%
    mutate(gene = strsplit(ifelse(is.na(gene_raw), "", gene_raw), ",")) %>%
    unnest(gene, keep_empty = TRUE) %>%
    mutate(gene = trimws(gene), gene = na_if(gene, "")) %>%
    dplyr::select(sample_raw, sample_core, gene, chrom_in, start_in, end_in)
}

coord_tbl <- if (coord_source == "metadata") {
  meta_expanded %>% transmute(sample_raw, sample_core, group, patient_id, gene, chrom_in, start_in, end_in)
} else {
  meta_expanded %>% dplyr::select(sample_raw, sample_core, group, patient_id) %>% distinct() %>% left_join(ins_tbl, by = c("sample_raw", "sample_core"))
}

if (coord_build == "hg19") {
  if (!file.exists(chain_path)) stop("Need chain file for hg19->hg38 liftover.")
  chain <- import.chain(chain_path)
  gr_hg19 <- GRanges(seqnames = paste0("chr", norm_chr(coord_tbl$chrom_in)),
                     ranges = IRanges(coord_tbl$start_in, coord_tbl$end_in),
                     sample_raw = coord_tbl$sample_raw, sample_core = coord_tbl$sample_core,
                     group = coord_tbl$group, patient_id = coord_tbl$patient_id, gene = coord_tbl$gene)
  lo <- liftOver(gr_hg19, chain)
  lo_first <- unlist(lo[lengths(lo) >= 1], use.names = FALSE)
  coord_tbl <- tibble(
    sample_raw = mcols(lo_first)$sample_raw,
    sample_core = mcols(lo_first)$sample_core,
    group = mcols(lo_first)$group,
    patient_id = mcols(lo_first)$patient_id,
    gene = mcols(lo_first)$gene,
    chrom = norm_chr(as.character(seqnames(lo_first))),
    start = start(lo_first), end = end(lo_first)
  )
} else {
  coord_tbl <- coord_tbl %>% transmute(sample_raw, sample_core, group, patient_id, gene, chrom = norm_chr(chrom_in), start = start_in, end = end_in)
}
coord_tbl <- coord_tbl %>%
  mutate(has_interval = !is.na(chrom) & is.finite(start) & is.finite(end)) %>%
  mutate(start2 = ifelse(has_interval, pmin(start, end), start),
         end2   = ifelse(has_interval, pmax(start, end), end),
         start = start2, end = end2) %>%
  dplyr::select(-start2, -end2)
write_tsv(coord_tbl, file.path(table_dir, "target_coordinates_used.tsv"))

msg("Importing GTF and building annotated intron catalog...")
gtf <- import(gtf_path)
genes <- gtf[gtf$type == "gene"]
if (!"gene_name" %in% names(mcols(genes))) stop("GTF lacks gene_name.")
genes <- genes[!is.na(mcols(genes)$gene_name)]
seqlevelsStyle(genes) <- "UCSC"

genes_df <- tibble(gene = mcols(genes)$gene_name,
                   chrom = norm_chr(as.character(seqnames(genes))),
                   gene_start = start(genes), gene_end = end(genes)) %>% distinct()
write_tsv(genes_df, file.path(table_dir, "genes_from_gtf.tsv"))

requested_genes <- coord_tbl %>%
  dplyr::filter(!is.na(gene), nzchar(gene)) %>%
  dplyr::pull(gene) %>%
  unique()

exon_records_for_introns <- gtf[gtf$type == "exon"]
if (!"transcript_id" %in% names(S4Vectors::mcols(exon_records_for_introns))) {
  stop("GTF lacks transcript_id on exon records.")
}
# Reference splice junctions
annot_introns_df <- tibble(
  gene = as.character(S4Vectors::mcols(exon_records_for_introns)$gene_name),
  transcript_id = as.character(S4Vectors::mcols(exon_records_for_introns)$transcript_id),
  chrom = norm_chr(as.character(GenomicRanges::seqnames(exon_records_for_introns))),
  exon_start = GenomicRanges::start(exon_records_for_introns),
  exon_end = GenomicRanges::end(exon_records_for_introns)
) %>%
  dplyr::filter(
    gene %in% requested_genes,
    !is.na(transcript_id),
    nzchar(transcript_id)
  ) %>%
  dplyr::distinct(gene, transcript_id, chrom, exon_start, exon_end) %>%
  dplyr::group_by(gene, transcript_id, chrom) %>%
  dplyr::arrange(exon_start, exon_end, .by_group = TRUE) %>%
  dplyr::mutate(next_exon_start = dplyr::lead(exon_start)) %>%
  dplyr::filter(!is.na(next_exon_start), next_exon_start > exon_end + 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    gene,
    junction_id = paste0(chrom, ":", exon_end, ":", next_exon_start - 1L),
    annotated_exact = TRUE
  ) %>%
  dplyr::distinct(gene, junction_id, .keep_all = TRUE)
write_tsv(annot_introns_df, file.path(table_dir, "annotated_introns_from_gtf.tsv"))

gene_gtf_audit <- tibble::tibble(gene = requested_genes) %>%
  dplyr::mutate(
    present_in_gtf = gene %in% unique(genes_df$gene),
    included_for_bam_analysis = present_in_gtf
  )
write_tsv(gene_gtf_audit, file.path(table_dir, "metadata_gene_gtf_audit.tsv"))
msg("\n=== Metadata gene / GTF audit ===")
print(gene_gtf_audit)

missing_gtf_genes <- gene_gtf_audit %>%
  dplyr::filter(!present_in_gtf) %>%
  dplyr::pull(gene)
if (length(missing_gtf_genes) > 0) {
  warning(
    "Skipping metadata genes absent from the GTF gene_name field: ",
    paste(missing_gtf_genes, collapse = ", ")
  )
}
selected_genes <- base::intersect(requested_genes, unique(genes_df$gene))
if (length(selected_genes) == 0) stop("No target genes from metadata matched the GTF gene_name field.")
selected_genes_df <- genes_df %>% dplyr::filter(gene %in% selected_genes)
write_tsv(selected_genes_df, file.path(table_dir, "selected_target_genes.tsv"))

msg("Discovering BAM files...")
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) stop("No BAM files found in: ", bam_dir)

bam_tbl <- tibble(
  bam_path = bam_files,
  bam_file = basename(bam_files),
  sample_core = canonical_sample_core(basename(bam_files))
) %>%
  left_join(meta_expanded %>% distinct(sample_raw, sample_core, group, patient_id), by = "sample_core") %>%
  mutate(match_status = ifelse(is.na(sample_raw), "unmatched", "matched"), sample = coalesce(sample_raw, sample_core))

write_tsv(bam_tbl, file.path(table_dir, "sample_matching_audit_bams.tsv"))
msg("\n=== BAM matching audit ===")
print(dplyr::count(bam_tbl, match_status))

if (only_matched_bams) bam_tbl <- bam_tbl %>% dplyr::filter(match_status == "matched")
if (nrow(bam_tbl) == 0) stop("No matched BAMs remain after filtering.")

msg("Extracting splice junctions from BAM files (target genes only)...")

get_gene_region <- function(g) {
  x <- selected_genes_df %>% dplyr::filter(gene == g) %>% slice_head(n = 1)
  GRanges(seqnames = paste0("chr", x$chrom[1]), ranges = IRanges(start = max(1, x$gene_start[1] - fetch_pad_bp), end = x$gene_end[1] + fetch_pad_bp), gene = g)
}

gene_regions <- do.call(c, lapply(selected_genes, get_gene_region))
seqlevelsStyle(gene_regions) <- "UCSC"

extract_junctions_from_bam <- function(bam_path, sample_name, sample_core, group, patient_id) {
  msg("Processing BAM: ", basename(bam_path))
  bam_seqlevels <- get_bam_seqlevels(bam_path)
  all_j <- list()
  for (i in seq_along(gene_regions)) {
    gr <- gene_regions[i]
    gene_name <- mcols(gr)$gene
    gr_bam <- convert_gr_to_bam_seqstyle(gr, bam_seqlevels)
    if (length(gr_bam) == 0) {
      msg(
        "  Skipping ", gene_name,
        ": chromosome ", as.character(GenomicRanges::seqnames(gr))[[1]],
        " is absent from the BAM header."
      )
      next
    }
    sbp <- ScanBamParam(which = gr_bam, flag = scanBamFlag(isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE))
    ga <- tryCatch(readGAlignments(bam_path, use.names = FALSE, param = sbp), error = function(e) NULL)
    if (is.null(ga) || length(ga) == 0) next

    jg <- tryCatch(GenomicAlignments::summarizeJunctions(ga), error = function(e) NULL)
    if (is.null(jg) || length(jg) == 0) next
    df <- as.data.frame(jg)
    if (!"score" %in% names(df)) {
      stop("summarizeJunctions() did not return its expected score column for ", basename(bam_path))
    }
    out <- tibble(
      sample = sample_name,
      sample_core = sample_core,
      group = group,
      patient_id = patient_id,
      gene = gene_name,
      chrom = norm_chr(as.character(df$seqnames)),
      start = df$start - 1,
      end = df$end,
      count = as.numeric(df$score),
      junction_id = paste0(norm_chr(as.character(df$seqnames)), ":", df$start - 1, ":", df$end)
    )
    all_j[[length(all_j)+1]] <- out
  }
  bind_rows(all_j)
}

bam_junctions <- pmap_dfr(
  list(bam_tbl$bam_path, bam_tbl$sample, bam_tbl$sample_core, bam_tbl$group, bam_tbl$patient_id),
  extract_junctions_from_bam
)

if (nrow(bam_junctions) == 0) stop("No splice junctions were extracted from BAMs for the selected target gene regions.")

bam_junctions <- bam_junctions %>%
  left_join(annot_introns_df, by = c("gene", "junction_id")) %>%
  mutate(annotation = ifelse(!is.na(annotated_exact) & annotated_exact, "annotated_exact", "novel_or_nonexact"))
write_tsv(bam_junctions, file.path(table_dir, "bam_extracted_junctions.tsv"))

indexed_junctions <- bam_junctions %>%
  mutate(.junction_row_id = row_number())

target_flags <- indexed_junctions %>%
  dplyr::select(.junction_row_id, sample_core, gene, chrom, start, end) %>%
  left_join(
    coord_tbl %>%
      dplyr::select(
        sample_core,
        gene,
        chrom_target = chrom,
        start_target = start,
        end_target = end,
        has_interval
      ),
    by = c("sample_core", "gene")
  ) %>%
  mutate(
    spans_target_interval = ifelse(has_interval & chrom == chrom_target, start <= end_target & end >= start_target, FALSE),
    donor_near_target = ifelse(has_interval & chrom == chrom_target, abs(start - start_target) <= local_pad_bp | abs(start - end_target) <= local_pad_bp, FALSE),
    acceptor_near_target = ifelse(has_interval & chrom == chrom_target, abs(end - start_target) <= local_pad_bp | abs(end - end_target) <= local_pad_bp, FALSE),
    near_target = spans_target_interval | donor_near_target | acceptor_near_target
  ) %>%
  group_by(.junction_row_id) %>%
  summarise(
    has_interval = any(has_interval %in% TRUE),
    spans_target_interval = any(spans_target_interval %in% TRUE),
    donor_near_target = any(donor_near_target %in% TRUE),
    acceptor_near_target = any(acceptor_near_target %in% TRUE),
    near_target = any(near_target %in% TRUE),
    .groups = "drop"
  )

analysis_rows <- indexed_junctions %>%
  left_join(target_flags, by = ".junction_row_id") %>%
  left_join(
    selected_genes_df %>% dplyr::select(gene, gene_start, gene_end),
    by = "gene"
  ) %>%
  dplyr::select(-.junction_row_id)

usage_rows <- analysis_rows %>%
  group_by(sample, sample_core, group, patient_id, gene) %>%
  mutate(gene_total_reads = sum(count, na.rm = TRUE), usage = count / pmax(gene_total_reads, 1)) %>%
  ungroup()

sample_gene <- usage_rows %>%
  group_by(sample, sample_core, group, patient_id, gene) %>%
  summarise(
    gene_total_reads = sum(count, na.rm = TRUE),
    n_junctions = n_distinct(junction_id),
    annotated_reads = sum(count[annotation == "annotated_exact"], na.rm = TRUE),
    novel_reads     = sum(count[annotation != "annotated_exact"], na.rm = TRUE),
    novel_fraction  = novel_reads / pmax(gene_total_reads, 1),
    local_reads     = sum(count[near_target], na.rm = TRUE),
    local_fraction  = local_reads / pmax(gene_total_reads, 1),
    local_novel_reads = sum(count[near_target & annotation != "annotated_exact"], na.rm = TRUE),
    local_novel_fraction = local_novel_reads / pmax(local_reads, 1),
    entropy = -sum(usage * log(usage + 1e-6), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(gene_total_reads >= min_gene_total_reads)
write_tsv(sample_gene, file.path(table_dir, "sample_gene_metrics_from_bam.tsv"))

msg("\n=== sample_gene summary from BAM ===")
print(sample_gene %>% summarise(
  n_rows = n(), n_genes = n_distinct(gene), n_samples = n_distinct(sample),
  median_gene_reads = median(gene_total_reads, na.rm = TRUE),
  median_novel_fraction = median(novel_fraction, na.rm = TRUE)
))

case_tbl <- coord_tbl %>% dplyr::filter(!is.na(gene), !toupper(group) %in% c("CONTROL")) %>% distinct(sample_core, sample_raw, gene, group)
control_meta <- meta_expanded %>% dplyr::filter(toupper(group) %in% c("CONTROL")) %>% distinct(sample_core, sample_raw, group)

single_case <- list(); groupwise <- list(); per_junction <- list()
for (g in unique(case_tbl$gene)) {
  sg <- sample_gene %>% dplyr::filter(gene == g)
  if (nrow(sg) == 0) next
  case_samples <- base::intersect(case_tbl$sample_core[case_tbl$gene == g], sg$sample_core)
  ctrl_samples <- base::intersect(control_meta$sample_core, sg$sample_core)
  if (length(case_samples) == 0 || length(ctrl_samples) < min_controls_per_gene) next

  for (s in case_samples) {
    one <- sg %>% dplyr::filter(sample_core == s)
    ref <- sg %>% dplyr::filter(sample_core %in% ctrl_samples)
    if (nrow(one) == 0 || nrow(ref) < min_controls_per_gene) next
    single_case[[length(single_case) + 1]] <- tibble(
      gene = g, sample = one$sample[1], group = one$group[1], n_controls = nrow(ref),
      gene_total_reads = one$gene_total_reads[1],
      entropy = one$entropy[1], entropy_z = make_z(one$entropy[1], ref$entropy), entropy_emp_p = empirical_p_2sided(one$entropy[1], ref$entropy),
      novel_fraction = one$novel_fraction[1], novel_fraction_z = make_z(one$novel_fraction[1], ref$novel_fraction), novel_fraction_emp_p = empirical_p_2sided(one$novel_fraction[1], ref$novel_fraction),
      local_fraction = one$local_fraction[1], local_fraction_z = make_z(one$local_fraction[1], ref$local_fraction), local_fraction_emp_p = empirical_p_2sided(one$local_fraction[1], ref$local_fraction),
      local_novel_fraction = one$local_novel_fraction[1], local_novel_fraction_z = make_z(one$local_novel_fraction[1], ref$local_novel_fraction), local_novel_fraction_emp_p = empirical_p_2sided(one$local_novel_fraction[1], ref$local_novel_fraction)
    )
  }

  if (length(case_samples) >= 2) {
    cases <- sg %>% dplyr::filter(sample_core %in% case_samples)
    ctrls <- sg %>% dplyr::filter(sample_core %in% ctrl_samples)
    gtest <- function(x, y) tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    groupwise[[length(groupwise) + 1]] <- tibble(
      gene = g, n_cases = nrow(cases), n_controls = nrow(ctrls),
      entropy_p = gtest(cases$entropy, ctrls$entropy), novel_fraction_p = gtest(cases$novel_fraction, ctrls$novel_fraction),
      local_fraction_p = gtest(cases$local_fraction, ctrls$local_fraction), local_novel_fraction_p = gtest(cases$local_novel_fraction, ctrls$local_novel_fraction)
    )
  }

  ug_observed <- usage_rows %>% dplyr::filter(gene == g)
  top_j <- ug_observed %>%
    group_by(junction_id) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = max(25, label_top_junc_n)) %>%
    pull(junction_id)

  comparison_samples <- sg %>%
    dplyr::filter(sample_core %in% c(case_samples, ctrl_samples)) %>%
    dplyr::select(sample_core, sample) %>%
    distinct(sample_core, .keep_all = TRUE)
  junction_meta <- ug_observed %>%
    dplyr::filter(junction_id %in% top_j) %>%
    group_by(junction_id) %>%
    summarise(
      annotation = dplyr::first(annotation),
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
        summarise(count = dplyr::first(count), usage = dplyr::first(usage), .groups = "drop"),
      by = c("sample_core", "junction_id")
    ) %>%
    mutate(
      count = replace_na(count, 0),
      usage = replace_na(usage, 0)
    )

  for (jid in unique(ug$junction_id)) {
    dd <- ug %>% dplyr::filter(junction_id == jid)
    ref_usage <- dd %>% dplyr::filter(sample_core %in% ctrl_samples) %>% pull(usage)
    for (s in case_samples) {
      obs <- dd %>% dplyr::filter(sample_core == s)
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
        annotation = obs$annotation[1],
        inserted_usage = obs$usage[1],
        inserted_count = obs$count[1],
        n_controls_total = length(ref_usage),
        n_controls_with_junction = sum(ref_usage > 0, na.rm = TRUE),
        control_mean_usage = mean(ref_usage, na.rm = TRUE),
        control_median_usage = median(ref_usage, na.rm = TRUE),
        usage_z = make_z(obs$usage[1], ref_usage),
        usage_emp_p = empirical_p_2sided(obs$usage[1], ref_usage),
        near_target = near_target,
        spans_target_interval = spans_target
      )
    }
  }
}

single_case_df <- bind_rows(single_case)
if (nrow(single_case_df) > 0) single_case_df <- single_case_df %>% mutate(
  entropy_fdr = p.adjust(entropy_emp_p, "BH"), novel_fraction_fdr = p.adjust(novel_fraction_emp_p, "BH"),
  local_fraction_fdr = p.adjust(local_fraction_emp_p, "BH"), local_novel_fraction_fdr = p.adjust(local_novel_fraction_emp_p, "BH")
)
write_tsv(single_case_df, file.path(table_dir, "single_case_disruption_results_from_bam.tsv"))

groupwise_df <- bind_rows(groupwise)
if (nrow(groupwise_df) > 0) groupwise_df <- groupwise_df %>% mutate(
  entropy_fdr = p.adjust(entropy_p, "BH"), novel_fraction_fdr = p.adjust(novel_fraction_p, "BH"),
  local_fraction_fdr = p.adjust(local_fraction_p, "BH"), local_novel_fraction_fdr = p.adjust(local_novel_fraction_p, "BH")
)
write_tsv(groupwise_df, file.path(table_dir, "groupwise_disruption_results_from_bam.tsv"))

per_junction_df <- bind_rows(per_junction)
if (nrow(per_junction_df) > 0) per_junction_df <- per_junction_df %>% mutate(usage_fdr = p.adjust(usage_emp_p, "BH"))
write_tsv(per_junction_df, file.path(table_dir, "per_junction_case_control_results_from_bam.tsv"))

msg("\n=== Top BAM-based single-case results ===")
if (nrow(single_case_df) > 0) print(single_case_df %>% arrange(local_fraction_fdr, entropy_fdr) %>% slice_head(n = 20))

ann_df <- usage_rows %>%
  mutate(annotation_simple = ifelse(annotation == "annotated_exact", "annotated_exact", "novel_or_nonexact")) %>%
  group_by(sample, annotation_simple) %>% summarise(reads = sum(count), .groups = "drop") %>%
  group_by(sample) %>% mutate(frac = reads / sum(reads)) %>% ungroup()

p_ann <- ggplot(ann_df, aes(x = reorder(sample, frac, FUN = sum), y = frac, fill = annotation_simple)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("annotated_exact" = "#4C6A87", "novel_or_nonexact" = "#9CB7CE")) +
  labs(title = "BAM-derived junction annotation composition by sample", x = "Sample", y = "Fraction of junction reads", fill = NULL)
save_plot_all(p_ann, "qc_bam_annotation_fraction_by_sample", width = 7.8, height = 6.8)

sample_gene_plot <- sample_gene %>% mutate(GroupClass = ifelse(toupper(group) == "CONTROL", "Control", "Inserted/Case"))
for (metric in c("novel_fraction", "entropy", "local_fraction", "local_novel_fraction")) {
  p <- ggplot(sample_gene_plot, aes(x = GroupClass, y = .data[[metric]], fill = GroupClass)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) + geom_jitter(width = 0.15, size = 1.7, alpha = 0.75) +
    scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
    labs(title = paste0("BAM-derived ", metric, " by group"), x = NULL, y = metric) + theme(legend.position = "none")
  save_plot_all(p, paste0("global_bam_", metric, "_by_group"), width = 5.8, height = 4.8)
}

top_single_genes <- if (nrow(single_case_df) > 0) {
  single_case_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
} else {
  character()
}
top_group_genes <- if (nrow(groupwise_df) > 0) {
  groupwise_df %>%
    arrange(local_fraction_fdr, entropy_fdr) %>%
    slice_head(n = 20) %>%
    pull(gene)
} else {
  character()
}
plot_genes <- unique(c(top_single_genes, top_group_genes))
plot_genes <- plot_genes[!is.na(plot_genes)]
for (g in plot_genes) {
  sg <- sample_gene_plot %>% dplyr::filter(gene == g)
  if (nrow(sg) == 0) next
  sg_long <- sg %>% dplyr::select(sample, GroupClass, entropy, novel_fraction, local_fraction, local_novel_fraction) %>%
    pivot_longer(cols = c(entropy, novel_fraction, local_fraction, local_novel_fraction), names_to = "metric", values_to = "value")
  p_metrics <- ggplot(sg_long, aes(x = GroupClass, y = value, fill = GroupClass)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95) + geom_jitter(width = 0.15, size = 1.6, alpha = 0.75) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = c("Control" = "#C0CEDD", "Inserted/Case" = "#4C6A87")) +
    labs(title = paste0(g, ": BAM-derived splice disruption metrics"), x = NULL, y = NULL) +
    theme(legend.position = "none")
  save_plot_all(p_metrics, paste0("bam_gene_metrics_", g), width = 8.2, height = 5.6)

  ug <- usage_rows %>% dplyr::filter(gene == g)
  top_j <- ug %>% group_by(junction_id) %>% summarise(total = sum(count), .groups = "drop") %>% arrange(desc(total)) %>% slice_head(n = label_top_junc_n) %>% pull(junction_id)
  hm_df <- ug %>% dplyr::filter(junction_id %in% top_j) %>% dplyr::select(junction_id, sample, usage) %>% pivot_wider(names_from = sample, values_from = usage, values_fill = 0)
  if (nrow(hm_df) >= 2 && ncol(hm_df) >= 3) {
    hm_mat <- as.data.frame(hm_df); rownames(hm_mat) <- hm_mat$junction_id; hm_mat$junction_id <- NULL; hm_mat <- as.matrix(hm_mat)
    ann_col <- sg %>% distinct(sample, GroupClass) %>% as.data.frame(); rownames(ann_col) <- ann_col$sample; ann_col$sample <- NULL
    ann_col <- ann_col[colnames(hm_mat), , drop = FALSE]
    pheat <- pheatmap(hm_mat, scale = "row", annotation_col = ann_col,
                      clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                      main = paste0(g, ": BAM-derived top junction usage"),
                      fontsize = heatmap_text_size,
                      fontsize_row = heatmap_text_size,
                      fontsize_col = heatmap_text_size,
                      angle_col = 45,
                      silent = TRUE)
    save_pheatmap_all(pheat, paste0("bam_gene_heatmap_", g), width = 8.2, height = 7.2)
  }
}

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
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
msg("\nDone. Outputs written to: ", normalizePath(out_dir))
