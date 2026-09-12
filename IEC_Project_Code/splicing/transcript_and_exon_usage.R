#!/usr/bin/env Rscript
# Transcript and exon usage
# Measure BAM-derived junction usage, transcript annotations, and exon coverage at insertion-associated genes.

args <- commandArgs(trailingOnly = TRUE)
default_base_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
base_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_base_dir
base_dir <- normalizePath(path.expand(base_dir), mustWork = TRUE)

bam_dir        <- file.path(base_dir, "bam")
meta_path <- file.path(base_dir, "splicing_metadata.tsv")
gtf_path <- file.path(base_dir, "GRCh38.gtf.gz")

# Gene selection
target_genes <- NULL

coord_source <- "metadata"
coord_build  <- "hg38"
only_matched_bams <- TRUE

# Coverage and reference parameters
fetch_pad_bp <- 1000
min_gene_total_junction_reads <- 5
min_gene_total_exon_reads <- 10
min_reference_samples <- 3
reference_mode <- "all_other"
local_pad_bp <- 1000
boundary_tol_bp <- 2
label_top_junc_n <- 20
min_mapq <- 0

altered_abs_usage_z <- 1.5
altered_abs_logcount_z <- 1.5

base_font_size    <- 14
plot_title_size   <- 16
axis_title_size   <- 14
axis_text_size    <- 12
legend_title_size <- 13
legend_text_size  <- 12
strip_text_size   <- 13
heatmap_text_size <- 12

out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "splicing", "transcript_and_exon_usage")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
log_dir   <- file.path(out_dir, "logs")
dirs <- c(out_dir, plot_dir, table_dir, log_dir)
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

meta_sample_col   <- "sample"
meta_paper_col    <- "paper_sample_id"
meta_patient_col  <- "patient_id"
meta_group_col    <- "group"
meta_gene_col     <- "gene"
meta_chr_col      <- "chr"
meta_start_col    <- "start"
meta_end_col      <- "end"

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
  theme_classic(base_size = base_font_size, base_family = base_font) +
    theme(
      plot.title = element_text(face = "bold", size = plot_title_size),
      axis.text = element_text(color = "gray20", size = axis_text_size),
      axis.title = element_text(color = "gray20", size = axis_title_size),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = strip_text_size),
      legend.position = "right"
    )
)

foggy_sf <- c("#4C6A87", "#7B99B6", "#9CB7CE", "#C0CEDD",
              "#9AA6B2", "#72808E", "#B8A9B4", "#C9BCC6", "#A6B8BE")
get_muted_palette <- function(n) {
  if (n <= length(foggy_sf)) return(foggy_sf[seq_len(n)])
  grDevices::colorRampPalette(foggy_sf)(n)
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
  candidates <- unique(c(
    s,
    gsub("^chr", "", s, ignore.case = TRUE),
    paste0("chr", gsub("^chr", "", s, ignore.case = TRUE)),
    ifelse(s %in% c("chrM", "M"), "MT", s),
    ifelse(s %in% c("chrMT", "MT"), "M", s),
    ifelse(s %in% c("M", "MT"), "chrM", s)
  ))
  hit <- candidates[candidates %in% bam_seqlevels]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

convert_gr_to_bam_seqstyle <- function(gr, bam_seqlevels) {
  if (length(gr) == 0) return(gr)
  new_seq <- vapply(as.character(GenomicRanges::seqnames(gr)), match_seqname_to_bam, character(1), bam_seqlevels = bam_seqlevels)
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

make_granges_for_bam <- function(chrom, start, end, bam_seqlevels) {
  seqn <- vapply(paste0("chr", norm_chr(chrom)), match_seqname_to_bam, character(1), bam_seqlevels = bam_seqlevels)
  keep <- !is.na(seqn)
  if (!any(keep)) return(GenomicRanges::GRanges())
  out <- GenomicRanges::GRanges(
    seqnames = seqn[keep],
    ranges = IRanges::IRanges(start = start[keep], end = end[keep])
  )

  S4Vectors::mcols(out) <- S4Vectors::DataFrame(
    .source_row = as.integer(which(keep))
  )
  out
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
  mu <- mean(ref, na.rm = TRUE)
  sdv <- stats::sd(ref, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(NA_real_)
  (x - mu) / sdv
}
collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}
safe_mcol <- function(gr, colname, default = NA_character_) {
  if (colname %in% names(S4Vectors::mcols(gr))) {
    as.character(S4Vectors::mcols(gr)[[colname]])
  } else {
    rep(default, length(gr))
  }
}

msg("Loading metadata...")
meta <- as_tibble(read_delim_flex(meta_path))
names(meta) <- norm_names(names(meta))
needed_meta <- c(meta_sample_col, meta_group_col, meta_gene_col)
missing_meta <- base::setdiff(needed_meta, names(meta))
if (length(missing_meta) > 0) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))

meta <- meta %>%
  dplyr::mutate(
    sample_raw = .data[[meta_sample_col]],
    sample_core = canonical_sample_core(.data[[meta_sample_col]]),
    group = as.character(.data[[meta_group_col]]),
    paper_sample_id = if (meta_paper_col %in% names(meta)) as.character(.data[[meta_paper_col]]) else NA_character_,
    patient_id = if (meta_patient_col %in% names(meta)) as.character(.data[[meta_patient_col]]) else NA_character_,
    gene_raw = as.character(.data[[meta_gene_col]]),
    chrom_in = if (meta_chr_col %in% names(meta)) norm_chr(.data[[meta_chr_col]]) else NA_character_,
    start_in = if (meta_start_col %in% names(meta)) safe_num(.data[[meta_start_col]]) else NA_real_,
    end_in   = if (meta_end_col %in% names(meta)) safe_num(.data[[meta_end_col]]) else NA_real_,
    display_label = dplyr::coalesce(
      if (meta_paper_col %in% names(meta)) as.character(.data[[meta_paper_col]]) else NA_character_,
      if (meta_patient_col %in% names(meta)) as.character(.data[[meta_patient_col]]) else NA_character_,
      as.character(.data[[meta_sample_col]])
    )
  ) %>%
  dplyr::mutate(gene_raw = ifelse(gene_raw %in% c("NA", "", "Na", "na"), NA_character_, gene_raw))

meta_expanded <- meta %>%
  dplyr::mutate(gene = strsplit(ifelse(is.na(gene_raw), "", gene_raw), ",")) %>%
  tidyr::unnest(gene, keep_empty = TRUE) %>%
  dplyr::mutate(gene = trimws(gene), gene = dplyr::na_if(gene, ""))

if (!is.null(target_genes)) {
  selected_genes <- target_genes
} else {
  selected_genes <- meta_expanded %>% dplyr::filter(!is.na(gene)) %>% dplyr::pull(gene) %>% unique()
}
selected_genes <- unique(selected_genes[!is.na(selected_genes)])
if (length(selected_genes) == 0) stop("No selected genes to analyze.")
msg("Selected genes: ", paste(selected_genes, collapse = ", "))

write_tsv(meta_expanded, file.path(table_dir, "metadata_expanded_gene_rows.tsv"))

coord_tbl <- meta_expanded %>%
  dplyr::filter(gene %in% selected_genes) %>%
  dplyr::transmute(
    sample_raw, sample_core, group, paper_sample_id, patient_id, display_label,
    gene, chrom = norm_chr(chrom_in), start = start_in, end = end_in,
    liftover_n = NA_integer_
  ) %>%
  dplyr::mutate(
    has_interval = !is.na(chrom) & is.finite(start) & is.finite(end),
    start2 = ifelse(has_interval, pmin(start, end), start),
    end2 = ifelse(has_interval, pmax(start, end), end),
    start = start2,
    end = end2
  ) %>%
  dplyr::select(-start2, -end2)
write_tsv(coord_tbl, file.path(table_dir, "target_coordinates_used.tsv"))

msg("Importing GTF and building gene/transcript models...")
gtf <- rtracklayer::import(gtf_path)
seqlevelsStyle(gtf) <- "UCSC"

genes <- gtf[gtf$type == "gene"]
if (!"gene_name" %in% names(S4Vectors::mcols(genes))) stop("GTF lacks gene_name.")
genes <- genes[!is.na(S4Vectors::mcols(genes)$gene_name)]

genes_df <- tibble(
  gene = as.character(S4Vectors::mcols(genes)$gene_name),
  gene_id = safe_mcol(genes, "gene_id"),
  chrom = norm_chr(as.character(GenomicRanges::seqnames(genes))),
  strand = as.character(GenomicRanges::strand(genes)),
  gene_start = GenomicRanges::start(genes),
  gene_end = GenomicRanges::end(genes)
) %>%
  dplyr::filter(gene %in% selected_genes) %>%
  dplyr::group_by(gene, chrom, strand) %>%
  dplyr::summarise(gene_id = dplyr::first(gene_id), gene_start = min(gene_start), gene_end = max(gene_end), .groups = "drop")
if (nrow(genes_df) == 0) stop("Selected genes did not match GTF gene_name values.")

requested_genes <- selected_genes
gene_gtf_audit <- tibble::tibble(gene = requested_genes) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    present_in_gtf = gene %in% unique(genes_df$gene),
    included_for_bam_analysis = present_in_gtf
  )
readr::write_tsv(gene_gtf_audit, file.path(table_dir, "metadata_gene_gtf_audit.tsv"))
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
if (length(selected_genes) == 0) {
  stop("No metadata target genes remain after matching the GTF.")
}
write_tsv(genes_df, file.path(table_dir, "selected_gene_bodies.tsv"))

exon_records <- gtf[gtf$type == "exon"]
cds_records <- gtf[gtf$type == "CDS"]
transcript_records <- gtf[gtf$type == "transcript"]

# Transcript models
exon_model <- tibble(
  gene = safe_mcol(exon_records, "gene_name"),
  gene_id = safe_mcol(exon_records, "gene_id"),
  transcript_id = safe_mcol(exon_records, "transcript_id"),
  transcript_name = safe_mcol(exon_records, "transcript_name"),
  transcript_biotype = dplyr::coalesce(safe_mcol(exon_records, "transcript_biotype"), safe_mcol(exon_records, "transcript_type")),
  exon_number_raw = safe_mcol(exon_records, "exon_number"),
  chrom = norm_chr(as.character(GenomicRanges::seqnames(exon_records))),
  strand = as.character(GenomicRanges::strand(exon_records)),
  exon_start = GenomicRanges::start(exon_records),
  exon_end = GenomicRanges::end(exon_records)
) %>%
  dplyr::filter(gene %in% selected_genes, !is.na(transcript_id), transcript_id != "") %>%
  dplyr::mutate(exon_number = dplyr::if_else(is.na(exon_number_raw) | exon_number_raw == "", NA_character_, exon_number_raw)) %>%
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
  dplyr::filter(gene %in% selected_genes, !is.na(transcript_id), transcript_id != "")

transcript_summary <- exon_model %>%
  dplyr::group_by(gene, transcript_id, transcript_name, transcript_biotype, chrom, strand) %>%
  dplyr::summarise(
    n_exons = dplyr::n(),
    tx_start = min(exon_start),
    tx_end = max(exon_end),
    cds_bases = {
      tx <- dplyr::first(transcript_id)
      sum((cds_model %>% dplyr::filter(transcript_id == tx) %>% dplyr::mutate(w = cds_end - cds_start + 1) %>% dplyr::pull(w)), na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(is_protein_coding = transcript_biotype %in% c("protein_coding", "protein coding")) %>%
  dplyr::arrange(gene, dplyr::desc(is_protein_coding), dplyr::desc(cds_bases), dplyr::desc(n_exons))

write_tsv(exon_model, file.path(table_dir, "transcript_exon_model.tsv"))
write_tsv(cds_model, file.path(table_dir, "transcript_cds_model.tsv"))
write_tsv(transcript_summary, file.path(table_dir, "transcript_summary.tsv"))

transcript_introns <- exon_model %>%
  dplyr::group_by(gene, transcript_id) %>%
  dplyr::arrange(exon_rank_genomic, .by_group = TRUE) %>%
  dplyr::mutate(
    next_exon_start = dplyr::lead(exon_start),
    next_exon_number = dplyr::lead(exon_number),
    next_exon_rank_genomic = dplyr::lead(exon_rank_genomic)
  ) %>%
  dplyr::filter(!is.na(next_exon_start), next_exon_start > exon_end + 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    gene, transcript_id, transcript_name, transcript_biotype, chrom, strand,
    left_exon_number = exon_number,
    right_exon_number = next_exon_number,
    left_exon_rank_genomic = exon_rank_genomic,
    right_exon_rank_genomic = next_exon_rank_genomic,
    left_exon_end = exon_end,
    right_exon_start = next_exon_start,
    intron_start = exon_end + 1,
    intron_end = next_exon_start - 1,
    junction_id = paste0(chrom, ":", left_exon_end, ":", intron_end)
  )
write_tsv(transcript_introns, file.path(table_dir, "transcript_introns_catalog.tsv"))

exon_bin_rows <- list()
for (g in selected_genes) {
  eg <- exon_records[safe_mcol(exon_records, "gene_name") == g]
  if (length(eg) == 0) next
  bins <- GenomicRanges::disjoin(eg, ignore.strand = TRUE)
  exon_bin_rows[[length(exon_bin_rows) + 1]] <- tibble(
    gene = g,
    chrom = norm_chr(as.character(GenomicRanges::seqnames(bins))),
    exon_bin_start = GenomicRanges::start(bins),
    exon_bin_end = GenomicRanges::end(bins)
  ) %>%
    dplyr::arrange(exon_bin_start, exon_bin_end) %>%
    dplyr::mutate(exon_bin_rank = dplyr::row_number(),
                  exon_bin_label = paste0(g, "_bin", exon_bin_rank, "_", exon_bin_start, "_", exon_bin_end))
}
exon_bins <- dplyr::bind_rows(exon_bin_rows)
write_tsv(exon_bins, file.path(table_dir, "exon_bin_model.tsv"))

msg("Discovering BAM files...")
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) stop("No BAM files found in: ", bam_dir)

bam_tbl <- tibble(
  bam_path = bam_files,
  bam_file = basename(bam_files),
  sample_core = canonical_sample_core(bam_files)
) %>%
  dplyr::left_join(
    meta_expanded %>% dplyr::distinct(sample_raw, sample_core, group, paper_sample_id, patient_id, display_label),
    by = "sample_core"
  ) %>%
  dplyr::mutate(match_status = ifelse(is.na(sample_raw), "unmatched", "matched"),
                sample = dplyr::coalesce(sample_raw, sample_core))
write_tsv(bam_tbl, file.path(table_dir, "sample_matching_audit_bams.tsv"))

bam_seqlevel_audit <- purrr::map_dfr(bam_tbl$bam_path, function(bp) {
  seqs <- tryCatch(get_bam_seqlevels(bp), error = function(e) character(0))
  tibble(
    bam_file = basename(bp),
    has_chr17 = "chr17" %in% seqs,
    has_17 = "17" %in% seqs,
    has_chr_prefix_any = any(grepl("^chr", seqs)),
    first_seqlevels = paste(utils::head(seqs, 12), collapse = ";")
  )
})
write_tsv(bam_seqlevel_audit, file.path(table_dir, "bam_seqlevel_style_audit.tsv"))

msg("\n=== BAM matching audit ===")
print(bam_tbl %>% dplyr::count(match_status))
msg("\n=== BAM seqlevel style audit ===")
print(bam_seqlevel_audit %>% dplyr::count(has_chr17, has_17, has_chr_prefix_any))
if (only_matched_bams) bam_tbl <- bam_tbl %>% dplyr::filter(match_status == "matched")
if (nrow(bam_tbl) == 0) stop("No matched BAMs remain after filtering.")

make_gene_region <- function(g) {
  x <- genes_df %>% dplyr::filter(gene == g) %>% dplyr::slice_head(n = 1)
  if (nrow(x) == 0) stop("Missing gene body for ", g)
  GenomicRanges::GRanges(
    seqnames = paste0("chr", x$chrom[1]),
    ranges = IRanges::IRanges(start = max(1, x$gene_start[1] - fetch_pad_bp), end = x$gene_end[1] + fetch_pad_bp),
    gene = g
  )
}
gene_regions <- do.call(c, lapply(selected_genes, make_gene_region))
seqlevelsStyle(gene_regions) <- "UCSC"

msg("Extracting BAM junctions and exon coverage over selected gene regions...")

extract_bam_gene <- function(bam_path, sample, sample_core, group, paper_sample_id, patient_id, display_label) {
  msg("Processing BAM: ", basename(bam_path))
  bam_seqlevels <- get_bam_seqlevels(bam_path)
  junction_list <- list()
  exon_bin_cov_list <- list()
  tx_exon_cov_list <- list()
  extraction_audit_list <- list()

  for (i in seq_along(gene_regions)) {
    gr <- gene_regions[i]
    gene_name <- as.character(S4Vectors::mcols(gr)$gene)
    gr_bam <- convert_gr_to_bam_seqstyle(gr, bam_seqlevels)
    if (length(gr_bam) == 0) {
      msg("  Skipping ", gene_name, ": no matching chromosome name in BAM header for requested seqlevel ",
          as.character(GenomicRanges::seqnames(gr))[1])
      extraction_audit_list[[length(extraction_audit_list) + 1]] <- tibble::tibble(
        sample = sample,
        sample_core = sample_core,
        gene = gene_name,
        chromosome_in_bam = FALSE,
        bam_read_status = "chromosome_absent_from_bam",
        n_gene_region_alignments = NA_integer_,
        n_exon_bins_in_bam = 0L,
        n_transcript_exons_in_bam = 0L
      )
      next
    }

    sbp <- Rsamtools::ScanBamParam(
      which = gr_bam,
      flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE,
                                    isSupplementaryAlignment = FALSE,
                                    isUnmappedQuery = FALSE),
      mapqFilter = min_mapq
    )

    ga <- tryCatch(
      GenomicAlignments::readGAlignments(bam_path, use.names = FALSE, param = sbp),
      error = function(e) {
        msg("  BAM read error in ", basename(bam_path), " for ", gene_name, ": ", e$message)
        NULL
      }
    )
    if (is.null(ga)) {
      extraction_audit_list[[length(extraction_audit_list) + 1]] <- tibble::tibble(
        sample = sample,
        sample_core = sample_core,
        gene = gene_name,
        chromosome_in_bam = TRUE,
        bam_read_status = "bam_read_error",
        n_gene_region_alignments = NA_integer_,
        n_exon_bins_in_bam = NA_integer_,
        n_transcript_exons_in_bam = NA_integer_
      )
      next
    }

    jg <- if (length(ga) > 0) {
      tryCatch(
        GenomicAlignments::summarizeJunctions(ga),
        error = function(e) {
          msg(
            "  Junction summarization error in ", basename(bam_path),
            " for ", gene_name, ": ", e$message
          )
          NULL
        }
      )
    } else {
      NULL
    }
    if (!is.null(jg) && length(jg) > 0) {
      jdf <- as.data.frame(jg)
      if (nrow(jdf) == 0) {
        msg("  No spliced junctions detected for ", gene_name, ".")
      } else {
        if (!"score" %in% names(jdf)) {
          jdf <- dplyr::mutate(
            tibble::as_tibble(jdf),
            score = rep(1, nrow(jdf))
          )
        }
        jout <- tibble(
          sample = sample,
          sample_core = sample_core,
          group = group,
          paper_sample_id = paper_sample_id,
          patient_id = patient_id,
          display_label = display_label,
          gene = gene_name,
          chrom = norm_chr(as.character(jdf$seqnames)),
          start = jdf$start - 1,
          end = jdf$end,
          count = jdf$score,
          junction_id = paste0(
            norm_chr(as.character(jdf$seqnames)),
            ":", jdf$start - 1, ":", jdf$end
          )
        ) %>%
          dplyr::group_by(
            sample, sample_core, group, paper_sample_id, patient_id,
            display_label, gene, chrom, start, end, junction_id
          ) %>%
          dplyr::summarise(
            count = sum(count, na.rm = TRUE),
            .groups = "drop"
          )
        junction_list[[length(junction_list) + 1]] <- jout
      }
    } else if (length(ga) > 0) {
      msg("  Reads present but no spliced junctions detected for ", gene_name, ".")
    }

    n_exon_bins_in_bam <- 0L
    bins_g <- exon_bins %>% dplyr::filter(gene == gene_name)
    if (nrow(bins_g) > 0) {
      bins_gr <- make_granges_for_bam(bins_g$chrom, bins_g$exon_bin_start, bins_g$exon_bin_end, bam_seqlevels)
      if (length(bins_gr) > 0) {
        source_rows <- as.integer(S4Vectors::mcols(bins_gr)$.source_row)
        bins_g <- bins_g[source_rows, , drop = FALSE]
        n_exon_bins_in_bam <- length(bins_gr)
        cov_counts <- GenomicRanges::countOverlaps(bins_gr, ga, ignore.strand = TRUE)
        stopifnot(length(cov_counts) == nrow(bins_g))
        exon_bin_cov_list[[length(exon_bin_cov_list) + 1]] <- bins_g %>%
          dplyr::mutate(
            sample = sample,
            sample_core = sample_core,
            group = group,
            paper_sample_id = paper_sample_id,
            patient_id = patient_id,
            display_label = display_label,
            count = as.numeric(cov_counts)
          )
      } else {
        msg("  No exon-bin intervals for ", gene_name, " matched this BAM header.")
      }
    }

    n_transcript_exons_in_bam <- 0L
    ex_tx_g <- exon_model %>% dplyr::filter(gene == gene_name)
    if (nrow(ex_tx_g) > 0) {
      ex_gr <- make_granges_for_bam(ex_tx_g$chrom, ex_tx_g$exon_start, ex_tx_g$exon_end, bam_seqlevels)
      if (length(ex_gr) > 0) {
        source_rows <- as.integer(S4Vectors::mcols(ex_gr)$.source_row)
        ex_tx_g <- ex_tx_g[source_rows, , drop = FALSE]
        n_transcript_exons_in_bam <- length(ex_gr)
        tx_counts <- GenomicRanges::countOverlaps(ex_gr, ga, ignore.strand = TRUE)
        stopifnot(length(tx_counts) == nrow(ex_tx_g))
        tx_exon_cov_list[[length(tx_exon_cov_list) + 1]] <- ex_tx_g %>%
          dplyr::mutate(
            sample = sample,
            sample_core = sample_core,
            group = group,
            paper_sample_id = paper_sample_id,
            patient_id = patient_id,
            display_label = display_label,
            count = as.numeric(tx_counts)
          )
      } else {
        msg("  No transcript-exon intervals for ", gene_name, " matched this BAM header.")
      }
    }

    extraction_audit_list[[length(extraction_audit_list) + 1]] <- tibble::tibble(
      sample = sample,
      sample_core = sample_core,
      gene = gene_name,
      chromosome_in_bam = TRUE,
      bam_read_status = ifelse(length(ga) == 0, "no_reads_in_gene_region", "reads_present"),
      n_gene_region_alignments = length(ga),
      n_exon_bins_in_bam = n_exon_bins_in_bam,
      n_transcript_exons_in_bam = n_transcript_exons_in_bam
    )
  }

  list(
    junctions = dplyr::bind_rows(junction_list),
    exon_bins = dplyr::bind_rows(exon_bin_cov_list),
    transcript_exons = dplyr::bind_rows(tx_exon_cov_list),
    extraction_audit = dplyr::bind_rows(extraction_audit_list)
  )
}

bam_results <- pmap(
  list(bam_tbl$bam_path, bam_tbl$sample, bam_tbl$sample_core, bam_tbl$group,
       bam_tbl$paper_sample_id, bam_tbl$patient_id, bam_tbl$display_label),
  extract_bam_gene
)

bam_junctions <- dplyr::bind_rows(lapply(bam_results, `[[`, "junctions"))
exon_bin_counts <- dplyr::bind_rows(lapply(bam_results, `[[`, "exon_bins"))
transcript_exon_counts <- dplyr::bind_rows(lapply(bam_results, `[[`, "transcript_exons"))
bam_extraction_audit <- dplyr::bind_rows(lapply(bam_results, `[[`, "extraction_audit"))

exon_gene_totals <- exon_bin_counts %>%
  dplyr::group_by(sample, sample_core, gene) %>%
  dplyr::summarise(
    total_exon_bin_overlaps = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

bam_target_gene_expression_audit <- bam_extraction_audit %>%
  dplyr::left_join(
    exon_gene_totals,
    by = c("sample", "sample_core", "gene")
  ) %>%
  dplyr::mutate(
    total_exon_bin_overlaps = dplyr::if_else(
      chromosome_in_bam & n_exon_bins_in_bam > 0L &
        is.na(total_exon_bin_overlaps),
      0,
      total_exon_bin_overlaps
    ),
    expression_status = dplyr::case_when(
      bam_read_status == "bam_read_error" ~ "not_assessable_bam_read_error",
      !chromosome_in_bam ~ "not_assessable_chromosome_absent",
      n_exon_bins_in_bam == 0L ~ "not_assessable_no_exon_intervals",
      total_exon_bin_overlaps > 0 ~ "detected",
      TRUE ~ "not_detected_zero_exon_coverage"
    )
  )

bam_target_gene_expression_summary <- bam_target_gene_expression_audit %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_bams = dplyr::n(),
    n_assessable_bams = sum(grepl("^(detected|not_detected)", expression_status)),
    n_bams_detected = sum(expression_status == "detected"),
    n_bams_not_detected = sum(
      expression_status == "not_detected_zero_exon_coverage"
    ),
    median_exon_bin_overlaps = if (any(is.finite(total_exon_bin_overlaps))) {
      median(total_exon_bin_overlaps, na.rm = TRUE)
    } else {
      NA_real_
    },
    overall_expression_status = dplyr::case_when(
      n_bams_detected > 0 ~ "detected_in_at_least_one_bam",
      n_assessable_bams > 0 ~ "not_detected_in_any_assessable_bam",
      TRUE ~ "not_assessable"
    ),
    .groups = "drop"
  )

write_tsv(
  bam_target_gene_expression_audit,
  file.path(table_dir, "bam_target_gene_expression_audit.tsv")
)
write_tsv(
  bam_target_gene_expression_summary,
  file.path(table_dir, "bam_target_gene_expression_summary.tsv")
)
msg("\n=== BAM target-gene expression audit ===")
print(bam_target_gene_expression_summary)

if (nrow(bam_junctions) == 0) stop("No splice junctions extracted from BAMs for selected genes.")
if (nrow(exon_bin_counts) == 0) stop("No exon-bin coverage extracted from BAMs for selected genes.")
if (nrow(transcript_exon_counts) == 0) msg("Warning: no transcript-exon coverage rows extracted.")

bam_junctions <- bam_junctions %>%
  dplyr::left_join(
    transcript_introns %>%
      dplyr::group_by(gene, junction_id) %>%
      dplyr::summarise(
        annotated_exact = TRUE,
        annotated_transcripts = collapse_unique(transcript_id),
        annotated_transcript_names = collapse_unique(transcript_name),
        annotated_left_exons = collapse_unique(left_exon_number),
        annotated_right_exons = collapse_unique(right_exon_number),
        .groups = "drop"
      ),
    by = c("gene", "junction_id")
  ) %>%
  dplyr::mutate(
    annotation = ifelse(!is.na(annotated_exact) & annotated_exact, "annotated_exact", "novel_or_nonexact"),
    annotation_simple = ifelse(annotation == "annotated_exact", "annotated", "novel")
  )
write_tsv(bam_junctions, file.path(table_dir, "bam_extracted_junctions.tsv"))

target_info <- coord_tbl %>%
  dplyr::transmute(
    sample_core,
    gene,
    chrom_target = chrom,
    start_target = start,
    end_target = end,
    target_has_interval = has_interval
  ) %>%
  dplyr::distinct(sample_core, gene, .keep_all = TRUE)

selected_genes_df <- genes_df %>% dplyr::filter(gene %in% selected_genes)

analysis_rows <- bam_junctions %>%
  dplyr::left_join(target_info, by = c("sample_core", "gene"), relationship = "many-to-one") %>%
  dplyr::left_join(selected_genes_df %>% dplyr::select(gene, gene_start, gene_end), by = "gene") %>%
  dplyr::mutate(
    chrom_target = ifelse(is.na(chrom_target), chrom, chrom_target),
    target_has_interval = ifelse(is.na(target_has_interval), FALSE, target_has_interval),
    spans_target_interval = ifelse(target_has_interval & chrom == chrom_target, start <= end_target & end >= start_target, FALSE),
    donor_near_target = ifelse(target_has_interval & chrom == chrom_target,
                               abs(start - start_target) <= local_pad_bp | abs(start - end_target) <= local_pad_bp, FALSE),
    acceptor_near_target = ifelse(target_has_interval & chrom == chrom_target,
                                  abs(end - start_target) <= local_pad_bp | abs(end - end_target) <= local_pad_bp, FALSE),
    near_target = spans_target_interval | donor_near_target | acceptor_near_target
  )
write_tsv(analysis_rows, file.path(table_dir, "bam_junction_analysis_rows.tsv"))

all_sample_gene <- analysis_rows %>%
  dplyr::distinct(sample, sample_core, group, paper_sample_id, patient_id, display_label, gene)

all_gene_junction <- analysis_rows %>%
  dplyr::distinct(gene, chrom, start, end, junction_id, annotation_simple, annotation)

junction_complete <- all_sample_gene %>%
  dplyr::inner_join(all_gene_junction, by = "gene") %>%
  dplyr::left_join(
    target_info,
    by = c("sample_core", "gene"),
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    target_has_interval = ifelse(is.na(target_has_interval), FALSE, target_has_interval),
    chrom_target = ifelse(is.na(chrom_target), chrom, chrom_target),
    spans_target_interval = ifelse(
      target_has_interval & chrom == chrom_target,
      start <= end_target & end >= start_target,
      FALSE
    ),
    donor_near_target = ifelse(
      target_has_interval & chrom == chrom_target,
      abs(start - start_target) <= local_pad_bp | abs(start - end_target) <= local_pad_bp,
      FALSE
    ),
    acceptor_near_target = ifelse(
      target_has_interval & chrom == chrom_target,
      abs(end - start_target) <= local_pad_bp | abs(end - end_target) <= local_pad_bp,
      FALSE
    ),
    near_target = spans_target_interval | donor_near_target | acceptor_near_target
  ) %>%
  dplyr::left_join(
    analysis_rows %>% dplyr::select(sample, sample_core, gene, junction_id, count),
    by = c("sample", "sample_core", "gene", "junction_id")
  ) %>%
  dplyr::mutate(count = ifelse(is.na(count), 0, count), detect = count > 0)

j_gene_totals <- junction_complete %>%
  dplyr::group_by(sample, sample_core, group, paper_sample_id, patient_id, display_label, gene) %>%
  dplyr::summarise(gene_total_junction_reads = sum(count, na.rm = TRUE), .groups = "drop")

junction_complete <- junction_complete %>%
  dplyr::left_join(j_gene_totals, by = c("sample", "sample_core", "group", "paper_sample_id", "patient_id", "display_label", "gene")) %>%
  dplyr::mutate(usage = ifelse(gene_total_junction_reads > 0, count / gene_total_junction_reads, 0),
                logcount = log1p(count))
write_tsv(junction_complete, file.path(table_dir, "bam_junction_complete_long.tsv"))

sample_gene_junction_metrics <- junction_complete %>%
  dplyr::group_by(sample, sample_core, group, paper_sample_id, patient_id, display_label, gene) %>%
  dplyr::summarise(
    gene_total_junction_reads = dplyr::first(gene_total_junction_reads),
    n_junctions_detected = sum(count > 0, na.rm = TRUE),
    n_junctions_universe = dplyr::n(),
    annotated_reads = sum(count[annotation_simple == "annotated"], na.rm = TRUE),
    novel_reads = sum(count[annotation_simple == "novel"], na.rm = TRUE),
    novel_fraction = novel_reads / pmax(gene_total_junction_reads, 1),
    entropy = -sum(usage * log(usage + 1e-6), na.rm = TRUE),
    local_reads = sum(count[near_target], na.rm = TRUE),
    local_fraction = local_reads / pmax(gene_total_junction_reads, 1),
    .groups = "drop"
  ) %>%
  dplyr::mutate(detected_fraction = n_junctions_detected / pmax(n_junctions_universe, 1)) %>%
  dplyr::filter(gene_total_junction_reads >= min_gene_total_junction_reads)
write_tsv(sample_gene_junction_metrics, file.path(table_dir, "bam_sample_gene_junction_metrics.tsv"))

msg("\n=== BAM junction sample-gene summary ===")
print(sample_gene_junction_metrics %>% dplyr::summarise(
  n_rows = dplyr::n(),
  n_genes = dplyr::n_distinct(gene),
  n_samples = dplyr::n_distinct(sample),
  median_junction_reads = median(gene_total_junction_reads, na.rm = TRUE),
  median_entropy = median(entropy, na.rm = TRUE),
  median_novel_fraction = median(novel_fraction, na.rm = TRUE)
))

case_tbl <- coord_tbl %>%
  dplyr::filter(gene %in% selected_genes, !is.na(gene), !toupper(group) %in% c("CONTROL")) %>%
  dplyr::distinct(sample_core, sample_raw, group, display_label, gene)

get_reference_samples <- function(g, case_sample_core, sample_gene_tbl, metadata_tbl, mode = reference_mode) {
  observed <- sample_gene_tbl %>% dplyr::filter(gene == g) %>% dplyr::pull(sample_core) %>% unique()
  refs <- base::setdiff(observed, case_sample_core)
  if (mode == "controls") {
    ctrl <- metadata_tbl %>% dplyr::filter(toupper(group) == "CONTROL") %>% dplyr::pull(sample_core) %>% unique()
    refs <- base::intersect(refs, ctrl)
  }
  refs
}

single_case_gene <- list()
for (i in seq_len(nrow(case_tbl))) {
  g <- case_tbl$gene[i]
  s_core <- case_tbl$sample_core[i]
  sg <- sample_gene_junction_metrics %>% dplyr::filter(gene == g)
  if (nrow(sg) == 0) next
  refs <- get_reference_samples(g, s_core, sample_gene_junction_metrics, meta_expanded, reference_mode)
  if (length(refs) < min_reference_samples) next
  one <- sg %>% dplyr::filter(sample_core == s_core)
  ref <- sg %>% dplyr::filter(sample_core %in% refs)
  if (nrow(one) == 0 || nrow(ref) < min_reference_samples) next
  single_case_gene[[length(single_case_gene) + 1]] <- tibble(
    gene = g,
    sample = dplyr::first(one$sample),
    display_label = dplyr::first(one$display_label),
    group = dplyr::first(one$group),
    reference_mode = reference_mode,
    n_reference = nrow(ref),
    gene_total_junction_reads = dplyr::first(one$gene_total_junction_reads),
    entropy = dplyr::first(one$entropy),
    entropy_z = make_z(dplyr::first(one$entropy), ref$entropy),
    entropy_emp_p = empirical_p_2sided(dplyr::first(one$entropy), ref$entropy),
    novel_fraction = dplyr::first(one$novel_fraction),
    novel_fraction_z = make_z(dplyr::first(one$novel_fraction), ref$novel_fraction),
    novel_fraction_emp_p = empirical_p_2sided(dplyr::first(one$novel_fraction), ref$novel_fraction),
    detected_fraction = dplyr::first(one$detected_fraction),
    detected_fraction_z = make_z(dplyr::first(one$detected_fraction), ref$detected_fraction),
    detected_fraction_emp_p = empirical_p_2sided(dplyr::first(one$detected_fraction), ref$detected_fraction)
  )
}
single_case_gene_df <- dplyr::bind_rows(single_case_gene)
if (nrow(single_case_gene_df) > 0) {
  single_case_gene_df <- single_case_gene_df %>%
    dplyr::mutate(
      entropy_fdr = p.adjust(entropy_emp_p, "BH"),
      novel_fraction_fdr = p.adjust(novel_fraction_emp_p, "BH"),
      detected_fraction_fdr = p.adjust(detected_fraction_emp_p, "BH")
    )
}
write_tsv(single_case_gene_df, file.path(table_dir, "bam_single_case_gene_junction_results.tsv"))

per_junction_stats <- list()
for (i in seq_len(nrow(case_tbl))) {
  g <- case_tbl$gene[i]
  s_core <- case_tbl$sample_core[i]
  refs <- get_reference_samples(g, s_core, sample_gene_junction_metrics, meta_expanded, reference_mode)
  if (length(refs) < min_reference_samples) next
  gmat <- junction_complete %>% dplyr::filter(gene == g)
  case_rows <- gmat %>% dplyr::filter(sample_core == s_core)
  ref_rows <- gmat %>% dplyr::filter(sample_core %in% refs)
  if (nrow(case_rows) == 0 || nrow(ref_rows) == 0) next

  ref_summary <- ref_rows %>%
    dplyr::group_by(gene, junction_id, chrom, start, end, annotation_simple, annotation, near_target, spans_target_interval) %>%
    dplyr::summarise(
      ref_n = dplyr::n(),
      ref_detect_n = sum(detect, na.rm = TRUE),
      ref_detect_rate = mean(detect, na.rm = TRUE),
      ref_mean_count = mean(count, na.rm = TRUE),
      ref_median_count = median(count, na.rm = TRUE),
      ref_mean_logcount = mean(logcount, na.rm = TRUE),
      ref_mean_usage = mean(usage, na.rm = TRUE),
      ref_median_usage = median(usage, na.rm = TRUE),
      .groups = "drop"
    )
  ref_vectors <- ref_rows %>%
    dplyr::group_by(gene, junction_id) %>%
    dplyr::summarise(usage_vector = list(usage), logcount_vector = list(logcount), count_vector = list(count), .groups = "drop")
  case_summary <- case_rows %>%
    dplyr::group_by(gene, sample, sample_core, display_label, junction_id, chrom, start, end, annotation_simple, annotation, near_target, spans_target_interval) %>%
    dplyr::summarise(
      case_count = dplyr::first(count),
      case_logcount = dplyr::first(logcount),
      case_usage = dplyr::first(usage),
      case_detect = dplyr::first(detect),
      case_gene_total_junction_reads = dplyr::first(gene_total_junction_reads),
      .groups = "drop"
    )
  joined <- case_summary %>%
    dplyr::left_join(ref_summary, by = c("gene", "junction_id", "chrom", "start", "end", "annotation_simple", "annotation", "near_target", "spans_target_interval")) %>%
    dplyr::left_join(ref_vectors, by = c("gene", "junction_id")) %>%
    dplyr::mutate(
      reference_mode = reference_mode,
      usage_z = purrr::map2_dbl(case_usage, usage_vector, make_z),
      usage_emp_p = purrr::map2_dbl(case_usage, usage_vector, empirical_p_2sided),
      logcount_z = purrr::map2_dbl(case_logcount, logcount_vector, make_z),
      logcount_emp_p = purrr::map2_dbl(case_logcount, logcount_vector, empirical_p_2sided),
      count_z = purrr::map2_dbl(case_count, count_vector, make_z),
      count_emp_p = purrr::map2_dbl(case_count, count_vector, empirical_p_2sided),
      delta_usage = case_usage - ref_median_usage,
      absent_in_case_supported_elsewhere = (!case_detect) & ref_detect_n >= 2 & ref_detect_rate >= 0.25,
      absent_loss_score = ifelse(absent_in_case_supported_elsewhere, ref_detect_rate * ref_mean_usage, 0),
      expressed_only_in_case = case_detect & ref_detect_n == 0,
      gained_novel_score = ifelse(expressed_only_in_case & annotation_simple == "novel", case_usage, 0),
      junction_label = paste0(chrom, ":", start, ":", end)
    )
  per_junction_stats[[length(per_junction_stats) + 1]] <- joined
}
per_junction_df <- dplyr::bind_rows(per_junction_stats)
if (nrow(per_junction_df) > 0) {
  per_junction_df <- per_junction_df %>%
    dplyr::mutate(
      usage_fdr = p.adjust(usage_emp_p, "BH"),
      logcount_fdr = p.adjust(logcount_emp_p, "BH"),
      count_fdr = p.adjust(count_emp_p, "BH")
    )
}
write_tsv(per_junction_df, file.path(table_dir, "bam_per_junction_case_vs_allother_results.tsv"))
write_tsv(per_junction_df %>% dplyr::filter(absent_in_case_supported_elsewhere) %>% dplyr::arrange(dplyr::desc(absent_loss_score)),
          file.path(table_dir, "bam_junctions_absent_in_case_but_present_elsewhere.tsv"))
write_tsv(per_junction_df %>% dplyr::filter(expressed_only_in_case) %>% dplyr::arrange(dplyr::desc(gained_novel_score)),
          file.path(table_dir, "bam_junctions_expressed_only_in_case.tsv"))

exon_bin_counts <- exon_bin_counts %>%
  dplyr::group_by(sample, sample_core, group, paper_sample_id, patient_id, display_label, gene) %>%
  dplyr::mutate(gene_total_exon_bin_reads = sum(count, na.rm = TRUE), exon_bin_fraction = count / pmax(gene_total_exon_bin_reads, 1)) %>%
  dplyr::ungroup()
write_tsv(exon_bin_counts, file.path(table_dir, "bam_exon_bin_coverage_long.tsv"))

case_lookup <- case_tbl %>% dplyr::select(gene, sample_core) %>% dplyr::mutate(is_case = TRUE)
exon_bin_counts <- exon_bin_counts %>%
  dplyr::left_join(case_lookup, by = c("gene", "sample_core")) %>%
  dplyr::mutate(is_case = ifelse(is.na(is_case), FALSE, is_case), Role = ifelse(is_case, "Case", "Reference"))

exon_bin_summary <- exon_bin_counts %>%
  dplyr::filter(gene_total_exon_bin_reads >= min_gene_total_exon_reads) %>%
  dplyr::group_by(gene, chrom, exon_bin_rank, exon_bin_start, exon_bin_end, exon_bin_label) %>%
  dplyr::summarise(
    case_count = dplyr::first(count[Role == "Case"]),
    case_fraction = dplyr::first(exon_bin_fraction[Role == "Case"]),
    case_total_exon_bin_reads = dplyr::first(gene_total_exon_bin_reads[Role == "Case"]),
    ref_mean_count = mean(count[Role == "Reference"], na.rm = TRUE),
    ref_median_count = median(count[Role == "Reference"], na.rm = TRUE),
    ref_mean_fraction = mean(exon_bin_fraction[Role == "Reference"], na.rm = TRUE),
    ref_median_fraction = median(exon_bin_fraction[Role == "Reference"], na.rm = TRUE),
    ref_sd_fraction = stats::sd(exon_bin_fraction[Role == "Reference"], na.rm = TRUE),
    n_reference = sum(Role == "Reference"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    fraction_delta = case_fraction - ref_median_fraction,
    fraction_z = ifelse(is.finite(ref_sd_fraction) & ref_sd_fraction > 0, (case_fraction - ref_mean_fraction) / ref_sd_fraction, NA_real_),
    usage_call = dplyr::case_when(
      is.finite(fraction_z) & fraction_z <= -1.5 & fraction_delta <= -0.01 ~ "underused_in_case",
      is.finite(fraction_z) & fraction_z >=  1.5 & fraction_delta >=  0.01 ~ "overused_in_case",
      TRUE ~ "no_strong_shift"
    )
  ) %>%
  dplyr::arrange(gene, exon_bin_rank)
write_tsv(exon_bin_summary, file.path(table_dir, "bam_exon_bin_usage_summary.tsv"))

annotate_one_exon_bin <- function(gene_value, chrom_value, bin_start, bin_end) {
  exon_keep <-
    exon_model$gene == gene_value &
    exon_model$chrom == chrom_value &
    exon_model$exon_start <= bin_end &
    exon_model$exon_end >= bin_start
  exon_keep[is.na(exon_keep)] <- FALSE
  exon_hits <- exon_model[exon_keep, , drop = FALSE]

  cds_keep <-
    cds_model$gene == gene_value &
    cds_model$chrom == chrom_value &
    cds_model$cds_start <= bin_end &
    cds_model$cds_end >= bin_start
  cds_keep[is.na(cds_keep)] <- FALSE
  cds_hits <- cds_model[cds_keep, , drop = FALSE]

  cds_bases <- if (nrow(cds_hits) == 0) {
    0L
  } else {
    overlap_bases <- pmax(
      0,
      pmin(cds_hits$cds_end, bin_end) -
        pmax(cds_hits$cds_start, bin_start) + 1
    )
    as.integer(max(overlap_bases, na.rm = TRUE))
  }

  tibble::tibble(
    overlapping_transcripts = collapse_unique(exon_hits$transcript_id),
    overlapping_transcript_names = collapse_unique(exon_hits$transcript_name),
    overlapping_exon_numbers = collapse_unique(
      paste0(exon_hits$transcript_id, ":exon", exon_hits$exon_number)
    ),
    cds_bases_in_bin = cds_bases
  )
}

bin_transcript_annotations <- if (nrow(exon_bin_summary) == 0) {
  tibble::tibble(
    overlapping_transcripts = character(),
    overlapping_transcript_names = character(),
    overlapping_exon_numbers = character(),
    cds_bases_in_bin = integer()
  )
} else {
  purrr::pmap_dfr(
    list(
      exon_bin_summary$gene,
      exon_bin_summary$chrom,
      exon_bin_summary$exon_bin_start,
      exon_bin_summary$exon_bin_end
    ),
    annotate_one_exon_bin
  )
}

map_bins_to_transcripts <- dplyr::bind_cols(
  exon_bin_summary,
  bin_transcript_annotations
) %>%
  dplyr::mutate(
    cds_overlap_class = ifelse(
      cds_bases_in_bin > 0,
      "CDS_overlapping_bin",
      "UTR_or_nonCDS_bin"
    )
  )
write_tsv(map_bins_to_transcripts, file.path(table_dir, "bam_exon_bin_usage_summary_transcript_mapped.tsv"))

if (nrow(transcript_exon_counts) > 0) {
  transcript_exon_counts <- transcript_exon_counts %>%
    dplyr::group_by(sample, sample_core, group, paper_sample_id, patient_id, display_label, gene, transcript_id) %>%
    dplyr::mutate(transcript_total_exon_reads = sum(count, na.rm = TRUE), transcript_exon_fraction = count / pmax(transcript_total_exon_reads, 1)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(case_lookup, by = c("gene", "sample_core")) %>%
    dplyr::mutate(is_case = ifelse(is.na(is_case), FALSE, is_case), Role = ifelse(is_case, "Case", "Reference"))
  write_tsv(transcript_exon_counts, file.path(table_dir, "bam_transcript_exon_coverage_long.tsv"))

  transcript_exon_summary <- transcript_exon_counts %>%
    dplyr::filter(transcript_total_exon_reads >= min_gene_total_exon_reads) %>%
    dplyr::group_by(gene, transcript_id, transcript_name, transcript_biotype, chrom, strand, exon_number, exon_rank_genomic, exon_start, exon_end) %>%
    dplyr::summarise(
      case_count = dplyr::first(count[Role == "Case"]),
      case_fraction = dplyr::first(transcript_exon_fraction[Role == "Case"]),
      case_total_transcript_exon_reads = dplyr::first(transcript_total_exon_reads[Role == "Case"]),
      ref_mean_count = mean(count[Role == "Reference"], na.rm = TRUE),
      ref_median_count = median(count[Role == "Reference"], na.rm = TRUE),
      ref_mean_fraction = mean(transcript_exon_fraction[Role == "Reference"], na.rm = TRUE),
      ref_median_fraction = median(transcript_exon_fraction[Role == "Reference"], na.rm = TRUE),
      ref_sd_fraction = stats::sd(transcript_exon_fraction[Role == "Reference"], na.rm = TRUE),
      n_reference = sum(Role == "Reference"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      fraction_delta = case_fraction - ref_median_fraction,
      fraction_z = ifelse(is.finite(ref_sd_fraction) & ref_sd_fraction > 0, (case_fraction - ref_mean_fraction) / ref_sd_fraction, NA_real_),
      usage_call = dplyr::case_when(
        is.finite(fraction_z) & fraction_z <= -1.5 & fraction_delta <= -0.01 ~ "underused_in_case",
        is.finite(fraction_z) & fraction_z >=  1.5 & fraction_delta >=  0.01 ~ "overused_in_case",
        TRUE ~ "no_strong_shift"
      )
    ) %>%
    dplyr::mutate(
      cds_overlap_class = ifelse(
        purrr::pmap_lgl(
          list(transcript_id, exon_start, exon_end),
          function(tx_value, exon_start_value, exon_end_value) {
            hit <-
              cds_model$transcript_id == tx_value &
              cds_model$cds_start <= exon_end_value &
              cds_model$cds_end >= exon_start_value
            any(hit %in% TRUE)
          }
        ),
        "CDS_overlapping_exon",
        "UTR_or_nonCDS_exon"
      )
    ) %>%
    dplyr::arrange(gene, transcript_id, exon_rank_genomic)
  write_tsv(transcript_exon_summary, file.path(table_dir, "bam_transcript_exon_usage_summary.tsv"))
}

classify_junction_to_transcript <- function(j_start, j_end, exon_tbl, tol = 2) {

  left_hits <- exon_tbl %>% dplyr::filter(abs(exon_end - j_start) <= tol)
  right_hits <- exon_tbl %>% dplyr::filter(abs(exon_start - (j_end + 1)) <= tol)

  if (nrow(left_hits) == 0 && nrow(right_hits) == 0) {
    return(tibble(event_type = "unmapped_to_transcript_exons", left_exon_number = NA_character_, right_exon_number = NA_character_,
                  left_exon_rank_genomic = NA_integer_, right_exon_rank_genomic = NA_integer_, skipped_exon_numbers = NA_character_,
                  skipped_exon_rank_genomic = NA_character_, n_skipped_exons = 0L))
  }
  if (nrow(left_hits) > 0 && nrow(right_hits) > 0) {
    out <- list()
    for (i in seq_len(nrow(left_hits))) {
      for (j in seq_len(nrow(right_hits))) {
        l <- left_hits[i, ]; r <- right_hits[j, ]
        if (r$exon_rank_genomic == l$exon_rank_genomic + 1) {
          event_type <- "canonical_adjacent_junction"; skipped <- exon_tbl[0, ]
        } else if (r$exon_rank_genomic > l$exon_rank_genomic + 1) {
          event_type <- "exon_skipping_candidate"
          skipped <- exon_tbl %>% dplyr::filter(exon_rank_genomic > l$exon_rank_genomic, exon_rank_genomic < r$exon_rank_genomic)
        } else {
          event_type <- "noncanonical_reverse_or_overlapping_boundary"; skipped <- exon_tbl[0, ]
        }
        out[[length(out) + 1]] <- tibble(
          event_type = event_type,
          left_exon_number = as.character(l$exon_number),
          right_exon_number = as.character(r$exon_number),
          left_exon_rank_genomic = as.integer(l$exon_rank_genomic),
          right_exon_rank_genomic = as.integer(r$exon_rank_genomic),
          skipped_exon_numbers = collapse_unique(skipped$exon_number),
          skipped_exon_rank_genomic = collapse_unique(skipped$exon_rank_genomic),
          n_skipped_exons = nrow(skipped)
        )
      }
    }
    return(dplyr::bind_rows(out))
  }
  if (nrow(left_hits) > 0 && nrow(right_hits) == 0) {
    return(tibble(event_type = "alternative_acceptor_candidate", left_exon_number = collapse_unique(left_hits$exon_number), right_exon_number = NA_character_,
                  left_exon_rank_genomic = suppressWarnings(min(left_hits$exon_rank_genomic, na.rm = TRUE)), right_exon_rank_genomic = NA_integer_,
                  skipped_exon_numbers = NA_character_, skipped_exon_rank_genomic = NA_character_, n_skipped_exons = 0L))
  }
  tibble(event_type = "alternative_donor_candidate", left_exon_number = NA_character_, right_exon_number = collapse_unique(right_hits$exon_number),
         left_exon_rank_genomic = NA_integer_, right_exon_rank_genomic = suppressWarnings(min(right_hits$exon_rank_genomic, na.rm = TRUE)),
         skipped_exon_numbers = NA_character_, skipped_exon_rank_genomic = NA_character_, n_skipped_exons = 0L)
}

calc_cds_impact <- function(tx_id, event_row, exon_tbl, cds_tbl) {
  tx_cds <- cds_tbl %>% dplyr::filter(transcript_id == tx_id)
  if (nrow(tx_cds) == 0) {
    return(tibble(affected_exon_numbers = NA_character_, cds_bases_affected = NA_integer_, cds_bases_affected_mod3 = NA_integer_,
                  cds_overlap_class = "no_CDS_annotation_for_transcript", predicted_functional_impact = "noncoding_or_unannotated_CDS"))
  }
  event_type <- event_row$event_type[1]
  affected_exons <- exon_tbl[0, ]
  if (event_type == "exon_skipping_candidate" && !is.na(event_row$skipped_exon_rank_genomic[1])) {
    ranks <- as.integer(strsplit(event_row$skipped_exon_rank_genomic[1], ";")[[1]])
    affected_exons <- exon_tbl %>% dplyr::filter(exon_rank_genomic %in% ranks)
  } else if (event_type %in% c("canonical_adjacent_junction", "alternative_acceptor_candidate", "alternative_donor_candidate")) {
    ranks <- c(event_row$left_exon_rank_genomic[1], event_row$right_exon_rank_genomic[1])
    ranks <- ranks[is.finite(ranks)]
    affected_exons <- exon_tbl %>% dplyr::filter(exon_rank_genomic %in% ranks)
  }
  if (nrow(affected_exons) == 0) {
    return(tibble(affected_exon_numbers = NA_character_, cds_bases_affected = 0L, cds_bases_affected_mod3 = NA_integer_,
                  cds_overlap_class = "no_mapped_affected_exon", predicted_functional_impact = "unmapped_event_functional_impact_uncertain"))
  }
  cds_bases <- 0L
  for (i in seq_len(nrow(affected_exons))) {
    ex_start <- affected_exons$exon_start[i]; ex_end <- affected_exons$exon_end[i]
    ov <- tx_cds %>% dplyr::mutate(ov_start = pmax(cds_start, ex_start), ov_end = pmin(cds_end, ex_end), ov_width = pmax(0, ov_end - ov_start + 1))
    cds_bases <- cds_bases + sum(ov$ov_width, na.rm = TRUE)
  }
  cds_mod <- ifelse(cds_bases > 0, cds_bases %% 3L, NA_integer_)
  impact <- dplyr::case_when(
    event_type == "exon_skipping_candidate" & cds_bases > 0 & cds_mod == 0 ~ "predicted_in_frame_CDS_skip",
    event_type == "exon_skipping_candidate" & cds_bases > 0 & cds_mod != 0 ~ "predicted_frameshift_CDS_skip_possible_NMD",
    event_type == "exon_skipping_candidate" & cds_bases == 0 ~ "predicted_UTR_or_nonCDS_skip",
    event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate", "canonical_adjacent_junction") & cds_bases > 0 ~ "splice_boundary_change_CDS_impact_uncertain",
    event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate", "canonical_adjacent_junction") & cds_bases == 0 ~ "splice_boundary_change_UTR_or_nonCDS",
    TRUE ~ "functional_impact_uncertain"
  )
  tibble(affected_exon_numbers = collapse_unique(affected_exons$exon_number),
         cds_bases_affected = as.integer(cds_bases),
         cds_bases_affected_mod3 = as.integer(cds_mod),
         cds_overlap_class = ifelse(cds_bases > 0, "CDS_overlapping_event", "UTR_or_nonCDS_event"),
         predicted_functional_impact = impact)
}

altered_junctions <- per_junction_df %>%
  dplyr::filter(gene %in% selected_genes) %>%
  dplyr::filter(absent_in_case_supported_elsewhere | expressed_only_in_case |
                  abs(usage_z) >= altered_abs_usage_z | abs(logcount_z) >= altered_abs_logcount_z) %>%
  dplyr::mutate(alteration_class = dplyr::case_when(
    absent_in_case_supported_elsewhere ~ "lost_in_case_supported_elsewhere",
    expressed_only_in_case ~ "gained_case_only",
    usage_z <= -altered_abs_usage_z | logcount_z <= -altered_abs_logcount_z ~ "reduced_in_case",
    usage_z >= altered_abs_usage_z | logcount_z >= altered_abs_logcount_z ~ "increased_in_case",
    TRUE ~ "altered"
  ))
write_tsv(altered_junctions, file.path(table_dir, "bam_altered_junctions_for_transcript_annotation.tsv"))

annotation_rows <- list()
for (i in seq_len(nrow(altered_junctions))) {
  jr <- altered_junctions[i, ]
  tx_ids <- exon_model %>% dplyr::filter(gene == jr$gene) %>% dplyr::pull(transcript_id) %>% unique()
  for (tx in tx_ids) {
    ex_tx <- exon_model %>% dplyr::filter(gene == jr$gene, transcript_id == tx) %>% dplyr::arrange(exon_rank_genomic)
    class_tbl <- classify_junction_to_transcript(jr$start, jr$end, ex_tx, boundary_tol_bp)
    for (k in seq_len(nrow(class_tbl))) {
      ct <- class_tbl[k, ]
      impact <- calc_cds_impact(tx, ct, ex_tx, cds_model)
      annotation_rows[[length(annotation_rows) + 1]] <- dplyr::bind_cols(
        jr,
        tibble(transcript_id = tx,
               transcript_name = dplyr::first(ex_tx$transcript_name),
               transcript_biotype = dplyr::first(ex_tx$transcript_biotype),
               strand = dplyr::first(ex_tx$strand),
               n_exons_transcript = nrow(ex_tx)),
        ct,
        impact
      )
    }
  }
}
transcript_aware_junction_annotation <- dplyr::bind_rows(annotation_rows)
if (nrow(transcript_aware_junction_annotation) > 0) {
  transcript_aware_junction_annotation <- transcript_aware_junction_annotation %>%
    dplyr::mutate(
      transcript_event_score = dplyr::case_when(
        event_type == "exon_skipping_candidate" ~ 100,
        event_type %in% c("alternative_acceptor_candidate", "alternative_donor_candidate") ~ 50,
        event_type == "canonical_adjacent_junction" ~ 25,
        TRUE ~ 0
      ) + dplyr::coalesce(abs(usage_z), 0) + dplyr::coalesce(abs(logcount_z), 0) +
        dplyr::coalesce(absent_loss_score, 0) + dplyr::coalesce(gained_novel_score, 0)
    ) %>%
    dplyr::arrange(gene, dplyr::desc(transcript_event_score))
}
write_tsv(transcript_aware_junction_annotation, file.path(table_dir, "bam_transcript_aware_junction_annotation.tsv"))

top_transcript_events <- transcript_aware_junction_annotation %>%
  dplyr::filter(event_type != "unmapped_to_transcript_exons") %>%
  dplyr::group_by(gene, junction_id, alteration_class) %>%
  dplyr::slice_max(order_by = transcript_event_score, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(gene, dplyr::desc(transcript_event_score))
write_tsv(top_transcript_events, file.path(table_dir, "bam_top_transcript_aware_events.tsv"))

transcript_event_summary <- transcript_aware_junction_annotation %>%
  dplyr::group_by(gene, transcript_id, transcript_name, transcript_biotype, event_type, predicted_functional_impact) %>%
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
  dplyr::mutate(dplyr::across(c(max_abs_usage_z, max_abs_logcount_z, max_absent_loss_score, max_gained_novel_score, max_cds_bases_affected),
                              ~ ifelse(is.infinite(.x), NA_real_, .x))) %>%
  dplyr::arrange(gene, dplyr::desc(n_lost + n_gained), dplyr::desc(max_absent_loss_score), dplyr::desc(max_gained_novel_score))
write_tsv(transcript_event_summary, file.path(table_dir, "bam_transcript_event_summary.tsv"))

msg("Generating plots...")

case_lookup_entropy <- case_tbl %>% dplyr::select(gene, sample_core) %>% dplyr::distinct() %>% dplyr::mutate(is_case_gene = TRUE)
entropy_plot_df <- sample_gene_junction_metrics %>%
  dplyr::left_join(case_lookup_entropy, by = c("gene", "sample_core")) %>%
  dplyr::mutate(is_case_gene = ifelse(is.na(is_case_gene), FALSE, is_case_gene),
                Role = ifelse(is_case_gene, "Inserted case", "All other samples")) %>%
  dplyr::filter(is.finite(entropy))

gene_order_entropy <- entropy_plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(case_entropy = dplyr::first(entropy[Role == "Inserted case"]),
                   ref_median_entropy = median(entropy[Role == "All other samples"], na.rm = TRUE),
                   entropy_delta = case_entropy - ref_median_entropy,
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(abs(entropy_delta))) %>% dplyr::pull(gene)
entropy_plot_df$gene <- factor(entropy_plot_df$gene, levels = gene_order_entropy)

p_entropy <- ggplot(entropy_plot_df, aes(x = gene, y = entropy)) +
  geom_point(data = entropy_plot_df %>% dplyr::filter(Role == "All other samples"),
             color = "black", size = 2, position = position_jitter(width = 0.15, height = 0)) +
  geom_point(data = entropy_plot_df %>% dplyr::filter(Role == "Inserted case"),
             color = "red", size = 3) +
  coord_flip() +
  labs(title = "BAM-derived splice-junction entropy by target gene",
       subtitle = "Black = all other samples; red = vector-inserted case sample",
       x = "Target gene", y = "Junction usage entropy")
save_plot_all(p_entropy, "bam_entropy_by_gene_case_vs_allother", width = 7.5, height = 4.8)

entropy_summary <- entropy_plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(case_sample = dplyr::first(sample[Role == "Inserted case"]),
                   case_label = dplyr::first(display_label[Role == "Inserted case"]),
                   case_entropy = dplyr::first(entropy[Role == "Inserted case"]),
                   ref_median_entropy = median(entropy[Role == "All other samples"], na.rm = TRUE),
                   entropy_delta = case_entropy - ref_median_entropy,
                   n_reference = sum(Role == "All other samples"),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(abs(entropy_delta)))
write_tsv(entropy_summary, file.path(table_dir, "bam_entropy_by_gene_case_vs_allother_summary.tsv"))

for (g in selected_genes) {
  ug <- junction_complete %>% dplyr::filter(gene == g)
  sg <- sample_gene_junction_metrics %>% dplyr::filter(gene == g)
  if (nrow(ug) == 0 || nrow(sg) == 0) next
  case_samples_g <- case_tbl %>% dplyr::filter(gene == g) %>% dplyr::pull(sample_raw) %>% unique()
  case_cores_g <- case_tbl %>% dplyr::filter(gene == g) %>% dplyr::pull(sample_core) %>% unique()
  refs_g <- sg %>% dplyr::filter(!sample_core %in% case_cores_g) %>% dplyr::pull(sample) %>% unique()

  stats_g <- per_junction_df %>% dplyr::filter(gene == g) %>% dplyr::mutate(score = dplyr::coalesce(abs(usage_z), 0) + absent_loss_score + gained_novel_score)
  top_j <- stats_g %>% dplyr::arrange(dplyr::desc(score), dplyr::desc(ref_mean_usage), dplyr::desc(case_usage)) %>%
    dplyr::slice_head(n = label_top_junc_n) %>% dplyr::pull(junction_id) %>% unique()
  if (length(top_j) == 0) next

  plot_df <- ug %>%
    dplyr::filter(junction_id %in% top_j) %>%
    dplyr::mutate(Role = ifelse(sample_core %in% case_cores_g, "Case", "Reference"),
                  JunctionLabel = ifelse(near_target | spans_target_interval, paste0(junction_id, " *"), junction_id))

  p_j <- ggplot(plot_df, aes(x = JunctionLabel, y = usage)) +
    geom_point(data = plot_df %>% dplyr::filter(Role == "Reference"), color = "black", size = 2,
               position = position_jitter(width = 0.12, height = 0)) +
    geom_point(data = plot_df %>% dplyr::filter(Role == "Case"), color = "red", size = 3) +
    coord_flip() +
    labs(title = paste0(g, ": BAM-derived altered junction usage"),
         subtitle = "Black = all other samples; red = case; * near/spans target interval",
         x = "Junction", y = "Within-gene junction usage")
  save_plot_all(p_j, paste0("bam_junction_usage_case_vs_allother_", g), width = 9.5, height = 6.2)

  hm_df <- ug %>%
    dplyr::filter(junction_id %in% top_j) %>%
    dplyr::mutate(HeatmapSample = ifelse(sample_core %in% case_cores_g, paste0("[CASE] ", display_label), display_label)) %>%
    dplyr::select(junction_id, HeatmapSample, usage) %>%
    dplyr::distinct(junction_id, HeatmapSample, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = HeatmapSample, values_from = usage, values_fill = 0)
  if (nrow(hm_df) >= 2 && ncol(hm_df) >= 3) {
    hm_mat <- as.data.frame(hm_df); rownames(hm_mat) <- hm_mat$junction_id; hm_mat$junction_id <- NULL; hm_mat <- as.matrix(hm_mat)
    ann_col <- data.frame(Role = ifelse(grepl("^\\[CASE\\]", colnames(hm_mat)), "Case", "Reference"), row.names = colnames(hm_mat))
    pheat <- pheatmap::pheatmap(hm_mat, scale = "row", annotation_col = ann_col,
                                main = paste0(g, ": BAM-derived top junction usage"),
                                fontsize = heatmap_text_size,
                                fontsize_row = heatmap_text_size,
                                fontsize_col = heatmap_text_size,
                                angle_col = 45,
                                silent = TRUE)
    save_pheatmap_all(pheat, paste0("bam_junction_heatmap_", g), width = 9.0, height = 7.2)
  }
}

for (g in selected_genes) {
  plot_df <- exon_bin_counts %>% dplyr::filter(gene == g)
  sum_df <- exon_bin_summary %>% dplyr::filter(gene == g)
  if (nrow(plot_df) == 0 || nrow(sum_df) == 0) next

  p_bin <- ggplot(plot_df, aes(x = factor(exon_bin_rank), y = exon_bin_fraction)) +
    geom_point(data = plot_df %>% dplyr::filter(Role == "Reference"), color = "black", size = 1.7,
               position = position_jitter(width = 0.12, height = 0)) +
    geom_point(data = plot_df %>% dplyr::filter(Role == "Case"), color = "red", size = 3) +
    labs(title = paste0(g, ": BAM exon-bin usage fractions"),
         subtitle = "Black = all other samples; red = inserted case sample",
         x = "Exon-bin rank", y = "Fraction of gene exonic coverage")
  save_plot_all(p_bin, paste0("bam_exon_bin_usage_fraction_", g), width = 9.0, height = 5.6)

  p_delta <- ggplot(sum_df, aes(x = factor(exon_bin_rank), y = fraction_delta, fill = usage_call)) +
    geom_col() +
    labs(title = paste0(g, ": BAM exon-bin usage delta"),
         subtitle = "Case fraction minus reference median fraction",
         x = "Exon-bin rank", y = "Coverage fraction delta", fill = NULL)
  save_plot_all(p_delta, paste0("bam_exon_bin_usage_delta_", g), width = 9.2, height = 5.6)
}

if (exists("transcript_exon_summary") && nrow(transcript_exon_summary) > 0) {
  top_tx <- transcript_summary %>%
    dplyr::filter(gene %in% selected_genes) %>%
    dplyr::group_by(gene) %>%
    dplyr::arrange(dplyr::desc(is_protein_coding), dplyr::desc(cds_bases), dplyr::desc(n_exons), .by_group = TRUE) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup()
  for (i in seq_len(nrow(top_tx))) {
    tx <- top_tx$transcript_id[i]
    g <- top_tx$gene[i]
    plot_df <- transcript_exon_counts %>% dplyr::filter(gene == g, transcript_id == tx)
    sum_df <- transcript_exon_summary %>% dplyr::filter(gene == g, transcript_id == tx)
    if (nrow(plot_df) == 0 || nrow(sum_df) == 0) next
    p_tx <- ggplot(plot_df, aes(x = factor(exon_rank_genomic), y = transcript_exon_fraction)) +
      geom_point(data = plot_df %>% dplyr::filter(Role == "Reference"), color = "black", size = 1.7,
                 position = position_jitter(width = 0.12, height = 0)) +
      geom_point(data = plot_df %>% dplyr::filter(Role == "Case"), color = "red", size = 3) +
      labs(title = paste0(g, " ", tx, ": BAM transcript-exon usage"),
           subtitle = "Black = all other samples; red = inserted case sample",
           x = "Transcript exon rank", y = "Fraction of transcript exon coverage")
    save_plot_all(p_tx, paste0("bam_transcript_exon_usage_", g, "_", gsub("[^A-Za-z0-9]+", "_", tx)), width = 9.2, height = 5.6)
  }
}

msg("\n=== Top BAM single-case gene junction results ===")
if (nrow(single_case_gene_df) > 0) print(single_case_gene_df)

msg("\n=== Top BAM lost junctions ===")
if (nrow(per_junction_df) > 0) {
  print(per_junction_df %>% dplyr::filter(absent_in_case_supported_elsewhere) %>%
          dplyr::arrange(dplyr::desc(absent_loss_score)) %>%
          dplyr::select(gene, display_label, junction_label, annotation_simple, case_count, case_usage, ref_detect_n, ref_detect_rate, ref_mean_usage, usage_z, logcount_z, absent_loss_score) %>%
          dplyr::slice_head(n = 20))
}

msg("\n=== Top BAM gained junctions ===")
if (nrow(per_junction_df) > 0) {
  print(per_junction_df %>% dplyr::filter(expressed_only_in_case) %>%
          dplyr::arrange(dplyr::desc(gained_novel_score)) %>%
          dplyr::select(gene, display_label, junction_label, annotation_simple, case_count, case_usage, ref_detect_n, ref_detect_rate, usage_z, logcount_z, gained_novel_score) %>%
          dplyr::slice_head(n = 20))
}

msg("\n=== Top transcript-aware BAM events ===")
if (exists("top_transcript_events") && nrow(top_transcript_events) > 0) {
  print(top_transcript_events %>%
          dplyr::select(gene, display_label, alteration_class, junction_label, annotation_simple, transcript_id, transcript_name,
                        event_type, left_exon_number, right_exon_number, skipped_exon_numbers,
                        affected_exon_numbers, cds_bases_affected, cds_bases_affected_mod3,
                        predicted_functional_impact, transcript_event_score) %>%
          dplyr::slice_head(n = 30))
}

msg("\n=== Top BAM exon-bin usage shifts ===")
print(map_bins_to_transcripts %>%
        dplyr::arrange(gene, dplyr::desc(abs(fraction_z))) %>%
        dplyr::select(gene, exon_bin_rank, exon_bin_start, exon_bin_end, case_fraction, ref_median_fraction, fraction_delta, fraction_z, usage_call, cds_overlap_class, overlapping_transcript_names, overlapping_exon_numbers) %>%
        dplyr::slice_head(n = 40))

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
msg("Key tables:")
msg("  - ", file.path(table_dir, "bam_per_junction_case_vs_allother_results.tsv"))
msg("  - ", file.path(table_dir, "bam_junctions_absent_in_case_but_present_elsewhere.tsv"))
msg("  - ", file.path(table_dir, "bam_junctions_expressed_only_in_case.tsv"))
msg("  - ", file.path(table_dir, "bam_transcript_aware_junction_annotation.tsv"))
msg("  - ", file.path(table_dir, "bam_top_transcript_aware_events.tsv"))
msg("  - ", file.path(table_dir, "bam_exon_bin_usage_summary_transcript_mapped.tsv"))
msg("  - ", file.path(table_dir, "bam_transcript_exon_usage_summary.tsv"))
msg("Key plots: see ", normalizePath(plot_dir))
