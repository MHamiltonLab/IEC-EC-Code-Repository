#!/usr/bin/env Rscript
# TCR repertoire diversity
# Collapse MiXCR TRB clones by CDR3 amino-acid sequence and compare clonality and diversity.

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(RColorBrewer)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
mixcr_out <- file.path(data_dir, "mixcr")
chain     <- "TRB"

out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "immune_repertoire", "tcr_diversity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
collapsed_dir <- file.path(out_dir, "collapsed_tables")
dir.create(collapsed_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

BASE_SIZE   <- 20
base_family <- "sans"
use_showtext <- FALSE

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE) &&
    requireNamespace("systemfonts", quietly = TRUE)) {

  ar <- systemfonts::match_font("Arial")
  if (is.list(ar) && !is.null(ar$path) && nzchar(ar$path)) {
    sysfonts::font_add("Arial", regular = ar$path)
    showtext::showtext_auto(enable = TRUE)
    base_family <- "Arial"
    use_showtext <- TRUE
  } else {
    message("[WARN] Arial not found on this node. Using 'sans' (will not crash).")
  }
} else {
  message("[WARN] showtext/sysfonts/systemfonts not available. Using 'sans' (will not crash).")
}

theme_set(
  theme_classic(base_size = BASE_SIZE, base_family = base_family) +
    theme(
      plot.title    = element_text(size = BASE_SIZE + 4, face = "bold"),
      plot.subtitle = element_text(size = BASE_SIZE),
      axis.title    = element_text(size = BASE_SIZE),
      axis.text     = element_text(size = BASE_SIZE - 2),
      legend.title  = element_text(size = BASE_SIZE),
      legend.text   = element_text(size = BASE_SIZE - 3),
      strip.text    = element_text(size = BASE_SIZE, face = "bold"),
      plot.margin   = margin(10, 10, 10, 10)
    )
)

save_both <- function(p, stem, width = 12, height = 8) {
  pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
  eps_file <- file.path(fig_dir, paste0(stem, ".eps"))

  ggplot2::ggsave(
    filename = pdf_file, plot = p, width = width, height = height,
    device = grDevices::cairo_pdf, dpi = 300
  )

  ggplot2::ggsave(
    filename = eps_file, plot = p, width = width, height = height,
    device = grDevices::cairo_ps, dpi = 300
  )

  invisible(list(pdf = pdf_file, eps = eps_file))
}

reverse_sample_factor <- function(p) {
  p + scale_y_discrete(limits = rev)
}

col_other_salmon <- "#E6A3A3"
col_top_aqua     <- "#73C6C6"

pal_iec <- c(
  "No camelid infiltrate"     = "#4C6A87",
  "Camelid infiltrate present" = "#B22222"
)

pal_car <- c("CAR−" = "#4C6A87", "CAR+" = "#B22222")

muted_palette <- function(n) {
  base <- brewer.pal(8, "Set2")
  colorRampPalette(base)(n)
}

tissue_normalize <- function(x) {
  x <- str_squish(x)
  str_replace_all(x, c(
    "^Terminal-ileum$" = "Terminal Ileum",
    "^Terminal ileum$" = "Terminal Ileum",
    "^TI$"             = "Terminal Ileum",
    "^Co$"             = "Colon",
    "^Du$"             = "Duodenum"
  ))
}
tissue_order <- c("Marrow", "Duodenum", "Terminal Ileum", "Colon", "LN")

# Repertoire metadata
meta <- readr::read_tsv(file.path(data_dir, "repertoire_metadata.tsv"), show_col_types = FALSE) %>%
  transmute(
    Sample    = Sample_ID,
    Paper_ID  = Paper_ID,
    Day       = as.numeric(Day),
    Tissue    = tissue_normalize(Type),
    Timepoint = Timepoint,
    Pathology = str_replace(Pathology, "\\.$", ""),
    Group     = case_when(
      Tissue == "Marrow" ~ "Marrow",
      Tissue == "LN"     ~ "LN",
      TRUE               ~ "GI"
    )
  ) %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Label  = paste0(Paper_ID, " | ", as.character(Tissue)),
    CamelidGroup = case_when(
      Pathology == "IEC"     ~ "Camelid infiltrate present",
      Pathology == "Control" ~ "No camelid infiltrate",
      TRUE                   ~ Pathology
    )
  )

time_label <- function(sample_vec) {
  meta %>%
    filter(Sample %in% sample_vec) %>%
    arrange(Day) %>%
    transmute(Sample, TimeLabel = paste0("Day ", Day, " | ", as.character(Tissue))) %>%
    deframe()
}

pattern <- paste0("\\.clones_", chain, "\\.tsv$")
files_all <- list.files(mixcr_out, pattern = pattern, recursive = TRUE, full.names = TRUE)
if (length(files_all) == 0) stop("No MiXCR clones files found matching: ", pattern)

sample_from_path <- function(p) str_replace(basename(p), pattern, "")
mixcr_index <- tibble(
  Sample = vapply(files_all, sample_from_path, character(1)),
  file   = files_all
)

missing_files <- setdiff(meta$Sample, mixcr_index$Sample)
if (length(missing_files) > 0) message("[WARN] Metadata samples without MiXCR files: ", paste(missing_files, collapse = ", "))

mixcr_index <- mixcr_index %>% filter(Sample %in% meta$Sample)
meta <- meta %>% filter(Sample %in% mixcr_index$Sample)

guess_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}
clean_gene <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  first <- str_split(x, ",", simplify = TRUE)[1]
  str_replace(first, "\\(.*\\)$", "")
}

# Clone-table import and CDR3 collapse
read_and_collapse_mixcr <- function(file, sample_name) {
  df_raw <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE))

  aa_col   <- guess_col(df_raw, c("CDR3.aa", "aa", "aaSeqCDR3", "cdr3aa"))
  cnt_col  <- guess_col(df_raw, c("Clones", "readCount", "cloneCount"))
  frac_col <- guess_col(df_raw, c("Proportion", "readFraction", "cloneFraction"))

  if (is.na(aa_col)) stop("No AA column found in ", sample_name, ". Columns: ", paste(colnames(df_raw), collapse = ", "))

  if (!is.na(cnt_col)) {
    df <- df_raw %>%
      transmute(
        cdr3aa = .data[[aa_col]],
        count  = suppressWarnings(as.numeric(.data[[cnt_col]])),
        v_raw  = if ("allVHitsWithScore" %in% colnames(df_raw)) df_raw$allVHitsWithScore else NA_character_,
        j_raw  = if ("allJHitsWithScore" %in% colnames(df_raw)) df_raw$allJHitsWithScore else NA_character_
      )
  } else if (!is.na(frac_col)) {
    df <- df_raw %>%
      transmute(
        cdr3aa = .data[[aa_col]],
        count  = suppressWarnings(as.numeric(.data[[frac_col]])),
        v_raw  = if ("allVHitsWithScore" %in% colnames(df_raw)) df_raw$allVHitsWithScore else NA_character_,
        j_raw  = if ("allJHitsWithScore" %in% colnames(df_raw)) df_raw$allJHitsWithScore else NA_character_
      )
  } else {
    stop("No count/fraction column found in ", sample_name)
  }

  df <- df %>% filter(!is.na(cdr3aa), cdr3aa != "", is.finite(count), count > 0)

  collapsed <- df %>%
    group_by(cdr3aa) %>%
    summarise(
      Clones = sum(count, na.rm = TRUE),
      V.name = clean_gene(v_raw[which.max(count)]),
      J.name = clean_gene(j_raw[which.max(count)]),
      .groups = "drop"
    ) %>%
    arrange(desc(Clones)) %>%
    mutate(
      Proportion = Clones / sum(Clones),
      `CDR3.aa`  = cdr3aa
    ) %>%
    select(Clones, Proportion, `CDR3.aa`, V.name, J.name)

  list(raw_n = nrow(df_raw), collapsed_n = nrow(collapsed), table = collapsed)
}

parsed   <- pmap(mixcr_index, function(Sample, file) read_and_collapse_mixcr(file, Sample))
names(parsed) <- mixcr_index$Sample
rep_list <- map(parsed, "table")

iwalk(rep_list, ~ readr::write_tsv(.x, file.path(collapsed_dir, paste0(.y, ".collapsed_by_cdr3aa.tsv"))))

collapse_qc <- tibble(
  Sample = names(parsed),
  n_rows_raw = map_int(parsed, "raw_n"),
  n_rows_collapsed = map_int(parsed, "collapsed_n"),
  frac_remaining = n_rows_collapsed / pmax(n_rows_raw, 1)
) %>%
  left_join(meta, by = "Sample") %>%
  arrange(Paper_ID, Tissue, Day)

readr::write_tsv(collapse_qc, file.path(out_dir, paste0("qc_collapse_by_cdr3aa_", chain, ".tsv")))

# Repertoire metrics
calc_metrics_one <- function(df, sample_name) {
  p <- df$Proportion
  p <- p[is.finite(p) & !is.na(p) & p > 0]
  if (length(p) == 0) {
    return(tibble(
      Sample = sample_name,
      n_clonotypes = 0,
      top1_prop = NA_real_,
      shannon = NA_real_,
      invsimpson = NA_real_,
      pielou = NA_real_,
      clonality_1mPielou = NA_real_
    ))
  }
  p <- p / sum(p)
  p_sorted <- sort(p, decreasing = TRUE)

  conc <- sum(p^2)
  sh   <- -sum(p * log(p))
  S    <- length(p_sorted)

  J <- if (S > 1) sh / log(S) else 0
  J <- pmin(pmax(J, 0), 1)

  tibble(
    Sample            = sample_name,
    n_clonotypes      = S,
    top1_prop         = p_sorted[1],
    shannon           = sh,
    invsimpson        = 1 / conc,
    pielou            = J,
    clonality_1mPielou = 1 - J
  )
}

metrics <- imap_dfr(rep_list, ~ calc_metrics_one(.x, .y)) %>%
  left_join(meta, by = "Sample") %>%
  arrange(Paper_ID, Tissue, Day)

readr::write_tsv(metrics, file.path(out_dir, paste0("metrics_", chain, "_cdr3aa.tsv")))

gi_iec_control <- metrics %>%
  filter(Group == "GI", Pathology %in% c("IEC", "Control"),
         !(Sample %in% exclude_for_iec_control)) %>%
  mutate(CamelidGroup = factor(
    CamelidGroup,
    levels = c("No camelid infiltrate", "Camelid infiltrate present"))
  )

w <- suppressWarnings(wilcox.test(top1_prop ~ CamelidGroup, data = gi_iec_control, exact = FALSE))
stats_out <- tibble(metric = "top1_prop", p_value = w$p.value,
                    n_no  = sum(gi_iec_control$CamelidGroup == "No camelid infiltrate"),
                    n_yes = sum(gi_iec_control$CamelidGroup == "Camelid infiltrate present"))
readr::write_tsv(stats_out, file.path(out_dir, paste0("wilcox_top1prop_camelid_", chain, ".tsv")))

geom_col_outlined <- function(...) geom_col(..., color = "black", linewidth = 0.25)

box_blackdots <- function(df, y, ylab, title) {
  ggplot(df, aes(x = CamelidGroup, y = .data[[y]], fill = CamelidGroup)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65, color = "black", linewidth = 0.25) +
    geom_jitter(width = 0.10, height = 0, size = 3.0, alpha = 0.95, color = "black") +
    scale_fill_manual(values = pal_iec) +
    labs(title = title, x = NULL, y = ylab) +
    theme(legend.position = "none")
}

plot_gi_top1_vs_other <- function(df_metrics) {
  d <- df_metrics %>%
    filter(is.finite(top1_prop)) %>%
    mutate(
      Top = top1_prop,
      Other = pmax(0, 1 - top1_prop),
      SampleLab = paste0(Paper_ID, " | ", as.character(Tissue))
    ) %>%
    arrange(CamelidGroup, Paper_ID, Tissue, Day) %>%
    mutate(SampleLab = factor(SampleLab, levels = unique(SampleLab)))

  long <- d %>%
    select(CamelidGroup, SampleLab, Top, Other) %>%
    pivot_longer(c(Top, Other), names_to = "Component", values_to = "Fraction") %>%
    mutate(Component = recode(Component,
                              Top = "Top clone (CDR3aa)",
                              Other = "All other clones"))

  ggplot(long, aes(x = Fraction, y = SampleLab, fill = Component)) +
    geom_col_outlined(width = 0.85) +
    facet_grid(CamelidGroup ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Top clone (CDR3aa)" = col_top_aqua,
                                 "All other clones"   = col_other_salmon)) +
    labs(
      title = "TRB GI clonality (collapsed by CDR3aa)",
      subtitle = "Stacked: Top clone vs All other clones. Excludes Marrow + LN.",
      x = "Fraction of repertoire", y = NULL, fill = "Component"
    ) +
    theme(strip.background = element_rect(fill = "grey95", color = NA))
}

plot_pair_top1_vs_other <- function(samples, title) {
  samples <- intersect(samples, metrics$Sample)
  if (length(samples) == 0) return(NULL)

  d <- metrics %>% filter(Sample %in% samples)

  if (!all(c("Day", "Tissue", "Paper_ID") %in% colnames(d))) {
    d <- d %>%
      left_join(meta %>% select(Sample, Day, Tissue, Paper_ID), by = "Sample",
                suffix = c("", ".meta")) %>%
      mutate(
        Day      = dplyr::coalesce(.data$Day, .data$Day.meta),
        Tissue   = dplyr::coalesce(.data$Tissue, .data$Tissue.meta),
        Paper_ID = dplyr::coalesce(.data$Paper_ID, .data$Paper_ID.meta)
      ) %>%
      select(-ends_with(".meta"))
  }

  d <- d %>%
    mutate(Day = as.numeric(Day)) %>%
    arrange(Day) %>%
    mutate(
      Top = top1_prop,
      Other = pmax(0, 1 - top1_prop),
      SampleLab = paste0(Paper_ID, " | ", "Day ", Day, " | ", as.character(Tissue))
    ) %>%
    mutate(SampleLab = factor(SampleLab, levels = SampleLab))

  long <- d %>%
    select(SampleLab, Top, Other) %>%
    pivot_longer(c(Top, Other), names_to = "Component", values_to = "Fraction") %>%
    mutate(Component = recode(Component,
                              Top = "Top clone (CDR3aa)",
                              Other = "All other clones"))

  ggplot(long, aes(x = Fraction, y = SampleLab, fill = Component)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = c("Top clone (CDR3aa)" = col_top_aqua,
                                 "All other clones"   = col_other_salmon)) +
    labs(title = title, x = "Fraction of repertoire", y = NULL, fill = "Component")
}

plot_pair_topN <- function(samples, topN = 10, title) {
  samples <- intersect(samples, names(rep_list))
  if (length(samples) == 0) return(NULL)

  ord    <- meta %>% filter(Sample %in% samples) %>% arrange(Day) %>% pull(Sample)
  lab_map <- time_label(ord)

  long <- map_dfr(ord, function(s) {
    rep_list[[s]] %>% select(`CDR3.aa`, Proportion) %>% mutate(Sample = s)
  })

  top_clones <- long %>%
    group_by(`CDR3.aa`) %>%
    summarise(max_p = max(Proportion, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_p)) %>%
    slice_head(n = topN) %>%
    pull(`CDR3.aa`)

  long2 <- long %>%
    mutate(Clone = if_else(`CDR3.aa` %in% top_clones, `CDR3.aa`, "Other")) %>%
    group_by(Sample, Clone) %>%
    summarise(Proportion = sum(Proportion), .groups = "drop") %>%
    mutate(Sample = factor(Sample, levels = ord))

  clones <- setdiff(unique(long2$Clone), "Other")
  cols   <- muted_palette(max(1, length(clones)))
  names(cols) <- clones
  cols <- c(cols, Other = col_other_salmon)

  ggplot(long2, aes(x = Sample, y = Proportion, fill = Clone)) +
    geom_col_outlined(width = 0.85) +
    coord_flip() +
    scale_x_discrete(labels = lab_map) +
    scale_fill_manual(values = cols) +
    labs(title = title, x = NULL, y = "Fraction of repertoire", fill = "Clone") +
    theme(legend.position = "right")
}

p_gi_top1 <- plot_gi_top1_vs_other(gi_iec_control)

p_box_top1 <- box_blackdots(
  gi_iec_control, "top1_prop", "Top clone fraction",
  "GI clonality: Top clone fraction"
) + labs(subtitle = paste0("Wilcoxon p = ", signif(stats_out$p_value[1], 3)))

p_box_clonality <- box_blackdots(
  gi_iec_control, "clonality_1mPielou", "Clonality (1 − Pielou)",
  "GI clonality: 1 − Pielou's evenness"
)

p_box_inv <- box_blackdots(
  gi_iec_control, "invsimpson", "Inverse Simpson (higher = more diverse)",
  "GI diversity: Inverse Simpson"
)

# Optional within-subject comparisons
pair_path <- file.path(data_dir, "repertoire_pairs.tsv")
pair_metadata <- if (file.exists(pair_path)) {
  readr::read_tsv(pair_path, show_col_types = FALSE)
} else {
  tibble::tibble(comparison_id = character(), Sample = character(), title = character())
}
for (comparison in unique(pair_metadata$comparison_id)) {
  pair <- dplyr::filter(pair_metadata, comparison_id == comparison)
  pair_title <- pair$title[[1]]
  plot_top10 <- plot_pair_topN(pair$Sample, topN = 10, title = pair_title)
  plot_top1 <- plot_pair_top1_vs_other(pair$Sample, title = pair_title)
  stem <- gsub("[^A-Za-z0-9_-]", "_", comparison)
  if (!is.null(plot_top10)) {
    save_both(reverse_sample_factor(plot_top10), paste0(stem, "_top10"), width = 12, height = 8)
  }
  if (!is.null(plot_top1)) {
    save_both(reverse_sample_factor(plot_top1), paste0(stem, "_top1"), width = 12, height = 6)
  }
}

save_both(p_gi_top1,        "Fig_GI_topclone_vs_other",   width = 14, height = 10)
save_both(p_box_top1,       "Fig_Box_topclone_fraction",  width = 9,  height = 7)
save_both(p_box_clonality,  "Fig_Box_clonality_1mPielou", width = 9,  height = 7)
save_both(p_box_inv,        "Fig_Box_invsimpson",         width = 9,  height = 7)

strip_long <- gi_iec_control %>%
  select(Sample, CamelidGroup, shannon, invsimpson, clonality_1mPielou) %>%
  pivot_longer(
    cols = c(shannon, invsimpson, clonality_1mPielou),
    names_to = "Metric", values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      shannon            = "Shannon",
      invsimpson         = "Inverse Simpson",
      clonality_1mPielou = "Clonality (1 − Pielou)"
    ),

    Metric = factor(Metric, levels = c("Shannon", "Inverse Simpson", "Clonality (1 − Pielou)"))
  )

w_strip <- strip_long %>%
  group_by(Metric) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ CamelidGroup, data = cur_data(), exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(p_label = ifelse(is.na(p_value), "p = NA",
                          ifelse(p_value < 1e-4, "p < 1e-4",
                                 paste0("p = ", formatC(p_value, format = "e", digits = 2)))))

lab_pos <- strip_long %>%
  group_by(Metric) %>%
  summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop") %>%
  left_join(w_strip, by = "Metric")

p_strip_metrics <- ggplot(strip_long, aes(x = CamelidGroup, y = Value, fill = CamelidGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65, color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.10, height = 0, size = 3.0, alpha = 0.95, color = "black") +
  scale_fill_manual(values = pal_iec) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
  labs(
    title = "Repertoire diversity/clonality (GI only)",
    subtitle = "Shannon, Inverse Simpson, and Clonality (1 − Pielou)",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme(legend.position = "none")

p_strip_metrics <- p_strip_metrics +
  geom_text(
    data = lab_pos,
    aes(x = 1.5, y = ymax * 1.05, label = p_label),
    inherit.aes = FALSE, vjust = 0, size = 5
  )

print(p_strip_metrics)
save_both(p_strip_metrics, "Fig_Strip_Shannon_InverseSimpson_1minusPielou", width = 16, height = 6)

message("[DONE] Outputs:")
message("  - Tables: ", out_dir)
message("  - Figures (PDF+EPS): ", fig_dir)
if (use_showtext) message("  - Font rendering: showtext enabled (Arial registered)")
