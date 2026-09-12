#!/usr/bin/env Rscript
# Delayed neurotoxicity associations
# Relate delayed neurotoxicity to enterocolitis, lymphocyte counts, CAR exposure, and persistence.

suppressPackageStartupMessages({
  library(readr);  library(dplyr);  library(tidyr);  library(stringr)
  library(ggplot2); library(scales); library(forcats); library(tibble)
  library(broom); library(purrr)
  library(splines); library(lme4); library(lmerTest); library(emmeans)
  library(grid); library(gridExtra)
})

emmeans::emm_options(pbkrtest.limit = 50000, lmerTest.limit = 50000)

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")

path_df3 <- file.path(data_dir, "clinical_metadata.tsv")
path_df4 <- file.path(data_dir, "car_measurements.tsv")
path_alc_long <- file.path(data_dir, "lymphocyte_measurements.tsv")

out_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "car_kinetics", "delayed_neurotoxicity"), "figures")
stats_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "car_kinetics", "delayed_neurotoxicity"), "stats")
tables_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "car_kinetics", "delayed_neurotoxicity"), "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

.is_mac <- identical(Sys.info()[["sysname"]], "Darwin")
.has_cairo_device <- isTRUE(capabilities("cairo"))
.has_Cairo_pkg <- requireNamespace("Cairo", quietly = TRUE)

if (.is_mac) {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

save_pdf <- function(plot, filename, width = 9, height = 7) {
  f <- file.path(out_dir, filename)
  device_fun <- NULL

  if (.has_Cairo_pkg) {
    device_fun <- function(file, width, height, ...) {
      Cairo::CairoPDF(file = file, width = width, height = height, family = "Arial", ...)
    }
  } else if (.is_mac) {
    device_fun <- function(file, width, height, ...) {
      grDevices::quartz(file = file, type = "pdf", width = width, height = height, family = "Arial", ...)
    }
  } else {
    device_fun <- function(file, width, height, ...) {
      grDevices::pdf(file = file, width = width, height = height, family = "Helvetica", ...)
    }
  }

  ggplot2::ggsave(filename = f, plot = plot, device = device_fun,
                  width = width, height = height, units = "in")
  message("Saved: ", f)
}

save_eps <- function(plot, filename, width = 9, height = 7) {
  f <- file.path(out_dir, filename)
  device_fun <- NULL

  if (.has_cairo_device) {

    device_fun <- function(file, width, height, ...) {
      grDevices::cairo_ps(
        filename = file,
        width = width,
        height = height,
        onefile = FALSE,
        fallback_resolution = 600,
        family = "Arial",
        ...
      )
    }
  } else {

    device_fun <- function(file, width, height, ...) {
      grDevices::postscript(
        file = file,
        onefile = FALSE,
        paper = "special",
        width = width,
        height = height,
        family = "Helvetica",
        horizontal = FALSE,
        ...
      )
    }
  }

  ggplot2::ggsave(filename = f, plot = plot, device = device_fun,
                  width = width, height = height, units = "in")
  message("Saved: ", f)
}

save_csv <- function(df, filename) {
  f <- file.path(stats_dir, filename)
  readr::write_csv(df, f)
  message("Saved: ", f)
}

save_tsv <- function(df, filename) {
  f <- file.path(tables_dir, filename)
  readr::write_tsv(df, f)
  message("Saved: ", f)
}

safe_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(readr::parse_number(as.character(x)))
}

safe_chr <- function(x) {
  y <- as.character(x)
  y <- stringr::str_squish(y)
  y[y %in% c("", "NA", "Na", "na", "N/A", "n/a")] <- NA_character_
  y
}

pick_col_num <- function(df, aliases) {
  nm <- aliases[aliases %in% names(df)][1]
  if (length(nm) == 0 || is.na(nm)) return(rep(NA_real_, nrow(df)))
  safe_num(df[[nm]])
}

pick_col_chr <- function(df, aliases) {
  nm <- aliases[aliases %in% names(df)][1]
  if (length(nm) == 0 || is.na(nm)) return(rep(NA_character_, nrow(df)))
  safe_chr(df[[nm]])
}

theme_arial <- function(base_size = 22) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = base_size),
      axis.text  = element_text(size = base_size - 4)
    )
}

std_ci_cols <- function(df) {
  if (all(c("lower.CL", "upper.CL") %in% names(df))) return(df)
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(df))) {
    return(dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL))
  }
  if (all(c("LCL", "UCL") %in% names(df))) {
    return(dplyr::rename(df, lower.CL = LCL, upper.CL = UCL))
  }
  stop("Could not locate CI columns. Found: ", paste(names(df), collapse = ", "))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 3))
}

fisher_tbl <- function(tab, comparison) {
  ft <- fisher.test(tab)
  tibble(
    comparison = comparison,
    test = "Fisher exact",
    odds_ratio = unname(ft$estimate),
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2],
    p.value = ft$p.value
  )
}

summ_cont <- function(df, group, value) {
  df %>%
    filter(!is.na(.data[[group]]), is.finite(.data[[value]])) %>%
    group_by(.data[[group]]) %>%
    summarise(
      n = n(),
      median = median(.data[[value]], na.rm = TRUE),
      q1 = quantile(.data[[value]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[value]], 0.75, na.rm = TRUE),
      min = min(.data[[value]], na.rm = TRUE),
      max = max(.data[[value]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = 1)
}

boxplot_log_by_dnt <- function(df, yvar, ylab, title, subtitle, filename_base, width = 7.5, height = 7.5,
                               hline = NULL, hline_label = NULL) {
  dd <- df %>% filter(!is.na(DNT_flag), is.finite(.data[[yvar]]), .data[[yvar]] > 0)
  if (nrow(dd) < 3 || nlevels(droplevels(dd$DNT_flag)) < 2) {
    message("Skipping ", filename_base, ": insufficient data.")
    return(NULL)
  }

  p <- ggplot(dd, aes(x = DNT_flag, y = .data[[yvar]], fill = DNT_flag)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.9) +
    geom_jitter(width = 0.10, alpha = 0.55, size = 2, color = "gray20") +
    scale_fill_manual(values = pal_DNT, guide = "none") +
    scale_y_log10(ylab, labels = scales::label_number(accuracy = 0.01, big.mark = ",")) +
    xlab("Delayed neurotoxicity") +
    labs(title = title, subtitle = subtitle) +
    theme_arial(22)

  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey35", linewidth = 0.7)
    if (!is.null(hline_label)) {
      p <- p + annotate("text", x = 1.5, y = hline, label = hline_label,
                        vjust = -0.5, size = 4.5, family = "Arial")
    }
  }

  print(p)
  save_pdf(p, paste0(filename_base, ".pdf"), width = width, height = height)
  save_eps(p, paste0(filename_base, ".eps"), width = width, height = height)
  p
}

bar_percent_yes <- function(summary_df, xvar, fill_values, title, subtitle, xlab, filename_base,
                            width = 7.5, height = 7) {
  dd <- summary_df %>% filter(DNT_flag == "Yes")
  p <- ggplot(dd, aes(x = .data[[xvar]], y = percent_DNT, fill = .data[[xvar]])) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_text(
      aes(label = sprintf("%.1f%%\n%d/%d", percent_DNT, n, group_n)),
      vjust = -0.25,
      size = 5,
      family = "Arial"
    ) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_y_continuous(
      "Delayed neurotoxicity (%)",
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.18))
    ) +
    xlab(xlab) +
    labs(title = title, subtitle = subtitle) +
    theme_arial(22)

  print(p)
  save_pdf(p, paste0(filename_base, ".pdf"), width = width, height = height)
  save_eps(p, paste0(filename_base, ".eps"), width = width, height = height)
  p
}

pal_IEC <- c("No IEC" = "#4C6A87", "IEC" = "#B22222")
pal_DNT <- c("No" = "#4C6A87", "Yes" = "#D95F02")
pal_ALC <- c("\u22643 K/uL" = "#4C6A87", ">3 K/uL" = "#B22222")
LOD_ND <- 0.1

# Clinical metadata and delayed neurotoxicity
df3_raw <- readr::read_delim(path_df3, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
if (!"Study_ID" %in% names(df3_raw) && "ID" %in% names(df3_raw)) {
  df3_raw <- dplyr::rename(df3_raw, Study_ID = ID)
}

if ("Day" %in% names(df3_raw)) {
  df3_raw <- df3_raw %>% dplyr::rename(DNT_day = Day)
}

cat("\n[DNT columns present]\n")
print(intersect(c("Parkinsonian_Neurotoxicity", "Other_Neurotoxicity",
                  "Any_Delayed_Neurotoxicity", "What_Neurotoxicity", "DNT_day"),
                names(df3_raw)))

df3 <- df3_raw %>%
  mutate(
    IEC_enteritis = safe_num(IEC_enteritis),
    IEC_flag = factor(if_else(IEC_enteritis == 1, "IEC", "No IEC"), levels = c("No IEC", "IEC")),
    Any_Delayed_Neurotoxicity_chr = pick_col_chr(cur_data(), c("Any_Delayed_Neurotoxicity")),
    DNT_flag = factor(
      if_else(str_to_lower(str_squish(Any_Delayed_Neurotoxicity_chr)) == "yes", "Yes", "No"),
      levels = c("No", "Yes")
    ),
    Parkinsonian_DNT = factor(
      if_else(str_to_lower(str_squish(pick_col_chr(cur_data(), c("Parkinsonian_Neurotoxicity")))) == "yes", "Yes", "No"),
      levels = c("No", "Yes")
    ),
    Other_DNT = factor(
      if_else(str_to_lower(str_squish(pick_col_chr(cur_data(), c("Other_Neurotoxicity")))) == "yes", "Yes", "No"),
      levels = c("No", "Yes")
    ),
    DNT_day = pick_col_num(cur_data(), c("DNT_day", "Delayed_Neurotoxicity_Day", "Neurotoxicity_Day")),
    What_Neurotoxicity = pick_col_chr(cur_data(), c("What_Neurotoxicity"))
  )

dnt_meta <- df3 %>%
  select(Study_ID, IEC_flag, DNT_flag, Parkinsonian_DNT, Other_DNT, What_Neurotoxicity, DNT_day)

cat("\n================ DNT prevalence ================\n")
dnt_prev <- dnt_meta %>%
  count(DNT_flag, name = "n") %>%
  mutate(percent = 100 * n / sum(n))
print(dnt_prev)
cat(sprintf("\nDNT prevalence: %.1f%% (%d/%d)\n",
            100 * sum(dnt_meta$DNT_flag == "Yes", na.rm = TRUE) / nrow(dnt_meta),
            sum(dnt_meta$DNT_flag == "Yes", na.rm = TRUE), nrow(dnt_meta)))
save_csv(dnt_prev, "DNT_prevalence.csv")

cat("\n================ DNT details among cases ================\n")
dnt_cases <- dnt_meta %>% filter(DNT_flag == "Yes") %>% arrange(DNT_day, Study_ID)
print(dnt_cases, n = Inf)
save_csv(dnt_cases, "DNT_case_details.csv")

cat("\n================ IEC-EC vs delayed neurotoxicity ================\n")
tab_iec_dnt <- table(IEC = dnt_meta$IEC_flag, DNT = dnt_meta$DNT_flag)
print(tab_iec_dnt)
fisher_iec_dnt <- fisher.test(tab_iec_dnt)
print(fisher_iec_dnt)

iec_dnt_summary <- as.data.frame(tab_iec_dnt) %>%
  as_tibble() %>%
  group_by(IEC) %>%
  mutate(group_n = sum(Freq), percent_within_IEC = 100 * Freq / group_n) %>%
  ungroup()
save_csv(iec_dnt_summary, "DNT_by_IEC_counts.csv")

iec_dnt_test_tbl <- tibble(
  comparison = "IEC-EC vs Any delayed neurotoxicity",
  test = "Fisher exact",
  odds_ratio = unname(fisher_iec_dnt$estimate),
  conf_low = fisher_iec_dnt$conf.int[1],
  conf_high = fisher_iec_dnt$conf.int[2],
  p.value = fisher_iec_dnt$p.value
)
save_csv(iec_dnt_test_tbl, "DNT_by_IEC_Fisher.csv")

p_iec_dnt_bar <- iec_dnt_summary %>%
  filter(DNT == "Yes") %>%
  ggplot(aes(x = IEC, y = percent_within_IEC, fill = IEC)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%\n%d/%d", percent_within_IEC, Freq, group_n)),
            vjust = -0.25, size = 5, family = "Arial") +
  scale_fill_manual(values = pal_IEC, guide = "none") +
  scale_y_continuous("Delayed neurotoxicity (%)",
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.18))) +
  xlab("") +
  labs(title = "Delayed neurotoxicity by IEC-EC status",
       subtitle = paste0("Fisher exact p=", fmt_p(fisher_iec_dnt$p.value))) +
  theme_arial(22)
print(p_iec_dnt_bar)
save_pdf(p_iec_dnt_bar, "DNT_by_IEC_bar.pdf", width = 7.5, height = 7)
save_eps(p_iec_dnt_bar, "DNT_by_IEC_bar.eps", width = 7.5, height = 7)

# Longitudinal lymphocyte counts
read_alc_kperul <- function(path, src_label) {
  raw <- readr::read_delim(path, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
  req <- c("Study_ID", "ALC", "Day")
  if (!all(req %in% names(raw))) {
    stop("ALC file lacks required columns Study_ID, ALC, Day. Found: ", paste(names(raw), collapse = ", "))
  }
  df <- raw %>%
    transmute(Study_ID, Day = safe_num(Day), ALC = safe_num(ALC)) %>%
    filter(!is.na(Study_ID), is.finite(Day), Day >= 0, is.finite(ALC)) %>%
    mutate(ALC_before = ALC, ALC = if_else(ALC <= 0, 0.01, ALC))
  audit <- df %>%
    filter(ALC_before <= 0) %>%
    transmute(source = src_label, Study_ID, Day,
              ALC_original_KperuL = ALC_before,
              ALC_used_KperuL = ALC,
              reason = "zero_floored_to_0.01")
  if (nrow(audit) > 0) save_csv(audit, "QC_DNT_ALC_zero_floor_audit.csv")
  df %>% select(Study_ID, Day, ALC)
}

alc_long <- read_alc_kperul(path_alc_long, basename(path_alc_long))
alc_dnt_long <- alc_long %>%
  left_join(dnt_meta, by = "Study_ID") %>%
  mutate(DNT_flag = factor(DNT_flag, levels = c("No", "Yes")))

alc_cutoff_0_100 <- alc_dnt_long %>%
  filter(Day >= 0, Day <= 100, is.finite(ALC)) %>%
  group_by(Study_ID) %>%
  summarise(
    n_ALC_0_100 = n(),
    Peak_ALC_0_100 = max(ALC, na.rm = TRUE),
    Ever_ALC_gt3_0_100 = any(ALC > 3, na.rm = TRUE),
    First_Day_ALC_gt3_0_100 = ifelse(any(ALC > 3, na.rm = TRUE), min(Day[ALC > 3], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(ALC_gt3_0_100 = factor(if_else(Ever_ALC_gt3_0_100, ">3 K/uL", "\u22643 K/uL"),
                                levels = c("\u22643 K/uL", ">3 K/uL"))) %>%
  left_join(dnt_meta, by = "Study_ID")

cat("\n================ ALC >3 K/uL D0-100 vs DNT ================\n")
tab_alc_gt3 <- table(ALC_gt3_D0_100 = alc_cutoff_0_100$ALC_gt3_0_100,
                     DNT = alc_cutoff_0_100$DNT_flag)
print(tab_alc_gt3)
fisher_alc_gt3 <- fisher.test(tab_alc_gt3)
print(fisher_alc_gt3)

alc_gt3_test_tbl <- tibble(
  comparison = "Any ALC >3 K/uL from D0-100 vs Any delayed neurotoxicity",
  cutoff = ">3 K/uL, equivalent to >3000/uL",
  test = "Fisher exact",
  odds_ratio = unname(fisher_alc_gt3$estimate),
  conf_low = fisher_alc_gt3$conf.int[1],
  conf_high = fisher_alc_gt3$conf.int[2],
  p.value = fisher_alc_gt3$p.value
)
save_csv(alc_gt3_test_tbl, "DNT_ALC_gt3_D0_100_Fisher.csv")

alc_gt3_summary <- alc_cutoff_0_100 %>%
  count(ALC_gt3_0_100, DNT_flag, name = "n") %>%
  group_by(ALC_gt3_0_100) %>%
  mutate(group_n = sum(n), percent_DNT = 100 * n / group_n) %>%
  ungroup()
save_csv(alc_gt3_summary, "DNT_ALC_gt3_D0_100_counts.csv")

p_alc_gt3_bar <- bar_percent_yes(
  summary_df = alc_gt3_summary,
  xvar = "ALC_gt3_0_100",
  fill_values = pal_ALC,
  title = "Delayed neurotoxicity by ALC >3 K/uL",
  subtitle = paste0("Fisher exact p=", fmt_p(fisher_alc_gt3$p.value)),
  xlab = "Peak ALC threshold from D0-100",
  filename_base = "DNT_by_ALC_gt3_D0_100_bar",
  width = 8,
  height = 7
)

cat("\n================ Peak ALC D0-100 by DNT ================\n")
wilcox_peak_alc <- wilcox.test(Peak_ALC_0_100 ~ DNT_flag, data = alc_cutoff_0_100, exact = FALSE)
print(wilcox_peak_alc)
peak_alc_summary <- summ_cont(alc_cutoff_0_100, "DNT_flag", "Peak_ALC_0_100")
save_csv(peak_alc_summary, "DNT_Peak_ALC_D0_100_summary.csv")

p_peak_alc <- boxplot_log_by_dnt(
  df = alc_cutoff_0_100,
  yvar = "Peak_ALC_0_100",
  ylab = "Peak ALC D0-100 (K/uL)",
  title = "Peak ALC by delayed neurotoxicity status",
  subtitle = paste0("Wilcoxon p=", fmt_p(wilcox_peak_alc$p.value)),
  filename_base = "Peak_ALC_D0_100_by_DNT",
  hline = 3,
  hline_label = "3 K/uL cutoff"
)

alc_window_summary <- function(df, start_day, end_day, label) {
  df %>%
    filter(Day >= start_day, Day <= end_day, is.finite(ALC)) %>%
    group_by(Study_ID) %>%
    summarise(window = label, start_day = start_day, end_day = end_day,
              n_ALC = n(), Peak_ALC = max(ALC, na.rm = TRUE),
              Ever_ALC_gt3 = any(ALC > 3, na.rm = TRUE), .groups = "drop") %>%
    mutate(ALC_gt3 = factor(if_else(Ever_ALC_gt3, ">3 K/uL", "\u22643 K/uL"),
                            levels = c("\u22643 K/uL", ">3 K/uL"))) %>%
    left_join(dnt_meta, by = "Study_ID")
}

alc_windows <- bind_rows(
  alc_window_summary(alc_long, 0, 30, "D0-30"),
  alc_window_summary(alc_long, 31, 100, "D31-100"),
  alc_window_summary(alc_long, 0, 200, "D0-200")
)

alc_window_tests <- alc_windows %>%
  group_by(window) %>%
  group_modify(~{
    tab <- table(.x$ALC_gt3, .x$DNT_flag)
    ft <- fisher.test(tab)
    tibble(n_patients = n_distinct(.x$Study_ID),
           odds_ratio = unname(ft$estimate),
           conf_low = ft$conf.int[1], conf_high = ft$conf.int[2], p.value = ft$p.value)
  }) %>%
  ungroup()
cat("\n================ ALC >3 K/uL window sensitivity tests ================\n")
print(alc_window_tests)
save_csv(alc_window_tests, "DNT_ALC_gt3_window_sensitivity_Fisher.csv")
save_csv(alc_windows, "DNT_ALC_window_patient_level_metrics.csv")

alc_pre_event <- alc_dnt_long %>%
  mutate(eligible_pre_event = case_when(
    DNT_flag == "Yes" & is.finite(DNT_day) ~ Day <= DNT_day,
    DNT_flag == "Yes" & !is.finite(DNT_day) ~ Day <= 100,
    DNT_flag == "No" ~ Day <= 100,
    TRUE ~ FALSE
  )) %>%
  filter(Day >= 0, eligible_pre_event, is.finite(ALC)) %>%
  group_by(Study_ID) %>%
  summarise(n_ALC_pre_event = n(),
            Peak_ALC_pre_event = max(ALC, na.rm = TRUE),
            Ever_ALC_gt3_pre_event = any(ALC > 3, na.rm = TRUE), .groups = "drop") %>%
  mutate(ALC_gt3_pre_event = factor(if_else(Ever_ALC_gt3_pre_event, ">3 K/uL", "\u22643 K/uL"),
                                    levels = c("\u22643 K/uL", ">3 K/uL"))) %>%
  left_join(dnt_meta, by = "Study_ID")

cat("\n================ Pre-event ALC >3 K/uL sensitivity vs DNT ================\n")
tab_alc_pre <- table(alc_pre_event$ALC_gt3_pre_event, alc_pre_event$DNT_flag)
print(tab_alc_pre)
fisher_alc_pre <- fisher.test(tab_alc_pre)
print(fisher_alc_pre)
save_csv(tibble(comparison = "Pre-event ALC >3 K/uL sensitivity",
                test = "Fisher exact",
                odds_ratio = unname(fisher_alc_pre$estimate),
                conf_low = fisher_alc_pre$conf.int[1],
                conf_high = fisher_alc_pre$conf.int[2],
                p.value = fisher_alc_pre$p.value),
         "DNT_ALC_gt3_pre_event_sensitivity_Fisher.csv")

alc_plot_df <- alc_dnt_long %>%
  filter(Day >= 0, Day <= 200, is.finite(ALC), ALC > 0) %>%
  mutate(DNT_flag = factor(DNT_flag, levels = c("No", "Yes")))

p_alc_dnt_spaghetti <- ggplot(alc_plot_df, aes(x = Day, y = ALC, group = Study_ID, color = DNT_flag)) +
  geom_line(alpha = 0.25, linewidth = 0.55) +
  geom_smooth(aes(group = DNT_flag, color = DNT_flag), method = "loess", se = TRUE, span = 0.55, linewidth = 1.25) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "grey35", linewidth = 0.65) +
  annotate("text", x = 190, y = 3, label = "3 K/uL", hjust = 1, vjust = -0.5, family = "Arial", size = 4) +
  scale_color_manual(values = pal_DNT, name = "Delayed neurotoxicity") +
  scale_y_log10("ALC (K/uL)", labels = scales::label_number(accuracy = 0.01)) +
  scale_x_continuous("Days post infusion", limits = c(0, 200)) +
  labs(title = "ALC over time by delayed neurotoxicity status",
       subtitle = "Spaghetti plot with LOESS group trends; dashed line indicates 3 K/uL") +
  theme_arial(22)
print(p_alc_dnt_spaghetti)
save_pdf(p_alc_dnt_spaghetti, "LONG_ALC_by_DNT_spaghetti_0_200.pdf", width = 11, height = 6.2)
save_eps(p_alc_dnt_spaghetti, "LONG_ALC_by_DNT_spaghetti_0_200.eps", width = 11, height = 6.2)

alc_lmm_df <- alc_dnt_long %>%
  filter(Day >= 0, Day <= 100, is.finite(ALC), ALC > 0) %>%
  mutate(DNT_flag = factor(DNT_flag, levels = c("No", "Yes")))

if (nlevels(droplevels(alc_lmm_df$DNT_flag)) == 2 &&
    n_distinct(alc_lmm_df$Study_ID[alc_lmm_df$DNT_flag == "Yes"]) >= 3) {
  fit_alc_dnt <- lmerTest::lmer(log10(ALC) ~ DNT_flag * splines::ns(Day, df = 4) + (1 | Study_ID),
                                data = alc_lmm_df, REML = TRUE)
  cat("\n================ ALC LMM by DNT status, D0-100 ================\n")
  print(anova(fit_alc_dnt, type = 3))
  anova_alc_dnt <- as.data.frame(anova(fit_alc_dnt, type = 3)) %>% rownames_to_column("term")
  save_csv(anova_alc_dnt, "DNT_ALC_LMM_type3_ANOVA.csv")

  days_grid <- c(0, 7, 14, 21, 28, 42, 60, 90, 100)
  emm_link <- emmeans::emmeans(fit_alc_dnt, ~ DNT_flag | Day, at = list(Day = days_grid))
  contr_dnt_alc <- summary(emmeans::contrast(emm_link, method = "revpairwise", by = "Day"),
                           infer = c(TRUE, TRUE), level = 0.95) %>%
    as.data.frame() %>% std_ci_cols() %>%
    mutate(ratio_DNTyes_No = 10^estimate,
           lower_ratio = 10^lower.CL,
           upper_ratio = 10^upper.CL)
  save_csv(contr_dnt_alc, "DNT_ALC_LMM_pairwise_contrasts_by_Day.csv")

  p_alc_lmm_ratio <- ggplot(contr_dnt_alc, aes(x = Day, y = ratio_DNTyes_No)) +
    geom_ribbon(aes(ymin = lower_ratio, ymax = upper_ratio), fill = "#D95F02", alpha = 0.18) +
    geom_line(color = "#D95F02", linewidth = 1.1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_y_log10("DNT Yes / No geometric-mean ratio of ALC",
                  breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5),
                  labels = scales::label_number(accuracy = 0.01)) +
    scale_x_continuous("Days post infusion", breaks = days_grid) +
    labs(title = "ALC ratio over time by delayed neurotoxicity status",
         subtitle = "LMM: log10(ALC) ~ DNT * ns(Day, df=4) + (1|Study_ID)") +
    theme_arial(22)
  print(p_alc_lmm_ratio)
  save_pdf(p_alc_lmm_ratio, "DNT_ALC_LMM_Ratio_over_Time.pdf", width = 10.5, height = 6.2)
  save_eps(p_alc_lmm_ratio, "DNT_ALC_LMM_Ratio_over_Time.eps", width = 10.5, height = 6.2)
} else {
  message("Skipping ALC LMM by DNT: insufficient DNT-positive patients or only one DNT level.")
}

# CAR exposure windows
auc0_30_const_strict <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  x <- as.numeric(x[ok]); y <- pmax(as.numeric(y[ok]), 0)
  keep <- x >= 0
  x <- x[keep]; y <- y[keep]
  if (!length(x)) return(NA_real_)
  if (any(abs(x - 0) < 1e-9)) {
    i0 <- which.min(abs(x - 0)); y[i0] <- 0
  } else {
    x <- c(0, x); y <- c(0, y)
  }
  ord <- order(x, y); x <- x[ord]; y <- y[ord]
  d <- data.frame(x = x, y = y)
  d <- d[!duplicated(d$x), , drop = FALSE]
  n_nz_0_30 <- sum(d$x > 0 & d$x <= 30 & d$y > 0)
  if (n_nz_0_30 < 2) return(NA_real_)
  y30 <- approx(d$x, d$y, xout = 30, rule = 2)$y
  in_rng <- d$x >= 0 & d$x <= 30
  xi <- c(0, d$x[in_rng], 30)
  yi <- c(0, d$y[in_rng], y30)
  dd <- data.frame(x = xi, y = yi)
  dd <- dd[order(dd$x, dd$y), , drop = FALSE]
  dd <- dd[!duplicated(dd$x), , drop = FALSE]
  if (nrow(dd) < 2) return(NA_real_)
  sum(diff(dd$x) * (head(dd$y, -1) + tail(dd$y, -1)) / 2)
}

auc_window <- function(x, y, start = 0, end = 100, min_nonzero = 2) {

  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)

  x <- as.numeric(x[ok])
  y <- pmax(as.numeric(y[ok]), 0)

  keep <- x >= 0
  x <- x[keep]
  y <- y[keep]
  if (!length(x)) return(NA_real_)

  ord <- order(x, y)
  x <- x[ord]
  y <- y[ord]

  d <- data.frame(x = x, y = y)
  d <- d[!duplicated(d$x), , drop = FALSE]

  if (nrow(d) < 2) return(NA_real_)

  n_nz <- sum(d$x >= start & d$x <= end & d$y > 0)
  if (n_nz < min_nonzero) return(NA_real_)

  y_start <- approx(d$x, d$y, xout = start, rule = 2)$y
  y_end   <- approx(d$x, d$y, xout = end,   rule = 2)$y

  in_rng <- d$x >= start & d$x <= end
  xi <- c(start, d$x[in_rng], end)
  yi <- c(y_start, d$y[in_rng], y_end)

  dd <- data.frame(x = xi, y = yi)
  dd <- dd[order(dd$x, dd$y), , drop = FALSE]
  dd <- dd[!duplicated(dd$x), , drop = FALSE]

  if (nrow(dd) < 2) return(NA_real_)

  sum(diff(dd$x) * (head(dd$y, -1) + tail(dd$y, -1)) / 2)
}

df4 <- readr::read_delim(path_df4, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
if (!"Day" %in% names(df4) && "Day_rel_infusion" %in% names(df4)) {
  df4 <- dplyr::rename(df4, Day = Day_rel_infusion)
}
if (!"CAR_abs" %in% names(df4)) {
  stop("Expansion file must contain CAR_abs. Found: ", paste(names(df4), collapse = ", "))
}

exp_all <- df4 %>%
  transmute(Study_ID, Day = safe_num(Day), CAR_abs = safe_num(CAR_abs)) %>%
  filter(!is.na(Study_ID), is.finite(Day), Day >= 0) %>%
  mutate(CAR_abs_imp = if_else(!is.finite(CAR_abs) | CAR_abs <= 0, LOD_ND, CAR_abs)) %>%
  left_join(dnt_meta, by = "Study_ID") %>%
  mutate(DNT_flag = factor(DNT_flag, levels = c("No", "Yes")))

aucs_car <- exp_all %>%
  group_by(Study_ID) %>%
  summarise(AUC_0_30 = auc0_30_const_strict(Day, CAR_abs),
            AUC_0_100 = auc_window(Day, CAR_abs, start = 0, end = 100, min_nonzero = 2),
            n_tp_any = sum(is.finite(Day)),
            n_tp_0_30 = sum(is.finite(Day) & Day <= 30),
            n_nonzero_0_30 = sum(is.finite(Day) & Day > 0 & Day <= 30 & is.finite(CAR_abs) & CAR_abs > 0),
            .groups = "drop") %>%
  left_join(dnt_meta, by = "Study_ID")

cat("\n================ CAR AUC availability ================\n")
auc_avail <- tibble(metric = c("Unique with any expansion rows", "Unique with valid AUC0-30", "Unique with valid AUC0-100"),
                    N = c(n_distinct(exp_all$Study_ID), sum(is.finite(aucs_car$AUC_0_30)), sum(is.finite(aucs_car$AUC_0_100))))
print(auc_avail)
save_tsv(auc_avail, "DNT_CAR_AUC_availability.tsv")

auc_dnt <- aucs_car %>% filter(is.finite(AUC_0_30), !is.na(DNT_flag))
cat("\n================ CAR AUC0-30 by DNT ================\n")
print(table(auc_dnt$DNT_flag))
wilcox_auc_0_30 <- wilcox.test(AUC_0_30 ~ DNT_flag, data = auc_dnt, exact = FALSE)
print(wilcox_auc_0_30)
save_csv(summ_cont(auc_dnt, "DNT_flag", "AUC_0_30"), "DNT_CAR_AUC0_30_summary.csv")

p_auc_dnt <- boxplot_log_by_dnt(
  df = auc_dnt,
  yvar = "AUC_0_30",
  ylab = "CAR AUC0-30 (cells*day/uL)",
  title = "CAR expansion AUC0-30 by delayed neurotoxicity status",
  subtitle = paste0("Wilcoxon p=", fmt_p(wilcox_auc_0_30$p.value)),
  filename_base = "CAR_AUC0_30_by_DNT"
)

auc100_dnt <- aucs_car %>% filter(is.finite(AUC_0_100), !is.na(DNT_flag))
if (nlevels(droplevels(auc100_dnt$DNT_flag)) == 2 && nrow(auc100_dnt) >= 3) {
  cat("\n================ CAR AUC0-100 by DNT ================\n")
  wilcox_auc_0_100 <- wilcox.test(AUC_0_100 ~ DNT_flag, data = auc100_dnt, exact = FALSE)
  print(wilcox_auc_0_100)
  save_csv(summ_cont(auc100_dnt, "DNT_flag", "AUC_0_100"), "DNT_CAR_AUC0_100_summary.csv")
  p_auc100_dnt <- boxplot_log_by_dnt(
    df = auc100_dnt,
    yvar = "AUC_0_100",
    ylab = "CAR AUC0-100 (cells*day/uL)",
    title = "CAR expansion AUC0-100 by delayed neurotoxicity status",
    subtitle = paste0("Wilcoxon p=", fmt_p(wilcox_auc_0_100$p.value)),
    filename_base = "CAR_AUC0_100_by_DNT"
  )
}

persist_car <- exp_all %>%
  filter(Day > 50) %>%
  group_by(Study_ID) %>%
  summarise(
    n_CAR_post50 = n(),
    n_unique_CAR_days_post50 = n_distinct(Day[is.finite(Day)]),
    Max_CAR_post50 = max(CAR_abs_imp, na.rm = TRUE),
    Last_CAR_day_post50 = suppressWarnings(max(Day[is.finite(CAR_abs_imp) & CAR_abs_imp > LOD_ND], na.rm = TRUE)),
    Any_CAR_detected_post90 = any(Day >= 90 & CAR_abs_imp > LOD_ND, na.rm = TRUE),
    Any_CAR_detected_post100 = any(Day >= 100 & CAR_abs_imp > LOD_ND, na.rm = TRUE),

    AUC_50_200 = auc_window(Day, CAR_abs, start = 50, end = 200, min_nonzero = 1),
    .groups = "drop"
  ) %>%
  mutate(
    Last_CAR_day_post50 = ifelse(is.infinite(Last_CAR_day_post50), NA_real_, Last_CAR_day_post50),
    CAR_post90_flag = factor(if_else(Any_CAR_detected_post90, "Yes", "No"), levels = c("No", "Yes")),
    CAR_post100_flag = factor(if_else(Any_CAR_detected_post100, "Yes", "No"), levels = c("No", "Yes"))
  ) %>%
  left_join(dnt_meta, by = "Study_ID")

cat("
================ CAR persistence AUC availability ================
")
print(tibble(
  metric = c("Patients with post-D50 CAR rows", "Patients with >=2 unique post-D50 CAR days", "Patients with valid AUC50-200"),
  N = c(nrow(persist_car), sum(persist_car$n_unique_CAR_days_post50 >= 2), sum(is.finite(persist_car$AUC_50_200)))
))
save_csv(persist_car, "DNT_CAR_persistence_patient_level_metrics.csv")

# Late CAR persistence
persist_max <- persist_car %>% filter(is.finite(Max_CAR_post50), !is.na(DNT_flag))
cat("\n================ Max CAR post-D50 by DNT ================\n")
print(table(persist_max$DNT_flag))
if (nlevels(droplevels(persist_max$DNT_flag)) == 2 && nrow(persist_max) >= 3) {
  wilcox_max_car_post50 <- wilcox.test(Max_CAR_post50 ~ DNT_flag, data = persist_max, exact = FALSE)
  print(wilcox_max_car_post50)
  save_csv(summ_cont(persist_max, "DNT_flag", "Max_CAR_post50"), "DNT_Max_CAR_post50_summary.csv")
  p_max_car_post50 <- boxplot_log_by_dnt(
    df = persist_max,
    yvar = "Max_CAR_post50",
    ylab = "Maximum CAR after day 50 (cells/uL)",
    title = "CAR persistence after day 50 by delayed neurotoxicity status",
    subtitle = paste0("Wilcoxon p=", fmt_p(wilcox_max_car_post50$p.value)),
    filename_base = "CAR_Max_postD50_by_DNT",
    hline = LOD_ND,
    hline_label = "LOD"
  )
}

persist90 <- persist_car %>% filter(!is.na(CAR_post90_flag), !is.na(DNT_flag))
if (nlevels(droplevels(persist90$CAR_post90_flag)) == 2 && nlevels(droplevels(persist90$DNT_flag)) == 2) {
  cat("\n================ Detectable CAR after D90 vs DNT ================\n")
  tab_post90 <- table(CAR_post90 = persist90$CAR_post90_flag, DNT = persist90$DNT_flag)
  print(tab_post90)
  fisher_post90 <- fisher.test(tab_post90)
  print(fisher_post90)
  save_csv(tibble(comparison = "Detectable CAR after D90 vs Any delayed neurotoxicity",
                  test = "Fisher exact",
                  odds_ratio = unname(fisher_post90$estimate),
                  conf_low = fisher_post90$conf.int[1],
                  conf_high = fisher_post90$conf.int[2],
                  p.value = fisher_post90$p.value),
           "DNT_CAR_detectable_postD90_Fisher.csv")
}

exp_plot_df <- exp_all %>% filter(Day >= 0, Day <= 200, is.finite(CAR_abs_imp), !is.na(DNT_flag))
p_car_dnt_spaghetti <- ggplot(exp_plot_df, aes(x = Day, y = CAR_abs_imp, group = Study_ID, color = DNT_flag)) +
  geom_line(alpha = 0.25, linewidth = 0.55) +
  geom_smooth(aes(group = DNT_flag, color = DNT_flag), method = "loess", se = TRUE, span = 0.55, linewidth = 1.25) +
  geom_hline(yintercept = LOD_ND, linetype = "dashed", color = "grey35") +
  scale_color_manual(values = pal_DNT, name = "Delayed neurotoxicity") +
  scale_y_log10("Absolute CAR (cells/uL)", labels = label_number(accuracy = 1, big.mark = ",")) +
  scale_x_continuous("Days post infusion", limits = c(0, 200)) +
  labs(title = "CAR expansion and persistence by delayed neurotoxicity status",
       subtitle = "Spaghetti plot with LOESS group trends") +
  theme_arial(22)
print(p_car_dnt_spaghetti)
save_pdf(p_car_dnt_spaghetti, "CAR_by_DNT_spaghetti_0_200.pdf", width = 11, height = 6.2)
save_eps(p_car_dnt_spaghetti, "CAR_by_DNT_spaghetti_0_200.eps", width = 11, height = 6.2)

binary_results_tbl <- bind_rows(
  iec_dnt_test_tbl %>% transmute(marker = "IEC-EC", test, odds_ratio, conf_low, conf_high, p.value),
  alc_gt3_test_tbl %>% transmute(marker = "Any ALC >3 K/uL D0-100", test, odds_ratio, conf_low, conf_high, p.value),
  tibble(marker = "Pre-event ALC >3 K/uL sensitivity", test = "Fisher exact",
         odds_ratio = unname(fisher_alc_pre$estimate), conf_low = fisher_alc_pre$conf.int[1],
         conf_high = fisher_alc_pre$conf.int[2], p.value = fisher_alc_pre$p.value)
)

if (exists("fisher_post90")) {
  binary_results_tbl <- bind_rows(binary_results_tbl,
                                  tibble(marker = "Detectable CAR after D90", test = "Fisher exact",
                                         odds_ratio = unname(fisher_post90$estimate),
                                         conf_low = fisher_post90$conf.int[1], conf_high = fisher_post90$conf.int[2],
                                         p.value = fisher_post90$p.value))
}

continuous_results_tbl <- bind_rows(
  tibble(marker = "Peak ALC D0-100, K/uL", test = "Wilcoxon rank-sum", p.value = wilcox_peak_alc$p.value),
  tibble(marker = "CAR AUC0-30, cells*day/uL", test = "Wilcoxon rank-sum", p.value = wilcox_auc_0_30$p.value)
)

if (exists("wilcox_auc_0_100")) {
  continuous_results_tbl <- bind_rows(continuous_results_tbl,
                                      tibble(marker = "CAR AUC0-100, cells*day/uL", test = "Wilcoxon rank-sum", p.value = wilcox_auc_0_100$p.value))
}
if (exists("wilcox_max_car_post50")) {
  continuous_results_tbl <- bind_rows(continuous_results_tbl,
                                      tibble(marker = "Max CAR after D50, cells/uL", test = "Wilcoxon rank-sum", p.value = wilcox_max_car_post50$p.value))
}

save_csv(binary_results_tbl, "DNT_binary_marker_tests.csv")
save_csv(continuous_results_tbl, "DNT_continuous_marker_tests.csv")

cat("\n================ Exploratory binary marker tests ================\n")
print(binary_results_tbl)
cat("\n================ Exploratory continuous marker tests ================\n")
print(continuous_results_tbl)

summary_lines <- c(
  "===== Delayed Neurotoxicity Analysis Summary =====",
  "",
  sprintf("Delayed neurotoxicity prevalence: %.1f%% (%d/%d)",
          100 * sum(dnt_meta$DNT_flag == "Yes", na.rm = TRUE) / nrow(dnt_meta),
          sum(dnt_meta$DNT_flag == "Yes", na.rm = TRUE), nrow(dnt_meta)),
  "",
  "Primary comparisons:",
  sprintf("1) IEC-EC vs delayed neurotoxicity: Fisher exact p=%s; OR=%.3f (95%% CI %.3f-%.3f)",
          fmt_p(fisher_iec_dnt$p.value), unname(fisher_iec_dnt$estimate),
          fisher_iec_dnt$conf.int[1], fisher_iec_dnt$conf.int[2]),
  sprintf("2) ALC >3 K/uL D0-100 vs delayed neurotoxicity: Fisher exact p=%s; OR=%.3f (95%% CI %.3f-%.3f)",
          fmt_p(fisher_alc_gt3$p.value), unname(fisher_alc_gt3$estimate),
          fisher_alc_gt3$conf.int[1], fisher_alc_gt3$conf.int[2]),
  "",
  "Secondary exploratory analyses:",
  sprintf("Peak ALC D0-100 by DNT: Wilcoxon p=%s", fmt_p(wilcox_peak_alc$p.value)),
  sprintf("CAR AUC0-30 by DNT: Wilcoxon p=%s", fmt_p(wilcox_auc_0_30$p.value)),
  sprintf("Max CAR after D50 by DNT: Wilcoxon p=%s",
          ifelse(exists("wilcox_max_car_post50"), fmt_p(wilcox_max_car_post50$p.value), "NA")),
  "",
  "Note: ALC is analyzed in K/uL. Therefore, ALC >3 K/uL corresponds to >3000/uL."
)

summary_file <- file.path(stats_dir, "DNT_analysis_summary.txt")
writeLines(summary_lines, summary_file)
message("Saved: ", summary_file)
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
