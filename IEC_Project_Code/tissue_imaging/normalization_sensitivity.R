#!/usr/bin/env Rscript
# Tissue normalization sensitivity
# Compare CD3 and total-cell normalization with tissue models and patient-cluster bootstrap summaries.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
})

args <- commandArgs(trailingOnly = TRUE)

default_data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
default_input_name <- "tissue_cell_counts.tsv"

input_arg <- if (length(args) >= 1 && nzchar(args[[1]])) {
  path.expand(args[[1]])
} else {
  default_data_dir
}

if (dir.exists(input_arg)) {
  data_dir <- normalizePath(input_arg, mustWork = TRUE)
  input_path <- file.path(data_dir, default_input_name)
} else {
  input_path <- normalizePath(input_arg, mustWork = TRUE)
  data_dir <- dirname(input_path)
}

output_dir <- if (length(args) >= 2 && nzchar(args[[2]])) {
  path.expand(args[[2]])
} else {
  file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "tissue_imaging", "normalization_sensitivity")
}

figure_dir <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
log_dir    <- file.path(output_dir, "logs")
invisible(lapply(
  c(output_dir, figure_dir, table_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# Eligibility and bootstrap parameters
min_tissue_n <- 3L

min_patient_tissues <- 2L

bootstrap_B <- 2000L
bootstrap_seed <- 20260726L

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}

BASE_FAMILY <- "Helvetica"

ggplot2::theme_set(
  ggplot2::theme_classic(base_family = BASE_FAMILY, base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 11, color = "grey20"),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey94", color = NA)
    )
)

pal_tissue <- c(
  "Colon" = "#8DA0CB",
  "Duodenum" = "#66C2A5",
  "Terminal Ileum" = "#FC8D62",
  "Other" = "#BDBDBD"
)

pal_status <- c(
  "Camelid Negative" = "#4C6A87",
  "Camelid Positive" = "#B22222",
  "Untreated Controls" = "#9E9E9E"
)

pal_method <- c(
  "CD3-normalized" = "#4C6A87",
  "Total-cell-normalized" = "#B22222"
)

status_label <- function(cohort) {
  z <- as.character(cohort)
  dplyr::case_when(
    z %in% c("IEC_EC", "IEC-EC") ~ "IEC-EC",
    z %in% c("CONTROL", "Control") ~ "Untreated Controls",
    z %in% c(
      "NOT_IEC", "NOT_IEC_EC", "Not_IEC_EC", "NOT_IEC-EC",
      "Indeterminate", "INDETERMINATE"
    ) ~ "Not IEC-EC",
    TRUE ~ z
  )
}

status_display <- function(x) {
  dplyr::case_when(
    x == "IEC-EC" ~ "Camelid Positive",
    x == "Not IEC-EC" ~ "Camelid Negative",
    x == "Untreated Controls" ~ "Untreated Controls",
    TRUE ~ as.character(x)
  )
}

tissue_fine <- function(tt) {
  z <- tolower(gsub("\\s+", "_", as.character(tt)))
  dplyr::case_when(
    grepl("^duodenum$", z) ~ "Duodenum",
    grepl("terminal_ileum", z) ~ "Terminal Ileum",
    grepl("ileocecal", z) ~ "Ileocecal Valve",
    grepl("^colon_left$", z) ~ "Colon (left)",
    grepl("^colon_right$", z) ~ "Colon (right)",
    grepl("^colon_random$", z) ~ "Colon (random)",
    grepl("^colon$", z) ~ "Colon",
    grepl("^rectum$", z) ~ "Rectum",
    grepl("^esophagus$", z) ~ "Esophagus",
    grepl("^stomach_antrum$", z) ~ "Stomach (antrum)",
    grepl("^stomach_polyp$", z) ~ "Stomach (polyp)",
    grepl("^stomach$", z) ~ "Stomach",
    TRUE ~ stringr::str_to_title(gsub("_", " ", as.character(tt)))
  )
}

tissue_primary_from_fine <- function(f) {
  dplyr::case_when(
    f == "Duodenum" ~ "Duodenum",
    f == "Terminal Ileum" ~ "Terminal Ileum",
    f %in% c(
      "Colon", "Colon (left)", "Colon (right)", "Colon (random)",
      "Ileocecal Valve"
    ) ~ "Colon",
    TRUE ~ "Other"
  )
}

frac_safe <- function(num, den) {
  num <- as.numeric(num)
  den <- as.numeric(den)
  ifelse(
    is.finite(num) & is.finite(den) & den > 0 & num >= 0 & num <= den,
    num / den,
    NA_real_
  )
}

log_prop_corrected <- function(num, den) {
  num <- as.numeric(num)
  den <- as.numeric(den)
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & den > 0 & num >= 0 & num <= den
  out[ok] <- log((num[ok] + 0.5) / (den[ok] + 1))
  out
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  sx <- stats::sd(x[ok])
  sy <- stats::sd(y[ok])
  if (
    sum(ok) < 3 ||
      !is.finite(sx) || !is.finite(sy) ||
      sx == 0 || sy == 0
  ) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  test <- suppressWarnings(
    stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  )
  c(
    rho = unname(test$estimate),
    p = test$p.value,
    n = sum(ok)
  )
}

safe_paired_wilcox <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) {
    return(c(n = sum(ok), statistic = NA_real_, p = NA_real_))
  }
  test <- tryCatch(
    stats::wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(c(n = sum(ok), statistic = NA_real_, p = NA_real_))
  }
  c(
    n = sum(ok),
    statistic = unname(test$statistic),
    p = test$p.value
  )
}

cluster_resample <- function(dat, id_col = "ID") {
  id_values <- unique(as.character(dat[[id_col]]))
  id_values <- id_values[!is.na(id_values) & nzchar(id_values)]
  if (length(id_values) == 0) return(dat[0, , drop = FALSE])
  row_index <- split(seq_len(nrow(dat)), as.character(dat[[id_col]]))
  sampled_ids <- sample(id_values, size = length(id_values), replace = TRUE)
  dat[unlist(row_index[sampled_ids], use.names = FALSE), , drop = FALSE]
}

partial_r2_one <- function(dat, outcome_col, tissue_col) {
  d <- dat %>%
    dplyr::transmute(
      .y = .data[[outcome_col]],
      .tissue = as.factor(.data[[tissue_col]]),
      .status = as.factor(Status_main)
    ) %>%
    dplyr::filter(is.finite(.y), !is.na(.tissue), !is.na(.status)) %>%
    droplevels()

  if (nrow(d) < 6 || nlevels(d$.tissue) < 2) {
    return(c(partial_r2 = NA_real_, p = NA_real_, n = nrow(d)))
  }

  reduced_formula <- if (nlevels(d$.status) >= 2) {
    stats::as.formula(".y ~ .status")
  } else {
    stats::as.formula(".y ~ 1")
  }
  full_formula <- if (nlevels(d$.status) >= 2) {
    stats::as.formula(".y ~ .status + .tissue")
  } else {
    stats::as.formula(".y ~ .tissue")
  }

  reduced <- tryCatch(stats::lm(reduced_formula, data = d), error = function(e) NULL)
  full <- tryCatch(stats::lm(full_formula, data = d), error = function(e) NULL)
  if (is.null(reduced) || is.null(full)) {
    return(c(partial_r2 = NA_real_, p = NA_real_, n = nrow(d)))
  }

  sse_reduced <- sum(stats::residuals(reduced)^2)
  sse_full <- sum(stats::residuals(full)^2)
  partial_r2 <- if (is.finite(sse_reduced) && sse_reduced > 0) {
    max(0, min(1, (sse_reduced - sse_full) / sse_reduced))
  } else {
    NA_real_
  }

  comparison <- tryCatch(
    as.data.frame(stats::anova(reduced, full)),
    error = function(e) NULL
  )
  p <- if (!is.null(comparison) && "Pr(>F)" %in% names(comparison)) {
    comparison[["Pr(>F)"]][2]
  } else {
    NA_real_
  }

  c(partial_r2 = partial_r2, p = p, n = nrow(d))
}

prepare_tissue_data <- function(dat, tissue_col, min_n = min_tissue_n) {
  keep_levels <- dat %>%
    dplyr::filter(!is.na(.data[[tissue_col]])) %>%
    dplyr::count(.data[[tissue_col]], name = "n") %>%
    dplyr::filter(n >= min_n) %>%
    dplyr::pull(1)

  dat %>%
    dplyr::filter(.data[[tissue_col]] %in% keep_levels) %>%
    dplyr::mutate(
      "{tissue_col}" := droplevels(as.factor(.data[[tissue_col]]))
    )
}

fit_beta_binomial_tissue <- function(
    dat,
    success_col,
    total_col,
    tissue_col,
    outcome_label,
    scope_label
) {
  d <- prepare_tissue_data(dat, tissue_col = tissue_col) %>%
    dplyr::transmute(
      ID = as.factor(ID),
      Status_main = as.factor(Status_main),
      Tissue = as.factor(.data[[tissue_col]]),
      success_raw = as.numeric(.data[[success_col]]),
      total_raw = as.numeric(.data[[total_col]])
    ) %>%
    dplyr::filter(
      is.finite(success_raw),
      is.finite(total_raw),
      total_raw > 0,
      success_raw >= 0,
      success_raw <= total_raw,
      !is.na(ID),
      !is.na(Status_main),
      !is.na(Tissue)
    ) %>%
    dplyr::mutate(
      total = as.integer(round(total_raw)),
      success = pmax(0L, pmin(as.integer(round(success_raw)), total)),
      failure = total - success
    ) %>%
    dplyr::filter(total > 0) %>%
    droplevels()

  empty_result <- tibble::tibble(
    scope = scope_label,
    tissue_definition = tissue_col,
    outcome = outcome_label,
    model = NA_character_,
    n = nrow(d),
    n_patients = dplyr::n_distinct(d$ID),
    n_tissues = dplyr::n_distinct(d$Tissue),
    likelihood_ratio_chisq = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    AIC_reduced = NA_real_,
    AIC_full = NA_real_
  )

  if (nrow(d) < 6 || dplyr::n_distinct(d$Tissue) < 2) {
    return(empty_result)
  }

  status_term <- if (dplyr::n_distinct(d$Status_main) >= 2) {
    "Status_main"
  } else {
    "1"
  }
  use_random_id <- any(table(d$ID) > 1)
  random_term <- if (use_random_id) " + (1 | ID)" else ""

  reduced_formula <- stats::as.formula(
    paste0("cbind(success, failure) ~ ", status_term, random_term)
  )
  full_formula <- stats::as.formula(
    paste0("cbind(success, failure) ~ ", status_term, " + Tissue", random_term)
  )

  fit_pair <- function(f0, f1, model_name) {
    reduced <- suppressWarnings(tryCatch(
      glmmTMB::glmmTMB(
        f0,
        family = glmmTMB::betabinomial(link = "logit"),
        data = d
      ),
      error = function(e) NULL
    ))
    full <- suppressWarnings(tryCatch(
      glmmTMB::glmmTMB(
        f1,
        family = glmmTMB::betabinomial(link = "logit"),
        data = d
      ),
      error = function(e) NULL
    ))
    if (is.null(reduced) || is.null(full)) return(NULL)
    list(reduced = reduced, full = full, model_name = model_name)
  }

  fits <- fit_pair(
    reduced_formula,
    full_formula,
    if (use_random_id) "beta-binomial GLMM" else "beta-binomial model"
  )

  if (is.null(fits) && use_random_id) {
    reduced_formula_fixed <- stats::as.formula(
      paste0("cbind(success, failure) ~ ", status_term)
    )
    full_formula_fixed <- stats::as.formula(
      paste0("cbind(success, failure) ~ ", status_term, " + Tissue")
    )
    fits <- fit_pair(
      reduced_formula_fixed,
      full_formula_fixed,
      "beta-binomial model (random-intercept fallback)"
    )
  }

  if (is.null(fits)) return(empty_result)

  comparison <- tryCatch(
    as.data.frame(stats::anova(fits$reduced, fits$full)),
    error = function(e) NULL
  )
  if (is.null(comparison) || nrow(comparison) < 2) {
    return(empty_result %>% dplyr::mutate(model = fits$model_name))
  }

  p_col <- grep("^Pr\\(", names(comparison), value = TRUE)
  chisq_col <- base::intersect(c("Chisq", "LRT"), names(comparison))
  df_col <- base::intersect(c("Chi Df", "Df"), names(comparison))

  empty_result %>%
    dplyr::mutate(
      model = fits$model_name,
      likelihood_ratio_chisq = if (length(chisq_col)) {
        as.numeric(comparison[[chisq_col[[1]]]][2])
      } else {
        NA_real_
      },
      df = if (length(df_col)) {
        as.numeric(comparison[[df_col[[1]]]][2])
      } else {
        NA_real_
      },
      p_value = if (length(p_col)) {
        as.numeric(comparison[[p_col[[1]]]][2])
      } else {
        NA_real_
      },
      AIC_reduced = stats::AIC(fits$reduced),
      AIC_full = stats::AIC(fits$full)
    )
}

save_plot_all <- function(plot, file_stub, width = 7, height = 5, dpi = 300) {
  png_file <- file.path(figure_dir, paste0(file_stub, ".png"))
  pdf_file <- file.path(figure_dir, paste0(file_stub, ".pdf"))
  eps_file <- file.path(figure_dir, paste0(file_stub, ".eps"))

  ggplot2::ggsave(
    png_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = BASE_FAMILY,
    useDingbats = FALSE,
    version = "1.4",
    colormodel = "srgb"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::postscript(
    eps_file,
    width = width,
    height = height,
    family = BASE_FAMILY,
    onefile = FALSE,
    horizontal = FALSE,
    paper = "special",
    colormodel = "srgb"
  )
  print(plot)
  grDevices::dev.off()

  invisible(c(png_file, pdf_file, eps_file))
}

show_and_save <- function(plot, file_stub, width = 7, height = 5) {
  print(plot)
  save_plot_all(plot, file_stub, width = width, height = height)
  invisible(plot)
}

if (!file.exists(input_path)) {
  stop("COMET input file not found: ", input_path)
}

df_raw <- readr::read_delim(
  input_path,
  delim = "\t",
  trim_ws = TRUE,
  show_col_types = FALSE
)

required_columns <- c(
  "Patient_ID", "Cohort", "Tissue_Type",
  "total_cells", "CD3", "CD3_Camelid"
)
missing_columns <- base::setdiff(required_columns, names(df_raw))
if (length(missing_columns) > 0) {
  stop(
    "COMET input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

df <- df_raw %>%
  dplyr::transmute(
    ID = as.character(Patient_ID),
    Cohort = as.character(Cohort),
    Tissue_Type = as.character(Tissue_Type),
    total_cells = suppressWarnings(as.numeric(total_cells)),
    CD3 = suppressWarnings(as.numeric(CD3)),
    CD3_Camelid = suppressWarnings(as.numeric(CD3_Camelid))
  ) %>%
  dplyr::mutate(
    Status_main = status_label(Cohort),
    Status_plot = factor(
      status_display(Status_main),
      levels = c(
        "Untreated Controls",
        "Camelid Negative",
        "Camelid Positive"
      )
    ),
    Tissue_fine = tissue_fine(Tissue_Type),
    Tissue_primary = factor(
      tissue_primary_from_fine(Tissue_fine),
      levels = c("Duodenum", "Terminal Ileum", "Colon", "Other")
    ),
    valid_total = is.finite(total_cells) & total_cells > 0,
    valid_cd3 = is.finite(CD3) & CD3 > 0,
    valid_camelid = is.finite(CD3_Camelid) & CD3_Camelid >= 0,
    cd3_le_total = valid_total & is.finite(CD3) & CD3 >= 0 & CD3 <= total_cells,
    camelid_le_cd3 = valid_cd3 & valid_camelid & CD3_Camelid <= CD3,
    valid_for_comparison =
      valid_total & valid_cd3 & valid_camelid &
      cd3_le_total & camelid_le_cd3 &
      !is.na(ID) & nzchar(ID) &
      !is.na(Status_plot) &
      !is.na(Tissue_primary),
    prop_cd3_total = frac_safe(CD3, total_cells),
    prop_camelid_cd3 = frac_safe(CD3_Camelid, CD3),
    prop_camelid_total = frac_safe(CD3_Camelid, total_cells),
    pct_cd3_total = 100 * prop_cd3_total,
    pct_camelid_cd3 = 100 * prop_camelid_cd3,
    pct_camelid_total = 100 * prop_camelid_total,
    log_prop_cd3_total = log_prop_corrected(CD3, total_cells),
    log_prop_camelid_cd3 = log_prop_corrected(CD3_Camelid, CD3),
    log_prop_camelid_total = log_prop_corrected(CD3_Camelid, total_cells)
  )

count_audit <- df %>%
  dplyr::summarise(
    n_input_rows = dplyr::n(),
    n_valid_for_comparison = sum(valid_for_comparison, na.rm = TRUE),
    n_missing_or_nonpositive_total = sum(!valid_total, na.rm = TRUE),
    n_missing_or_nonpositive_cd3 = sum(!valid_cd3, na.rm = TRUE),
    n_invalid_camelid = sum(!valid_camelid, na.rm = TRUE),
    n_cd3_gt_total = sum(valid_total & is.finite(CD3) & CD3 > total_cells, na.rm = TRUE),
    n_camelid_gt_cd3 = sum(valid_cd3 & valid_camelid & CD3_Camelid > CD3, na.rm = TRUE)
  )

analysis_df <- df %>%
  dplyr::filter(valid_for_comparison) %>%
  droplevels()

if (nrow(analysis_df) < 6) {
  stop("Fewer than six valid COMET biopsies remain after count validation.")
}

identity_audit <- analysis_df %>%
  dplyr::summarise(
    max_absolute_identity_error = max(
      abs(
        prop_camelid_total -
          prop_camelid_cd3 * prop_cd3_total
      ),
      na.rm = TRUE
    ),
    median_absolute_identity_error = median(
      abs(
        prop_camelid_total -
          prop_camelid_cd3 * prop_cd3_total
      ),
      na.rm = TRUE
    )
  )

tissue_summary <- analysis_df %>%
  dplyr::group_by(Tissue_primary) %>%
  dplyr::summarise(
    n_biopsies = dplyr::n(),
    n_patients = dplyr::n_distinct(ID),
    median_pct_CD3_of_total = median(pct_cd3_total, na.rm = TRUE),
    IQR_pct_CD3_of_total = stats::IQR(pct_cd3_total, na.rm = TRUE),
    median_pct_Camelid_of_CD3 = median(pct_camelid_cd3, na.rm = TRUE),
    IQR_pct_Camelid_of_CD3 = stats::IQR(pct_camelid_cd3, na.rm = TRUE),
    median_pct_Camelid_of_total = median(pct_camelid_total, na.rm = TRUE),
    IQR_pct_Camelid_of_total = stats::IQR(pct_camelid_total, na.rm = TRUE),
    .groups = "drop"
  )

formal_primary <- prepare_tissue_data(
  analysis_df,
  tissue_col = "Tissue_primary"
)

formal_fine <- prepare_tissue_data(
  analysis_df,
  tissue_col = "Tissue_fine"
)

primary_counts <- formal_primary %>%
  dplyr::count(Tissue_primary, name = "n")
primary_labels <- setNames(
  paste0(as.character(primary_counts$Tissue_primary), "\nn = ", primary_counts$n),
  as.character(primary_counts$Tissue_primary)
)

scope_data <- list(
  "All cohorts (adjusted for status)" = analysis_df,
  "IEC-EC only" = analysis_df %>% dplyr::filter(Status_main == "IEC-EC")
)

# Normalization definitions
model_outcomes <- tribble(
  ~success_col, ~total_col, ~outcome,
  "CD3", "total_cells", "CD3 / total cells",
  "CD3_Camelid", "CD3", "Camelid+ / CD3",
  "CD3_Camelid", "total_cells", "Camelid+ / total cells"
)

bb_results <- list()
for (scope_name in names(scope_data)) {
  scope_df <- scope_data[[scope_name]]
  for (tissue_name in c("Tissue_primary", "Tissue_fine")) {
    for (i in seq_len(nrow(model_outcomes))) {
      bb_results[[length(bb_results) + 1L]] <-
        fit_beta_binomial_tissue(
          dat = scope_df,
          success_col = model_outcomes$success_col[[i]],
          total_col = model_outcomes$total_col[[i]],
          tissue_col = tissue_name,
          outcome_label = model_outcomes$outcome[[i]],
          scope_label = scope_name
        )
    }
  }
}
bb_results <- dplyr::bind_rows(bb_results)

correlation_results <- list()
correlation_bootstrap <- list()

set.seed(bootstrap_seed)
for (scope_name in names(scope_data)) {
  scope_df <- scope_data[[scope_name]]

  cor_cd3 <- safe_spearman(
    scope_df$prop_camelid_cd3,
    scope_df$prop_cd3_total
  )
  cor_total <- safe_spearman(
    scope_df$prop_camelid_total,
    scope_df$prop_cd3_total
  )

  correlation_results[[length(correlation_results) + 1L]] <- tibble::tibble(
    scope = scope_name,
    normalization = c("CD3-normalized", "Total-cell-normalized"),
    outcome = c("Camelid+ / CD3", "Camelid+ / total cells"),
    n = c(cor_cd3[["n"]], cor_total[["n"]]),
    spearman_rho = c(cor_cd3[["rho"]], cor_total[["rho"]]),
    p_value = c(cor_cd3[["p"]], cor_total[["p"]])
  )

  boot_delta <- rep(NA_real_, bootstrap_B)
  boot_rho_cd3 <- rep(NA_real_, bootstrap_B)
  boot_rho_total <- rep(NA_real_, bootstrap_B)

  if (dplyr::n_distinct(scope_df$ID) >= 2) {
    for (b in seq_len(bootstrap_B)) {
      boot_df <- cluster_resample(scope_df)
      rho_cd3 <- safe_spearman(
        boot_df$prop_camelid_cd3,
        boot_df$prop_cd3_total
      )[["rho"]]
      rho_total <- safe_spearman(
        boot_df$prop_camelid_total,
        boot_df$prop_cd3_total
      )[["rho"]]
      boot_rho_cd3[[b]] <- rho_cd3
      boot_rho_total[[b]] <- rho_total
      boot_delta[[b]] <- abs(rho_total) - abs(rho_cd3)
    }
  }

  observed_delta <- abs(cor_total[["rho"]]) - abs(cor_cd3[["rho"]])
  correlation_bootstrap[[length(correlation_bootstrap) + 1L]] <-
    tibble::tibble(
      scope = scope_name,
      statistic = "|rho(total-normalized)| - |rho(CD3-normalized)|",
      estimate = observed_delta,
      bootstrap_median = safe_median(boot_delta),
      CI_low = safe_quantile(boot_delta, 0.025),
      CI_high = safe_quantile(boot_delta, 0.975),
      bootstrap_B = bootstrap_B,
      interpretation = paste(
        "Positive values mean total-cell normalization remains more strongly",
        "coupled to lymphocyte density."
      )
    )
}

correlation_results <- dplyr::bind_rows(correlation_results)
correlation_bootstrap <- dplyr::bind_rows(correlation_bootstrap)

r2_results <- list()
r2_bootstrap_results <- list()

set.seed(bootstrap_seed + 1L)
for (scope_name in names(scope_data)) {
  scope_df <- scope_data[[scope_name]]

  for (tissue_name in c("Tissue_primary", "Tissue_fine")) {
    d_scope <- prepare_tissue_data(scope_df, tissue_col = tissue_name)

    r_cd3 <- partial_r2_one(
      d_scope,
      outcome_col = "log_prop_camelid_cd3",
      tissue_col = tissue_name
    )
    r_total <- partial_r2_one(
      d_scope,
      outcome_col = "log_prop_camelid_total",
      tissue_col = tissue_name
    )
    r_lymph <- partial_r2_one(
      d_scope,
      outcome_col = "log_prop_cd3_total",
      tissue_col = tissue_name
    )

    r2_results[[length(r2_results) + 1L]] <- tibble::tibble(
      scope = scope_name,
      tissue_definition = tissue_name,
      normalization = c(
        "Lymphocyte density",
        "CD3-normalized",
        "Total-cell-normalized"
      ),
      outcome = c(
        "CD3 / total cells",
        "Camelid+ / CD3",
        "Camelid+ / total cells"
      ),
      n = c(r_lymph[["n"]], r_cd3[["n"]], r_total[["n"]]),
      tissue_partial_R2 = c(
        r_lymph[["partial_r2"]],
        r_cd3[["partial_r2"]],
        r_total[["partial_r2"]]
      ),
      p_value = c(r_lymph[["p"]], r_cd3[["p"]], r_total[["p"]])
    )

    boot_r_cd3 <- rep(NA_real_, bootstrap_B)
    boot_r_total <- rep(NA_real_, bootstrap_B)
    boot_delta <- rep(NA_real_, bootstrap_B)

    if (dplyr::n_distinct(d_scope$ID) >= 2) {
      for (b in seq_len(bootstrap_B)) {
        boot_df <- cluster_resample(d_scope)
        b_cd3 <- partial_r2_one(
          boot_df,
          outcome_col = "log_prop_camelid_cd3",
          tissue_col = tissue_name
        )[["partial_r2"]]
        b_total <- partial_r2_one(
          boot_df,
          outcome_col = "log_prop_camelid_total",
          tissue_col = tissue_name
        )[["partial_r2"]]
        boot_r_cd3[[b]] <- b_cd3
        boot_r_total[[b]] <- b_total
        boot_delta[[b]] <- b_total - b_cd3
      }
    }

    r2_bootstrap_results[[length(r2_bootstrap_results) + 1L]] <-
      dplyr::bind_rows(
        tibble::tibble(
          scope = scope_name,
          tissue_definition = tissue_name,
          normalization = "CD3-normalized",
          estimate = r_cd3[["partial_r2"]],
          CI_low = safe_quantile(boot_r_cd3, 0.025),
          CI_high = safe_quantile(boot_r_cd3, 0.975),
          bootstrap_B = bootstrap_B
        ),
        tibble::tibble(
          scope = scope_name,
          tissue_definition = tissue_name,
          normalization = "Total-cell-normalized",
          estimate = r_total[["partial_r2"]],
          CI_low = safe_quantile(boot_r_total, 0.025),
          CI_high = safe_quantile(boot_r_total, 0.975),
          bootstrap_B = bootstrap_B
        ),
        tibble::tibble(
          scope = scope_name,
          tissue_definition = tissue_name,
          normalization =
            "Difference: total-cell minus CD3-normalized",
          estimate =
            r_total[["partial_r2"]] - r_cd3[["partial_r2"]],
          CI_low = safe_quantile(boot_delta, 0.025),
          CI_high = safe_quantile(boot_delta, 0.975),
          bootstrap_B = bootstrap_B
        )
      )
  }
}

r2_results <- dplyr::bind_rows(r2_results)
r2_bootstrap_results <- dplyr::bind_rows(r2_bootstrap_results)

patient_tissue_values <- analysis_df %>%
  dplyr::group_by(ID, Status_main, Status_plot, Tissue_primary) %>%
  dplyr::summarise(
    n_biopsies = dplyr::n(),
    median_log_camelid_cd3 = median(log_prop_camelid_cd3, na.rm = TRUE),
    median_log_camelid_total = median(log_prop_camelid_total, na.rm = TRUE),
    .groups = "drop"
  )

patient_spread <- patient_tissue_values %>%
  dplyr::group_by(ID, Status_main, Status_plot) %>%
  dplyr::summarise(
    n_tissues = dplyr::n_distinct(Tissue_primary),
    n_biopsies = sum(n_biopsies),
    SD_log_CD3_normalized = if (n_tissues >= 2) {
      stats::sd(median_log_camelid_cd3, na.rm = TRUE)
    } else {
      NA_real_
    },
    SD_log_total_normalized = if (n_tissues >= 2) {
      stats::sd(median_log_camelid_total, na.rm = TRUE)
    } else {
      NA_real_
    },
    IQR_log_CD3_normalized = if (n_tissues >= 2) {
      stats::IQR(median_log_camelid_cd3, na.rm = TRUE)
    } else {
      NA_real_
    },
    IQR_log_total_normalized = if (n_tissues >= 2) {
      stats::IQR(median_log_camelid_total, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  dplyr::filter(n_tissues >= min_patient_tissues)

spread_test_results <- list()
for (scope_name in names(scope_data)) {
  spread_scope <- if (scope_name == "IEC-EC only") {
    patient_spread %>% dplyr::filter(Status_main == "IEC-EC")
  } else {
    patient_spread
  }

  test_sd <- safe_paired_wilcox(
    spread_scope$SD_log_total_normalized,
    spread_scope$SD_log_CD3_normalized
  )
  test_iqr <- safe_paired_wilcox(
    spread_scope$IQR_log_total_normalized,
    spread_scope$IQR_log_CD3_normalized
  )

  spread_test_results[[length(spread_test_results) + 1L]] <-
    dplyr::bind_rows(
      tibble::tibble(
        scope = scope_name,
        variability_metric = "SD of log proportion across tissue medians",
        n_patients = test_sd[["n"]],
        median_total_cell_normalized = safe_median(
          spread_scope$SD_log_total_normalized
        ),
        median_CD3_normalized = safe_median(
          spread_scope$SD_log_CD3_normalized
        ),
        median_paired_difference_total_minus_CD3 = safe_median(
          spread_scope$SD_log_total_normalized -
            spread_scope$SD_log_CD3_normalized
        ),
        wilcoxon_statistic = test_sd[["statistic"]],
        p_value = test_sd[["p"]]
      ),
      tibble::tibble(
        scope = scope_name,
        variability_metric = "IQR of log proportion across tissue medians",
        n_patients = test_iqr[["n"]],
        median_total_cell_normalized = safe_median(
          spread_scope$IQR_log_total_normalized
        ),
        median_CD3_normalized = safe_median(
          spread_scope$IQR_log_CD3_normalized
        ),
        median_paired_difference_total_minus_CD3 = safe_median(
          spread_scope$IQR_log_total_normalized -
            spread_scope$IQR_log_CD3_normalized
        ),
        wilcoxon_statistic = test_iqr[["statistic"]],
        p_value = test_iqr[["p"]]
      )
    )
}
spread_test_results <- dplyr::bind_rows(spread_test_results)

p1 <- ggplot2::ggplot(
  formal_primary,
  ggplot2::aes(x = Tissue_primary, y = pct_cd3_total)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = Tissue_primary),
    width = 0.68,
    alpha = 0.28,
    outlier.shape = NA,
    color = "grey35"
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Status_plot),
    width = 0.14,
    height = 0,
    size = 2.3,
    alpha = 0.88
  ) +
  ggplot2::scale_fill_manual(values = pal_tissue, guide = "none") +
  ggplot2::scale_color_manual(values = pal_status, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = primary_labels) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(suffix = "%", accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Lymphocyte density differs across gastrointestinal tissues",
    subtitle = "CD3+ cells as a percentage of all counted cells",
    x = NULL,
    y = "CD3+ cells / total cells",
    color = "COMET status"
  ) +
  ggplot2::theme(legend.position = "bottom")

abundance_long <- formal_primary %>%
  dplyr::select(
    ID, Status_plot, Tissue_primary,
    pct_camelid_cd3, pct_camelid_total
  ) %>%
  tidyr::pivot_longer(
    cols = c(pct_camelid_cd3, pct_camelid_total),
    names_to = "normalization",
    values_to = "pct_camelid"
  ) %>%
  dplyr::mutate(
    normalization = dplyr::recode(
      normalization,
      pct_camelid_cd3 = "CD3-normalized",
      pct_camelid_total = "Total-cell-normalized"
    ),
    normalization = factor(
      normalization,
      levels = c("CD3-normalized", "Total-cell-normalized")
    )
  )

p2 <- ggplot2::ggplot(
  abundance_long,
  ggplot2::aes(x = Tissue_primary, y = pct_camelid)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = Tissue_primary),
    width = 0.68,
    alpha = 0.28,
    outlier.shape = NA,
    color = "grey35"
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Status_plot),
    width = 0.14,
    height = 0,
    size = 2.1,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(normalization),
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::scale_fill_manual(values = pal_tissue, guide = "none") +
  ggplot2::scale_color_manual(values = pal_status, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = primary_labels) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(suffix = "%", accuracy = 0.1)
  ) +
  ggplot2::labs(
    title = "Camelid+ abundance by tissue and denominator",
    subtitle = paste(
      "Facets use independent y-axis ranges; formal comparisons use",
      "count-based and log-proportion models."
    ),
    x = NULL,
    y = "Camelid+ cells",
    color = "COMET status"
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )

coupling_long <- formal_primary %>%
  dplyr::select(
    ID, Status_plot, Tissue_primary, pct_cd3_total,
    pct_camelid_cd3, pct_camelid_total
  ) %>%
  tidyr::pivot_longer(
    cols = c(pct_camelid_cd3, pct_camelid_total),
    names_to = "normalization",
    values_to = "pct_camelid"
  ) %>%
  dplyr::mutate(
    normalization = dplyr::recode(
      normalization,
      pct_camelid_cd3 = "CD3-normalized",
      pct_camelid_total = "Total-cell-normalized"
    ),
    normalization = factor(
      normalization,
      levels = c("CD3-normalized", "Total-cell-normalized")
    )
  )

p3 <- ggplot2::ggplot(
  coupling_long,
  ggplot2::aes(
    x = pct_cd3_total,
    y = pct_camelid,
    color = Status_plot
  )
) +
  ggplot2::geom_point(size = 2.4, alpha = 0.88) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey30",
    fill = "grey80",
    linewidth = 0.65
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(normalization),
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::scale_color_manual(values = pal_status, drop = FALSE) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(suffix = "%", accuracy = 1)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(suffix = "%", accuracy = 0.1)
  ) +
  ggplot2::labs(
    title = "Dependence of Camelid+ abundance on lymphocyte density",
    subtitle = paste(
      "A weaker relationship after CD3 normalization supports removal of",
      "stromal/non-lymphocyte composition effects."
    ),
    x = "CD3+ cells / total cells",
    y = "Camelid+ abundance",
    color = "COMET status"
  ) +
  ggplot2::theme(legend.position = "bottom")

spread_plot_df <- patient_spread %>%
  dplyr::select(
    ID, Status_plot, n_tissues,
    SD_log_CD3_normalized, SD_log_total_normalized
  ) %>%
  tidyr::pivot_longer(
    cols = c(SD_log_CD3_normalized, SD_log_total_normalized),
    names_to = "normalization",
    values_to = "SD_log_proportion"
  ) %>%
  dplyr::mutate(
    normalization = dplyr::recode(
      normalization,
      SD_log_CD3_normalized = "CD3-normalized",
      SD_log_total_normalized = "Total-cell-normalized"
    ),
    normalization = factor(
      normalization,
      levels = c("CD3-normalized", "Total-cell-normalized")
    )
  )

p4 <- ggplot2::ggplot(
  spread_plot_df,
  ggplot2::aes(
    x = normalization,
    y = SD_log_proportion,
    group = ID,
    color = Status_plot
  )
) +
  ggplot2::geom_line(alpha = 0.45, linewidth = 0.55) +
  ggplot2::geom_point(size = 2.7, alpha = 0.92) +
  ggplot2::scale_color_manual(values = pal_status, drop = FALSE) +
  ggplot2::labs(
    title = "Within-patient variability across tissue sites",
    subtitle = paste(
      "SD of log abundance across patient-specific tissue medians;",
      "lower values indicate less multiplicative variability."
    ),
    x = NULL,
    y = "Within-patient SD of log proportion",
    color = "COMET status"
  ) +
  ggplot2::theme(legend.position = "bottom")

r2_plot_df <- r2_bootstrap_results %>%
  dplyr::filter(
    scope == "All cohorts (adjusted for status)",
    tissue_definition == "Tissue_primary",
    normalization %in% c("CD3-normalized", "Total-cell-normalized")
  ) %>%
  dplyr::mutate(
    normalization = factor(
      normalization,
      levels = c("CD3-normalized", "Total-cell-normalized")
    )
  )

p5 <- ggplot2::ggplot(
  r2_plot_df,
  ggplot2::aes(
    x = normalization,
    y = estimate,
    color = normalization
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "grey70",
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = CI_low, ymax = CI_high),
    width = 0.12,
    linewidth = 0.8
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::scale_color_manual(values = pal_method, guide = "none") +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(
      0,
      max(c(r2_plot_df$CI_high, r2_plot_df$estimate, 0.05), na.rm = TRUE) * 1.12
    )
  ) +
  ggplot2::labs(
    title = "Fraction of abundance variability attributable to tissue",
    subtitle = paste(
      "Partial R-squared after adjustment for COMET status;",
      "95% patient-cluster bootstrap confidence intervals."
    ),
    x = NULL,
    y = "Tissue partial R-squared"
  )

message("\nPrinting Figure 1 of 5...")
show_and_save(
  p1,
  "fig1_CD3_fraction_of_total_by_tissue",
  width = 7.6,
  height = 5.4
)

message("\nPrinting Figure 2 of 5...")
show_and_save(
  p2,
  "fig2_Camelid_abundance_by_denominator_and_tissue",
  width = 11.0,
  height = 5.8
)

message("\nPrinting Figure 3 of 5...")
show_and_save(
  p3,
  "fig3_Camelid_abundance_vs_lymphocyte_density",
  width = 10.5,
  height = 5.6
)

message("\nPrinting Figure 4 of 5...")
show_and_save(
  p4,
  "fig4_within_patient_cross_tissue_variability",
  width = 7.6,
  height = 5.5
)

message("\nPrinting Figure 5 of 5...")
show_and_save(
  p5,
  "fig5_tissue_partial_R2_by_normalization",
  width = 7.2,
  height = 5.4
)

primary_r2_delta <- r2_bootstrap_results %>%
  dplyr::filter(
    scope == "All cohorts (adjusted for status)",
    tissue_definition == "Tissue_primary",
    normalization ==
      "Difference: total-cell minus CD3-normalized"
  ) %>%
  dplyr::slice_head(n = 1)

primary_cor_delta <- correlation_bootstrap %>%
  dplyr::filter(scope == "All cohorts (adjusted for status)") %>%
  dplyr::slice_head(n = 1)

primary_spread <- spread_test_results %>%
  dplyr::filter(
    scope == "All cohorts (adjusted for status)",
    variability_metric ==
      "SD of log proportion across tissue medians"
  ) %>%
  dplyr::slice_head(n = 1)

evidence_summary <- dplyr::bind_rows(
  tibble::tibble(
    test = "Tissue partial R2",
    estimate_supporting_CD3_normalization =
      primary_r2_delta$estimate,
    CI_low = primary_r2_delta$CI_low,
    CI_high = primary_r2_delta$CI_high,
    p_value = NA_real_,
    support = is.finite(primary_r2_delta$CI_low) &&
      primary_r2_delta$CI_low > 0,
    interpretation = paste(
      "Positive difference means tissue explains more variability under",
      "total-cell normalization."
    )
  ),
  tibble::tibble(
    test = "Coupling to lymphocyte density",
    estimate_supporting_CD3_normalization =
      primary_cor_delta$estimate,
    CI_low = primary_cor_delta$CI_low,
    CI_high = primary_cor_delta$CI_high,
    p_value = NA_real_,
    support = is.finite(primary_cor_delta$CI_low) &&
      primary_cor_delta$CI_low > 0,
    interpretation = paste(
      "Positive difference means total-cell normalization is more strongly",
      "coupled to CD3/total."
    )
  ),
  tibble::tibble(
    test = "Within-patient cross-tissue SD",
    estimate_supporting_CD3_normalization =
      primary_spread$median_paired_difference_total_minus_CD3,
    CI_low = NA_real_,
    CI_high = NA_real_,
    p_value = primary_spread$p_value,
    support =
      is.finite(primary_spread$median_paired_difference_total_minus_CD3) &&
      primary_spread$median_paired_difference_total_minus_CD3 > 0 &&
      is.finite(primary_spread$p_value) &&
      primary_spread$p_value < 0.05,
    interpretation = paste(
      "Positive difference means greater within-patient variability under",
      "total-cell normalization."
    )
  )
)

n_support <- sum(evidence_summary$support, na.rm = TRUE)
overall_interpretation <- dplyr::case_when(
  n_support >= 2 ~ paste(
    "Overall result: the analyses support CD3/lymphocyte normalization as",
    "reducing tissue-composition-associated variability."
  ),
  n_support == 1 ~ paste(
    "Overall result: evidence is mixed; CD3 normalization is supported by",
    "one of the three primary criteria."
  ),
  TRUE ~ paste(
    "Overall result: these data do not provide consistent evidence that CD3",
    "normalization reduces tissue-composition-associated variability."
  )
)

readr::write_tsv(
  count_audit,
  file.path(table_dir, "count_and_denominator_audit.tsv")
)
readr::write_tsv(
  identity_audit,
  file.path(table_dir, "normalization_algebra_identity_audit.tsv")
)
readr::write_tsv(
  tissue_summary,
  file.path(table_dir, "tissue_descriptive_summary.tsv")
)
readr::write_tsv(
  bb_results,
  file.path(table_dir, "beta_binomial_tissue_tests.tsv")
)
readr::write_tsv(
  correlation_results,
  file.path(table_dir, "lymphocyte_density_spearman_correlations.tsv")
)
readr::write_tsv(
  correlation_bootstrap,
  file.path(table_dir, "lymphocyte_density_correlation_bootstrap.tsv")
)
readr::write_tsv(
  r2_results,
  file.path(table_dir, "tissue_partial_R2_results.tsv")
)
readr::write_tsv(
  r2_bootstrap_results,
  file.path(table_dir, "tissue_partial_R2_cluster_bootstrap.tsv")
)
readr::write_tsv(
  patient_tissue_values,
  file.path(table_dir, "patient_tissue_median_log_abundance.tsv")
)
readr::write_tsv(
  patient_spread,
  file.path(table_dir, "within_patient_cross_tissue_variability.tsv")
)
readr::write_tsv(
  spread_test_results,
  file.path(table_dir, "within_patient_variability_tests.tsv")
)
readr::write_tsv(
  evidence_summary,
  file.path(table_dir, "normalization_evidence_summary.tsv")
)

writeLines(
  c(
    "COMET lymphocyte-normalization validation",
    paste0("Run date: ", Sys.Date()),
    paste0("Input: ", input_path),
    paste0("Output: ", output_dir),
    paste0("Minimum tissue n: ", min_tissue_n),
    paste0("Minimum patient tissue groups: ", min_patient_tissues),
    paste0("Patient-cluster bootstrap B: ", bootstrap_B),
    paste0("Bootstrap seed: ", bootstrap_seed),
    "Figures are printed to the active R graphics device.",
    "Saved PDF/EPS text remains live and selectable."
  ),
  con = file.path(log_dir, "analysis_settings.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

cat("\n\n")
cat("============================================================\n")
cat("COMET NORMALIZATION VALIDATION: STATISTICAL RESULTS\n")
cat("============================================================\n")

cat("\n1) Count and denominator audit\n")
print(count_audit, n = Inf, width = Inf)

cat("\n2) Algebraic identity audit\n")
cat(
  "Camelid/total = (Camelid/CD3) x (CD3/total).\n",
  "Numerical errors below should be approximately zero.\n",
  sep = ""
)
print(identity_audit, n = Inf, width = Inf)

cat("\n3) Tissue descriptive summary\n")
print(tissue_summary, n = Inf, width = Inf)

cat("\n4) Adjusted beta-binomial tissue-effect tests\n")
cat(
  "These likelihood-ratio tests compare models with versus without tissue,\n",
  "while adjusting for COMET status and using a patient random intercept when possible.\n",
  sep = ""
)
print(
  bb_results %>%
    dplyr::mutate(p_label = vapply(p_value, fmt_p, character(1))),
  n = Inf,
  width = Inf
)

cat("\n5) Spearman coupling to lymphocyte density (CD3/total)\n")
print(
  correlation_results %>%
    dplyr::mutate(p_label = vapply(p_value, fmt_p, character(1))),
  n = Inf,
  width = Inf
)

cat("\n6) Patient-cluster bootstrap comparison of correlation strength\n")
print(correlation_bootstrap, n = Inf, width = Inf)

cat("\n7) Tissue partial R-squared on the log-proportion scale\n")
cat(
  "Partial R2 is the additional fraction of variability explained by tissue\n",
  "after accounting for COMET status.\n",
  sep = ""
)
print(
  r2_results %>%
    dplyr::mutate(p_label = vapply(p_value, fmt_p, character(1))),
  n = Inf,
  width = Inf
)

cat("\n8) Patient-cluster bootstrap confidence intervals for tissue partial R2\n")
print(r2_bootstrap_results, n = Inf, width = Inf)

cat("\n9) Within-patient cross-tissue variability tests\n")
print(
  spread_test_results %>%
    dplyr::mutate(p_label = vapply(p_value, fmt_p, character(1))),
  n = Inf,
  width = Inf
)

cat("\n10) Primary evidence summary\n")
print(
  evidence_summary %>%
    dplyr::mutate(p_label = vapply(p_value, fmt_p, character(1))),
  n = Inf,
  width = Inf
)

cat("\n", overall_interpretation, "\n", sep = "")
cat(
  "\nInterpretive caution: stronger coupling of Camelid/total to CD3/total is\n",
  "partly mathematical because the total-cell-normalized measure contains the\n",
  "lymphocyte-density term by definition. The tissue partial-R2 and paired\n",
  "within-patient analyses quantify whether that coupling creates materially\n",
  "greater tissue-associated variability in this dataset.\n",
  sep = ""
)

cat("\nAnalysis complete. Outputs written to:\n", output_dir, "\n", sep = "")
