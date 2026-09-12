#!/usr/bin/env Rscript
# CAR expansion and persistence
# Summarize early exposure and late persistence with clinical covariate associations.

suppressPackageStartupMessages({
  library(readr);  library(dplyr);  library(tidyr);  library(stringr)
  library(ggplot2); library(scales); library(forcats); library(broom)
  library(mgcv);   library(purrr)
  library(tableone)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
path_df3 <- file.path(data_dir, "clinical_metadata.tsv")
path_df4 <- file.path(data_dir, "car_measurements.tsv")

theme_arial <- function(base_size = 22){
  theme_classic(base_size = base_size) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text  = element_text(size = base_size - 4))
}

safe_num <- function(x){
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(readr::parse_number(as.character(x)))
}

clean_ethnicity <- function(x){
  raw <- ifelse(is.na(x), "Unknown", as.character(x))
  raw <- stringr::str_replace_all(raw, "_", " ")
  raw <- stringr::str_squish(raw)
  low <- stringr::str_to_lower(raw)
  out <- dplyr::case_when(
    stringr::str_detect(low, "^not\\s*hispanic") |
      stringr::str_detect(low, "^non\\s*-?\\s*hispanic") ~ "Not Hispanic/Latino",
    stringr::str_detect(low, "hispanic|latino|latinx")   ~ "Hispanic/Latino",
    low %in% c("", "unknown", "na", "n/a", "declined")   ~ "Unknown",
    TRUE ~ "Unknown"
  )
  factor(out, levels = c("Hispanic/Latino","Not Hispanic/Latino","Unknown"))
}

clean_race3 <- function(x){
  low <- tolower(ifelse(is.na(x), "unknown", as.character(x)))
  lev <- dplyr::case_when(
    stringr::str_detect(low, "black|african")        ~ "Black",
    stringr::str_detect(low, "white|caucasian")      ~ "White",
    low %in% c("", "unknown", "na", "n/a", "declined") ~ "Other/Unknown",
    TRUE ~ "Other/Unknown"
  )
  factor(lev, levels = c("White","Black","Other/Unknown"))
}

pick_col <- function(nms, choices){ hit <- choices[choices %in% nms]; if (length(hit)) hit[[1]] else NA_character_ }

trapz_auc_d0_d30 <- function(x_day, y_car, xmin = 0, xmax = 30, lod = 0.1){
  df_raw <- tibble(x = as.numeric(x_day), y = as.numeric(y_car)) %>%
    filter(is.finite(x), is.finite(y)) %>%
    arrange(x)
  nz_in_0_30 <- sum(df_raw$y[df_raw$x >= xmin & df_raw$x <= xmax] > 0, na.rm = TRUE)
  if (nz_in_0_30 < 2) return(NA_real_)

  df <- df_raw %>%
    mutate(y = if_else(y <= 0, lod, y)) %>%
    mutate(x = pmax(x, xmin)) %>%
    distinct(x, .keep_all = TRUE)
  if (nrow(df) < 2) return(NA_real_)

  add_at <- function(x0) if (!any(abs(df$x - x0) < 1e-9))
    tibble(x = x0, y = approx(df$x, df$y, xout = x0, rule = 2)$y) else NULL

  df2 <- bind_rows(df, add_at(xmin), add_at(xmax)) %>% arrange(x) %>% filter(x >= xmin, x <= xmax)
  df2$y[df2$x == xmin] <- 0
  sum(diff(df2$x) * (head(df2$y, -1) + tail(df2$y, -1)) / 2)
}

suppressPackageStartupMessages({ library(grid); library(gridExtra) })

out_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "car_kinetics", "expansion_and_persistence"), "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (Sys.info()[["sysname"]] == "Darwin") {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

save_pdf <- function(plot, filename, width = 9, height = 7) {
  f <- file.path(out_dir, filename)
  device_fun <- function(file, width, height, ...) {
    grDevices::quartz(file = file, type = "pdf",
                      width = width, height = height, family = "Arial", ...)
  }
  ggplot2::ggsave(filename = f, plot = plot, device = device_fun,
                  width = width, height = height, units = "in")
  message("Saved: ", f)
}

save_grob_pdf <- function(grob, filename, width = 9, height = 7) {
  f <- file.path(out_dir, filename)
  grDevices::quartz(file = f, type = "pdf", width = width, height = height, family = "Arial")
  grid::grid.newpage(); grid::grid.draw(grob)
  grDevices::dev.off()
  message("Saved: ", f)
}

save_csv <- function(df, filename) {
  f <- file.path(out_dir, filename)
  readr::write_csv(df, f)
  message("Saved: ", f)
}

# Clinical covariates
df3 <- read_delim(path_df3, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
if (!"Study_ID" %in% names(df3) && "ID" %in% names(df3)) df3 <- rename(df3, Study_ID = ID)

meta <- df3 %>%
  mutate(
    Age_at_infusion   = safe_num(Age_at_infusion),
    Sex               = factor(Sex, levels = c("Male","Female")),
    Prior_Lines       = safe_num(Prior_Lines),
    Race3             = clean_race3(Race),
    Ethnicity3        = clean_ethnicity(Ethnicity),
    PreLD_ALC         = safe_num(ALC_PreLD),
    Baseline_LDH      = safe_num(Baseline_pre_LD_LDH),
    Baseline_CRP      = safe_num(`Baseline_pre-LD_CRP`),
    Baseline_Ferritin = safe_num(`Baseline_pre-LD_Ferritin`),
    M_spike           = safe_num(`preLD_M-spike`),
    IEC_flag          = factor(if_else(safe_num(IEC_enteritis) == 1, "IEC", "No IEC"),
                               levels = c("No IEC","IEC"))
  ) %>%
  select(Study_ID, IEC_flag, Age_at_infusion, Sex, Prior_Lines, Race3, Ethnicity3,
         PreLD_ALC, Baseline_LDH, Baseline_CRP, Baseline_Ferritin, M_spike)

cat("\n== IEC counts ==\n"); print(table(meta$IEC_flag, useNA = "ifany"))

# CAR measurements
df4 <- read_delim(path_df4, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)

day_col <- pick_col(names(df4), c("Day","day","Day_rel_infusion","DayRelative","Day_from_infusion","Day_rel"))
car_col <- pick_col(names(df4), c("CAR_abs","CAR_Abs","Absolute_CAR","CARcells_per_uL","CAR_cells_per_uL"))
if (is.na(day_col) || is.na(car_col)) {
  stop("Could not detect Day and/or CAR columns in df4. Found: ",
       paste(names(df4), collapse = ", "))
}

LOD_ND <- 0.1
exp_all <- df4 %>%
  transmute(
    Study_ID,
    Day     = safe_num(.data[[day_col]]),
    CAR_abs = safe_num(.data[[car_col]])
  ) %>%
  filter(!is.na(Study_ID), is.finite(Day), Day >= 0) %>%
  mutate(CAR_abs_imp = if_else(!is.finite(CAR_abs) | CAR_abs <= 0, LOD_ND, CAR_abs)) %>%
  inner_join(meta %>% select(Study_ID, IEC_flag), by = "Study_ID")

# Sampling coverage
exp_counts <- exp_all %>%
  mutate(win_0_30 = Day >= 0  & Day <= 30,
         win_gt50 = Day > 50 & Day <= 365) %>%
  group_by(Study_ID, IEC_flag) %>%
  summarise(n_0_30 = sum(win_0_30),
            n_gt50 = sum(win_gt50),
            .groups = "drop")

cat("\n== Per-patient expansion counts by IEC group ==\n")
print(exp_counts %>%
        group_by(IEC_flag) %>%
        summarise(across(c(n_0_30, n_gt50),
                         list(median = median, IQR = IQR)),
                  .groups = "drop"))

p_hist_0_30 <- ggplot(exp_counts, aes(n_0_30, fill = IEC_flag)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "white",
                 position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("No IEC" = "#6baed6", "IEC" = "#74c476"), name = "IEC") +
  labs(title = "Counts of expansion measurements per patient (Days 0–30)",
       x = "Count in 0–30 days", y = "Patients") +
  theme_arial(22)

p_hist_gt50 <- ggplot(exp_counts, aes(n_gt50, fill = IEC_flag)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "white",
                 position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("No IEC" = "#6baed6", "IEC" = "#74c476"), name = "IEC") +
  labs(title = "Counts of expansion measurements per patient (Days >50–365)",
       x = "Count in >50–365 days", y = "Patients") +
  theme_arial(22)

print(p_hist_0_30); print(p_hist_gt50)

cat("\n== Duplicate (Study_ID, Day) rows ==\n")
dup_tbl <- exp_all %>% count(Study_ID, Day) %>% filter(n > 1)
print(head(dup_tbl, 10))
cat(sprintf("Total exact duplicate day rows: %d\n", nrow(dup_tbl)))

cat("\n== Day heaping (integer vs non-integer) ==\n")
int_prop <- mean(abs(exp_all$Day - round(exp_all$Day)) < 1e-9)
cat(sprintf("Proportion integer days: %.1f%%\n", 100*int_prop))

cat("\n== Non-detect fraction (raw CAR_abs ≤ 0 or NA) ==\n")
nd_frac <- mean(exp_all$CAR_abs <= 0 | is.na(exp_all$CAR_abs))
cat(sprintf("Rows imputed to LOD (0.1): %.1f%%\n", 100*nd_frac))

vars_pat <- c("Age_at_infusion","Sex","Prior_Lines","Race3","Ethnicity3",
              "PreLD_ALC","Baseline_LDH","Baseline_CRP","Baseline_Ferritin","M_spike","IEC_flag")
miss_pat <- meta %>%
  select(all_of(vars_pat)) %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "miss_frac") %>%
  arrange(desc(miss_frac))
cat("\n== Patient-level missingness (selected covariates) ==\n"); print(miss_pat)

p_miss_pat <- ggplot(miss_pat, aes(x = reorder(variable, miss_frac), y = miss_frac)) +
  geom_col() + coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Missingness by variable (patient data)", x = NULL, y = "Missing (%)") +
  theme_arial(22)
print(p_miss_pat)

aucs <- exp_all %>%
  filter(Day <= 45) %>%
  group_by(Study_ID) %>%
  summarise(
    AUC_0_30        = trapz_auc_d0_d30(Day, CAR_abs, xmin = 0, xmax = 30, lod = LOD_ND),
    n_nonzero_0_30  = sum(Day >= 0 & Day <= 30 & CAR_abs > 0, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n== AUC build summary ==\n")
print(aucs %>% summarise(
  n_total = n(),
  n_with_AUC = sum(is.finite(AUC_0_30)),
  pct_with_AUC = 100*mean(is.finite(AUC_0_30)),
  median_nonzero = median(n_nonzero_0_30, na.rm = TRUE)
))

dat0 <- meta %>%
  inner_join(aucs, by = "Study_ID") %>%
  filter(is.finite(AUC_0_30)) %>%
  mutate(
    logAUC    = log10(AUC_0_30 + 1),
    z_age     = as.numeric(scale(Age_at_infusion)),
    z_lines   = as.numeric(scale(Prior_Lines)),
    z_log_ALC = as.numeric(scale(log10(PreLD_ALC + 1))),
    z_log_LDH = as.numeric(scale(log10(Baseline_LDH + 1))),
    z_log_CRP = as.numeric(scale(log10(Baseline_CRP + 1))),
    z_log_Fer = as.numeric(scale(log10(Baseline_Ferritin + 1))),
    z_Mspike  = as.numeric(scale(M_spike))
  )

cat("\n== Analysis dataset missingness (rows used in OLS) ==\n")
print(summarise(dat0, across(everything(), ~mean(is.na(.)))))

num_preds <- c("z_age","z_lines","z_log_ALC","z_log_LDH","z_log_CRP","z_log_Fer","z_Mspike")
linearity_tbl <- map_dfr(num_preds, function(v){
  f_lin <- reformulate(v, response = "logAUC")
  f_gam <- as.formula(paste0("logAUC ~ s(", v, ", k=4)"))
  m_lin <- lm(f_lin, data = dat0)
  m_gam <- mgcv::gam(f_gam, data = dat0, method = "REML")
  a <- try(anova(m_lin, m_gam, test = "F"), silent = TRUE)
  p_nl <- if (inherits(a, "try-error")) NA_real_ else a$`Pr(>F)`[2]
  tibble::tibble(
    predictor = v,
    AIC_linear = AIC(m_lin),
    AIC_gam    = AIC(m_gam),
    dAIC       = AIC(m_lin) - AIC(m_gam),
    p_nonlinearity = p_nl
  )
})
cat("\n== Linearity screen (positive dAIC favors GAM; small p suggests nonlinearity) ==\n")
print(linearity_tbl)

f_ols <- logAUC ~ z_age + z_lines + z_log_ALC + z_log_LDH + z_log_CRP + z_log_Fer + z_Mspike +
  Sex + Race3 + Ethnicity3
fit_ols <- lm(f_ols, data = dat0)
cat("\n=== Multivariable OLS: log10(AUC+1) ~ pre-treatment predictors ===\n")
print(summary(fit_ols))

infl <- broom::augment(fit_ols) %>%
  mutate(cooks = .cooksd, hat = .hat, std_resid = .std.resid) %>%
  arrange(desc(cooks))
cat("\n== Top 5 Cook's distance observations ==\n")
print(head(infl %>% select(.rownames, cooks, hat, std_resid), 5))

if (requireNamespace("car", quietly = TRUE)) {
  cat("\n== Collinearity (VIF) ==\n"); print(car::vif(fit_ols))
} else {
  cat("\n== 'car' not installed; pairwise correlations among numeric predictors ==\n")
  num_mat <- dat0 %>% select(all_of(c("z_age","z_lines","z_log_ALC","z_log_LDH","z_log_CRP","z_log_Fer","z_Mspike"))) %>% as.matrix()
  print(round(cor(num_mat, use = "pairwise.complete.obs"), 2))
}

t1_vars <- c("Age_at_infusion","Sex","Race3","Ethnicity3","Prior_Lines",
             "PreLD_ALC","Baseline_LDH","Baseline_CRP","Baseline_Ferritin","M_spike")
# Cohort balance summaries
t1_df   <- meta %>% select(IEC_flag, all_of(t1_vars))
factorVars <- c("Sex","Race3","Ethnicity3")
nonnormal  <- c("Age_at_infusion","Prior_Lines","PreLD_ALC","Baseline_LDH","Baseline_CRP","Baseline_Ferritin","M_spike")

cat("\n=== TableOne (IEC vs No IEC) — SMDs reflect covariate imbalance ===\n")
t1 <- CreateTableOne(
  data = t1_df, vars = t1_vars, strata = "IEC_flag",
  factorVars = intersect(factorVars, names(t1_df)), includeNA = TRUE
)
print(t1, smd = TRUE, nonnormal = intersect(nonnormal, names(t1_df)), quote = TRUE, noSpaces = TRUE)

first_post50 <- exp_all %>%
  filter(Day > 50, Day <= 365) %>%
  arrange(Study_ID, Day) %>%
  group_by(Study_ID) %>%
  slice(1) %>%
  ungroup()

cat("\n== First post-50 sample day summary ==\n")
print(summary(first_post50$Day))

p_first_post50 <- ggplot(first_post50, aes(Day)) +
  geom_histogram(binwidth = 7, boundary = 0, color = "white") +
  labs(title = "Distribution of first post-50 expansion timepoint", x = "Day", y = "Patients") +
  theme_arial(22)
print(p_first_post50)

if (!exists("meta") || !exists("exp_all")) {
  library(readr); library(dplyr); library(tidyr)
  df3 <- read_delim(path_df3,
                    delim="\t", trim_ws=TRUE, show_col_types=FALSE)
  df4 <- read_delim(path_df4,
                    delim="\t", trim_ws=TRUE, show_col_types=FALSE)

  meta <- df3 %>%
    transmute(Study_ID,
              IEC_flag = factor(if_else(IEC_enteritis == 1, "IEC", "No IEC"),
                                levels = c("No IEC","IEC")))

  day_col <- dplyr::coalesce(
    names(df4)[match(c("Day","Day_rel_infusion"), names(df4))][1],
    NA_character_)
  stopifnot(!is.na(day_col))
  exp_all <- df4 %>%
    transmute(Study_ID,
              Day = as.numeric(.data[[day_col]]),
              CAR_abs = suppressWarnings(as.numeric(CAR_abs))) %>%
    inner_join(meta, by="Study_ID") %>%
    filter(is.finite(Day), Day >= 0)
}

day_low  <- 50
day_high <- Inf

iec_post50 <- exp_all %>%
  filter(Day > day_low, Day <= day_high, IEC_flag == "IEC")

iec_post50_ids <- iec_post50 %>% distinct(Study_ID)
n_iec_post50   <- nrow(iec_post50_ids)

cat(sprintf("IEC patients with any sample after Day %d (≤%s): %d\n",
            day_low, ifelse(is.infinite(day_high), "Inf", as.character(day_high)), n_iec_post50))

first_post50_IEC <- iec_post50 %>%
  arrange(Study_ID, Day) %>%
  group_by(Study_ID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(bin = cut(Day, breaks = c(50, 100, 200, Inf),
                   labels = c("51–100", "101–200", "201+"), right = TRUE))

cat("\nBreakdown by day of FIRST post-50 sample (IEC only):\n")
print(first_post50_IEC %>% count(bin))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(scales)
})

trapz_auc_d0_d30 <- function(x, y, xmin = 0, xmax = 30, lod = 0.1) {
  df <- tibble(x = as.numeric(x), y = as.numeric(y)) |>
    filter(is.finite(x), is.finite(y)) |>
    arrange(x) |>
    mutate(y = if_else(y <= 0, lod, y),
           x = pmin(pmax(x, xmin), xmax)) |>
    distinct(x, .keep_all = TRUE)
  if (nrow(df) < 2) return(NA_real_)
  add_at <- function(x0) if (!any(abs(df$x - x0) < 1e-9))
    tibble(x = x0, y = approx(df$x, df$y, xout = x0, rule = 2)$y) else NULL
  df2 <- bind_rows(df, add_at(xmin), add_at(xmax)) |> arrange(x)
  sum(diff(df2$x) * (head(df2$y, -1) + tail(df2$y, -1)) / 2)
}

LOD_ND <- 0.1

aucs <- exp_all %>%
  filter(Day <= 45) %>%
  group_by(Study_ID) %>%
  summarise(
    AUC_0_30       = trapz_auc_d0_d30(Day, CAR_abs, xmin = 0, xmax = 30, lod = LOD_ND),
    n_meas_0_30    = sum(Day >= 0 & Day <= 30, na.rm = TRUE),
    n_nonzero_0_30 = sum(Day >= 0 & Day <= 30 & CAR_abs > 0, na.rm = TRUE),
    .groups = "drop"
  )

evaluable <- aucs %>% filter(is.finite(AUC_0_30)) %>% select(Study_ID)
cat(sprintf("Patients in meta: %d | AUC-evaluable: %d (%.1f%%)\n",
            nrow(meta), nrow(evaluable),
            100 * nrow(evaluable) / nrow(meta)))

vars_pat <- c("Age_at_infusion","Sex","Prior_Lines","Race3","Ethnicity3",
              "PreLD_ALC","Baseline_LDH","Baseline_CRP","Baseline_Ferritin","M_spike","IEC_flag")

miss_all <- meta %>%
  select(all_of(vars_pat)) %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "miss_frac") %>%
  mutate(cohort = "All patients")

miss_eval <- meta %>%
  semi_join(evaluable, by = "Study_ID") %>%
  select(all_of(vars_pat)) %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "miss_frac") %>%
  mutate(cohort = "AUC-evaluable")

miss_both <- bind_rows(miss_all, miss_eval) %>%
  arrange(variable, cohort)

cat("\n== Missingness summary (fraction NA) ==\n")
print(miss_both)

p_miss <- miss_both %>%
  mutate(variable = reorder(variable, miss_frac)) %>%
  ggplot(aes(variable, miss_frac, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.7)) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Missingness by variable: All vs AUC-evaluable",
       x = NULL, y = "Missing (%)", fill = NULL) +
  theme_classic(base_size = 16)
print(p_miss)

excluded <- meta %>% anti_join(evaluable, by = "Study_ID") %>%
  left_join(aucs, by = "Study_ID") %>%
  mutate(reason = case_when(
    is.na(n_meas_0_30) | n_meas_0_30 == 0 ~ "No measurements in 0–30",
    n_nonzero_0_30 == 0                   ~ "All ≤ LOD (0–30)",
    n_nonzero_0_30 == 1                   ~ "Only 1 > LOD in 0–30 (needs ≥2)",
    !is.finite(AUC_0_30)                  ~ "AUC not computable",
    TRUE                                  ~ "Other"
  )) %>%
  select(Study_ID, n_meas_0_30, n_nonzero_0_30, AUC_0_30, reason)

cat("\n== Patients NOT AUC-evaluable (with reason) ==\n")
print(excluded)

counts <- tibble(
  cohort = c("All patients", "AUC-evaluable"),
  n      = c(nrow(meta), nrow(evaluable))
)
print(counts)

save_pdf(p_hist_0_30, "EXP_counts_hist_0_30.pdf",  width = 9,  height = 7)
save_pdf(p_hist_gt50, "EXP_counts_hist_50_365.pdf", width = 9,  height = 7)

save_pdf(p_miss_pat, "Missingness_patient_covariates.pdf", width = 9, height = 7)
save_csv(miss_pat,    "Missingness_patient_covariates.csv")

save_csv(aucs, "AUC0_30_per_patient.csv")

aucs_summary <- aucs %>% summarise(
  n_total = n(),
  n_with_AUC = sum(is.finite(AUC_0_30)),
  pct_with_AUC = 100*mean(is.finite(AUC_0_30)),
  median_nonzero = median(n_nonzero_0_30, na.rm = TRUE)
)
save_csv(aucs_summary, "AUC0_30_summary.csv")

save_csv(linearity_tbl, "Linearity_screen_LM_vs_GAM.csv")

save_csv(infl %>% select(.rownames, cooks, hat, std_resid),
         "OLS_influence_top.csv")

if (requireNamespace("car", quietly = TRUE)) {
  vif_tbl <- tibble::tibble(variable = names(car::vif(fit_ols)),
                            VIF      = as.numeric(car::vif(fit_ols)))
  save_csv(vif_tbl, "OLS_VIF.csv")
} else {

  num_mat <- dat0 %>%
    select(all_of(c("z_age","z_lines","z_log_ALC","z_log_LDH","z_log_CRP","z_log_Fer","z_Mspike"))) %>%
    as.matrix()
  save_csv(as.data.frame(round(cor(num_mat, use = "pairwise.complete.obs"), 3)),
           "OLS_numeric_cor_matrix.csv")
}

t1_df_out <- print(t1,
                   smd = TRUE,
                   nonnormal = intersect(nonnormal, names(t1_df)),
                   quote = TRUE, noSpaces = TRUE,
                   printToggle = FALSE)

save_csv(as.data.frame(t1_df_out), "Table1_IEC_vs_NoIEC.csv")

t1_mat <- as.matrix(t1_df_out)
tg <- gridExtra::tableGrob(
  t1_mat, rows = NULL,
  theme = gridExtra::ttheme_minimal(
    core    = list(fg_params = list(fontfamily = "Arial", cex = 0.60)),
    colhead = list(fg_params = list(fontfamily = "Arial", fontface = 2, cex = 0.70)),
    rowhead = list(fg_params = list(fontfamily = "Arial", fontface = 2, cex = 0.70))
  )
)
title_grob <- grid::textGrob("Table 1. Baseline Characteristics (IEC vs No IEC)",
                             gp = grid::gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"))
tbl_page <- gridExtra::arrangeGrob(title_grob, tg,
                                   heights = grid::unit.c(grid::unit(0.5, "in"),
                                                          grid::unit(1, "npc") - grid::unit(0.5, "in")))
save_grob_pdf(tbl_page, "Table1_IEC_vs_NoIEC.pdf", width = 11, height = 8.5)
