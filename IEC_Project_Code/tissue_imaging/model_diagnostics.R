#!/usr/bin/env Rscript
# Tissue-model diagnostics
# Check count validity, random-effect structures, convergence, and simulated residuals.

suppressPackageStartupMessages({
  library(tidyverse)
  library(forcats)
  library(glmmTMB)
  library(DHARMa)
  library(splines)
  library(performance)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
path_df7 <- file.path(data_dir, "tissue_cell_counts.tsv")
stopifnot(file.exists(path_df7))
df7_raw <- readr::read_delim(path_df7, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)

status_label <- function(cohort) {
  z <- as.character(cohort)
  dplyr::case_when(
    z %in% c("IEC_EC", "IEC-EC") ~ "IEC-EC",
    z %in% c("CONTROL", "Control") ~ "Untreated Controls",
    z %in% c("NOT_IEC","NOT_IEC_EC","Not_IEC_EC","NOT_IEC-EC","Indeterminate","INDETERMINATE") ~ "Not IEC-EC",
    TRUE ~ z
  )
}

tissue_fine <- function(tt) {
  z <- tolower(gsub("\\s+","_", as.character(tt)))
  dplyr::case_when(
    grepl("^duodenum$", z)               ~ "Duodenum",
    grepl("terminal_ileum", z)           ~ "Terminal Ileum",
    grepl("ileocecal", z)                ~ "Ileocecal Valve",
    grepl("^colon_left$", z)             ~ "Colon (left)",
    grepl("^colon_right$", z)            ~ "Colon (right)",
    grepl("^colon_random$", z)           ~ "Colon (random)",
    grepl("^colon$", z)                  ~ "Colon",
    grepl("^rectum$", z)                 ~ "Rectum",
    grepl("^esophagus$", z)              ~ "Esophagus",
    grepl("^stomach_antrum$", z)         ~ "Stomach (antrum)",
    grepl("^stomach_polyp$", z)          ~ "Stomach (polyp)",
    grepl("^stomach$", z)                ~ "Stomach",
    TRUE                                 ~ stringr::str_to_title(gsub("_"," ",tt))
  )
}

tissue_primary_from_fine <- function(f) {
  dplyr::case_when(
    f %in% c("Duodenum") ~ "Duodenum",
    f %in% c("Terminal Ileum") ~ "Terminal Ileum",
    f %in% c("Colon","Colon (left)","Colon (right)","Colon (random)","Ileocecal Valve") ~ "Colon",
    TRUE ~ "Other"
  )
}

gi_group3 <- function(f) {
  dplyr::case_when(
    f %in% c("Colon","Colon (left)","Colon (right)","Colon (random)","Ileocecal Valve") ~ "Colon",
    f %in% "Duodenum"       ~ "Duodenum",
    f %in% "Terminal Ileum" ~ "Terminal Ileum",
    TRUE ~ NA_character_
  )
}

frac_safe <- function(num, den) {
  num <- as.numeric(num); den <- as.numeric(den)
  ifelse(is.finite(num) & is.finite(den) & den > 0, num/den, NA_real_)
}

df7 <- df7_raw %>%
  transmute(
    ID          = Patient_ID,
    Cohort,
    Tissue_Type,
    total_cells,
    CD68,
    CD3,
    CD3_CD103             = `CD3-CD103`,
    CD3_Camelid           = CD3_Camelid,

    cam_CD4               = `CD3-Camelid-CD4`,
    cam_CD8               = `CD3-Camelid-CD8`,
    cam_DP                = `CD3-Camelid-DP`,
    cam_DN                = `CD3-Camelid-DN`,
    cam_GZMB              = `CD3-Camelid-GZMB`,
    cam_PD1               = `CD3-Camelid-PD1`,

    CD3_nonCamelid        = `CD3-Camelid-negative`,
    non_CD4               = `CD3-Camelid-negative-CD4`,
    non_CD8               = `CD3-Camelid-negative-CD8`,
    non_DP                = `CD3-Camelid-negative-DP`,
    non_DN                = `CD3-Camelid-negative-DN`,
    non_GZMB              = `CD3-Camelid-negative-GZMB`,
    non_PD1               = `CD3-Camelid-negative-PD1`,

    CD103_fraction_tcells,
    Status_main    = status_label(Cohort),
    Tissue_fine    = tissue_fine(Tissue_Type),
    Tissue_primary = tissue_primary_from_fine(tissue_fine(Tissue_Type)),
    gi_frac_cd3_camelid   = frac_safe(CD3_Camelid, CD3),
    frac_CD103            = CD103_fraction_tcells,
    frac_CD68_total       = frac_safe(CD68, total_cells)
  )

iec_ec    <- df7 %>% filter(Status_main == "IEC-EC")
not_iecec <- df7 %>% filter(Status_main == "Not IEC-EC")

cat("\n--- df7 summary (Status_main x Tissue_primary) ---\n")
print(df7 %>% count(Status_main, Tissue_primary))

# Beta-binomial model diagnostics
run_bb_checks <- function(
  dat, succ, total, fixed, random_full, random_reduced = NULL,
  family = betabinomial(link="logit"),
  label = "Model", id_col = NULL, extra_info = list()
) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("[DIAGNOSTICS] ", label, "\n", sep = "")
  cat(strrep("-", 80), "\n", sep = "")

  s <- dat[[succ]]; t <- dat[[total]]
  n0 <- nrow(dat)
  bad_nonfinite <- sum(!is.finite(s) | !is.finite(t))
  bad_total     <- sum(!(t > 0), na.rm = TRUE)
  bad_bounds    <- sum(!(s >= 0 & s <= t), na.rm = TRUE)

  cat(sprintf("Rows: %d | Non-finite succ/total: %d | total<=0: %d | succ out of bounds: %d\n",
              n0, bad_nonfinite, bad_total, bad_bounds))

  dat_ok <- dat %>%
    filter(is.finite(.data[[succ]]), is.finite(.data[[total]]),
           .data[[total]] > 0, .data[[succ]] >= 0, .data[[succ]] <= .data[[total]])
  cat(sprintf("Rows retained after basic count checks: %d\n", nrow(dat_ok)))

  f_full    <- as.formula(paste0("cbind(", succ, ", ", total, " - ", succ, ") ~ ", fixed, " + ", random_full))
  f_reduced <- if (!is.null(random_reduced))
    as.formula(paste0("cbind(", succ, ", ", total, " - ", succ, ") ~ ", fixed, " + ", random_reduced)) else NULL

  options(glmmTMB.verbose = FALSE)
  fit_full_M <- try(glmmTMB::glmmTMB(f_full, family = family, data = dat_ok), silent = TRUE)
  if (inherits(fit_full_M, "try-error")) {
    cat("[WARN] Full model failed to fit. Aborting diagnostics for this model.\n")
    return(invisible(NULL))
  }

  fit_re_M <- NULL
  if (!is.null(f_reduced)) {
    fit_re_M <- try(glmmTMB::glmmTMB(f_reduced, family = family, data = dat_ok), silent = TRUE)
    if (inherits(fit_re_M, "try-error")) fit_re_M <- NULL
  }

  cat("\n[AIC / Likelihood comparisons (ML)]\n")
  if (!is.null(fit_re_M)) {
    print(try(anova(fit_re_M, fit_full_M), silent = TRUE))
  } else {
    print(try(AIC(fit_full_M), silent = TRUE))
  }

  fit_final <- fit_full_M

  cat("\n[Convergence / Random-effects variances]\n")
  summ <- summary(fit_final)
  print(summ$fit)
  vc <- VarCorr(fit_final)
  print(vc)

  re_vars <- try(as.data.frame(vc), silent = TRUE)
  if (!inherits(re_vars, "try-error")) {
    tiny <- re_vars %>% filter(vcov < 1e-6)
    if (nrow(tiny)) {
      cat("\n[NOTE] Near-zero RE variances (possible overparameterization):\n")
      print(tiny)
    }
  }

  phi <- try(sigma(fit_final), silent = TRUE)
  if (!inherits(phi, "try-error") && is.finite(phi)) {
    rho <- 1 / (as.numeric(phi) + 1)
    cat(sprintf("\n[Overdispersion] Beta-binomial precision (phi) = %.3f; ICC ≈ %.3f\n",
                as.numeric(phi), rho))
  }

  cat("\n[R2: marginal/conditional]\n")
  r2out <- try(performance::r2(fit_final), silent = TRUE)
  if (!inherits(r2out, "try-error")) print(r2out) else cat("R2 not available.\n")

  cat("\n[Fitted probabilities / design conditioning]\n")
  eta_hat <- as.numeric(predict(fit_final, type = "link"))
  p_hat   <- plogis(eta_hat)
  cat(sprintf("Fitted p range: [%.4f, %.4f]; median=%.4f\n",
              min(p_hat, na.rm=TRUE), max(p_hat, na.rm=TRUE), median(p_hat, na.rm=TRUE)))

  mm <- try(model.matrix(formula(f_final <- update.formula(f_full, . ~ . - (1|dummy))), data = dat_ok), silent = TRUE)
  if (!inherits(mm, "try-error")) {
    kappa_val <- try(kappa(mm), silent = TRUE)
    if (!inherits(kappa_val, "try-error")) cat(sprintf("Fixed-effect design condition number κ ≈ %.1f\n", kappa_val))
  }

  set.seed(123)
  cat("\n[DHARMa residual diagnostics]\n")
  resDH <- try(DHARMa::simulateResiduals(fit_final, plot = FALSE, n = 1000), silent = TRUE)
  if (inherits(resDH, "try-error")) {
    cat("DHARMa simulation failed.\n")
  } else {
    cat("\n- Uniformity (KS):\n"); print(try(testUniformity(resDH), silent = TRUE))
    cat("\n- Dispersion:\n");       print(try(testDispersion(resDH),  silent = TRUE))
    cat("\n- Zero-inflation:\n");   print(try(testZeroInflation(resDH), silent = TRUE))
    cat("\n- Outliers:\n");         print(try(testOutliers(resDH),    silent = TRUE))
  }

  if (length(extra_info)) {
    cat("\n[Context]\n")
    print(extra_info)
  }

  invisible(fit_final)
}

colon_set   <- c("Colon","Colon (left)","Colon (right)","Colon (random)","Ileocecal Valve")
stomach_set <- c("Stomach","Stomach (antrum)","Stomach (polyp)")

iec_sites_glmm <- iec_ec %>%
  filter(is.finite(CD3_Camelid), is.finite(CD3), CD3 > 0) %>%
  mutate(
    Tissue_fine = tissue_fine(Tissue_Type),
    Site_simple = case_when(
      Tissue_fine %in% colon_set      ~ "Colon",
      Tissue_fine %in% stomach_set    ~ "Stomach",
      Tissue_fine == "Duodenum"       ~ "Duodenum",
      Tissue_fine == "Terminal Ileum" ~ "Terminal Ileum",
      Tissue_fine == "Rectum"         ~ "Rectum",
      Tissue_fine == "Esophagus"      ~ "Esophagus",
      TRUE                            ~ "Other"
    ),
    Site_simple = factor(Site_simple,
      levels = c("Esophagus","Stomach","Duodenum","Terminal Ileum","Colon","Rectum","Other")),
    succ = CD3_Camelid,
    total = CD3,
    fail = pmax(total - succ, 0)
  ) %>% filter(!is.na(Site_simple))

iec_sites_ref <- iec_sites_glmm %>%
  mutate(Site_simple = stats::relevel(Site_simple, ref = "Duodenum"))

run_bb_checks(
  dat = iec_sites_ref,
  succ = "succ", total = "total",
  fixed = "Site_simple",
  random_full = "(1|ID)",
  random_reduced = NULL,
  label = "IEC-EC: Site model (Duodenum ref)"
)

iec_sites_bin <- iec_sites_glmm %>%
  mutate(Site_duo = factor(if_else(Site_simple=="Duodenum","Duodenum","Other"),
                           levels=c("Other","Duodenum")))
run_bb_checks(
  dat = iec_sites_bin %>% mutate(succ = CD3_Camelid, total = CD3),
  succ = "succ", total = "total",
  fixed = "Site_duo",
  random_full = "(1|ID)",
  label = "IEC-EC: Duodenum vs Other"
)

cdsub_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary,
            CD3_Camelid, CD3_nonCamelid,
            cam_CD4, cam_CD8, cam_DP, cam_DN,
            non_CD4, non_CD8, non_DP, non_DN) %>%
  pivot_longer(cols = c(cam_CD4, cam_CD8, cam_DP, cam_DN, non_CD4, non_CD8, non_DP, non_DN),
               names_to = "key", values_to = "succ") %>%
  mutate(group  = if_else(startsWith(key,"cam_"), "Camelid+", "Camelid-"),
         subset = sub("^(cam_|non_)","", key),
         total  = if_else(group=="Camelid+", CD3_Camelid, CD3_nonCamelid)) %>%
  filter(is.finite(total), total > 0, is.finite(succ), succ >= 0, succ <= total) %>%
  mutate(group = factor(group, levels=c("Camelid-","Camelid+")),
         subset= factor(subset,levels=c("CD4","CD8","DP","DN")))

for (s in levels(cdsub_glmm_dat$subset)) {
  dat_s <- cdsub_glmm_dat %>% filter(subset == s)
  if (nrow(dat_s) == 0 || n_distinct(dat_s$group) < 2) next
  run_bb_checks(
    dat = dat_s,
    succ = "succ", total = "total",
    fixed = "group + Tissue_primary",
    random_full   = "(1|ID) + (1|ID:Tissue_primary)",
    random_reduced= "(1|ID)",
    label = paste0("IEC-EC: CD subset ", s, " (Camelid+ vs Camelid−)")
  )
}

pd1_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid, cam_PD1, non_PD1) %>%
  pivot_longer(c(cam_PD1, non_PD1), names_to="key", values_to="succ") %>%
  mutate(group = if_else(grepl("^cam_",key),"Camelid+","Camelid-"),
         total = if_else(group=="Camelid+", CD3_Camelid, CD3_nonCamelid)) %>%
  filter(is.finite(total), total>0, is.finite(succ), succ>=0, succ<=total) %>%
  mutate(group=factor(group, levels=c("Camelid-","Camelid+")))
if (nrow(pd1_glmm_dat)>0 && n_distinct(pd1_glmm_dat$group)==2) {
  run_bb_checks(
    dat = pd1_glmm_dat,
    succ = "succ", total = "total",
    fixed = "group + Tissue_primary",
    random_full   = "(1|ID) + (1|ID:Tissue_primary)",
    random_reduced= "(1|ID)",
    label = "IEC-EC: PD1 (Camelid+ vs Camelid−)"
  )
}

gzmb_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid, cam_GZMB, non_GZMB) %>%
  pivot_longer(c(cam_GZMB, non_GZMB), names_to="key", values_to="succ") %>%
  mutate(group=if_else(grepl("^cam_",key),"Camelid+","Camelid-"),
         total=if_else(group=="Camelid+",CD3_Camelid,CD3_nonCamelid)) %>%
  filter(is.finite(total), total>0, is.finite(succ), succ>=0, succ<=total) %>%
  mutate(group=factor(group,levels=c("Camelid-","Camelid+")))
if (nrow(gzmb_glmm_dat)>0 && n_distinct(gzmb_glmm_dat$group)==2) {
  run_bb_checks(
    dat = gzmb_glmm_dat,
    succ = "succ", total = "total",
    fixed = "group + Tissue_primary",
    random_full   = "(1|ID) + (1|ID:Tissue_primary)",
    random_reduced= "(1|ID)",
    label = "IEC-EC: GZMB (Camelid+ vs Camelid−)"
  )
}

cd103_glmm_bb <- df7 %>%
  filter(Status_main %in% c("IEC-EC","Not IEC-EC"),
         is.finite(CD3_CD103), is.finite(CD3), CD3 > 0) %>%
  mutate(Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
         succ = as.integer(round(CD3_CD103)),
         total= as.integer(round(CD3))) %>%
  mutate(succ = pmax(0L, pmin(succ, total)))
if (nrow(cd103_glmm_bb)>0 && n_distinct(cd103_glmm_bb$Status_two)==2) {
  run_bb_checks(
    dat = cd103_glmm_bb,
    succ = "succ", total = "total",
    fixed = "Status_two + Tissue_primary",
    random_full   = "(1|ID) + (1|ID:Tissue_primary)",
    random_reduced= "(1|ID)",
    label = "CD103: IEC-EC vs Not IEC-EC"
  )
}

cd68_glmm_bb <- df7 %>%
  filter(Status_main %in% c("IEC-EC","Not IEC-EC"),
         is.finite(CD68), is.finite(total_cells), total_cells > 0) %>%
  mutate(Status_two = factor(Status_main, levels = c("Not IEC-EC","IEC-EC")),
         succ = as.integer(round(CD68)),
         total= as.integer(round(total_cells))) %>%
  mutate(succ = pmax(0L, pmin(succ, total)))
if (nrow(cd68_glmm_bb)>0 && n_distinct(cd68_glmm_bb$Status_two)==2) {
  run_bb_checks(
    dat = cd68_glmm_bb,
    succ = "succ", total = "total",
    fixed = "Status_two + Tissue_primary",
    random_full   = "(1|ID) + (1|ID:Tissue_primary)",
    random_reduced= "(1|ID)",
    label = "CD68: IEC-EC vs Not IEC-EC"
  )
}

cat("\n=== Diagnostics complete. ===\n")
