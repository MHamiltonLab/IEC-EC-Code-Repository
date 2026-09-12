#!/usr/bin/env Rscript
# Cohort characteristics and follow-up
# Summarize baseline characteristics, survival, competing risks, and follow-up.

suppressPackageStartupMessages({
  library(readr);  library(dplyr);  library(tidyr);  library(stringr)
  library(ggplot2); library(scales); library(forcats); library(broom)
  library(mgcv);   library(purrr)
  library(tableone); library(survminer); library(survival)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
path_df3 <- file.path(data_dir, "clinical_metadata.tsv")

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
    stringr::str_detect(low, "hispanic|latino|latinx") ~ "Hispanic/Latino",
    low %in% c("", "unknown", "na", "n/a") ~ "Unknown",
    TRUE ~ "Unknown"
  )
  factor(out, levels = c("Hispanic/Latino","Not Hispanic/Latino","Unknown"))
}

clean_race <- function(x){
  raw <- ifelse(is.na(x), "Unknown", as.character(x))
  raw <- str_replace_all(raw, "_", " ")
  raw <- str_squish(raw)

  raw <- str_to_title(raw)
  factor(raw)
}

clean_heavy_chain <- function(x){
  s <- as.character(x)
  s <- str_replace_all(s, "[^A-Za-z]", "")
  s <- toupper(s)
  out <- dplyr::case_when(
    s %in% c("IGG")   ~ "IgG",
    s %in% c("IGA")   ~ "IgA",
    s %in% c("IGM")   ~ "IgM",
    s %in% c("IGD")   ~ "IgD",
    s %in% c("IGE")   ~ "IgE",
    is.na(s) | s==""  ~ "Unknown",
    TRUE              ~ "Unknown"
  )
  factor(out, levels = c("IgG","IgA","IgM","IgD","IgE","Unknown"))
}

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra) })

out_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "clinical", "cohort_characteristics"), "figures")
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

save_survminer_pdf <- function(ggs, filename_base, width = 9, height = 9) {
  arranged_ok <- FALSE
  if (requireNamespace("survminer", quietly = TRUE)) {
    arranged <- try(survminer::arrange_ggsurvplots(list(ggs), print = FALSE), silent = TRUE)
    if (!inherits(arranged, "try-error")) {

      save_grob_pdf(arranged, paste0(filename_base, ".pdf"), width = width, height = height)
      arranged_ok <- TRUE
    }
  }
  if (!arranged_ok) {

    if (!is.null(ggs$plot))  save_pdf(ggs$plot,  paste0(filename_base, "_plot.pdf"),      width = width, height = height * 0.7)
    if (!is.null(ggs$table)) save_grob_pdf(ggs$table, paste0(filename_base, "_risktable.pdf"), width = width, height = height * 0.4)
  }
}

# Clinical metadata
df3 <- read_delim(path_df3, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
if (!"Study_ID" %in% names(df3) && "ID" %in% names(df3)) df3 <- dplyr::rename(df3, Study_ID = ID)

df3_h <- df3 %>%
  mutate(

    IEC_enteritis     = safe_num(IEC_enteritis),
    IEC_flag          = factor(if_else(IEC_enteritis == 1, "IEC", "No IEC"),
                               levels = c("No IEC","IEC")),

    Age_at_infusion   = safe_num(Age_at_infusion),
    Sex               = factor(Sex, levels = c("Male","Female")),
    ECOG_at_apheresis = safe_num(ECOG_at_apheresis),
    ECOG_group        = factor(if_else(is.finite(ECOG_at_apheresis) & ECOG_at_apheresis > 1, ">1", "0–1"),
                               levels = c("0–1",">1")),
    Prior_Lines       = safe_num(Prior_Lines),

    Race              = clean_race(Race),
    Ethnicity         = clean_ethnicity(Ethnicity),

    Heavy_Chain       = clean_heavy_chain(Heavy_Chain),
    Light_Chain       = factor(str_to_title(coalesce(Light_Chain, "Unknown")),
                               levels = c("Kappa","Lambda","Unknown")),
    RISS              = factor(str_to_upper(coalesce(RISS, "Unknown")),
                               levels = c("I","II","III","UNKNOWN")),

    PreLD_ALC         = safe_num(ALC_PreLD),
    M_spike           = safe_num(`preLD_M-spike`),
    Baseline_LDH      = safe_num(Baseline_pre_LD_LDH),
    Baseline_CRP      = safe_num(`Baseline_pre-LD_CRP`),
    Baseline_Ferritin = safe_num(`Baseline_pre-LD_Ferritin`),

    Max_ICANS         = safe_num(Max_ICANS),
    ICANS_gt0         = factor(if_else(is.finite(Max_ICANS) & Max_ICANS > 0, "Yes", "No"),
                               levels = c("No","Yes")),
    Max_CRS_grade     = safe_num(Max_CRS_grade),
    CRS_grp           = factor(if_else(is.finite(Max_CRS_grade) & Max_CRS_grade > 1, ">1", "0–1"),
                               levels = c("0–1",">1")),

    PD_day            = safe_num(PD_day),
    Death_day         = safe_num(Death_day),
    FU_day            = safe_num(Last_follow_up_day)
  )

cat("\nEthnicity counts (after harmonization):\n")
print(table(df3_h$Ethnicity, useNA = "ifany"))

cat("\nHeavy_Chain counts (normalized):\n")
print(table(df3_h$Heavy_Chain, useNA = "ifany"))

# Baseline characteristics
tbl_vars <- c(
  "Age_at_infusion","Sex","Race","Ethnicity",
  "ECOG_group","Prior_Lines","Heavy_Chain","Light_Chain","RISS",
  "PreLD_ALC","M_spike","Baseline_LDH","Baseline_CRP","Baseline_Ferritin",
  "ICANS_gt0","CRS_grp"
)

table_df   <- df3_h %>% dplyr::select(IEC_flag, dplyr::any_of(tbl_vars))

factorVars <- intersect(
  c("Sex","Race","Ethnicity","ECOG_group","Heavy_Chain","Light_Chain","RISS","ICANS_gt0","CRS_grp"),
  names(table_df)
)

nonnormal  <- intersect(
  c("Age_at_infusion","Prior_Lines","PreLD_ALC","M_spike",
    "Baseline_LDH","Baseline_CRP","Baseline_Ferritin"),
  names(table_df)
)

cat("\n=== Table 1 (baseline + demographics with corrected Ethnicity) ===\n")
tbl1 <- CreateTableOne(
  data = table_df,
  vars = setdiff(names(table_df), "IEC_flag"),
  strata = "IEC_flag",
  factorVars = factorVars,
  includeNA = TRUE
)

print(tbl1, smd = TRUE, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE)

# Progression-free survival
pfs_df <- df3_h %>%
  transmute(
    Study_ID, IEC_flag,
    pd = PD_day, death = Death_day, fu = FU_day
  ) %>%
  mutate(
    t_event = pmin(coalesce(pd, Inf), coalesce(death, Inf)),
    time    = pmin(t_event, coalesce(fu, Inf)),
    event   = as.integer(is.finite(t_event) & (t_event <= coalesce(fu, Inf)))
  ) %>%
  filter(is.finite(time), time >= 0)

sf_pfs_overall <- survfit(Surv(time, event) ~ 1, data = pfs_df)
sf_pfs_byIEC   <- survfit(Surv(time, event) ~ IEC_flag, data = pfs_df)

pal_IEC <- c("No IEC" = "#4C6A87", "IEC" = "#B22222")
iec_levels <- levels(droplevels(pfs_df$IEC_flag))
iec_cols   <- unname(pal_IEC[iec_levels])

p_km_overall <- survminer::ggsurvplot(
  fit = sf_pfs_overall,
  data = pfs_df,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.title = "Number at risk",
  surv.scale = "percent",
  xlab = "Days since infusion",
  ylab = "PFS (Kaplan–Meier)",
  ggtheme = theme_arial(22),
  legend = "none",
  break.time.by = 60
)
print(p_km_overall)

p_km_byIEC <- survminer::ggsurvplot(
  fit = sf_pfs_byIEC,
  data = pfs_df,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.28,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.title = "Number at risk",
  surv.scale = "percent",
  palette = iec_cols,
  legend.title = "IEC status",
  legend.labs = iec_levels,
  pval = TRUE,
  xlab = "Days since infusion",
  ylab = "PFS (Kaplan–Meier)",
  ggtheme = theme_arial(22),
  break.time.by = 60
)
print(p_km_byIEC)

if (requireNamespace("cmprsk", quietly = TRUE)) {
  ci <- cmprsk::cuminc(ftime = cr_df$time, fstatus = cr_df$status, cencode = 0)
  keys <- names(ci)
  key1 <- keys[grepl("^1(\\b|\\s|$)", keys)]
  if (length(key1) > 0) {
    nrm <- ci[[key1[1]]]
    nrm_df <- data.frame(time = nrm$time, cif = nrm$est) |>
      dplyr::distinct(time, .keep_all = TRUE)

    est_at <- function(tt) if (nrow(nrm_df)) approx(nrm_df$time, nrm_df$cif, xout = tt, rule = 2)$y else NA_real_
    cat(sprintf("\nNRM cumulative incidence (CRR): ~%.1f%% at 180 d; ~%.1f%% at 365 d\n",
                100*est_at(180), 100*est_at(365)))

    p_nrm <- ggplot(nrm_df, aes(time, cif)) +
      geom_step(linewidth = 1.2, color = "#b22222") +
      scale_y_continuous("NRM cumulative incidence", labels = scales::percent_format(accuracy = 0.1), limits = c(0,1)) +
      scale_x_continuous("Days since infusion") +
      ggtitle("Non-relapse mortality (competing risks: PD as competing event)") +
      theme_arial(22)
    print(p_nrm)
  } else {
    message("No NRM (type-1) deaths observed; skipping NRM curve.")
  }
} else {
  message("Package 'cmprsk' not installed; skipping competing-risk NRM.")
}

suppressPackageStartupMessages({ library(survival); library(survminer) })

# Overall survival
os_df <- df3_h %>%
  dplyr::transmute(
    Study_ID, IEC_flag,
    death = Death_day,
    fu    = FU_day
  ) %>%
  dplyr::mutate(
    time  = pmin(coalesce(death, Inf), coalesce(fu, Inf)),
    event = as.integer(is.finite(death) & (death <= coalesce(fu, Inf)))
  ) %>%
  dplyr::filter(is.finite(time), time >= 0) %>%
  droplevels()

iec_levels <- levels(droplevels(os_df$IEC_flag))
iec_cols   <- unname(pal_IEC[iec_levels])

sf_os_overall <- survfit(Surv(time, event) ~ 1, data = os_df)
sf_os_byIEC   <- survfit(Surv(time, event) ~ IEC_flag, data = os_df)

p_os_overall <- ggsurvplot(
  fit = sf_os_overall, data = os_df,
  conf.int = TRUE, risk.table = TRUE,
  risk.table.height = 0.25, risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE, risk.table.title = "Number at risk",
  surv.scale = "percent",
  xlab = "Days since infusion", ylab = "Overall survival (KM)",
  ggtheme = theme_arial(22), legend = "none", break.time.by = 60
)
print(p_os_overall)

p_os_byIEC <- ggsurvplot(
  fit = sf_os_byIEC, data = os_df,
  conf.int = F, risk.table = TRUE, risk.table.col = "strata",
  risk.table.height = 0.28, risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE, risk.table.title = "Number at risk",
  surv.scale = "percent",
  palette = iec_cols, legend.title = "IEC status", legend.labs = iec_levels,
  pval = TRUE,
  xlab = "Days since infusion", ylab = "Overall survival (KM)",
  ggtheme = theme_arial(22), break.time.by = 60
)
print(p_os_byIEC)

km_at <- function(fit, t_days = 730) {
  s <- summary(fit, times = t_days)
  tibble::tibble(
    strata = if (is.null(s$strata)) "(overall)" else sub("^.*IEC_flag=", "", s$strata),
    surv   = s$surv,
    lower  = s$lower,
    upper  = s$upper
  )
}
cat("\n=== OS @ 2 years (730 d) ===\n")
os2_overall <- km_at(sf_os_overall, 730)
os2_byIEC   <- km_at(sf_os_byIEC,   730)
print(dplyr::bind_rows(os2_overall, os2_byIEC) %>%
        dplyr::mutate(
          OS_2y = scales::percent(surv, accuracy = 0.1),
          LCL   = scales::percent(lower, accuracy = 0.1),
          UCL   = scales::percent(upper, accuracy = 0.1)
        ) %>%
        dplyr::select(strata, OS_2y, LCL, UCL))

# Competing-risk outcomes
cr_df <- df3_h %>%
  dplyr::transmute(
    Study_ID, IEC_flag,
    pd    = PD_day,
    death = Death_day,
    fu    = FU_day
  ) %>%
  dplyr::mutate(
    first_event = pmin(coalesce(pd, Inf), coalesce(death, Inf), coalesce(fu, Inf)),
    status = dplyr::case_when(
      is.finite(death) & (death < coalesce(pd, Inf)) & (death <= coalesce(fu, Inf)) ~ 1L,
      is.finite(pd)    & (pd <= coalesce(death, Inf)) & (pd <= coalesce(fu, Inf))   ~ 2L,
      TRUE ~ 0L
    ),
    time = first_event
  ) %>%
  dplyr::filter(is.finite(time), time >= 0) %>%
  droplevels()

if (requireNamespace("cmprsk", quietly = TRUE)) {
  library(cmprsk)

  ci_overall <- cuminc(ftime = cr_df$time, fstatus = cr_df$status, cencode = 0)

  ci_byIEC   <- cuminc(ftime = cr_df$time, fstatus = cr_df$status,
                       group = cr_df$IEC_flag, cencode = 0)

  cif_at <- function(ci_obj, t_days = 730, event_code = 1L) {
    nms <- names(ci_obj)
    nms <- nms[!grepl("^Tests", nms)]
    keep <- grep(paste0("\\b", event_code, "$"), nms, value = TRUE)
    dplyr::bind_rows(lapply(keep, function(k) {
      est <- ci_obj[[k]]

      df <- dplyr::distinct(data.frame(time = est$time, cif = est$est), time, .keep_all = TRUE)
      tibble::tibble(
        curve = k,
        cif   = if (nrow(df)) approx(df$time, df$cif, xout = t_days, rule = 2)$y else NA_real_
      )
    }))
  }

  nrm2_overall <- cif_at(ci_overall, 730) %>%
    dplyr::mutate(strata = "(overall)") %>%
    dplyr::select(strata, cif)
  nrm2_byIEC   <- cif_at(ci_byIEC, 730) %>%
    dplyr::mutate(strata = sub("\\s+1$", "", curve),
                  strata = sub("^.*=", "", strata)) %>%
    dplyr::select(strata, cif)

  cat("\n=== NRM cumulative incidence @ 2 years (730 d) ===\n")
  print(dplyr::bind_rows(nrm2_overall, nrm2_byIEC) %>%
          dplyr::mutate(NRM_2y = scales::percent(cif, accuracy = 0.1)) %>%
          dplyr::select(strata, NRM_2y))

  k_overall <- names(ci_overall)[grep("\\b1$", names(ci_overall))][1]
  if (!is.na(k_overall)) {
    nrm <- ci_overall[[k_overall]]
    nrm_df <- dplyr::distinct(data.frame(time = nrm$time, cif = nrm$est), time, .keep_all = TRUE)
    p_nrm_overall <- ggplot(nrm_df, aes(time, cif)) +
      geom_step(linewidth = 1.2, color = "#b22222") +
      scale_y_continuous("NRM cumulative incidence",
                         labels = scales::percent_format(accuracy = 0.1), limits = c(0, 1)) +
      scale_x_continuous("Days since infusion") +
      ggtitle("NRM (overall, PD competing)") +
      theme_arial(22)
    print(p_nrm_overall)
  }

  nms <- names(ci_byIEC); nms <- nms[grep("\\b1$", nms)]
  if (length(nms)) {
    nrm_by <- dplyr::bind_rows(lapply(nms, function(k) {
      est <- ci_byIEC[[k]]
      df  <- dplyr::distinct(data.frame(time = est$time, cif = est$est), time, .keep_all = TRUE)
      df$IEC_flag <- sub("\\s+1$", "", k)
      df$IEC_flag <- sub("^.*=", "", df$IEC_flag)
      df
    }))
    p_nrm_byIEC <- ggplot(nrm_by, aes(time, cif, color = IEC_flag)) +
      geom_step(linewidth = 1.2) +
      scale_color_manual(values = pal_IEC, guide = "none") +
      scale_y_continuous("NRM cumulative incidence",
                         labels = scales::percent_format(accuracy = 0.1), limits = c(0, 1)) +
      scale_x_continuous("Days since infusion") +
      ggtitle("NRM by IEC status (PD competing)") +
      facet_wrap(~ IEC_flag) +
      theme_arial(22)
    print(p_nrm_byIEC)
  }
} else {
  message("Package 'cmprsk' not installed; skipping competing-risk NRM.")
}

t6 <- 183

cr6 <- df3 %>%
  transmute(
    Study_ID,
    IEC_day   = safe_num(IEC_Day),
    PD_day    = safe_num(PD_day),
    Death_day = safe_num(Death_day),
    FU_day    = safe_num(Last_follow_up_day)
  ) %>%
  mutate(

    c6 = if_else(!is.na(FU_day), pmin(FU_day, t6), t6),

    e1 = if_else(!is.na(IEC_day)   & IEC_day   <= c6, IEC_day,   Inf),
    e2 = if_else(!is.na(PD_day)    & PD_day    <= c6, PD_day,    Inf),
    e3 = if_else(!is.na(Death_day) & Death_day <= c6, Death_day, Inf),

    time   = pmin(e1, e2, e3, c6, na.rm = TRUE),
    status = dplyr::case_when(
      is.finite(e1) & time == e1 ~ 1L,
      is.finite(e2) & time == e2 ~ 2L,
      is.finite(e3) & time == e3 ~ 3L,
      TRUE                       ~ 0L
    )
  ) %>%
  filter(!is.na(time), !is.na(status)) %>%
  mutate(
    time   = pmax(0, as.numeric(time)),
    status = as.integer(status)
  )

cat("\nEvent counts within 6 months (1=IEC, 2=PD, 3=Death, 0=Censor):\n")
print(table(cr6$status))

ci <- with(cr6, cmprsk::cuminc(ftime = time, fstatus = status, cencode = 0))
nm <- names(ci)
ix1 <- if (any(nm == "1")) which(nm == "1") else grep("^1", nm)
c1 <- ci[[ix1[1]]]

iec_df <- tibble::tibble(time = c1$time, est = c1$est) %>%
  mutate(
    se  = if (!is.null(c1$var)) sqrt(c1$var) else 0,
    lcl = pmax(0, est - 1.96*se),
    ucl = pmin(1, est + 1.96*se)
  )

ix6  <- max(which(iec_df$time <= t6), na.rm = TRUE)
cif6 <- if (is.finite(ix6)) iec_df$est[ix6] else 0
lcl6 <- if (is.finite(ix6)) iec_df$lcl[ix6] else 0
ucl6 <- if (is.finite(ix6)) iec_df$ucl[ix6] else 0

cat(sprintf("\nIEC cumulative incidence at 6 months: %.1f%% (95%% CI %.1f–%.1f%%)\n",
            100*cif6, 100*lcl6, 100*ucl6))

p_cif <- ggplot(iec_df, aes(x = time, y = est)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "#B22222", alpha = 0.12) +
  geom_step(linewidth = 1.1, color = "#B22222") +
  geom_vline(xintercept = t6, linetype = "dashed") +
  annotate("point", x = t6, y = cif6, color = "#B22222", size = 3) +
  annotate(
    "label", x = t6, y = cif6,
    label = sprintf("6-mo IEC: %.1f%% (95%% CI %.1f–%.1f%%)", 100*cif6, 100*lcl6, 100*ucl6),
    hjust = 1, vjust = -0.25, size = 5, family = "Arial", label.size = 0
  ) +
  scale_x_continuous("Days from infusion (capped at 6 months)", limits = c(0, t6)) +
  scale_y_continuous("IEC cumulative incidence", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  ggtitle("Cumulative incidence of IEC with competing risks (PD, Death) at 6 months") +
  theme_classic(base_size = 22) +
  theme(
    text       = element_text(family = "Arial"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20)
  )

print(p_cif)

suppressPackageStartupMessages({ library(survival); library(survminer); library(dplyr); library(tibble) })

if (!exists("os_df")) {
  os_df <- df3_h %>%
    transmute(
      Study_ID, IEC_flag,
      death = Death_day,
      fu    = FU_day
    ) %>%
    mutate(
      time  = pmin(coalesce(death, Inf), coalesce(fu, Inf)),
      event = as.integer(is.finite(death) & (death <= coalesce(fu, Inf)))
    ) %>%
    filter(is.finite(time), time >= 0) %>%
    droplevels()
}

rkm_overall <- survfit(Surv(time, 1 - event) ~ 1, data = os_df)
rkm_byIEC   <- survfit(Surv(time, 1 - event) ~ IEC_flag, data = os_df)

.rkm_extract <- function(fit){
  tb <- summary(fit)$table
  if (is.matrix(tb) || is.data.frame(tb)) {
    df <- as.data.frame(tb)
    df$strata <- rownames(df)
    df %>%
      mutate(
        strata = sub("^.*=", "", strata),
        median_days = as.numeric(df[,"median"]),
        lcl_days    = as.numeric(df[,grep("LCL", colnames(df), value = TRUE)[1]]),
        ucl_days    = as.numeric(df[,grep("UCL", colnames(df), value = TRUE)[1]])
      ) %>%
      select(strata, median_days, lcl_days, ucl_days)
  } else {
    tibble(
      strata = "(overall)",
      median_days = as.numeric(tb["median"]),
      lcl_days    = as.numeric(tb["0.95LCL"]),
      ucl_days    = as.numeric(tb["0.95UCL"])
    )
  }
}

rkm_overall_tbl <- .rkm_extract(rkm_overall)
rkm_byIEC_tbl   <- .rkm_extract(rkm_byIEC)

if (nrow(rkm_overall_tbl) == 1) {
  md  <- rkm_overall_tbl$median_days[1]
  lcl <- rkm_overall_tbl$lcl_days[1]
  ucl <- rkm_overall_tbl$ucl_days[1]
  cat(sprintf("\n=== Median follow-up (Reverse KM, overall) ===\n%.1f days (95%% CI %.1f–%.1f)\n",
              md, lcl, ucl))
}

cat("\n=== Median follow-up (Reverse KM) by IEC status ===\n")
print(
  bind_rows(
    rkm_overall_tbl %>% mutate(strata = "(overall)"),
    rkm_byIEC_tbl
  ) %>%
    mutate(
      median_months = median_days/30.44,
      lcl_months    = lcl_days/30.44,
      ucl_months    = ucl_days/30.44
    ) %>%
    transmute(
      strata,
      median_days   = round(median_days, 1),
      `95%CI_days`  = sprintf("%.1f–%.1f", lcl_days, ucl_days),
      median_months = sprintf("%.1f", median_months),
      `95%CI_months`= sprintf("%.1f–%.1f", lcl_months, ucl_months)
    )
)

p_rkm_overall <- ggsurvplot(
  fit = rkm_overall,
  conf.int = TRUE,
  risk.table = FALSE,
  surv.scale = "percent",
  xlab = "Days since infusion",
  ylab = "Proportion with ≥ this follow-up (Reverse KM)",
  ggtheme = ggplot2::theme_classic(base_size = 22),
  break.time.by = 60
)
print(p_rkm_overall)

tbl1_df <- print(tbl1, smd = TRUE, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)
tbl1_path_csv <- file.path(out_dir, "Table1_IEC_demographics_baseline.csv")
readr::write_csv(as.data.frame(tbl1_df), tbl1_path_csv)
message("Saved: ", tbl1_path_csv)

if (exists("p_nrm") && inherits(p_nrm, "ggplot")) {
  save_pdf(p_nrm, "NRM_overall_quick.pdf", width = 9, height = 7)
}

save_survminer_pdf(p_os_overall, "KM_OS_overall")
save_survminer_pdf(p_os_byIEC,   "KM_OS_byIEC")

if (exists("p_nrm_overall") && inherits(p_nrm_overall, "ggplot")) {
  save_pdf(p_nrm_overall, "CIF_NRM_overall.pdf", width = 9, height = 7)
}

if (exists("p_nrm_byIEC") && inherits(p_nrm_byIEC, "ggplot")) {
  save_pdf(p_nrm_byIEC, "CIF_NRM_byIEC.pdf", width = 11, height = 8.5)
}

save_pdf(p_cif, "CIF_IEC_6mo.pdf", width = 11, height = 8.5)

save_survminer_pdf(p_rkm_overall, "ReverseKM_followup_overall")

suppressPackageStartupMessages({ library(survival); library(dplyr); library(readr); library(tibble) })

stats_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "clinical", "cohort_characteristics"), "stats")
if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

.km_median_extract <- function(fit) {
  tb <- summary(fit)$table
  if (is.matrix(tb) || is.data.frame(tb)) {
    df <- as.data.frame(tb)
    df$strata <- rownames(df)
    df %>%
      mutate(
        strata      = sub("^.*=", "", strata),
        median_days = suppressWarnings(as.numeric(df[, "median"])),
        lcl_days    = suppressWarnings(as.numeric(df[, grep("LCL", colnames(df), value = TRUE)[1]])),
        ucl_days    = suppressWarnings(as.numeric(df[, grep("UCL", colnames(df), value = TRUE)[1]]))
      ) %>%
      select(strata, median_days, lcl_days, ucl_days)
  } else {
    tibble(
      strata      = "(overall)",
      median_days = suppressWarnings(as.numeric(tb["median"])),
      lcl_days    = suppressWarnings(as.numeric(tb[grep("LCL", names(tb))[1]])),
      ucl_days    = suppressWarnings(as.numeric(tb[grep("UCL", names(tb))[1]]))
    )
  }
}

.km_at <- function(fit, t_days = 730) {
  s <- summary(fit, times = t_days)
  tibble(
    strata = if (is.null(s$strata)) "(overall)" else sub("^.*IEC_flag=", "", s$strata),
    surv   = s$surv,
    lower  = s$lower,
    upper  = s$upper
  )
}

if (!exists("rkm_overall") || !exists("rkm_byIEC")) {
  rkm_overall <- survfit(Surv(time, 1 - event) ~ 1, data = os_df)
  rkm_byIEC   <- survfit(Surv(time, 1 - event) ~ IEC_flag, data = os_df)
}

rkm_overall_tbl <- .km_median_extract(rkm_overall) %>% mutate(strata = "(overall)")
rkm_byIEC_tbl   <- .km_median_extract(rkm_byIEC)

followup_tbl <- bind_rows(rkm_overall_tbl, rkm_byIEC_tbl) %>%
  mutate(
    median_months = median_days/30.44,
    lcl_months    = lcl_days/30.44,
    ucl_months    = ucl_days/30.44,
    median_fmt    = ifelse(is.na(median_days), "NR", sprintf("%.1f d (%.1f–%.1f)", median_days, lcl_days, ucl_days)),
    median_mo     = ifelse(is.na(median_days), "NR", sprintf("%.1f mo (%.1f–%.1f)", median_months, lcl_months, ucl_months))
  )

cat("\n=== Median follow-up (Reverse KM) ===\n")
followup_tbl %>%
  transmute(strata, `Median follow-up` = median_fmt, `Median follow-up (months)` = median_mo) %>%
  print(n = Inf)

write_tsv(followup_tbl, file.path(stats_dir, "followup_median_reverseKM.tsv"))

pfs_med_overall <- .km_median_extract(sf_pfs_overall) %>% mutate(strata = "(overall)")
pfs_med_byIEC   <- .km_median_extract(sf_pfs_byIEC)

pfs_median_tbl <- bind_rows(pfs_med_overall, pfs_med_byIEC) %>%
  mutate(
    median_months = median_days/30.44,
    median_fmt    = ifelse(is.na(median_days), "NR", sprintf("%.1f d (%.1f–%.1f)", median_days, lcl_days, ucl_days)),
    median_mo     = ifelse(is.na(median_days), "NR", sprintf("%.1f mo (%.1f–%.1f)", median_months, lcl_days/30.44, ucl_days/30.44))
  )

cat("\n=== Median PFS (KM) ===\n")
pfs_median_tbl %>%
  transmute(strata, `Median PFS` = median_fmt, `Median PFS (months)` = median_mo) %>%
  print(n = Inf)

write_tsv(pfs_median_tbl, file.path(stats_dir, "pfs_median_time.tsv"))

pfs2_overall <- .km_at(sf_pfs_overall, t_days = 730)
pfs2_byIEC   <- .km_at(sf_pfs_byIEC,   t_days = 730)

pfs2_tbl <- bind_rows(pfs2_overall, pfs2_byIEC) %>%
  mutate(
    PFS_2y = scales::percent(surv,  accuracy = 0.1),
    LCL    = scales::percent(lower, accuracy = 0.1),
    UCL    = scales::percent(upper, accuracy = 0.1)
  ) %>%
  transmute(strata, PFS_2y, `95% CI` = paste0(LCL, "–", UCL))

cat("\n=== PFS at 2 years (730 d) ===\n")
print(pfs2_tbl, n = Inf)

write_tsv(pfs2_tbl, file.path(stats_dir, "pfs_at_2years.tsv"))
