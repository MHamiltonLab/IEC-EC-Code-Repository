#!/usr/bin/env Rscript
# CD3-normalized tissue abundance
# Quantify marker abundance and compartment phenotypes using beta-binomial mixed models.

suppressPackageStartupMessages({
  library(tidyverse)
  library(forcats)
  library(glmmTMB)
  library(emmeans)
  library(lme4)
  library(broom.mixed)
})

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto(enable = FALSE)
}
BASE_FAMILY <- "Helvetica"
ggplot2::theme_set(
  ggplot2::theme_classic(base_family = BASE_FAMILY, base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14, color = "grey20"),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13),
      strip.text = ggplot2::element_text(size = 14, face = "bold")
    )
)

args <- commandArgs(trailingOnly = TRUE)
default_data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
data_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_data_dir
data_dir <- normalizePath(path.expand(data_dir), mustWork = TRUE)
path_df7 <- file.path(data_dir, "tissue_cell_counts.tsv")

# Abundance denominator
ABUNDANCE_DENOMINATOR <- "CD3"
ABUNDANCE_DENOMINATOR_LABEL <- "CD3+ T cells"
ABUNDANCE_Y_LABEL <- "% Camelid+ among CD3+ T cells"
ABUNDANCE_REFERENCE_LINES <- c(5, 10, 20, 50)
if (!ABUNDANCE_DENOMINATOR %in% c("CD3", "total_cells")) {
  stop("ABUNDANCE_DENOMINATOR must be 'CD3' or 'total_cells'.")
}

output_root <- Sys.getenv("METHODS_OUTPUT_DIR", unset = "results")
output_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "tissue_imaging", "cd3_normalized_abundance")
fig_dir     <- file.path(output_dir, "figures")
stats_dir   <- file.path(output_dir, "stats")
log_dir     <- file.path(output_dir, "logs")
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,   recursive = TRUE, showWarnings = FALSE)

if (!file.exists(path_df7)) {
  stop("COMET input file not found: ", path_df7)
}
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

.std_ci <- function(df) {
  nc <- names(df)
  if (all(c("lower.CL","upper.CL") %in% nc)) return(df)
  if (all(c("asymp.LCL","asymp.UCL") %in% nc)) return(dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL))
  if (all(c("LCL","UCL") %in% nc)) return(dplyr::rename(df, lower.CL = LCL, upper.CL = UCL))
  df
}

status_display <- function(x) {
  dplyr::case_when(
    x == "IEC-EC"             ~ "Camelid Positive",
    x == "Not IEC-EC"         ~ "Camelid Negative",
    x == "Untreated Controls" ~ "Untreated Controls",
    TRUE ~ x
  )
}
TITLE_IECEC    <- "Camelid Positive"

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

xlab_with_n <- function(dat, xvar, nvar = NULL) {
  if (is.null(nvar)) {
    counts <- dat %>% dplyr::count(.data[[xvar]], name = "n")
  } else {
    counts <- dat %>%
      dplyr::group_by(.data[[xvar]]) %>%
      dplyr::summarise(n = sum(.data[[nvar]], na.rm = TRUE), .groups = "drop")
  }
  function(x) {
    sapply(x, function(k) {
      n <- counts$n[counts[[xvar]] == k]
      if (length(n)) paste0(k,"\n","n = ",n) else k
    })
  }
}

pal_primary_full <- c("Colon"="#8da0cb","Duodenum"="#66c2a5","Terminal Ileum"="#fc8d62","Other"="#bdbdbd")
scale_color_tissue_primary <- function(name=NULL) scale_color_manual(values = pal_primary_full, drop = TRUE, name = name)
scale_fill_tissue_primary  <- function(name=NULL) scale_fill_manual(values  = pal_primary_full, drop = TRUE, name = name)

render_to_device <- function(open_device, draw) {
  open_device()
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(NULL)
}

save_plot_all <- function(filename, plot, width = 7, height = 5, dpi = 300, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  stem <- sub("\\.(png|pdf|eps)$", "", filename, ignore.case = TRUE)
  file_png <- paste0(stem, ".png")
  file_pdf <- paste0(stem, ".pdf")
  file_eps <- paste0(stem, ".eps")

  ggplot2::ggsave(
    file_png, plot = plot, width = width, height = height, dpi = dpi,
    bg = "white", ...
  )
  render_to_device(
    function() grDevices::pdf(
      file_pdf, width = width, height = height, family = BASE_FAMILY,
      useDingbats = FALSE, version = "1.4", colormodel = "srgb"
    ),
    function() print(plot)
  )
  render_to_device(
    function() grDevices::postscript(
      file_eps, width = width, height = height, family = BASE_FAMILY,
      onefile = FALSE, horizontal = FALSE, paper = "special",
      colormodel = "srgb"
    ),
    function() print(plot)
  )
  invisible(c(file_png, file_pdf, file_eps))
}

# Cell-count harmonization
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
    non_KI67              = `CD3-Camelid-negative-KI67`,
    non_FOXP3             = `CD3-Camelid-negative-FOXP3`,
    non_PD1               = `CD3-Camelid-negative-PD1`,

    pct_cam_CD4           = `pct_CD3-Camelid-CD4`         / 100,
    pct_cam_CD8           = `pct-CD3-Camelid-CD8`         / 100,
    pct_cam_DP            = `pct-CD3-Camelid-DP`          / 100,
    pct_cam_DN            = `pct-CD3-Camelid-DN`          / 100,
    pct_cam_GZMB          = `pct-CD3-Camelid-GZMB`        / 100,
    pct_cam_KI67          = `pct-CD3-Camelid-KI67`        / 100,
    pct_cam_PD1           = `pct-CD3-Camelid-PD1`         / 100,

    frac_non_camelid_cd3  = `pct-CD3-Camelid-negative`,

    pct_non_CD4           = `pct_CD3-Camelid-negative-CD4`   / 100,
    pct_non_CD8           = `pct_CD3-Camelid-negative-CD8`   / 100,
    pct_non_DP            = `pct_CD3-Camelid-negative_DP`    / 100,
    pct_non_DN            = `pct_CD3-Camelid-negative_DN`    / 100,
    pct_non_GZMB          = `pct_CD3-Camelid-negative-GZMB`  / 100,
    pct_non_KI67          = `pct_CD3-Camelid-negative-KI67`  / 100,
    pct_non_FOXP3         = `pct_CD3-Camelid-negative-FOXP3` / 100,
    pct_non_PD1           = `pct_CD3-Camelid-negative-PD1`   / 100,

    CD3_fraction_camelid,
    CD4_CD8_ratio         = `CD4:CD8`,
    CD103_fraction_tcells,
    Total_cell_fraction_camelid
  ) %>%
  mutate(
    Status_main    = status_label(Cohort),
    Tissue_fine    = tissue_fine(Tissue_Type),
    Tissue_primary = tissue_primary_from_fine(Tissue_fine),
    abundance_total = if (identical(ABUNDANCE_DENOMINATOR, "CD3")) {
      as.numeric(CD3)
    } else {
      as.numeric(total_cells)
    },

    gi_frac_cd3_camelid   = frac_safe(CD3_Camelid, CD3),
    gi_frac_total_camelid = frac_safe(CD3_Camelid, total_cells),
    frac_camelid_abundance = frac_safe(CD3_Camelid, abundance_total),
    frac_cd103_abundance   = frac_safe(CD3_CD103, abundance_total),
    frac_CD68_total       = frac_safe(CD68, total_cells),
    frac_non_camelid_cd3_from_counts = frac_safe(CD3_nonCamelid, CD3),

    tot_CD4   = cam_CD4 + non_CD4,
    tot_CD8   = cam_CD8 + non_CD8,
    tot_DP    = cam_DP  + non_DP,
    tot_DN    = cam_DN  + non_DN,
    tot_GZMB  = cam_GZMB + non_GZMB,

    cam_KI67  = pct_cam_KI67 * CD3_Camelid,
    tot_KI67  = cam_KI67 + non_KI67,

    tot_FOXP3 = non_FOXP3,

    tot_PD1   = cam_PD1 + non_PD1,

    frac_tot_CD4   = frac_safe(tot_CD4,   CD3),
    frac_tot_CD8   = frac_safe(tot_CD8,   CD3),
    frac_tot_DP    = frac_safe(tot_DP,    CD3),
    frac_tot_DN    = frac_safe(tot_DN,    CD3),
    frac_gzmb_abundance = frac_safe(tot_GZMB, abundance_total),
    frac_tot_KI67  = frac_safe(tot_KI67,  CD3),
    frac_foxp3_abundance = frac_safe(tot_FOXP3, abundance_total),
    frac_tot_PD1   = frac_safe(tot_PD1,   CD3),

    frac_non_FOXP3 = frac_safe(non_FOXP3, CD3_nonCamelid),
    frac_non_KI67  = frac_safe(non_KI67,  CD3_nonCamelid)
  )

# Count and denominator audit
denominator_audit <- df7 %>%
  dplyr::summarise(
    denominator = ABUNDANCE_DENOMINATOR,
    n_rows = dplyr::n(),
    n_valid_denominator = sum(is.finite(abundance_total) & abundance_total > 0),
    n_missing_denominator = sum(!is.finite(abundance_total)),
    n_nonpositive_denominator = sum(is.finite(abundance_total) & abundance_total <= 0),
    n_camelid_gt_denominator = sum(CD3_Camelid > abundance_total, na.rm = TRUE),
    n_cd103_gt_denominator = sum(CD3_CD103 > abundance_total, na.rm = TRUE),
    n_gzmb_gt_denominator = sum(tot_GZMB > abundance_total, na.rm = TRUE),
    n_foxp3_gt_denominator = sum(tot_FOXP3 > abundance_total, na.rm = TRUE)
  )
readr::write_tsv(
  denominator_audit,
  file.path(stats_dir, "abundance_denominator_audit.tsv")
)
if (denominator_audit$n_valid_denominator == 0) {
  stop("No rows have a valid positive abundance denominator.")
}
if (sum(unlist(denominator_audit[c(
  "n_camelid_gt_denominator", "n_cd103_gt_denominator",
  "n_gzmb_gt_denominator", "n_foxp3_gt_denominator"
)], use.names = FALSE), na.rm = TRUE) > 0) {
  warning(
    "One or more marker counts exceed the configured denominator. ",
    "Inspect abundance_denominator_audit.tsv; model counts are bounded to the denominator."
  )
}

iec_ec    <- df7 %>% dplyr::filter(Status_main == "IEC-EC")

cat("\n--- df7 summary (Status_main x Tissue_primary) ---\n")
print(df7 %>% dplyr::count(Status_main, Tissue_primary))

counts_status <- df7 %>%
  dplyr::filter(Status_main %in% c("Untreated Controls","Not IEC-EC","IEC-EC")) %>%
  mutate(Status_plot = factor(status_display(Status_main),
                              levels=c("Untreated Controls","Camelid Negative","Camelid Positive"))) %>%
  dplyr::count(Status_plot, name="n")

p_counts <- ggplot(counts_status, aes(x=Status_plot, y=n, fill=Status_plot)) +
  geom_col(width=0.7) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  scale_fill_manual(values=c("Untreated Controls"="#bdbdbd","Camelid Negative"="#8da0cb","Camelid Positive"="#66c2a5"), guide="none") +
  labs(x=NULL, y="Samples", title="Samples by status (COMET)") +
  scale_x_discrete(labels = xlab_with_n(counts_status, "Status_plot", nvar = "n"))
save_plot_all(file.path(fig_dir,"fig_samples_by_status.pdf"), p_counts, width=7.2, height=4.6, dpi=300)
readr::write_tsv(counts_status, file.path(stats_dir,"samples_by_status.tsv"))

camelid_status <- df7 %>%
  dplyr::filter(Status_main %in% c("Untreated Controls","Not IEC-EC","IEC-EC"),
         is.finite(frac_camelid_abundance)) %>%
  mutate(Status_plot = factor(status_display(Status_main),
                              levels=c("Untreated Controls","Camelid Negative","Camelid Positive")),
         pct_camelid_abundance = 100*frac_camelid_abundance)

p_camelid_status <- ggplot(camelid_status,
                           aes(x=Status_plot, y=pct_camelid_abundance, color=Tissue_primary)) +
  geom_jitter(width=0.15, height=0, alpha=0.9, size=2.6) +
  stat_summary(fun=median, geom="crossbar", width=0.6, linewidth=0.55, color="black", fill=NA) +
  scale_color_tissue_primary("Tissue") +
  labs(x=NULL, y=ABUNDANCE_Y_LABEL, title="Camelid infiltration by status") +
  scale_x_discrete(labels = xlab_with_n(camelid_status,"Status_plot"))
save_plot_all(file.path(fig_dir,"fig_camelid_status_IECEC_NotIEC_Controls.pdf"),
       p_camelid_status, width=8.8, height=5.8, dpi=300)

camelid_summary <- camelid_status %>%
  group_by(Status_plot) %>%
  summarise(n=n(), median_pct=median(pct_camelid_abundance), IQR_low=quantile(pct_camelid_abundance,0.25),
            IQR_high=quantile(pct_camelid_abundance,0.75), .groups="drop")
readr::write_tsv(camelid_summary, file.path(stats_dir,"stats_camelid_status.tsv"))

colon_set   <- c("Colon","Colon (left)","Colon (right)","Colon (random)","Ileocecal Valve")
stomach_set <- c("Stomach","Stomach (antrum)","Stomach (polyp)")

iec_sites_lite <- iec_ec %>%
  dplyr::filter(is.finite(frac_camelid_abundance)) %>%
  mutate(
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
                         levels=c("Esophagus","Stomach","Duodenum","Terminal Ileum","Colon","Rectum","Other")),
    pct_camelid_abundance  = 100 * frac_camelid_abundance,
    Site_color_group = tissue_primary_from_fine(Site_simple)
  )

p_sites_lite <- ggplot(iec_sites_lite,
                       aes(x=Site_simple, y=pct_camelid_abundance, color=Site_color_group)) +
  geom_point(size=2.2, alpha=0.95, position=position_jitter(width=0.15, height=0)) +
  stat_summary(fun=median, geom="crossbar", width=0.6, linewidth=0.55, color="black", fill=NA) +
  scale_color_tissue_primary("Tissue") +
  labs(x=NULL, y=ABUNDANCE_Y_LABEL, title=paste0("Camelid infiltration by GI site (", TITLE_IECEC, ")")) +
  scale_x_discrete(labels = xlab_with_n(iec_sites_lite,"Site_simple"))
save_plot_all(file.path(fig_dir,"fig_camelid_by_site_IECEC_COLON_COLLAPSED.pdf"),
       p_sites_lite, width=8.8, height=5.8, dpi=300)

p_sites_lite_box <- ggplot(iec_sites_lite,
                           aes(x=Site_simple, y=pct_camelid_abundance, fill=Site_color_group, color=Site_color_group)) +
  geom_boxplot(width=0.65, alpha=0.25, outlier.shape=NA) +
  geom_jitter(width=0.15, height=0, alpha=0.95, size=2.2) +
  scale_color_tissue_primary("Tissue") + scale_fill_tissue_primary("Tissue") +
  labs(x=NULL, y=ABUNDANCE_Y_LABEL,
       title=paste0("Camelid infiltration by GI site (", TITLE_IECEC, ") — boxplot")) +
  scale_x_discrete(labels = xlab_with_n(iec_sites_lite,"Site_simple"))
save_plot_all(file.path(fig_dir,"fig_camelid_by_site_IECEC_COLON_COLLAPSED_boxplot.pdf"),
       p_sites_lite_box, width=8.8, height=5.8, dpi=300)

site_summary <- iec_sites_lite %>%
  group_by(Site_simple) %>%
  summarise(n=n(), median_pct=median(pct_camelid_abundance), IQR_low=quantile(pct_camelid_abundance,0.25),
            IQR_high=quantile(pct_camelid_abundance,0.75), .groups="drop")
readr::write_tsv(site_summary, file.path(stats_dir,"stats_camelid_by_site_IECEC.tsv"))

iec_sites_glmm <- iec_ec %>%
  dplyr::filter(
    is.finite(CD3_Camelid),
    is.finite(abundance_total),
    abundance_total > 0
  ) %>%
  mutate(
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
                         levels=c("Esophagus","Stomach","Duodenum","Terminal Ileum","Colon","Rectum","Other")),
    total_count = as.integer(round(abundance_total)),
    succ = pmax(0L, pmin(as.integer(round(CD3_Camelid)), total_count)),
    fail = pmax(0L, total_count - succ)
  ) %>%
  dplyr::filter(!is.na(Site_simple))

iec_sites_glmm_ref <- iec_sites_glmm %>%
  mutate(Site_simple = stats::relevel(Site_simple, ref = "Duodenum"))

# Tissue-site comparison
fit_site_bb <- glmmTMB(cbind(succ, fail) ~ Site_simple + (1 | ID),
                       family = betabinomial(link = "logit"),
                       data   = iec_sites_glmm_ref)

emm_site_bb_link <- emmeans(fit_site_bb, ~ Site_simple)
emm_site_bb_resp <- summary(emm_site_bb_link, type="response") %>%
  as.data.frame() %>% .std_ci() %>%
  dplyr::rename(pct_hat = prob, pct_LCL = lower.CL, pct_UCL = upper.CL) %>%
  mutate(across(c(pct_hat,pct_LCL,pct_UCL), ~ .x*100))
readr::write_tsv(emm_site_bb_resp,
                 file.path(stats_dir,"stats_glmmBB_camelid_by_site_IECEC_emmeans_percent.tsv"))

site_contrasts_bb <- summary(contrast(emm_site_bb_link, method="trt.vs.ctrl", ref="Duodenum"),
                             infer=TRUE, adjust="dunnett") %>%
  as.data.frame() %>% .std_ci() %>%
  mutate(
    OR_site_vs_duo = exp(estimate),
    OR_site_LCL    = exp(lower.CL),
    OR_site_UCL    = exp(upper.CL),
    OR_duo_vs_site = exp(-estimate),
    OR_duo_LCL     = exp(-upper.CL),
    OR_duo_UCL     = exp(-lower.CL),
    p_adj_BH       = p.adjust(p.value,"BH")
  )
readr::write_tsv(site_contrasts_bb,
                 file.path(stats_dir,"stats_glmmBB_camelid_by_site_IECEC_Duodenum_vs_each_OR.tsv"))
cat("\n[GLMM-BB: site] Duodenum vs each site (ORs, Dunnett-adjusted):\n")
print(site_contrasts_bb)

iec_sites_glmm_bin <- iec_sites_glmm %>%
  mutate(Site_duo = factor(if_else(Site_simple=="Duodenum","Duodenum","Other"),
                           levels=c("Other","Duodenum")))
fit_duo_other_bb <- glmmTMB(cbind(succ, fail) ~ Site_duo + (1 | ID),
                            family = betabinomial(link="logit"),
                            data   = iec_sites_glmm_bin)
emm_duo_bb_link <- emmeans(fit_duo_other_bb, ~ Site_duo)
duo_con_bb <- summary(contrast(emm_duo_bb_link, method="trt.vs.ctrl", ref="Other"),
                      infer=TRUE, adjust="none") %>%
  as.data.frame() %>% .std_ci() %>%
  mutate(OR = exp(estimate), OR_LCL = exp(lower.CL), OR_UCL = exp(upper.CL),
         p_adj_BH = p.adjust(p.value,"BH"))
duo_resp_bb <- summary(emm_duo_bb_link, type="response") %>%
  as.data.frame() %>% .std_ci() %>%
  dplyr::rename(pct_hat = prob, pct_LCL = lower.CL, pct_UCL = upper.CL) %>%
  mutate(across(c(pct_hat,pct_LCL,pct_UCL), ~ .x*100))

readr::write_tsv(duo_con_bb,  file.path(stats_dir,"stats_glmmBB_camelid_Duodenum_vs_Other_OR.tsv"))
readr::write_tsv(duo_resp_bb, file.path(stats_dir,"stats_glmmBB_camelid_Duodenum_vs_Other_percent.tsv"))
cat("\n[GLMM-BB: Duodenum vs Other] ORs and % scale:\n")
print(duo_con_bb); print(duo_resp_bb)

cdsub_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid,
            cam_CD4, cam_CD8, cam_DP, cam_DN, non_CD4, non_CD8, non_DP, non_DN) %>%
  pivot_longer(cols = c(cam_CD4, cam_CD8, cam_DP, cam_DN, non_CD4, non_CD8, non_DP, non_DN),
               names_to = "key", values_to = "succ") %>%
  mutate(group  = if_else(startsWith(key,"cam_"), "Camelid+", "Camelid-"),
         subset = sub("^(cam_|non_)","", key),
         total  = if_else(group=="Camelid+", CD3_Camelid, CD3_nonCamelid)) %>%
  dplyr::filter(is.finite(total), total > 0, is.finite(succ), succ >= 0, succ <= total) %>%
  mutate(group = factor(group, levels=c("Camelid-","Camelid+")),
         subset= factor(subset,levels=c("CD4","CD8","DP","DN")))

subset_levels <- c("CD4","CD8","DP","DN")
bb_emm_list <- vector("list", length(subset_levels))
bb_con_list <- vector("list", length(subset_levels))

for (i in seq_along(subset_levels)) {
  s <- subset_levels[i]
  dat_s <- cdsub_glmm_dat %>% dplyr::filter(subset == s)
  if (nrow(dat_s)==0 || dplyr::n_distinct(dat_s$group) < 2) next

  fit_bb <- tryCatch(
    glmmTMB(cbind(succ, total - succ) ~ group + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
            family=betabinomial(link="logit"), data=dat_s),
    error=function(e) glmmTMB(cbind(succ, total - succ) ~ group + Tissue_primary + (1 | ID),
                              family=betabinomial(link="logit"), data=dat_s)
  )
  emm_link_bb <- emmeans(fit_bb, ~ group)
  emm_resp_bb <- summary(emm_link_bb, type="response") %>%
    as.data.frame() %>% .std_ci() %>%
    transmute(subset=s, group=as.character(group),
              pct_hat = 100*prob, pct_LCL = 100*lower.CL, pct_UCL = 100*upper.CL)
  con_bb <- summary(contrast(emm_link_bb, "pairwise"), infer=TRUE, adjust="none") %>%
    as.data.frame() %>% .std_ci() %>%
    mutate(subset=s,
           OR_pos_vs_neg = exp(-estimate),
           OR_pos_LCL    = exp(-upper.CL),
           OR_pos_UCL    = exp(-lower.CL),
           p_adj_BH      = p.adjust(p.value,"BH"))
  bb_emm_list[[i]] <- emm_resp_bb
  bb_con_list[[i]] <- con_bb
}
bb_emm_tbl <- bind_rows(bb_emm_list)
bb_con_tbl <- bind_rows(bb_con_list)
if (nrow(bb_emm_tbl)>0) readr::write_tsv(bb_emm_tbl, file.path(stats_dir,"stats_glmmBB_CDsub_emmeans_percent.tsv"))
if (nrow(bb_con_tbl)>0) readr::write_tsv(bb_con_tbl, file.path(stats_dir,"stats_glmmBB_CDsub_contrasts_OR.tsv"))
if (nrow(bb_con_tbl)>0) {
  cat("\n[CD subsets] Beta-binomial OR (Camelid+ vs Camelid−):\n")
  print(bb_con_tbl %>% transmute(subset, OR_CamPos_vs_CamNeg = OR_pos_vs_neg,
                                 LCL=OR_pos_LCL, UCL=OR_pos_UCL, p=p.value, q_BH=p_adj_BH),
        digits=4)
}

pd1_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid, cam_PD1, non_PD1) %>%
  pivot_longer(c(cam_PD1, non_PD1), names_to="key", values_to="succ") %>%
  mutate(group = if_else(grepl("^cam_",key),"Camelid+","Camelid-"),
         total = if_else(group=="Camelid+", CD3_Camelid, CD3_nonCamelid)) %>%
  dplyr::filter(is.finite(total), total>0, is.finite(succ), succ>=0, succ<=total) %>%
  mutate(group=factor(group, levels=c("Camelid-","Camelid+")))
if (nrow(pd1_glmm_dat)>0 && dplyr::n_distinct(pd1_glmm_dat$group)==2) {
  fit_pd1_bb <- tryCatch(
    glmmTMB(cbind(succ, total - succ) ~ group + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
            family=betabinomial(link="logit"), data=pd1_glmm_dat),
    error=function(e) glmmTMB(cbind(succ, total - succ) ~ group + Tissue_primary + (1 | ID),
                              family=betabinomial(link="logit"), data=pd1_glmm_dat)
  )
  emm_pd1_link <- emmeans(fit_pd1_bb, ~ group)
  emm_pd1_resp <- summary(emm_pd1_link, type="response") %>%
    as.data.frame() %>% .std_ci() %>%
    transmute(group=as.character(group), pct_hat=100*prob, pct_LCL=100*lower.CL, pct_UCL=100*upper.CL)
  pd1_con_bb <- summary(contrast(emm_pd1_link, "pairwise"), infer=TRUE, adjust="none") %>%
    as.data.frame() %>% .std_ci() %>%
    mutate(OR_pos_vs_neg = exp(-estimate), OR_pos_LCL=exp(-upper.CL), OR_pos_UCL=exp(-lower.CL),
           p_adj_BH = p.adjust(p.value,"BH"))
  readr::write_tsv(emm_pd1_resp, file.path(stats_dir,"stats_glmmBB_PD1_emmeans_percent.tsv"))
  readr::write_tsv(pd1_con_bb,  file.path(stats_dir,"stats_glmmBB_PD1_contrast_OR.tsv"))
  cat("\n[GLMM-BB] PD1 EMMeans (%):\n"); print(emm_pd1_resp, digits=4, row.names=FALSE)
  cat("\n[GLMM-BB] PD1 OR (Camelid+ / Camelid−):\n")
  print(pd1_con_bb %>% transmute(contrast, OR=OR_pos_vs_neg, LCL=OR_pos_LCL, UCL=OR_pos_UCL,
                                 p.value, q_BH=p_adj_BH), digits=4, row.names=FALSE)
}

gzmb_long <- iec_ec %>%
  transmute(ID, Tissue_fine, GI_group=gi_group3(Tissue_fine),
            CD3_Camelid, CD3_nonCamelid, cam_GZMB, non_GZMB) %>%
  mutate(pct_cam_GZMB=frac_safe(cam_GZMB,CD3_Camelid),
         pct_non_GZMB=frac_safe(non_GZMB,CD3_nonCamelid)) %>%
  pivot_longer(c(pct_cam_GZMB, pct_non_GZMB), names_to="key", values_to="value") %>%
  mutate(group=if_else(grepl("^pct_cam",key),"Camelid+","Camelid-"), marker="GZMB") %>%
  dplyr::filter(!is.na(GI_group), !is.na(value))

gzmb_pairs <- gzmb_long %>%
  group_by(ID, GI_group, marker, group) %>%
  summarise(val=mean(value), .groups="drop") %>%
  group_by(ID, GI_group, marker) %>%
  dplyr::filter(all(c("Camelid+","Camelid-") %in% group)) %>%
  ungroup() %>% mutate(pair_id = interaction(ID, GI_group, marker, drop=TRUE))

if (nrow(gzmb_pairs)>0) {
  p_gzmb <- gzmb_pairs %>%
    ggplot(aes(x=group, y=100*val, color=GI_group, group=pair_id)) +
    geom_line(alpha=0.35, linewidth=0.45) + geom_point(size=2.5) +
    scale_color_manual(values = pal_primary_full[c("Colon","Duodenum","Terminal Ileum")], name="GI group") +
    facet_grid(rows=vars(GI_group)) +
    labs(x=NULL, y="% GZMB+", title=paste0("GZMB in Camelid+ vs Camelid− — ", TITLE_IECEC, " (by GI group)"))
  save_plot_all(file.path(fig_dir,"fig_GZMB_camelid_vs_non_IECEC_by_GIgroup.pdf"),
         p_gzmb, width=8.4, height=5.6, dpi=300)

  gzmb_tests <- gzmb_pairs %>%
    dplyr::select(ID, GI_group, group, val) %>%
    pivot_wider(names_from=group, values_from=val) %>%
    group_by(GI_group) %>%
    summarise(n_pairs=sum(!is.na(`Camelid+`)&!is.na(`Camelid-`)),
              W=ifelse(n_pairs>0, wilcox.test(`Camelid+`,`Camelid-`,paired=TRUE,exact=FALSE)$statistic,NA_real_),
              p_value=ifelse(n_pairs>0, wilcox.test(`Camelid+`,`Camelid-`,paired=TRUE,exact=FALSE)$p.value,NA_real_),
              .groups="drop") %>%
    mutate(marker="GZMB", p_adj_BH=p.adjust(p_value,"BH"))
  readr::write_tsv(gzmb_tests, file.path(stats_dir,"stats_GZMB_camelid_vs_non_IECEC_by_GIgroup.tsv"))
  cat("\n[Stats] GZMB Camelid+/Camelid- paired by GI group:\n"); print(gzmb_tests)
}

gzmb_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid, cam_GZMB, non_GZMB) %>%
  pivot_longer(c(cam_GZMB,non_GZMB), names_to="key", values_to="succ") %>%
  mutate(group=if_else(grepl("^cam_",key),"Camelid+","Camelid-"),
         total=if_else(group=="Camelid+",CD3_Camelid,CD3_nonCamelid)) %>%
  dplyr::filter(is.finite(total), total>0, is.finite(succ), succ>=0, succ<=total) %>%
  mutate(group=factor(group,levels=c("Camelid-","Camelid+")))
if (nrow(gzmb_glmm_dat)>0 && dplyr::n_distinct(gzmb_glmm_dat$group)==2) {
  fit_gzmb_bb <- tryCatch(
    glmmTMB(cbind(succ,total - succ) ~ group + Tissue_primary + (1|ID) + (1|ID:Tissue_primary),
            family=betabinomial(link="logit"), data=gzmb_glmm_dat),
    error=function(e) glmmTMB(cbind(succ,total - succ) ~ group + Tissue_primary + (1|ID),
                              family=betabinomial(link="logit"), data=gzmb_glmm_dat)
  )
  emm_gzmb_link <- emmeans(fit_gzmb_bb, ~ group)
  emm_gzmb_resp <- summary(emm_gzmb_link, type="response") %>% as.data.frame() %>% .std_ci() %>%
    transmute(group=as.character(group), pct_hat=100*prob, pct_LCL=100*lower.CL, pct_UCL=100*upper.CL)
  gzmb_con_bb <- summary(contrast(emm_gzmb_link,"pairwise"), infer=TRUE, adjust="none") %>%
    as.data.frame() %>% .std_ci() %>%
    mutate(OR_pos_vs_neg=exp(-estimate), OR_pos_LCL=exp(-upper.CL), OR_pos_UCL=exp(-lower.CL),
           p_adj_BH=p.adjust(p.value,"BH"))
  readr::write_tsv(emm_gzmb_resp, file.path(stats_dir,"stats_glmmBB_GZMB_emmeans_percent.tsv"))
  readr::write_tsv(gzmb_con_bb,  file.path(stats_dir,"stats_glmmBB_GZMB_contrast_OR.tsv"))
  cat("\n[GLMM-BB] GZMB EMMeans (%):\n"); print(emm_gzmb_resp, digits=4, row.names=FALSE)
  cat("\n[GLMM-BB] GZMB OR (Camelid+ / Camelid−):\n")
  print(gzmb_con_bb %>% transmute(contrast, OR=OR_pos_vs_neg, LCL=OR_pos_LCL, UCL=OR_pos_UCL,
                                  p.value, q_BH=p_adj_BH), digits=4, row.names=FALSE)
}

ki67_long <- iec_ec %>%
  transmute(ID, Tissue_fine, GI_group=gi_group3(Tissue_fine),
            CD3_Camelid, CD3_nonCamelid, cam_KI67, non_KI67) %>%
  mutate(pct_cam_KI67=frac_safe(cam_KI67,CD3_Camelid), pct_non_KI67=frac_safe(non_KI67,CD3_nonCamelid)) %>%
  pivot_longer(c(pct_cam_KI67,pct_non_KI67), names_to="key", values_to="value") %>%
  mutate(group=if_else(grepl("^pct_cam",key),"Camelid+","Camelid-"), marker="KI67") %>%
  dplyr::filter(!is.na(GI_group), !is.na(value))

ki67_pairs <- ki67_long %>%
  group_by(ID, GI_group, marker, group) %>% summarise(val=mean(value), .groups="drop") %>%
  group_by(ID, GI_group, marker) %>% dplyr::filter(all(c("Camelid+","Camelid-") %in% group)) %>%
  ungroup() %>% mutate(pair_id = interaction(ID,GI_group,marker,drop=TRUE))

if (nrow(ki67_pairs)>0) {
  p_ki67 <- ki67_pairs %>%
    ggplot(aes(x=group, y=100*val, color=GI_group, group=pair_id)) +
    geom_line(alpha=0.35, linewidth=0.45) + geom_point(size=2.5) +
    scale_color_manual(values = pal_primary_full[c("Colon","Duodenum","Terminal Ileum")], name="GI group") +
    facet_grid(rows=vars(GI_group)) +
    labs(x=NULL, y="% KI67+", title=paste0("KI67 in Camelid+ vs Camelid− — ", TITLE_IECEC, " (by GI group)"))
  save_plot_all(file.path(fig_dir,"fig_KI67_camelid_vs_non_IECEC_by_GIgroup.pdf"),
         p_ki67, width=8.4, height=5.6, dpi=300)

  ki67_tests <- ki67_pairs %>%
    dplyr::select(ID, GI_group, group, val) %>%
    pivot_wider(names_from=group, values_from=val) %>%
    group_by(GI_group) %>%
    summarise(n_pairs=sum(!is.na(`Camelid+`)&!is.na(`Camelid-`)),
              W=ifelse(n_pairs>0, wilcox.test(`Camelid+`,`Camelid-`,paired=TRUE,exact=FALSE)$statistic,NA_real_),
              p_value=ifelse(n_pairs>0, wilcox.test(`Camelid+`,`Camelid-`,paired=TRUE,exact=FALSE)$p.value,NA_real_),
              .groups="drop") %>%
    mutate(marker="KI67", p_adj_BH=p.adjust(p_value,"BH"))
  readr::write_tsv(ki67_tests, file.path(stats_dir,"stats_KI67_camelid_vs_non_IECEC_by_GIgroup.tsv"))
  cat("\n[Stats] KI67 Camelid+/Camelid- paired by GI group:\n"); print(ki67_tests)
}

ki67_glmm_dat <- iec_ec %>%
  transmute(ID, Tissue_primary, CD3_Camelid, CD3_nonCamelid, cam_KI67, non_KI67) %>%
  pivot_longer(c(cam_KI67, non_KI67), names_to="key", values_to="succ_raw") %>%
  mutate(group=if_else(grepl("^cam_",key),"Camelid+","Camelid-"),
         total=if_else(group=="Camelid+",CD3_Camelid,CD3_nonCamelid)) %>%
  dplyr::filter(is.finite(total), total>0, is.finite(succ_raw), succ_raw>=0, succ_raw<=total) %>%
  mutate(succ_round = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total)))),
         fail_round = pmax(0L, as.integer(round(total)) - succ_round),
         was_adjusted = abs(succ_round - succ_raw) > 1e-6,
         group=factor(group,levels=c("Camelid-","Camelid+")))
if (nrow(ki67_glmm_dat)>0) {
  n_adj <- sum(ki67_glmm_dat$was_adjusted, na.rm=TRUE)
  cat(sprintf("\n[KI67] Rounded %d non-integer subset counts to nearest integer (out of %d rows).\n", n_adj, nrow(ki67_glmm_dat)))
}
if (nrow(ki67_glmm_dat)>0 && dplyr::n_distinct(ki67_glmm_dat$group)==2) {
  fit_ki67_bb <- tryCatch(
    glmmTMB(cbind(succ_round, fail_round) ~ group + Tissue_primary + (1|ID) + (1|ID:Tissue_primary),
            family=betabinomial(link="logit"), data=ki67_glmm_dat),
    error=function(e) glmmTMB(cbind(succ_round, fail_round) ~ group + Tissue_primary + (1|ID),
                              family=betabinomial(link="logit"), data=ki67_glmm_dat)
  )
  emm_ki67_link <- emmeans(fit_ki67_bb, ~ group)

  emm_sum_resp <- summary(emm_ki67_link, type="response")
  emm_df <- as.data.frame(emm_sum_resp)
  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_ki67_link, type="link")))
    emm_ki67_resp <- emm_df %>% transmute(group=as.character(group),
                                          pct_hat=100*plogis(emmean),
                                          pct_LCL=100*plogis(lower.CL),
                                          pct_UCL=100*plogis(upper.CL))
  } else {
    emm_df <- .std_ci(emm_df)
    emm_ki67_resp <- emm_df %>% transmute(group=as.character(group),
                                          pct_hat=100*prob, pct_LCL=100*lower.CL, pct_UCL=100*upper.CL)
  }
  ki67_con_bb <- summary(contrast(emm_ki67_link,"pairwise"), infer=TRUE, adjust="none") %>%
    as.data.frame() %>% .std_ci() %>%
    mutate(OR_pos_vs_neg=exp(-estimate), OR_pos_LCL=exp(-upper.CL), OR_pos_UCL=exp(-lower.CL),
           p_adj_BH=p.adjust(p.value,"BH"))
  readr::write_tsv(emm_ki67_resp, file.path(stats_dir,"stats_glmmBB_KI67_emmeans_percent.tsv"))
  readr::write_tsv(ki67_con_bb,  file.path(stats_dir,"stats_glmmBB_KI67_contrast_OR.tsv"))
  cat("\n[GLMM-BB] KI67 EMMeans (%):\n"); print(emm_ki67_resp, digits=4, row.names=FALSE)
  cat("\n[GLMM-BB] KI67 OR (Camelid+ / Camelid−):\n")
  print(ki67_con_bb %>% transmute(contrast, OR=OR_pos_vs_neg, LCL=OR_pos_LCL, UCL=OR_pos_UCL,
                                  p.value, q_BH=p_adj_BH), digits=4, row.names=FALSE)
}

cd103_df <- df7 %>%
  dplyr::filter(
    Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
    is.finite(frac_cd103_abundance)
  ) %>%
  dplyr::mutate(
    Status_plot = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct = 100 * frac_cd103_abundance
  )

p_cd103_status <- ggplot(cd103_df,
                         aes(x = Status_plot, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = NULL,
    y = paste0("% CD103+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "CD103 abundance by status"
  ) +
  scale_x_discrete(labels = xlab_with_n(cd103_df, "Status_plot"))

save_plot_all(file.path(fig_dir, "fig_cd103_status_IECEC_vs_NotIECEC_wControls.pdf"),
       p_cd103_status, width = 8.6, height = 5.2, dpi = 300)

dd_cd103_iecec <- iec_ec %>%
  dplyr::filter(is.finite(frac_cd103_abundance), is.finite(frac_camelid_abundance)) %>%
  dplyr::transmute(
    Tissue_primary,
    x_pct = 100 * frac_cd103_abundance,
    y_pct = 100 * frac_camelid_abundance
  )

if (nrow(dd_cd103_iecec) >= 3) {
  ct_sp  <- suppressWarnings(cor.test(dd_cd103_iecec$x_pct, dd_cd103_iecec$y_pct,
                                      method = "spearman", exact = FALSE))
  rho_val <- unname(ct_sp$estimate)
  p_val   <- ct_sp$p.value

  cat(sprintf("\n[Spearman] CD103 vs Camelid+ (IEC-EC): rho = %.3f, p = %s, N = %d\n",
              rho_val, formatC(p_val, format = "e", digits = 2), nrow(dd_cd103_iecec)))

  readr::write_tsv(
    tibble::tibble(
      metric      = "CD103_vs_Camelid_IECEC",
      N           = nrow(dd_cd103_iecec),
      spearman_r  = rho_val,
      spearman_p  = p_val
    ),
    file.path(stats_dir, "stats_cd103_vs_camelid_IECEC_spearman.tsv")
  )

  lab_txt <- paste0(
    "Spearman rho = ", formatC(rho_val, digits = 2, format = "f"),
    " (p = ", formatC(p_val, format = "e", digits = 2), ")"
  )
} else {
  lab_txt <- "Spearman: insufficient N"
  warning("[Spearman] Not enough data points (N < 3) for correlation in IEC-EC.")
}

p_cd103_corr <- ggplot(dd_cd103_iecec,
                       aes(x = x_pct, y = y_pct, color = Tissue_primary)) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "grey40") +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = paste0("% CD103+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    y = ABUNDANCE_Y_LABEL,
    title = paste0("TRM (CD103) vs Camelid infiltration — ", TITLE_IECEC, " only")
  ) +
  annotate("text", x = Inf, y = Inf, label = lab_txt,
           hjust = 1.02, vjust = 1.2, size = 3.8, color = "black")

save_plot_all(file.path(fig_dir, "fig_cd103_vs_camelid_IECEC.pdf"),
       p_cd103_corr, width = 8.2, height = 5.2, dpi = 300)

cd103_glmm_bb <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC"),
                is.finite(CD3_CD103),
                is.finite(abundance_total), abundance_total > 0) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = CD3_CD103,
    total_raw  = abundance_total
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(cd103_glmm_bb) > 0) {
  n_adj <- sum(cd103_glmm_bb$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[CD103] Rounded %d non-integer counts (out of %d rows) to satisfy count model.\n",
              n_adj, nrow(cd103_glmm_bb)))
}

if (nrow(cd103_glmm_bb) > 0 && dplyr::n_distinct(cd103_glmm_bb$Status_two) == 2) {
  fit_cd103_bb <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = cd103_glmm_bb
    ),
    error = function(e) {
      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID),
        family = glmmTMB::betabinomial(link = "logit"),
        data   = cd103_glmm_bb
      )
    }
  )

  emm_cd103_link <- emmeans::emmeans(fit_cd103_bb, ~ Status_two)

  emm_sum_resp <- summary(emm_cd103_link, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)

  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_cd103_link, type = "link")))
    emm_cd103_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_cd103_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  cd103_con_bb <- summary(contrast(emm_cd103_link, "pairwise"), infer = TRUE, adjust = "none")
  cd103_con_bb <- .std_ci(as.data.frame(cd103_con_bb)) %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_cd103_resp, file.path(stats_dir, "stats_glmmBB_CD103_emmeans_percent.tsv"))
  readr::write_tsv(cd103_con_bb,  file.path(stats_dir, "stats_glmmBB_CD103_contrast_OR.tsv"))

  cat("\n[GLMM-BB] CD103 EMMeans (% of ", ABUNDANCE_DENOMINATOR_LABEL, ") by status:\n", sep = "")
  print(emm_cd103_resp %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] CD103 OR (IEC-EC / Not IEC-EC):\n")
  print(cd103_con_bb %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

cd68_df <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC"),
                is.finite(frac_CD68_total)) %>%
  dplyr::mutate(
    Status_two = factor(status_display(Status_main),
                        levels = c("Camelid Negative", "Camelid Positive")),
    pct = 100 * frac_CD68_total
  )

p_cd68 <- ggplot(cd68_df,
                 aes(x = Status_two, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55,
               color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(x = NULL, y = "% CD68+ (of total cells)",
       title = "Macrophage infiltration (CD68): Camelid Positive vs Camelid Negative") +
  scale_x_discrete(labels = xlab_with_n(cd68_df, "Status_two"))

save_plot_all(file.path(fig_dir, "fig_cd68_IECEC_vs_NotIECEC.pdf"),
       p_cd68, width = 8.2, height = 5.0, dpi = 300)

cd68_glmm_bb <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC"),
                is.finite(CD68), is.finite(total_cells), total_cells > 0) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = CD68,
    total_raw  = total_cells
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(cd68_glmm_bb) > 0) {
  n_adj <- sum(cd68_glmm_bb$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[CD68] Rounded %d non-integer counts (out of %d rows) to satisfy count model.\n",
              n_adj, nrow(cd68_glmm_bb)))
}

if (nrow(cd68_glmm_bb) > 0 && dplyr::n_distinct(cd68_glmm_bb$Status_two) == 2) {
  fit_cd68_bb <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = cd68_glmm_bb
    ),
    error = function(e) {
      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID),
        family = glmmTMB::betabinomial(link = "logit"),
        data   = cd68_glmm_bb
      )
    }
  )

  emm_cd68_link <- emmeans::emmeans(fit_cd68_bb, ~ Status_two)

  emm_sum_resp <- summary(emm_cd68_link, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)

  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_cd68_link, type = "link")))
    emm_cd68_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_cd68_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  cd68_con_bb <- summary(contrast(emm_cd68_link, "pairwise"), infer = TRUE, adjust = "none")
  cd68_con_bb <- .std_ci(as.data.frame(cd68_con_bb)) %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_cd68_resp, file.path(stats_dir, "stats_glmmBB_CD68_emmeans_percent.tsv"))
  readr::write_tsv(cd68_con_bb,  file.path(stats_dir, "stats_glmmBB_CD68_contrast_OR.tsv"))

  cat("\n[GLMM-BB] CD68 EMMeans (% of total cells) by status:\n")
  print(emm_cd68_resp %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] CD68 OR (IEC-EC / Not IEC-EC):\n")
  print(cd68_con_bb %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

gzmb_total_df <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
                is.finite(frac_gzmb_abundance)) %>%
  dplyr::mutate(
    Status_plot = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct = 100 * frac_gzmb_abundance
  )

p_gzmb_total_status <- ggplot(gzmb_total_df,
                              aes(x = Status_plot, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = NULL,
    y = paste0("% GZMB+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "GZMB abundance by status"
  ) +
  scale_x_discrete(labels = xlab_with_n(gzmb_total_df, "Status_plot"))

save_plot_all(file.path(fig_dir, "fig_gzmb_total_status_IECEC_vs_NotIECEC_wControls.pdf"),
       p_gzmb_total_status, width = 8.6, height = 5.2, dpi = 300)

gzmb_total_glmm_bb <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC"),
                is.finite(tot_GZMB),
                is.finite(abundance_total), abundance_total > 0) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = tot_GZMB,
    total_raw  = abundance_total
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(gzmb_total_glmm_bb) > 0) {
  n_adj <- sum(gzmb_total_glmm_bb$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[GZMB-total] Rounded %d non-integer counts (out of %d rows) to satisfy count model.\n",
              n_adj, nrow(gzmb_total_glmm_bb)))
}

if (nrow(gzmb_total_glmm_bb) > 0 && dplyr::n_distinct(gzmb_total_glmm_bb$Status_two) == 2) {
  fit_gzmbtot_bb <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = gzmb_total_glmm_bb
    ),
    error = function(e) {
      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID),
        family = glmmTMB::betabinomial(link = "logit"),
        data   = gzmb_total_glmm_bb
      )
    }
  )

  emm_gzmbtot_link <- emmeans::emmeans(fit_gzmbtot_bb, ~ Status_two)

  emm_sum_resp <- summary(emm_gzmbtot_link, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)
  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_gzmbtot_link, type = "link")))
    emm_gzmbtot_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_gzmbtot_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  gzmbtot_con_bb <- summary(contrast(emm_gzmbtot_link, "pairwise"), infer = TRUE, adjust = "none") %>%
    as.data.frame() %>%
    .std_ci() %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_gzmbtot_resp, file.path(stats_dir, "stats_glmmBB_GZMB_total_emmeans_percent.tsv"))
  readr::write_tsv(gzmbtot_con_bb,  file.path(stats_dir, "stats_glmmBB_GZMB_total_contrast_OR.tsv"))

  cat("\n[GLMM-BB] GZMB (% of ", ABUNDANCE_DENOMINATOR_LABEL, ") EMMeans by status:\n", sep = "")
  print(emm_gzmbtot_resp %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] GZMB OR (IEC-EC / Not IEC-EC):\n")
  print(gzmbtot_con_bb %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

foxp3_total_df <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
                is.finite(frac_foxp3_abundance)) %>%
  dplyr::mutate(
    Status_plot = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct = 100 * frac_foxp3_abundance
  )

p_foxp3_total_status <- ggplot(foxp3_total_df,
                               aes(x = Status_plot, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = NULL,
    y = paste0("% FOXP3+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "FOXP3 abundance by status"
  ) +
  scale_x_discrete(labels = xlab_with_n(foxp3_total_df, "Status_plot"))

save_plot_all(file.path(fig_dir, "fig_foxp3_total_status_IECEC_vs_NotIECEC_wControls.pdf"),
       p_foxp3_total_status, width = 8.6, height = 5.2, dpi = 300)

foxp3_total_glmm_bb <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC"),
                is.finite(tot_FOXP3),
                is.finite(abundance_total), abundance_total > 0) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = tot_FOXP3,
    total_raw  = abundance_total
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(foxp3_total_glmm_bb) > 0) {
  n_adj <- sum(foxp3_total_glmm_bb$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[FOXP3-total] Rounded %d non-integer counts (out of %d rows) to satisfy count model.\n",
              n_adj, nrow(foxp3_total_glmm_bb)))
}

if (nrow(foxp3_total_glmm_bb) > 0 && dplyr::n_distinct(foxp3_total_glmm_bb$Status_two) == 2) {
  fit_foxp3tot_bb <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = foxp3_total_glmm_bb
    ),
    error = function(e) {
      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID),
        family = glmmTMB::betabinomial(link = "logit"),
        data   = foxp3_total_glmm_bb
      )
    }
  )

  emm_foxp3tot_link <- emmeans::emmeans(fit_foxp3tot_bb, ~ Status_two)

  emm_sum_resp <- summary(emm_foxp3tot_link, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)
  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_foxp3tot_link, type = "link")))
    emm_foxp3tot_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_foxp3tot_resp <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  foxp3tot_con_bb <- summary(contrast(emm_foxp3tot_link, "pairwise"), infer = TRUE, adjust = "none") %>%
    as.data.frame() %>%
    .std_ci() %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_foxp3tot_resp, file.path(stats_dir, "stats_glmmBB_FOXP3_total_emmeans_percent.tsv"))
  readr::write_tsv(foxp3tot_con_bb,  file.path(stats_dir, "stats_glmmBB_FOXP3_total_contrast_OR.tsv"))

  cat("\n[GLMM-BB] FOXP3 (% of ", ABUNDANCE_DENOMINATOR_LABEL, ") EMMeans by status:\n", sep = "")
  print(emm_foxp3tot_resp %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] FOXP3 OR (IEC-EC / Not IEC-EC):\n")
  print(foxp3tot_con_bb %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

gzmb_total_df_noOther <- df7 %>%
  dplyr::filter(
    Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
    is.finite(frac_gzmb_abundance),
    Tissue_primary != "Other"
  ) %>%
  dplyr::mutate(
    Status_plot = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct = 100 * frac_gzmb_abundance
  )

p_gzmb_total_status_noOther <- ggplot(gzmb_total_df_noOther,
                                      aes(x = Status_plot, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = NULL,
    y = paste0("% GZMB+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "GZMB abundance by status — excluding Other tissues"
  ) +
  scale_x_discrete(labels = xlab_with_n(gzmb_total_df_noOther, "Status_plot"))

save_plot_all(file.path(fig_dir, "fig_gzmb_total_status_IECEC_vs_NotIECEC_wControls_noOther.pdf"),
       p_gzmb_total_status_noOther, width = 8.6, height = 5.2, dpi = 300)

gzmb_total_glmm_bb_noOther <- df7 %>%
  dplyr::filter(
    Status_main %in% c("IEC-EC", "Not IEC-EC"),
    Tissue_primary != "Other",
    is.finite(tot_GZMB),
    is.finite(abundance_total), abundance_total > 0
  ) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = tot_GZMB,
    total_raw  = abundance_total
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(gzmb_total_glmm_bb_noOther) > 0) {
  n_adj <- sum(gzmb_total_glmm_bb_noOther$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[GZMB-total (noOther)] Rounded %d non-integer counts (out of %d rows).\n",
              n_adj, nrow(gzmb_total_glmm_bb_noOther)))
}

if (nrow(gzmb_total_glmm_bb_noOther) > 0 && dplyr::n_distinct(gzmb_total_glmm_bb_noOther$Status_two) == 2) {
  fit_gzmbtot_bb_noOther <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID) + (1 | ID:Tissue_primary),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = gzmb_total_glmm_bb_noOther
    ),
    error = function(e) {
      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two + Tissue_primary + (1 | ID),
        family = glmmTMB::betabinomial(link = "logit"),
        data   = gzmb_total_glmm_bb_noOther
      )
    }
  )

  emm_gzmbtot_link_noOther <- emmeans::emmeans(fit_gzmbtot_bb_noOther, ~ Status_two)
  emm_sum_resp <- summary(emm_gzmbtot_link_noOther, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)

  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_gzmbtot_link_noOther, type = "link")))
    emm_gzmbtot_resp_noOther <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_gzmbtot_resp_noOther <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  gzmbtot_con_bb_noOther <- summary(contrast(emm_gzmbtot_link_noOther, "pairwise"),
                                    infer = TRUE, adjust = "none") %>%
    as.data.frame() %>% .std_ci() %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_gzmbtot_resp_noOther,
                   file.path(stats_dir, "stats_glmmBB_GZMB_total_emmeans_percent_noOther.tsv"))
  readr::write_tsv(gzmbtot_con_bb_noOther,
                   file.path(stats_dir, "stats_glmmBB_GZMB_total_contrast_OR_noOther.tsv"))

  cat("\n[GLMM-BB] GZMB (% of ", ABUNDANCE_DENOMINATOR_LABEL,
      ") — excluding Other tissues — EMMeans by status:\n", sep = "")
  print(emm_gzmbtot_resp_noOther %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] GZMB — excluding Other — OR (IEC-EC / Not IEC-EC):\n")
  print(gzmbtot_con_bb_noOther %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

gzmb_total_df_duo <- df7 %>%
  dplyr::filter(
    Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
    is.finite(frac_gzmb_abundance),
    Tissue_primary == "Duodenum"
  ) %>%
  dplyr::mutate(
    Status_plot = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct = 100 * frac_gzmb_abundance
  )

p_gzmb_total_status_duo <- ggplot(gzmb_total_df_duo,
                                  aes(x = Status_plot, y = pct, color = Tissue_primary)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  labs(
    x = NULL,
    y = paste0("% GZMB+ among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "GZMB abundance by status — Duodenum only"
  ) +
  scale_x_discrete(labels = xlab_with_n(gzmb_total_df_duo, "Status_plot"))

save_plot_all(file.path(fig_dir, "fig_gzmb_total_status_IECEC_vs_NotIECEC_wControls_DuodenumOnly.pdf"),
       p_gzmb_total_status_duo, width = 8.6, height = 5.2, dpi = 300)

gzmb_total_glmm_bb_duo <- df7 %>%
  dplyr::filter(
    Status_main %in% c("IEC-EC", "Not IEC-EC"),
    Tissue_primary == "Duodenum",
    is.finite(tot_GZMB),
    is.finite(abundance_total), abundance_total > 0
  ) %>%
  dplyr::mutate(
    Status_two = factor(Status_main, levels = c("Not IEC-EC", "IEC-EC")),
    succ_raw   = tot_GZMB,
    total_raw  = abundance_total
  ) %>%
  dplyr::mutate(
    succ = pmax(0L, pmin(as.integer(round(succ_raw)), as.integer(round(total_raw)))),
    fail = pmax(0L, as.integer(round(total_raw)) - succ),
    was_adjusted = abs(succ - succ_raw) > 1e-6 | abs((succ + fail) - total_raw) > 1e-6
  )

if (nrow(gzmb_total_glmm_bb_duo) > 0) {
  n_adj <- sum(gzmb_total_glmm_bb_duo$was_adjusted, na.rm = TRUE)
  cat(sprintf("\n[GZMB-total (Duodenum)] Rounded %d non-integer counts (out of %d rows).\n",
              n_adj, nrow(gzmb_total_glmm_bb_duo)))
}

if (nrow(gzmb_total_glmm_bb_duo) > 0 && dplyr::n_distinct(gzmb_total_glmm_bb_duo$Status_two) == 2) {

  fit_gzmbtot_bb_duo <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(succ, fail) ~ Status_two + (1 | ID),
      family = glmmTMB::betabinomial(link = "logit"),
      data   = gzmb_total_glmm_bb_duo
    ),
    error = function(e) {

      glmmTMB::glmmTMB(
        cbind(succ, fail) ~ Status_two,
        family = glmmTMB::betabinomial(link = "logit"),
        data   = gzmb_total_glmm_bb_duo
      )
    }
  )

  emm_gzmbtot_link_duo <- emmeans::emmeans(fit_gzmbtot_bb_duo, ~ Status_two)
  emm_sum_resp <- summary(emm_gzmbtot_link_duo, type = "response")
  emm_df <- as.data.frame(emm_sum_resp)

  if (!("prob" %in% names(emm_df))) {
    emm_df <- .std_ci(as.data.frame(summary(emm_gzmbtot_link_duo, type = "link")))
    emm_gzmbtot_resp_duo <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * plogis(emmean),
        pct_LCL = 100 * plogis(lower.CL),
        pct_UCL = 100 * plogis(upper.CL)
      )
  } else {
    emm_df <- .std_ci(emm_df)
    emm_gzmbtot_resp_duo <- emm_df %>%
      dplyr::transmute(
        Status_two = as.character(Status_two),
        pct_hat = 100 * prob,
        pct_LCL = 100 * lower.CL,
        pct_UCL = 100 * upper.CL
      )
  }

  gzmbtot_con_bb_duo <- summary(contrast(emm_gzmbtot_link_duo, "pairwise"),
                                infer = TRUE, adjust = "none") %>%
    as.data.frame() %>% .std_ci() %>%
    dplyr::mutate(
      OR_IEC_vs_Not = exp(-estimate),
      OR_LCL        = exp(-upper.CL),
      OR_UCL        = exp(-lower.CL),
      p_adj_BH      = p.adjust(p.value, "BH")
    )

  readr::write_tsv(emm_gzmbtot_resp_duo,
                   file.path(stats_dir, "stats_glmmBB_GZMB_total_emmeans_percent_DuodenumOnly.tsv"))
  readr::write_tsv(gzmbtot_con_bb_duo,
                   file.path(stats_dir, "stats_glmmBB_GZMB_total_contrast_OR_DuodenumOnly.tsv"))

  cat("\n[GLMM-BB] GZMB (% of ", ABUNDANCE_DENOMINATOR_LABEL,
      ") — Duodenum only — EMMeans by status:\n", sep = "")
  print(emm_gzmbtot_resp_duo %>% dplyr::arrange(Status_two), digits = 4, row.names = FALSE)

  cat("\n[GLMM-BB] GZMB — Duodenum only — OR (IEC-EC / Not IEC-EC):\n")
  print(gzmbtot_con_bb_duo %>%
          dplyr::transmute(contrast, OR = OR_IEC_vs_Not, LCL = OR_LCL, UCL = OR_UCL,
                           p.value, q_BH = p_adj_BH),
        digits = 4, row.names = FALSE)
}

iec_markers_long <- iec_ec %>%
  dplyr::transmute(
    ID, Tissue_fine, Tissue_primary,
    CD3_Camelid, CD3_nonCamelid,
    cam_KI67 = cam_KI67,   non_KI67 = non_KI67,
    cam_GZMB = cam_GZMB,   non_GZMB = non_GZMB,
    cam_PD1  = cam_PD1,    non_PD1  = non_PD1
  ) %>%
  tidyr::pivot_longer(
    cols = c(cam_KI67, non_KI67, cam_GZMB, non_GZMB, cam_PD1, non_PD1),
    names_to = "key", values_to = "succ"
  ) %>%
  dplyr::mutate(
    group  = dplyr::if_else(grepl("^cam_", key), "Camelid+", "Camelid-"),
    marker = dplyr::case_when(
      grepl("KI67", key) ~ "KI67",
      grepl("GZMB", key) ~ "GZMB",
      grepl("PD1",  key) ~ "PD1",
      TRUE ~ NA_character_
    ),
    total  = dplyr::if_else(group == "Camelid+", CD3_Camelid, CD3_nonCamelid),
    pct    = 100 * frac_safe(succ, total)
  ) %>%
  dplyr::filter(!is.na(marker), is.finite(pct), is.finite(total), total > 0) %>%
  dplyr::mutate(
    group  = factor(group,  levels = c("Camelid-", "Camelid+")),
    marker = factor(marker, levels = c("KI67", "GZMB", "PD1"))
  )

iec_marker_pairs <- iec_markers_long %>%
  dplyr::group_by(ID, Tissue_fine, Tissue_primary, marker, group) %>%
  dplyr::summarise(val = mean(pct), .groups = "drop") %>%
  dplyr::group_by(ID, Tissue_fine, Tissue_primary, marker) %>%
  dplyr::filter(all(c("Camelid+", "Camelid-") %in% group)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pair_id = interaction(ID, Tissue_fine, marker, drop = TRUE))

p_camelid_facets <- iec_marker_pairs %>%
  ggplot2::ggplot(ggplot2::aes(x = group, y = val, color = Tissue_primary, group = pair_id)) +
  ggplot2::geom_line(alpha = 0.35, linewidth = 0.5) +
  ggplot2::geom_point(size = 2.4) +
  scale_color_tissue_primary("Tissue") +
  ggplot2::facet_wrap(~ marker, nrow = 1) +
  ggplot2::labs(
    x = NULL, y = "% positive",
    title = paste0("Camelid+ vs Camelid− — ", TITLE_IECEC, " (KI67, GZMB, PD1)")
  ) +
  ggplot2::theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

save_plot_all(file.path(fig_dir, "fig_camelid_pos_vs_neg_facets_KI67_GZMB_PD1_IECEC.pdf"),
       p_camelid_facets, width = 11.0, height = 4.8, dpi = 300)

status_levels3 <- c("Untreated Controls", "Camelid Negative", "Camelid Positive")

df_cd103_status <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
                is.finite(frac_cd103_abundance)) %>%
  dplyr::mutate(
    Status_plot = factor(status_display(Status_main), levels = status_levels3),
    pct         = 100 * frac_cd103_abundance,
    metric      = "CD103"
  )

df_gzmb_status <- df7 %>%
  dplyr::filter(Status_main %in% c("IEC-EC", "Not IEC-EC", "Untreated Controls"),
                is.finite(frac_gzmb_abundance)) %>%
  dplyr::mutate(
    Status_plot = factor(status_display(Status_main), levels = status_levels3),
    pct         = 100 * frac_gzmb_abundance,
    metric      = "GZMB"
  )

df_status_facets <- dplyr::bind_rows(df_cd103_status, df_gzmb_status) %>%
  dplyr::mutate(metric = factor(metric, levels = c("CD103", "GZMB")))

p_status_facets <- ggplot2::ggplot(
  df_status_facets,
  ggplot2::aes(x = Status_plot, y = pct, color = Tissue_primary)
) +
  ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.9, size = 2.6) +
  ggplot2::stat_summary(fun = stats::median, geom = "crossbar",
                        width = 0.6, linewidth = 0.55, color = "black", fill = NA) +
  scale_color_tissue_primary("Tissue") +
  ggplot2::facet_wrap(~ metric, nrow = 1) +
  ggplot2::labs(
    x = NULL,
    y = paste0("% positive among ", ABUNDANCE_DENOMINATOR_LABEL),
    title = "By status — CD103 and GZMB (controls included)"
  ) +
  ggplot2::theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

save_plot_all(file.path(fig_dir, "fig_status_facets_CD103_GZMB_wControls.pdf"),
       p_status_facets, width = 10.0, height = 4.8, dpi = 300)

# Optional display labels
label_path <- file.path(data_dir, "patient_labels.tsv")
map_id_to_paper <- if (file.exists(label_path)) {
  readr::read_tsv(label_path, show_col_types = FALSE) %>% dplyr::distinct()
} else {
  tibble::tibble(ID = character(), Paper_ID = character())
}

biopsy_base <- df7 %>%
  dplyr::filter(
    Status_main %in% c("Untreated Controls", "Not IEC-EC", "IEC-EC"),
    is.finite(frac_camelid_abundance),
    !is.na(ID), !is.na(Tissue_primary)
  ) %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::left_join(map_id_to_paper, by = "ID") %>%
  dplyr::mutate(
    Paper_ID = dplyr::coalesce(Paper_ID, ID),
    Status_panel = factor(
      status_display(Status_main),
      levels = c("Untreated Controls", "Camelid Negative", "Camelid Positive")
    ),
    pct_camelid_abundance = 100 * frac_camelid_abundance
  )

if (nrow(biopsy_base) == 0L) {
  warning("[Patient plot] No biopsies found after filtering. Skipping figure.")
} else {

  patient_counts <- biopsy_base %>%
    dplyr::count(Status_panel, Paper_ID, name = "n_biopsies")

  stats_by_patient <- biopsy_base %>%
    dplyr::group_by(Status_panel, Paper_ID) %>%
    dplyr::summarise(
      mean_pct = mean(pct_camelid_abundance, na.rm = TRUE),
      max_pct  = max(pct_camelid_abundance, na.rm = TRUE),
      .groups = "drop"
    )

  id_order <- stats_by_patient %>%
    dplyr::left_join(patient_counts, by = c("Status_panel","Paper_ID")) %>%
    dplyr::arrange(Status_panel, max_pct, dplyr::desc(n_biopsies), Paper_ID) %>%
    dplyr::group_by(Status_panel) %>%
    dplyr::mutate(
      Paper_ord = dplyr::row_number(),
      Paper_lab = paste0(Paper_ID, " (n=", n_biopsies, ")")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Status_panel, Paper_ID, Paper_ord, Paper_lab)

  biopsy_plot <- biopsy_base %>%
    dplyr::left_join(id_order, by = c("Status_panel", "Paper_ID")) %>%
    dplyr::group_by(Status_panel) %>%
    dplyr::mutate(Paper_lab = forcats::fct_reorder(Paper_lab, Paper_ord, .desc = FALSE)) %>%
    dplyr::ungroup()

  bars_df <- stats_by_patient %>%
    dplyr::left_join(id_order, by = c("Status_panel", "Paper_ID")) %>%
    dplyr::group_by(Status_panel) %>%
    dplyr::mutate(Paper_lab = forcats::fct_reorder(Paper_lab, Paper_ord, .desc = FALSE)) %>%
    dplyr::ungroup()

  p_pat_h <- ggplot2::ggplot(biopsy_plot, ggplot2::aes(y = Paper_lab)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Status_panel),
      scales = "free_y",
      space  = "free_y"
    ) +

    ggplot2::geom_vline(
      xintercept = ABUNDANCE_REFERENCE_LINES,
      linetype   = "dashed",
      color      = "grey75",
      linewidth  = 0.35,
      alpha      = 0.5
    ) +

    ggplot2::geom_col(
      data = bars_df,
      ggplot2::aes(x = max_pct, y = Paper_lab),
      fill = "grey95", color = "grey88", alpha = 0.55, width = 0.58, inherit.aes = FALSE
    ) +

    ggplot2::geom_point(
      ggplot2::aes(x = pct_camelid_abundance, color = Tissue_primary),
      size = 2.1, alpha = 0.95, position = ggplot2::position_identity()
    ) +
    scale_color_tissue_primary("Tissue") +
    ggplot2::scale_x_continuous(
      name = ABUNDANCE_Y_LABEL,
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      y = NULL,
      title = ""
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing.y = grid::unit(0.8, "lines"),
      axis.text.y = ggplot2::element_text(hjust = 1)
    )

  out_file <- file.path(fig_dir, "fig_patient_MAX_bar_and_camelid_dots_by_status_ALLSITES_PAPERID.pdf")
  save_plot_all(out_file, p_pat_h, dpi = 300, width = 10, height = 8)

  readr::write_tsv(
    patient_counts %>%
      dplyr::arrange(Status_panel, Paper_ID),
    file.path(stats_dir, "patient_biopsy_counts_by_status_ALLSITES_PAPERID.tsv")
  )

  readr::write_tsv(
    stats_by_patient %>%
      dplyr::arrange(Status_panel, Paper_ID),
    file.path(stats_dir, "per_patient_MEAN_MAX_camelid_abundance_ALLSITES_PAPERID.tsv")
  )

  readr::write_tsv(
    biopsy_plot %>%
      dplyr::transmute(
        Status_panel, Paper_ID, ID, Tissue_fine, Tissue_primary, pct_camelid_abundance
      ) %>%
      dplyr::arrange(Status_panel, Paper_ID, Tissue_primary, Tissue_fine),
    file.path(stats_dir, "per_biopsy_camelid_abundance_ALLSITES_PAPERID.tsv")
  )

  cat(sprintf("\n[Patient plot] MAX-ordered with reference guides saved (10 x 8 in): %s\n", out_file))
}

writeLines(
  c(
    "COMET analysis settings",
    paste0("Run date: ", Sys.Date()),
    paste0("Input: ", path_df7),
    paste0("Abundance denominator: ", ABUNDANCE_DENOMINATOR),
    paste0("Abundance label: ", ABUNDANCE_DENOMINATOR_LABEL),
    paste0("Output: ", output_dir),
    paste0("Vector font family: ", BASE_FAMILY),
    "PDF/EPS text: live and selectable"
  ),
  con = file.path(log_dir, "analysis_settings.txt")
)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
message("COMET analysis complete: ", output_dir)
