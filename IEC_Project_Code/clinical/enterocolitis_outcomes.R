#!/usr/bin/env Rscript
# Enterocolitis outcome comparisons
# Compare nutritional support, hospitalization, mortality, and infection outcomes by infiltrate status.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
out_dir <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "clinical", "enterocolitis_outcomes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Clinical input
df <- readr::read_tsv(file.path(data_dir, "clinical_outcomes.tsv"), show_col_types = FALSE) %>%
  mutate(
    CAR       = factor(CAR, levels = c("N","Y")),
    TPN       = factor(TPN, levels = c("N","Y")),
    Mortality = factor(Mortality, levels = c("N","Y"))
  )

df
table(df$CAR)

tab_tpn  <- table(df$CAR, df$TPN)
tab_mort <- table(df$CAR, df$Mortality)

fisher_tpn  <- fisher.test(tab_tpn,  alternative = "two.sided")
fisher_mort <- fisher.test(tab_mort, alternative = "two.sided")

or_ha <- function(tab_2x2) {

  a <- tab_2x2["Y","Y"]
  b <- tab_2x2["Y","N"]
  c <- tab_2x2["N","Y"]
  d <- tab_2x2["N","N"]

  a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
  OR <- (a2 * d2) / (b2 * c2)
  se <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
  ci <- exp(log(OR) + c(-1, 1) * 1.96 * se)
  tibble(OR_HA = OR, CI_low = ci[1], CI_high = ci[2])
}

or_tpn_ha  <- or_ha(tab_tpn)
or_mort_ha <- or_ha(tab_mort)

wilcox_hosp <- wilcox.test(HospitalDays ~ CAR, data = df,
                           alternative = "two.sided", conf.int = TRUE)

wilcox_inf <- wilcox.test(PostEndoscopyInfections ~ CAR, data = df,
                          alternative = "two.sided", conf.int = TRUE)

pois_inf <- glm(PostEndoscopyInfections ~ CAR, data = df, family = poisson())
irr_inf  <- exp(coef(pois_inf))
irr_ci   <- exp(confint(pois_inf))

desc <- df %>%
  group_by(CAR) %>%
  summarise(
    n = n(),
    TPN_n = sum(TPN == "Y"),
    TPN_prop = mean(TPN == "Y"),
    Mort_n = sum(Mortality == "Y"),
    Mort_prop = mean(Mortality == "Y"),
    Hosp_median = median(HospitalDays),
    Hosp_IQR = IQR(HospitalDays),
    Inf_median = median(PostEndoscopyInfections),
    Inf_IQR = IQR(PostEndoscopyInfections),
    .groups = "drop"
  )

tests <- tibble(
  outcome = c("TPN (Y/N)", "Mortality (Y/N)", "Hospital days", "Post-endoscopy infections"),
  test    = c("Fisher exact", "Fisher exact", "Wilcoxon rank-sum", "Wilcoxon rank-sum (plus Poisson IRR below)"),
  p_value = c(fisher_tpn$p.value, fisher_mort$p.value, wilcox_hosp$p.value, wilcox_inf$p.value)
)

list(
  descriptives = desc,
  tests = tests,
  fisher_TPN = fisher_tpn,
  fisher_Mortality = fisher_mort,
  OR_TPN_HA = or_tpn_ha,
  OR_Mortality_HA = or_mort_ha,
  wilcox_HospitalDays = wilcox_hosp,
  wilcox_Infections = wilcox_inf,
  poisson_Infections = summary(pois_inf),
  poisson_IRR = tibble(term = names(irr_inf), IRR = as.numeric(irr_inf),
                       CI_low = as.numeric(irr_ci[,1]), CI_high = as.numeric(irr_ci[,2]))
)

# Outcome summaries and figures
binom_ci_by_group <- function(data, outcome_col) {
  data %>%
    group_by(CAR) %>%
    summarise(
      n = n(),
      x = sum(.data[[outcome_col]] == "Y"),
      prop = x / n,
      ci_low  = binom.test(x, n)$conf.int[1],
      ci_high = binom.test(x, n)$conf.int[2],
      .groups = "drop"
    )
}

p_bin <- function(dat_ci, title) {
  ggplot(dat_ci, aes(x = CAR, y = prop)) +
    geom_col() +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion", title = title) +
    theme_classic(base_size = 13)
}

tpn_ci  <- binom_ci_by_group(df, "TPN")
mort_ci <- binom_ci_by_group(df, "Mortality")

p_tpn <- p_bin(tpn_ci,  "TPN use by CAR status")
p_mort <- p_bin(mort_ci, "Mortality by CAR status")

p_hosp <- ggplot(df, aes(x = CAR, y = HospitalDays)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  labs(x = NULL, y = "Hospital days", title = "Hospital days by CAR status") +
  theme_classic(base_size = 13)

p_inf <- ggplot(df, aes(x = CAR, y = PostEndoscopyInfections)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  labs(x = NULL, y = "Post-endoscopy infections", title = "Post-endoscopy infections by CAR status") +
  theme_classic(base_size = 13)

(p_hosp | p_inf)

ggsave(file.path(out_dir, "outcome_comparisons.pdf"), width = 10, height = 5)
