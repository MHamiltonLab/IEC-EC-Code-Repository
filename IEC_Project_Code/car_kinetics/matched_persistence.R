#!/usr/bin/env Rscript
# Time-matched CAR persistence
# Match late CAR measurements within a day caliper and compare patient-level persistence.

suppressPackageStartupMessages({
  library(readr);  library(dplyr);  library(tidyr)
  library(ggplot2); library(scales); library(stringr)
  library(grid);   library(gridExtra)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
path_df3 <- file.path(data_dir, "clinical_metadata.tsv")
path_df4 <- file.path(data_dir, "car_measurements.tsv")

# Detection limit and matching parameters
LOD_ND           <- 0.1
day_low          <- 50L
day_high         <- Inf
caliper_days     <- 14L
max_ratio        <- 2L
allow_replacement<- FALSE
pd_censor_pool   <- TRUE
pd_censor_post   <- TRUE
pal              <- c("No IEC"="#4C6A87", "IEC"="#B22222")

safe_num <- function(x){
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(readr::parse_number(as.character(x)))
}
theme_arial <- function(base_size = 22){
  theme_classic(base_size = base_size) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text  = element_text(size = base_size - 4))
}

df3 <- read_delim(path_df3, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
df4 <- read_delim(path_df4, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)

suppressPackageStartupMessages({ library(grid); library(gridExtra) })

out_dir <- file.path(file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "car_kinetics", "matched_persistence"), "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

.has_cairo <- requireNamespace("Cairo", quietly = TRUE)
.is_mac    <- identical(Sys.info()[["sysname"]], "Darwin")

if (.is_mac) {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

save_pdf <- function(plot, filename, width = 9, height = 7) {
  f <- file.path(out_dir, filename)
  if (.has_cairo) {
    device_fun <- function(file, width, height, ...) {
      Cairo::CairoPDF(file = file, width = width, height = height, family = "Arial", ...)
    }
  } else if (.is_mac) {
    device_fun <- function(file, width, height, ...) {
      grDevices::quartz(file = file, type = "pdf",
                        width = width, height = height, family = "Arial", ...)
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
  if (.has_cairo) {
    device_fun <- function(file, width, height, ...) {
      Cairo::CairoPS(file = file, width = width, height = height,
                     family = "Arial", onefile = FALSE)
    }
  } else {

    device_fun <- function(file, width, height, ...) {
      grDevices::postscript(file = file, onefile = FALSE, paper = "special",
                            width = width, height = height, family = "Helvetica",
                            horizontal = FALSE, ...)
    }
  }
  ggplot2::ggsave(filename = f, plot = plot, device = device_fun,
                  width = width, height = height, units = "in")
  message("Saved: ", f)
}

meta <- df3 %>%
  transmute(
    Study_ID,
    IEC_flag = factor(if_else(safe_num(IEC_enteritis) == 1, "IEC", "No IEC"),
                      levels = c("No IEC","IEC")),
    PD_day   = safe_num(PD_day)
  )

day_candidates <- c("Day","Day_rel_infusion","day","Day_from_infusion","Day_rel")
car_candidates <- c("CAR_abs","CAR_Abs","Absolute_CAR","CARcells_per_uL","CAR_cells_per_uL")

day_col <- day_candidates[day_candidates %in% names(df4)][1]
car_col <- car_candidates[car_candidates %in% names(df4)][1]
if (is.na(day_col) || is.na(car_col)) {
  stop("Could not detect day or CAR columns in df4.\n",
       "Looked for day in: ", paste(day_candidates, collapse=", "),
       "\nLooked for CAR in: ", paste(car_candidates, collapse=", "),
       "\nAvailable columns: ", paste(names(df4), collapse=", "))
}
message(sprintf("Using day column '%s' and CAR column '%s'.", day_col, car_col))

exp_all <- df4 %>%
  transmute(
    Study_ID,
    Day = as.numeric(.data[[day_col]]),
    CAR_abs = suppressWarnings(as.numeric(.data[[car_col]]))
  ) %>%
  inner_join(meta, by = "Study_ID") %>%
  filter(is.finite(Day), Day >= 0) %>%
  mutate(
    CAR_abs_imp = if_else(!is.finite(CAR_abs) | CAR_abs <= 0, LOD_ND, CAR_abs)
  )

exp_0_100 <- exp_all %>% filter(Day <= 100)

p_spaghetti_0_100 <-
  ggplot() +

  geom_line(
    data = dplyr::filter(exp_0_100, IEC_flag=="No IEC"),
    aes(Day, CAR_abs_imp, group = Study_ID, color = IEC_flag),
    linewidth = 0.7, alpha = 0.65
  ) +
  geom_line(
    data = dplyr::filter(exp_0_100, IEC_flag=="IEC"),
    aes(Day, CAR_abs_imp, group = Study_ID, color = IEC_flag),
    linewidth = 0.9, alpha = 0.95
  ) +
  geom_hline(yintercept = LOD_ND, linetype = "dashed") +
  annotate("text", x = 92, y = LOD_ND, label = "Not detected (0.1)",
           hjust = 1, vjust = -0.6, family = "Arial", size = 4) +
  scale_color_manual(values = pal, name = "Group") +
  scale_y_log10("Absolute CAR (cells/µL)",
                breaks = 10^(-1:5),
                labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  scale_x_continuous("Days post infusion", limits = c(0, 100)) +
  ggtitle("CAR-T expansion (0–100 d): IEC highlighted") +
  theme_arial(22)
print(p_spaghetti_0_100)

# Patient-level matching without replacement
match_iec_unique <- function(exp_all,
                             day_low=50, day_high=Inf,
                             max_ratio=2L,
                             caliper_fun=function(d) 14,
                             allow_replacement=FALSE,
                             pd_censor_pool=TRUE,
                             pd_censor_post=TRUE) {

  iec_pool0 <- exp_all %>%
    filter(IEC_flag=="IEC", Day > day_low, Day <= day_high) %>%
    arrange(Study_ID, Day)
  ctl_pool0 <- exp_all %>%
    filter(IEC_flag=="No IEC", Day > day_low, Day <= day_high) %>%
    arrange(Study_ID, Day)

  if (pd_censor_pool) {
    iec_pool0 <- iec_pool0 %>% filter(is.na(PD_day) | PD_day >= Day)
    ctl_pool0 <- ctl_pool0 %>% filter(is.na(PD_day) | PD_day >= Day)
  }

  if (nrow(iec_pool0)==0 || nrow(ctl_pool0)==0)
    stop("Empty pool after filters; check inputs.")

  iec_cands <- iec_pool0 %>%
    group_by(Study_ID) %>%
    summarise(cand_days = list(Day), .groups="drop")

  used_ctrl_ids <- character(0)
  controls_within <- function(day) {
    tol <- caliper_fun(day)
    avail <- if (allow_replacement) ctl_pool0 else ctl_pool0 %>% filter(! (Study_ID %in% used_ctrl_ids))
    if (!nrow(avail)) return(avail[0, ])
    avail %>%
      mutate(abs_diff = abs(Day - day)) %>%
      filter(abs_diff <= tol) %>%
      arrange(abs_diff) %>%
      group_by(Study_ID) %>%
      slice_min(abs_diff, n = 1, with_ties = FALSE) %>%
      ungroup()
  }

  scarcity_tbl <- iec_cands %>%
    rowwise() %>%
    mutate(
      min_avail = {
        if (!length(cand_days)) 0L else {
          vals <- vapply(cand_days, function(d) nrow(controls_within(d)), numeric(1))
          if (length(vals)==0) 0L else min(vals)
        }
      },
      feasible_days = {
        if (!length(cand_days)) 0L else {
          vals <- vapply(cand_days, function(d) nrow(controls_within(d)) > 0, logical(1))
          sum(vals)
        }
      },
      first_day = if (length(cand_days)) min(cand_days) else Inf
    ) %>%
    ungroup() %>%
    arrange(min_avail, first_day)

  matches <- list(); set_id <- 0L; iec_used <- character(0)

  for (i in seq_len(nrow(scarcity_tbl))) {
    sid  <- scarcity_tbl$Study_ID[i]
    if (sid %in% iec_used) next
    days <- sort(unique(unlist(iec_cands$cand_days[iec_cands$Study_ID==sid])))

    chosen_day <- NA_real_; chosen_ctls <- NULL
    for (d in days) {
      cand <- controls_within(d)
      if (nrow(cand) >= 1L) {
        chosen_day  <- d
        chosen_ctls <- cand %>% slice_head(n = max_ratio)
        break
      }
    }
    if (is.na(chosen_day)) next

    set_id <- set_id + 1L
    iec_row  <- iec_pool0 %>% filter(Study_ID==sid, Day==chosen_day) %>% slice(1) %>% mutate(set_id = set_id)
    ctl_rows <- chosen_ctls %>% mutate(set_id = set_id)

    matches[[length(matches)+1L]] <- list(iec = iec_row, ctrls = ctl_rows)

    if (!allow_replacement) used_ctrl_ids <- union(used_ctrl_ids, ctl_rows$Study_ID)
    iec_used <- c(iec_used, sid)
  }

  if (!length(matches)) stop("Matching failed: no sets created under current caliper.")

  iec_matched <- dplyr::bind_rows(lapply(matches, `[[`, "iec"))
  ctl_matched <- dplyr::bind_rows(lapply(matches, `[[`, "ctrls"))

  matched <- bind_rows(
    iec_matched %>% mutate(IEC_flag = factor("IEC",    levels=c("No IEC","IEC"))),
    ctl_matched %>% mutate(IEC_flag = factor("No IEC", levels=c("No IEC","IEC")))
  ) %>% select(set_id, Study_ID, IEC_flag, Day, PD_day, CAR_abs_imp)

  matched <- matched %>%
    group_by(set_id) %>%
    mutate(pd_prior_match = !is.na(PD_day) & PD_day < Day) %>%
    ungroup()

  if (pd_censor_post) {
    if (any(matched$pd_prior_match, na.rm = TRUE)) {
      matched <- matched %>% filter(!pd_prior_match)
      keep_sets <- matched %>%
        group_by(set_id) %>%
        summarise(n_iec = sum(IEC_flag=="IEC"),
                  n_ctl = sum(IEC_flag=="No IEC"), .groups="drop") %>%
        filter(n_iec==1, n_ctl>=1) %>% pull(set_id)
      matched <- matched %>% filter(set_id %in% keep_sets)
    }
  }

  if (!allow_replacement) {
    dup_ctrls <- matched %>% filter(IEC_flag=="No IEC") %>% count(Study_ID) %>% filter(n > 1)
    if (nrow(dup_ctrls) > 0) {
      stop(sprintf("Control reuse detected despite allow_replacement=FALSE: %s",
                   paste(dup_ctrls$Study_ID, collapse=", ")))
    }
  }

  set_summary <- matched %>%
    group_by(set_id) %>%
    summarise(
      IEC_val = CAR_abs_imp[IEC_flag=="IEC"][1],
      ctrl_gm = exp(mean(log(CAR_abs_imp[IEC_flag=="No IEC"]))),
      diff_log = log(IEC_val) - log(ctrl_gm),
      n_ctrl = sum(IEC_flag=="No IEC"),
      IEC_ID = Study_ID[IEC_flag=="IEC"][1],
      IEC_Day = Day[IEC_flag=="IEC"][1],
      .groups="drop"
    )

  wt_two    <- wilcox.test(CAR_abs_imp ~ IEC_flag, data = matched, exact = FALSE)
  wt_paired <- wilcox.test(set_summary$diff_log, mu = 0, exact = FALSE,
                           conf.int = TRUE, conf.level = 0.95)

  list(
    matched = matched,
    set_summary = set_summary,
    wt_two = wt_two,
    wt_paired = wt_paired,
    scarcity = scarcity_tbl %>% select(Study_ID, min_avail, feasible_days, first_day)
  )
}

# Matched comparison and diagnostics
caliper_fixed <- function(day) caliper_days

res <- match_iec_unique(
  exp_all,
  day_low = day_low, day_high = day_high,
  max_ratio = max_ratio,
  caliper_fun = caliper_fixed,
  allow_replacement = allow_replacement,
  pd_censor_pool = pd_censor_pool,
  pd_censor_post = pd_censor_post
)

matched     <- res$matched
set_summary <- res$set_summary
wt_two      <- res$wt_two
wt_paired   <- res$wt_paired

hl_est <- as.numeric(wt_paired$estimate)
hl_ci  <- as.numeric(wt_paired$conf.int)
fc_est <- exp(hl_est); fc_lcl <- exp(hl_ci[1]); fc_ucl <- exp(hl_ci[2])

cat(sprintf("\nMatched sets: %d   Controls used (unique patients): %d\n",
            dplyr::n_distinct(matched$set_id),
            dplyr::n_distinct(matched$Study_ID[matched$IEC_flag=='No IEC'])))
cat(sprintf("Two-sample Wilcoxon (matched rows): p = %.4g\n", wt_two$p.value))
cat(sprintf("Paired signed-rank (IEC vs GM controls): p = %.4g\n", wt_paired$p.value))
cat(sprintf("Median fold-change (IEC / controls): %.2f (95%% CI %.2f–%.2f) [HL]\n",
            fc_est, fc_lcl, fc_ucl))

by_set <- matched %>%
  group_by(set_id) %>%
  mutate(IEC_Day = Day[IEC_flag=="IEC"][1],
         IEC_ID  = Study_ID[IEC_flag=="IEC"][1],
         diff_days = if_else(IEC_flag=="No IEC", Day - IEC_Day, NA_real_)) %>%
  ungroup()

cat("\nPer-set day diffs (control − IEC):\n")
print(by_set %>% filter(IEC_flag=="No IEC") %>%
        select(set_id, IEC_ID, IEC_Day, Study_ID, Day, diff_days) %>%
        arrange(set_id, diff_days))

cat("\nPD prior to matched day (should be all FALSE after censoring):\n")
print(matched %>% select(set_id, IEC_flag, Study_ID, Day, PD_day, pd_prior_match) %>%
        arrange(set_id, IEC_flag))

iec_all_ids <- exp_all %>% filter(IEC_flag=="IEC", Day > day_low, Day <= day_high) %>% distinct(Study_ID) %>% pull()
iec_m_ids   <- set_summary$IEC_ID
iec_unmatched_ids <- setdiff(iec_all_ids, iec_m_ids)
cat(sprintf("\nIEC IDs unmatched (n=%d): %s\n",
            length(iec_unmatched_ids),
            ifelse(length(iec_unmatched_ids)==0, "None",
                   paste(sort(iec_unmatched_ids), collapse=", "))))
cat("\nIEC scarcity summary (min controls available across candidate days):\n")
print(res$scarcity %>% arrange(min_avail, first_day))

n_labs <- matched %>%
  count(IEC_flag, name="n") %>%
  tidyr::complete(IEC_flag=factor(c("No IEC","IEC"), levels=c("No IEC","IEC")),
                  fill=list(n=0)) %>%
  mutate(lbl = paste0(as.character(IEC_flag), "\n(n=", n, ")"))
x_labels <- setNames(n_labs$lbl, n_labs$IEC_flag)

# Persistence figures
p_box <- ggplot(matched, aes(IEC_flag, CAR_abs_imp, fill = IEC_flag)) +
  geom_boxplot(outlier.shape=NA, width=.62, alpha=.85) +
  geom_point(position = position_identity(),
             alpha=.80, size=2, color="gray20") +
  scale_fill_manual(values=pal, guide="none") +
  scale_x_discrete(labels = x_labels) +
  scale_y_log10("CAR cells/µL (first matched post-Day-50; ±14 d)",
                labels = label_number(accuracy=1, big.mark=",")) +
  geom_hline(yintercept = LOD_ND, linetype = "dashed") +
  ggtitle(sprintf("Time-matched CAR persistence (No-IEC:IEC ≤ %d:1, ±%d d)\nHL FC: %.2f (95%% CI %.2f–%.2f)",
                  max_ratio, caliper_days, fc_est, fc_lcl, fc_ucl)) +
  theme_arial(22)
print(p_box)

iec_points <- matched %>% filter(IEC_flag=="IEC") %>% transmute(set_id, IEC_y = CAR_abs_imp)
segments_df <- matched %>%
  filter(IEC_flag=="No IEC") %>%
  left_join(iec_points, by="set_id") %>%
  transmute(set_id,
            x_start = factor("No IEC", levels=c("No IEC","IEC")),
            x_end   = factor("IEC",    levels=c("No IEC","IEC")),
            y_start = CAR_abs_imp,
            y_end   = IEC_y)

p_box_linked <- ggplot() +
  geom_boxplot(data = matched, aes(IEC_flag, CAR_abs_imp, fill = IEC_flag),
               outlier.shape = NA, width = 0.62, alpha = 0.70) +
  geom_segment(data = segments_df,
               aes(x = x_start, xend = x_end, y = y_start, yend = y_end, group = set_id),
               linewidth = 0.6, alpha = 0.40, color = "gray35") +
  geom_point(data = matched, aes(IEC_flag, CAR_abs_imp),
             position = position_identity(),
             alpha = 0.80, size = 2, color = "gray20") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = x_labels) +
  scale_y_log10("CAR cells/µL (first matched post-Day-50; ±14 d)",
                labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  geom_hline(yintercept = LOD_ND, linetype = "dashed") +
  ggtitle("Matched persistence with IEC–control links") +
  theme_arial(22)
print(p_box_linked)

save_pdf(p_spaghetti_0_100, "Exp_spaghetti_0_100.pdf", width = 11, height = 6.2)

save_pdf(p_box, "Matched_box_nojitter.pdf", width = 7.5, height = 7.5)

save_pdf(p_box_linked, "Matched_box_linked_nojitter.pdf", width = 7.5, height = 7.5)
