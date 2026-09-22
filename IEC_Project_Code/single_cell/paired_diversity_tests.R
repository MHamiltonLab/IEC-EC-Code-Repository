#!/usr/bin/env Rscript
# Paired tests of patient-level repertoire summaries.
# Columns of a two-row matrix are matched patients, never subsampling iterations.

paired_t_from_matrix <- function(values, metric = "metric", conf_level = 0.95) {
  if (!is.matrix(values) || nrow(values) != 2L || !is.numeric(values)) {
    stop("values must be a numeric matrix with reference and comparison rows.")
  }
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) stop("conf_level must be between 0 and 1.")
  keep <- is.finite(values[1L, ]) & is.finite(values[2L, ])
  reference <- values[1L, keep]
  comparison <- values[2L, keep]
  n <- length(reference)
  out <- data.frame(
    metric = metric, n_pairs = n, alternative = "two.sided",
    adjustment = "none", confidence_level = conf_level,
    mean_reference = if (n) mean(reference) else NA_real_,
    mean_comparison = if (n) mean(comparison) else NA_real_,
    mean_difference = if (n) mean(comparison - reference) else NA_real_,
    t_statistic = NA_real_, df = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    status = "insufficient_pairs", stringsAsFactors = FALSE
  )
  if (n < 2L) return(out)
  fit <- tryCatch(
    stats::t.test(comparison, reference, paired = TRUE,
                  alternative = "two.sided", conf.level = conf_level),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out$status <- paste("not_estimable:", conditionMessage(fit))
    return(out)
  }
  out$t_statistic <- unname(fit$statistic)
  out$df <- unname(fit$parameter)
  out$ci_low <- unname(fit$conf.int[1L])
  out$ci_high <- unname(fit$conf.int[2L])
  out$p_value <- fit$p.value
  out$status <- "ok"
  out
}

paired_diversity_tests <- function(data, metrics = c("shannon", "inv_simpson", "clonality"),
                                   reference = "Blood", comparison = "Ileum",
                                   conf_level = 0.95) {
  required <- c("patient_id", "Tissue", metrics)
  if (!all(required %in% names(data))) {
    stop("Missing fields: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (!nrow(data)) stop("Supply patient-level median summaries.")
  if (identical(reference, comparison)) stop("Select two different tissues.")
  data$patient_id <- as.character(data$patient_id)
  data$Tissue <- sub("_downsampled$", "", as.character(data$Tissue))
  if (anyNA(data$patient_id) || any(!nzchar(trimws(data$patient_id))) || anyNA(data$Tissue)) {
    stop("Patient and tissue identifiers must be populated.")
  }
  if (!all(data$Tissue %in% c(reference, comparison))) stop("Unexpected tissue labels.")
  if (anyDuplicated(data[c("patient_id", "Tissue")])) {
    stop("Supply one summary per patient and tissue; repeated iterations are not independent pairs.")
  }
  first <- data[data$Tissue == reference, , drop = FALSE]
  second <- data[data$Tissue == comparison, , drop = FALSE]
  patients <- intersect(first$patient_id, second$patient_id)
  first <- first[match(patients, first$patient_id), , drop = FALSE]
  second <- second[match(patients, second$patient_id), , drop = FALSE]
  results <- lapply(metrics, function(metric) {
    if (!is.numeric(data[[metric]]) && !all(is.na(data[[metric]]))) {
      stop("Metric must be numeric: ", metric)
    }
    values <- rbind(as.numeric(first[[metric]]), as.numeric(second[[metric]]))
    colnames(values) <- patients
    out <- paired_t_from_matrix(values, metric, conf_level)
    out$reference <- reference
    out$comparison <- comparison
    out
  })
  do.call(rbind, results)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) {
    stop("Usage: Rscript single_cell/paired_diversity_tests.R [medians.tsv] [results.tsv]")
  }
  input_path <- if (length(args) >= 1L) args[[1L]] else file.path(
    Sys.getenv("METHODS_INPUT_DIR", unset = "inputs"), "paired_diversity_medians.tsv"
  )
  output_path <- if (length(args) >= 2L) args[[2L]] else file.path(
    Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "single_cell",
    "paired_diversity_tests", "paired_t_tests.tsv"
  )
  data <- utils::read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  result <- paired_diversity_tests(data)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(result, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  print(result, digits = 8)
}
