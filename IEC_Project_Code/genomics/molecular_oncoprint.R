#!/usr/bin/env Rscript
# Molecular feature comparison
# Build an annotated oncoprint and evaluate a supplied two-by-two comparison.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(wesanderson)
})

data_dir <- Sys.getenv("METHODS_INPUT_DIR", unset = "inputs")
# Molecular feature inputs
INFILE <- file.path(data_dir, "molecular_features.tsv")

OUTDIR <- file.path(Sys.getenv("METHODS_OUTPUT_DIR", unset = "results"), "genomics", "molecular_oncoprint")

GROUP_LEVELS <- c("IEC-EC", "PTCL")

PDF_WIDTH          <- 12
PDF_HEIGHT         <- 5.5
HEATMAP_BODY_HEIGHT <- 3.7

FONT_FAM <- "Arial"

MUT_COLOR    <- "#5C6D7E"
TET2_ORANGE  <- "#E69F00"
UNKNOWN_FILL <- "#D9D9D9"

group_colors <- c(
  "IEC-EC" = "#4C78A8",
  "PTCL"   = "#9E9E9E"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (capabilities("aqua")) {
  quartzFonts(
    Arial = quartzFont(
      c(
        "Arial",
        "Arial Bold",
        "Arial Italic",
        "Arial Bold Italic"
      )
    )
  )
}

save_heatmap_pdf <- function(
    ht,
    filename,
    width = PDF_WIDTH,
    height = PDF_HEIGHT
) {
  outfile <- file.path(OUTDIR, filename)

  quartz(
    file = outfile,
    type = "pdf",
    width = width,
    height = height,
    family = FONT_FAM
  )

  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )

  dev.off()

  message("Saved editable-text PDF: ", outfile)
}

meta_df <- read_delim(
  INFILE,
  delim = "\t",
  trim_ws = TRUE,
  show_col_types = FALSE
)

meta_df <- meta_df %>%
  rename_with(
    ~ str_replace_all(.x, "[\\s/-]+", "_")
  )

required_cols <- c(
  "sample_id",
  "Grouping",
  "TCL_Genotype"
)

missing_cols <- setdiff(
  required_cols,
  colnames(meta_df)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\n\nColumns present:\n",
    paste(colnames(meta_df), collapse = ", ")
  )
}

meta_df <- meta_df %>%
  mutate(
    across(
      where(is.character),
      ~ str_trim(.x)
    )
  ) %>%
  mutate(
    across(
      where(is.character),
      ~ replace(
        .x,
        is.na(.x) |
          .x == "" |
          tolower(.x) %in% c("na", "nan"),
        NA_character_
      )
    )
  )

if (!"Product" %in% names(meta_df)) {
  meta_df$Product <- "None"
}

if (!"CAR_Detected" %in% names(meta_df)) {
  meta_df$CAR_Detected <- "Unknown"
}

if (!"CD4_CD8" %in% names(meta_df)) {
  meta_df$CD4_CD8 <- "Unknown"
}

if (!"Location" %in% names(meta_df)) {
  meta_df$Location <- "Unknown"
}

meta_df <- meta_df %>%
  mutate(
    Product = coalesce(Product, "None"),
    CAR_Detected = coalesce(CAR_Detected, "Unknown"),
    CD4_CD8 = coalesce(CD4_CD8, "Unknown"),
    Location = coalesce(Location, "Unknown")
  )

meta_df <- meta_df %>%
  mutate(
    Grouping = factor(
      as.character(Grouping),
      levels = GROUP_LEVELS
    )
  ) %>%
  filter(!is.na(Grouping)) %>%
  arrange(Grouping, sample_id)

if (nrow(meta_df) == 0) {
  stop(
    "No cases remained after filtering for: ",
    paste(GROUP_LEVELS, collapse = ", ")
  )
}

if (anyDuplicated(meta_df$sample_id)) {
  duplicated_cases <- unique(
    meta_df$sample_id[
      duplicated(meta_df$sample_id)
    ]
  )

  stop(
    "sample_id must be unique. Duplicated cases: ",
    paste(duplicated_cases, collapse = ", ")
  )
}

split_genes <- function(x) {
  if (is.na(x)) {
    return(character(0))
  }

  x <- str_trim(as.character(x))

  if (
    x == "" ||
    tolower(x) %in% c("none", "na", "nan", "unknown")
  ) {
    return(character(0))
  }

  genes <- unlist(
    str_split(x, ",")
  )

  genes <- str_trim(genes)
  genes <- genes[genes != ""]
  genes <- genes[
    !tolower(genes) %in% c("none", "na", "nan", "unknown")
  ]

  unique(genes)
}

gene_list <- lapply(
  meta_df$TCL_Genotype,
  split_genes
)

all_genes <- sort(
  unique(
    unlist(gene_list)
  )
)

if (length(all_genes) == 0) {
  stop(
    paste0(
      "No genes were parsed from TCL_Genotype. ",
      "All entries may be missing, None, or NA."
    )
  )
}

mut_matrix <- matrix(
  "",
  nrow = length(all_genes),
  ncol = nrow(meta_df),
  dimnames = list(
    all_genes,
    meta_df$sample_id
  )
)

for (i in seq_along(gene_list)) {
  genes_i <- intersect(
    gene_list[[i]],
    rownames(mut_matrix)
  )

  if (length(genes_i) > 0) {
    mut_matrix[genes_i, i] <- "mutation"
  }
}

stopifnot(
  identical(
    colnames(mut_matrix),
    meta_df$sample_id
  )
)

iec_columns <- meta_df$Grouping == "IEC-EC"

if (any(iec_columns)) {
  iec_gene_counts <- rowSums(
    mut_matrix[
      ,
      iec_columns,
      drop = FALSE
    ] == "mutation"
  )
} else {
  iec_gene_counts <- rep(
    0,
    nrow(mut_matrix)
  )
}

total_gene_counts <- rowSums(
  mut_matrix == "mutation"
)

gene_order <- order(
  iec_gene_counts > 0,
  iec_gene_counts,
  total_gene_counts,
  decreasing = TRUE
)

mut_matrix <- mut_matrix[
  gene_order,
  ,
  drop = FALSE
]

if ("TET2" %in% rownames(mut_matrix)) {
  tet2_columns <- mut_matrix["TET2", ] == "mutation"

  mut_matrix[
    "TET2",
    tet2_columns
  ] <- "TET2_mut"
}

make_product_palette <- function(values) {
  levels_present <- unique(
    as.character(values)
  )

  n_levels <- length(levels_present)

  if (n_levels <= 5) {
    palette_values <- wes_palette(
      "GrandBudapest1",
      n = n_levels,
      type = "discrete"
    )
  } else {
    palette_values <- colorRampPalette(
      wes_palette(
        "GrandBudapest1",
        n = 5,
        type = "discrete"
      )
    )(n_levels)
  }

  palette_values <- setNames(
    palette_values,
    levels_present
  )

  gray_levels <- intersect(
    names(palette_values),
    c("None", "Unknown")
  )

  palette_values[gray_levels] <- UNKNOWN_FILL

  palette_values
}

complete_palette <- function(values, palette) {
  observed_levels <- unique(
    as.character(values)
  )

  missing_levels <- setdiff(
    observed_levels,
    names(palette)
  )

  if (length(missing_levels) > 0) {
    palette <- c(
      palette,
      setNames(
        rep(UNKNOWN_FILL, length(missing_levels)),
        missing_levels
      )
    )
  }

  palette[observed_levels]
}

product_colors <- make_product_palette(
  meta_df$Product
)

car_colors <- c(
  "Yes"     = "#4C7A6A",
  "No"      = "#7A7A7A",
  "Unknown" = UNKNOWN_FILL
)

car_colors <- complete_palette(
  meta_df$CAR_Detected,
  car_colors
)

location_colors <- c(
  "Institution_A"     = "#6F7C91",
  "Other"   = "#9A8F86",
  "Unknown" = UNKNOWN_FILL
)

location_colors <- complete_palette(
  meta_df$Location,
  location_colors
)

cd4_cd8_colors <- c(
  "CD4+"    = "#6B8F9C",
  "CD8+"    = "#8FA6A1",
  "DN"      = "#A08C7C",
  "Mixed"   = "#7E7D8A",
  "Unknown" = UNKNOWN_FILL,
  "NA"      = UNKNOWN_FILL
)

cd4_cd8_colors <- complete_palette(
  meta_df$CD4_CD8,
  cd4_cd8_colors
)

top_annot <- HeatmapAnnotation(
  Grouping = meta_df$Grouping,
  Product = meta_df$Product,
  CAR_Detected = meta_df$CAR_Detected,
  Location = meta_df$Location,

  col = list(
    Grouping = group_colors,
    Product = product_colors,
    CAR_Detected = car_colors,
    Location = location_colors

  ),

  simple_anno_size = unit(3.5, "mm"),
  gap = unit(0.8, "mm"),

  annotation_name_side = "left",

  annotation_name_gp = gpar(
    fontfamily = FONT_FAM,
    fontsize = 8
  ),

  annotation_legend_param = list(
    title_gp = gpar(
      fontfamily = FONT_FAM,
      fontface = "bold",
      fontsize = 9
    ),
    labels_gp = gpar(
      fontfamily = FONT_FAM,
      fontsize = 8
    ),
    grid_height = unit(3.5, "mm"),
    grid_width = unit(3.5, "mm")
  )
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(
      x = x,
      y = y,
      width = w,
      height = h,
      gp = gpar(
        fill = "white",
        col = "white"
      )
    )
  },

  mutation = function(x, y, w, h) {
    grid.rect(
      x = x,
      y = y,
      width = 0.9 * w,
      height = 0.9 * h,
      gp = gpar(
        fill = MUT_COLOR,
        col = "white",
        lwd = 0.5
      )
    )
  },

  TET2_mut = function(x, y, w, h) {
    grid.rect(
      x = x,
      y = y,
      width = 0.9 * w,
      height = 0.9 * h,
      gp = gpar(
        fill = TET2_ORANGE,
        col = "white",
        lwd = 0.5
      )
    )
  }
)

mutation_colors <- c(
  "mutation" = MUT_COLOR,
  "TET2_mut" = TET2_ORANGE
)

ht <- oncoPrint(
  mut_matrix,

  alter_fun = alter_fun,
  alter_fun_is_vectorized = FALSE,
  col = mutation_colors,

  height = unit(
    HEATMAP_BODY_HEIGHT,
    "in"
  ),

  top_annotation = top_annot,

  column_title = "IEC-EC vs PTCL: Genotype Oncoprint",

  column_title_gp = gpar(
    fontfamily = FONT_FAM,
    fontface = "bold",
    fontsize = 12
  ),

  column_split = factor(
    meta_df$Grouping,
    levels = GROUP_LEVELS
  ),

  column_gap = unit(2, "mm"),

  show_column_names = TRUE,

  column_names_gp = gpar(
    fontfamily = FONT_FAM,
    fontsize = 8
  ),

  show_row_names = TRUE,
  row_names_side = "left",

  row_names_gp = gpar(
    fontfamily = FONT_FAM,
    fontsize = 8
  ),

  pct_side = "right",

  remove_empty_columns = FALSE,
  remove_empty_rows = TRUE,

  right_annotation = rowAnnotation(
    barplot = anno_oncoprint_barplot(
      border = FALSE
    ),

    width = unit(0.9, "in"),

    annotation_name_side = "bottom",

    annotation_name_gp = gpar(
      fontfamily = FONT_FAM,
      fontsize = 8
    )
  ),

  heatmap_legend_param = list(
    title = "Genotype",
    at = c(
      "mutation",
      "TET2_mut"
    ),
    labels = c(
      "Mutated",
      "TET2 mutation"
    ),
    title_gp = gpar(
      fontfamily = FONT_FAM,
      fontface = "bold",
      fontsize = 9
    ),
    labels_gp = gpar(
      fontfamily = FONT_FAM,
      fontsize = 8
    ),
    grid_height = unit(3.5, "mm"),
    grid_width = unit(3.5, "mm")
  ),

  border = TRUE
)

save_heatmap_pdf(
  ht = ht,
  filename = "PTCL_vs_IEC_oncoprint.pdf",
  width = PDF_WIDTH,
  height = PDF_HEIGHT
)

# Supplied two-by-two comparison
comparison <- readr::read_tsv(
  file.path(data_dir, "molecular_comparison.tsv"), show_col_types = FALSE
)
stopifnot(all(c("Group", "Outcome", "n") %in% names(comparison)))
stopifnot(all(is.finite(comparison$n)), all(comparison$n >= 0),
          all(comparison$n == floor(comparison$n)))
comparison$Group <- factor(comparison$Group, levels = c("G1", "G2"))
comparison$Outcome <- factor(comparison$Outcome, levels = c("Success", "Failure"))
stopifnot(!anyNA(comparison[c("Group", "Outcome")]))
tab <- xtabs(n ~ Group + Outcome, data = comparison, drop.unused.levels = FALSE)
stopifnot(all(rowSums(tab) > 0), all(colSums(tab) > 0))

ft_two <- fisher.test(
  tab,
  alternative = "two.sided",
  conf.level = 0.95
)

ft_less <- fisher.test(
  tab,
  alternative = "less"
)

ft_greater <- fisher.test(
  tab,
  alternative = "greater"
)

odds_ratio <- unname(
  ft_two$estimate
)

ci_lower <- ft_two$conf.int[1]
ci_upper <- ft_two$conf.int[2]

p_two_sided <- ft_two$p.value
p_less <- ft_less$p.value
p_greater <- ft_greater$p.value

cat("\n============================================================\n")
cat("2 x 2 TABLE\n")
cat("============================================================\n\n")

print(tab)

cat("\n============================================================\n")
cat("FISHER EXACT TEST RESULTS\n")
cat("============================================================\n")

cat(
  sprintf(
    paste0(
      "Odds ratio: %.4f\n",
      "Exact 95%% CI: [%.6f, %.6f]\n",
      "Two-sided p-value: %.3g\n",
      "One-sided p-value, G1 < G2: %.3g\n",
      "One-sided p-value, G1 > G2: %.3g\n"
    ),
    odds_ratio,
    ci_lower,
    ci_upper,
    p_two_sided,
    p_less,
    p_greater
  )
)

cat("============================================================\n\n")

cat(
  sprintf(
    paste0(
      "Fisher exact test: OR = %.4f, 95%% CI [%.6f, %.6f], ",
      "two-sided p = %.3g; one-sided p for G1 < G2 = %.3g; ",
      "one-sided p for G1 > G2 = %.3g.\n"
    ),
    odds_ratio,
    ci_lower,
    ci_upper,
    p_two_sided,
    p_less,
    p_greater
  )
)
