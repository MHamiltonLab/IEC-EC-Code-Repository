# Empty input objects for the analysis interfaces.
# This file contains no observations and performs no file I/O.

empty_table <- function(character_fields = character(), numeric_fields = character(),
                        integer_fields = character()) {
  fields <- c(
    setNames(lapply(character_fields, function(x) character()), character_fields),
    setNames(lapply(numeric_fields, function(x) numeric()), numeric_fields),
    setNames(lapply(integer_fields, function(x) integer()), integer_fields)
  )
  as.data.frame(fields, check.names = FALSE, stringsAsFactors = FALSE)
}

input_schemas <- list(
  clinical_metadata = empty_table(
    character_fields = c(
      "Study_ID", "Sex", "Race", "Ethnicity", "Heavy_Chain", "Light_Chain", "RISS",
      "Any_Delayed_Neurotoxicity", "Parkinsonian_Neurotoxicity", "Other_Neurotoxicity",
      "What_Neurotoxicity"
    ),
    numeric_fields = c(
      "IEC_enteritis", "Age_at_infusion", "ECOG_at_apheresis", "Prior_Lines",
      "ALC_PreLD", "preLD_M-spike", "Baseline_pre_LD_LDH", "Baseline_pre-LD_CRP",
      "Baseline_pre-LD_Ferritin", "Max_ICANS", "Max_CRS_grade", "PD_day",
      "Death_day", "Last_follow_up_day", "IEC_Day", "DNT_day"
    )
  ),
  clinical_outcomes = empty_table(
    character_fields = c("ID", "CAR", "TPN", "Mortality"),
    numeric_fields = "HospitalDays",
    integer_fields = "PostEndoscopyInfections"
  ),
  car_measurements = empty_table(
    character_fields = "Study_ID", numeric_fields = c("Day", "CAR_abs")
  ),
  lymphocyte_measurements = empty_table(
    character_fields = "Study_ID", numeric_fields = c("Day", "ALC")
  ),
  tissue_cell_counts = empty_table(
    character_fields = c("Patient_ID", "Cohort", "Tissue_Type"),
    numeric_fields = c(
      "total_cells", "CD3", "CD68", "CD3-CD103", "CD3_Camelid",
      "CD3-Camelid-CD4", "CD3-Camelid-CD8", "CD3-Camelid-DP", "CD3-Camelid-DN",
      "CD3-Camelid-GZMB", "CD3-Camelid-PD1", "CD3-Camelid-negative",
      "CD3-Camelid-negative-CD4", "CD3-Camelid-negative-CD8",
      "CD3-Camelid-negative-DP", "CD3-Camelid-negative-DN",
      "CD3-Camelid-negative-GZMB", "CD3-Camelid-negative-KI67",
      "CD3-Camelid-negative-FOXP3", "CD3-Camelid-negative-PD1",
      "pct_CD3-Camelid-CD4", "pct-CD3-Camelid-CD8", "pct-CD3-Camelid-DP",
      "pct-CD3-Camelid-DN", "pct-CD3-Camelid-GZMB", "pct-CD3-Camelid-KI67",
      "pct-CD3-Camelid-PD1", "pct-CD3-Camelid-negative",
      "pct_CD3-Camelid-negative-CD4", "pct_CD3-Camelid-negative-CD8",
      "pct_CD3-Camelid-negative_DP", "pct_CD3-Camelid-negative_DN",
      "pct_CD3-Camelid-negative-GZMB", "pct_CD3-Camelid-negative-KI67",
      "pct_CD3-Camelid-negative-FOXP3", "pct_CD3-Camelid-negative-PD1",
      "CD3_fraction_camelid", "CD4:CD8", "CD103_fraction_tcells",
      "Total_cell_fraction_camelid"
    )
  ),
  patient_labels = empty_table(character_fields = c("ID", "Paper_ID")),
  repertoire_metadata = empty_table(
    character_fields = c("Sample_ID", "Paper_ID", "Type", "Timepoint", "Pathology"),
    numeric_fields = "Day"
  ),
  repertoire_pairs = empty_table(
    character_fields = c("comparison_id", "Sample", "title")
  ),
  single_cell_metadata = empty_table(
    character_fields = c("sample_id", "subject_id", "Tissue")
  ),
  expression_metadata = empty_table(
    character_fields = c(
      "Sample_RNAseq", "Cohort", "Tissue", "Timepoint", "Paper_ID", "source_id", "Insertion"
    ),
    numeric_fields = c("COMET_perc_CD3", "COMET_perc_total")
  ),
  expression_fpkm = empty_table(
    character_fields = c("GeneID", "Gene_Name", "Gene_Biotype")
  ),
  splicing_metadata = empty_table(
    character_fields = c("sample", "patient_id", "group", "gene", "chr", "paper_sample_id"),
    integer_fields = c("start", "end")
  ),
  insertion_coordinates = empty_table(
    character_fields = c("sample", "gene", "chrom"),
    integer_fields = c("start", "end", "position")
  ),
  molecular_features = empty_table(
    character_fields = c(
      "sample_id", "Grouping", "TCL_Genotype", "Product", "CAR_Detected", "CD4_CD8", "Location"
    )
  ),
  molecular_comparison = empty_table(
    character_fields = c("Group", "Outcome"), integer_fields = "n"
  ),
  insertion_clusters = empty_table(
    character_fields = c("sample", "insertion_id", "chrom", "vec_contig_mode"),
    integer_fields = c(
      "cluster_start_1based", "cluster_end_1based", "breakpoint_mid_1based",
      "support_total_unique_qname", "support_split_unique_qname",
      "support_discordant_unique_qname", "vec_pos_mode"
    )
  )
)

# Add sample-named numeric columns to expression_fpkm when supplying real inputs.
# Count fields remain numeric where the source analyses allow imported noninteger values.
