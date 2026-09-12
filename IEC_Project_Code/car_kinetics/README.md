# CAR kinetics

| Script | Function |
| --- | --- |
| [expansion_and_persistence.R](expansion_and_persistence.R) | Day 0–30 CAR exposure, sampling coverage, late measurements, covariate associations, and evaluability summaries. |
| [matched_persistence.R](matched_persistence.R) | Patient-level matching of post-day-50 measurements, with a ±14-day caliper, up to two controls per case, and no control reuse. |
| [delayed_neurotoxicity.R](delayed_neurotoxicity.R) | Delayed neurotoxicity associations with enterocolitis, lymphocyte counts, CAR exposure, and late persistence. |

**Inputs:** `clinical_metadata.tsv`, `car_measurements.tsv`, and, for neurotoxicity analyses, `lymphocyte_measurements.tsv`.

**Outputs:** Exposure and persistence summaries, covariate models, matched comparison plots, and neurotoxicity tables/figures under `results/car_kinetics/`.

CAR counts use cells/µL. Lymphocyte counts in the neurotoxicity script use thousands/µL, so the threshold of 3 corresponds to 3,000 cells/µL. The default CAR detection-limit imputation is 0.1 cells/µL. Exposure-window interpolation, minimum sampling requirements, and progression censoring are specified within each script; consult those definitions when comparing summaries.

Each script runs independently in a fresh R session. 