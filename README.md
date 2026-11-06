# Diagnostic delay in lumbar brucellar spondylitis: IPD meta-analysis of published cases

Analysis code for the manuscript:

> Factors associated with diagnostic delay in lumbar brucellar spondylitis: an individual patient data meta-analysis of published cases. 

Authors: Lin Yang, Yu Fu, Yang Yang, Mohamad Zaki bin Haji Mohd Amin, Jeffery anak Stephen, Yang Liu

Contact: Lin Yang (ORCID 0009-0006-5574-8106)

## Overview

We pooled individual patient data from published case reports and case series of lumbar brucellar spondylitis. The analysis has three parts:

- Diagnostic delay and the diagnoses made before brucellosis was confirmed (descriptive).
- Factors associated with diagnostic delay: a log-normal accelerated failure time model with cluster-robust (CR2) standard errors by report.
- The association of diagnostic delay with epidural abscess and neurological deficit: logistic regression.

The study is reported according to PRISMA-IPD.

## Repository structure

```
scripts/          R scripts (01–23)
scripts/python/   Python data cleaning script (02_clean_data.py)
outputs/          tables and figures produced by the scripts
```

## Software

- R version [TO CHECK: e.g. 4.4.1]
- Main R packages: [TO CHECK against sessionInfo(), e.g. survival, clubSandwich, logistf, rms, quantreg, lme4, lmerTest, glmnet, flextable, officer]
- Python [TO CHECK: version] for data cleaning (scripts/python)

Run `sessionInfo()` in R and replace this list with the exact package versions.

## How to run

1. Run `scripts/python/02_clean_data.py`. It reads `data_cleaned_final.xlsx` [TO CHECK] and writes `ipd_clean.csv`, which all R scripts use.
2. Run the R scripts in numerical order (01 to 23).

Script 23 combines the supplementary tables and figures into Online Resource 5. It reads the outputs of the earlier scripts, so it must be run last.

Paths are relative to the repository root. [TO CHECK: no setwd() with local paths]

## Scripts

### Python

| Script | Output |
|---|---|
| python/02_clean_data.py | Cleans the extraction file and derives the analysis variables, including the broad definition of neurological deficit; writes `ipd_clean.csv` |

### R

| Script | Output |
|---|---|
| 01_fig1_prisma.R | Fig. 1 (PRISMA-IPD flow diagram) |
| 02_fig3_jbi.R | Fig. 3 (quality assessment, JBI checklist) |
| 03_tableS1_fig2_map.R | Fig. 2 (map of countries); Table S1 (endemic classification) |
| 04_[TO CHECK].R | [TO CHECK, e.g. Table 1] |
| 05_[TO CHECK].R | [TO CHECK, e.g. main model, Table 2; VIF, Table S5] |
| 06_[TO CHECK].R | [TO CHECK, e.g. secondary analyses, Table 3] |
| 07_S1_report_type.R | Sensitivity analysis 1: report type (Table S9) |
| 08_S2_quantile.R | Sensitivity analysis 2: quantile regression (Table S10) |
| 09_S3_distributions.R | Sensitivity analysis 3: Weibull and log-logistic models (Table S11) |
| 10_S4_jbi.R | Sensitivity analysis 4: report quality (Table S12) |
| 11_S5_classification_year.R | Sensitivity analysis 5: endemic classification and year of publication (Table S13; Fig. S4) |
| 12_S6_no_fever.R | Sensitivity analysis 6: model without fever (Table S16) |
| 13_S7_fever_imputation.R | Sensitivity analysis 7: multiple imputation of fever (Table S15) |
| 14_S8_leave_one_out.R | Sensitivity analysis 8: removal of large case series, and post hoc removal of surgical cohorts (Table S14; Fig. S3) |
| 15_S9_variable_selection.R | Sensitivity analysis 9: LASSO and AIC-based selection (Table S17) |
| 16_figS2_diagnostics.R | Model diagnostics (Fig. S2; Table S7) |
| 17_S10_mixed_model.R | Unplanned: mixed-effects models (Table S18) |
| 18_S11_report_weighted.R | Unplanned: equal weighting of reports (Table S19) |
| 19_S12_period_interaction.R | Unplanned: interactions with publication period (Table S20) |
| 20_extraction_agreement.R | Extraction agreement (Table S2) |
| 21_tableS4_prior_diagnoses.R | Diagnoses before confirmation (Table S4) |
| 22_tableS21_neuro_unplanned.R | Unplanned: Firth regression and broad definition of neurological deficit (Table S21) |
| 23_supplement.R | Combines all supplementary tables and figures (Online Resource 5) |

## Data

The extracted individual patient dataset is not included in this repository. All data were extracted from the published reports listed in Online Resource 7 of the article. The data extraction manual is provided as Online Resource 4. The scripts therefore cannot be run without the input files below.

### Input files (not included)

| File | Content | Used by |
|---|---|---|
| `data_cleaned_final.xlsx` | Extracted individual patient data | `python/02_clean_data.py` [TO CHECK] |
| `endemic_classification.csv` | Endemic status of each reporting country (shown in Table S1) | [TO CHECK: e.g. 02_clean_data.py, 03_tableS1_fig2_map.R] |
| `jbi/` | JBI checklist ratings for each patient | [TO CHECK: e.g. 02_fig3_jbi.R, 10_S4_jbi.R] |
| `reextraction/Reextraction_form_filled.xlsx` | Independent re-extraction of 25 patients by the second reviewer | `20_extraction_agreement.R` |
| `reextraction/Neuro_agreement_form_25_filled.xlsx` | Independent coding of neurological deficit for the same 25 patients | `20_extraction_agreement.R` |
| `reextraction/Reextraction_KEY_original_values_DO_NOT_SEND.xlsx` | Original values for comparison with the re-extraction | `20_extraction_agreement.R` [TO CHECK] |

## License

MIT License (see the LICENSE file).

## Citation

If you use this code, please cite the article (citation to be added after publication) and the archived version of this repository on Zenodo: [DOI to be added].
