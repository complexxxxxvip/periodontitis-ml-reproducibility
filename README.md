# Reproducibility materials

## Study

**Development and External Validation of an Interpretable Machine-Learning Model for Severity Stratification of Current Periodontitis Using Routine Blood Biomarkers, with Exploratory Analyses of Age- and Sex-Specific Model Behavior**

This repository contains code and supporting metadata for the final model-development/refitting workflow, baseline-model comparison, external validation, SHAP analysis, and exploratory SHAP/RCS subgroup analyses.

## Recommended repository files

- `01_development_model.R` — final 16-predictor SVM + XGBoost stacking model, leakage-controlled multiple imputation, tuning, calibration, and full-development refitting.
- `02_baseline_models.R` — four unpenalized logistic-regression baseline models.
- `03_external_validation.R` — external validation using frozen development-derived preprocessing/model/calibration/threshold settings.
- `04_SHAP_analysis.R` — SHAP-based model interpretation.
- `05_SHAP_RCS.R` — exploratory sex- and age-stratified SHAP/RCS analyses.
- `06_functional_demo.R` — synthetic-data functional demonstration only.
- `data_dictionary.csv` — variable definitions, units, coding, and derivation rules.
- `hyperparameter_search_space.csv` — final-model tuning ranges and selected values.
- `synthetic_example_data.csv` — fully simulated example data; no real patient rows.
- `model_specification_public.txt` — compact final-model specification.
- `session_info_public.txt` — public software/package-version information with system-specific details removed.
- `external_validation_summary_public.txt` — aggregate external-validation summary.
- `baseline_model_summary_public.xlsx` — aggregate-only baseline-model results; patient-level OOF/fold sheets removed.
- `UPLOAD_CHECKLIST.md` — public-release safety checklist.

## Data availability and privacy

Patient-level development and external-validation data are not included because of privacy and ethical restrictions.

Fitted model objects are also not publicly released because serialized objects may retain patient-level feature information from the development cohort.

The public materials therefore provide the analysis code, model specification, final-model hyperparameter search settings, data dictionary, aggregate results, and a synthetic-data functional demonstration.

## Expected local folder structure

```text
project/
  01_development_model.R
  02_baseline_models.R
  03_external_validation.R
  04_SHAP_analysis.R
  05_SHAP_RCS.R
  06_functional_demo.R

  data_dictionary.csv
  hyperparameter_search_space.csv
  synthetic_example_data.csv
  model_specification_public.txt
  session_info_public.txt
  external_validation_summary_public.txt
  baseline_model_summary_public.xlsx
  UPLOAD_CHECKLIST.md

  private_data/          # not included in the public repository
  outputs/               # generated locally
```

Authorized users wishing to reproduce the study analyses with the original data should place the required private datasets in `private_data/` using the filenames expected by the public scripts.

## Functional demonstration

Run:

```r
source("06_functional_demo.R")
```

The demonstration uses completely simulated data to show the expected source-variable structure, multiple imputation, deterministic calculation of derived predictors, and the SVM + XGBoost stacking architecture.

**Important:** the synthetic dataset does not represent real participants and cannot reproduce the numerical AUC, calibration, or classification results reported in the manuscript. Demonstration predictions are not study results.

## Final model inputs

Final predictors:

`Sex`, `Smoke`, `Age`, `AG_ratio`, `eGFR`, `Fbg`, `PAR`, `TC`, `GGT`, `ALP`, `NPR`, `INR`, `HB`, `APTT`, `LDH`, and `TyG`.

MICE source variables:

`Sex`, `Smoke`, `Age`, `TP`, `ALB`, `CR`, `Fbg`, `PLT`, `TC`, `GGT`, `ALP`, `NE_C`, `PT`, `HB`, `APTT`, `LDH`, `GLU`, and `TG`.

Derived after each completed imputation and never directly imputed:

- `AG_ratio = ALB / (TP - ALB)`
- `PAR = PLT / ALB`
- `NPR = NE_C / PLT`
- `INR = (PT / 13.1)^1.31`
- `TyG = log(1594.26 * TG * GLU / 2)`
- `eGFR` is calculated from creatinine, age, and sex as specified in `data_dictionary.csv`.

## Outcome coding

- `Stage0`: Stage I-II periodontitis
- `Stage1`: Stage III-IV periodontitis

The model addresses current periodontal severity rather than future disease progression.

## Software

Analyses used R 4.5.1. Package versions recovered from the supplied analysis-session records are listed in `session_info_public.txt`.

## Important reproducibility scope

The supplied final-development script refits and evaluates the **prespecified final 16-predictor SVM + XGBoost stacking architecture** on the complete development cohort.

The earlier model-specification stage described in the manuscript — including screening of the initial 70 candidate variables, comparison of 12 individual algorithms, and selection among 57 stacking architectures — is not recreated by these five final-analysis scripts unless the corresponding earlier model-selection scripts are also deposited.

If those earlier scripts are available, adding them would provide the strongest response to a request for the **full modelling code**.

## Citation

Repository URL: `[ADD PUBLIC REPOSITORY URL]`

Article DOI / citation: `[ADD AFTER AVAILABLE]`
