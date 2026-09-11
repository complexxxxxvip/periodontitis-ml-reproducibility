# ==============================================================================
# SHAP Analysis
# Final SVM + XGBoost stacking model with multiple imputation
#
# Development cohort: n = 1321
#
# Outputs:
#   1. Mean absolute SHAP bar plot
#   2. SHAP beeswarm summary plot
#   3. SHAP dependence plots colored by Sex
#   4. SHAP dependence plots colored by Age
#   5. SHAP values / feature importance / zero-crossings / diagnostics
#
# IMPORTANT:
#   - Model variables remain AG_ratio and GGT
#   - Plot labels are displayed as A/G and γ-GT
#   - No external standardization is performed
#   - Five saved MICE datasets are reused
#   - Five SVM + five XGBoost fitted models are pooled exactly as in final model
# ==============================================================================


rm(list = ls())
gc()

# Project-relative paths. Run this script from the repository root.
PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PRIVATE_DIR <- file.path(PROJECT_DIR, "private")
DATA_DIR <- file.path(PRIVATE_DIR, "data")
MODEL_DIR <- file.path(PRIVATE_DIR, "model_objects")
OUTPUT_ROOT <- file.path(PRIVATE_DIR, "analysis_outputs")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# The private/ directory is intended for local use only and should not be
# committed to the public repository.

# ==============================================================================
# 0. USER-ADJUSTABLE PARAMETERS
# ==============================================================================

SEED <- 2025


# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

RAW_DATA_PATH <- file.path(DATA_DIR, "development_data.xlsx")

BUNDLE_PATH <- file.path(MODEL_DIR, "final_stacking_model_bundle.rds")

RESULTS_PATH <- file.path(MODEL_DIR, "development_model_results.xlsx")


# Output folder
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "shap_analysis")

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------------------------
# Run mode
#
# Optional run modes:
#   "quick" = reduced Monte-Carlo repetitions for a short code check
#   "full"  = settings used for the reported analysis
# ------------------------------------------------------------------------------

RUN_MODE <- "full"


# FastSHAP Monte-Carlo repetitions
#
# quick = 10
# full = 500
#
# NOTE:
# The final model contains five SVM and five XGBoost fits; the full SHAP run is computationally intensive.
# ------------------------------------------------------------------------------

NSIM_QUICK <- 10

NSIM_FULL <- 500


NSIM <- ifelse(
  RUN_MODE == "quick",
  NSIM_QUICK,
  NSIM_FULL
)


# ------------------------------------------------------------------------------
# LOESS parameter for dependence plots
# ------------------------------------------------------------------------------

LOESS_SPAN <- 0.80


# ------------------------------------------------------------------------------
# Search range for zero-crossings
#
# Restrict to central 95% of observed feature range
# to reduce unstable crossings at extreme tails
# ------------------------------------------------------------------------------

ZERO_QUANTILES <- c(
  0.025,
  0.975
)


# ------------------------------------------------------------------------------
# If multiple zero-crossings exist:
#
# "closest_to_median"
#     = choose the crossing nearest the feature median for figure annotation
#
# ALL crossings will still be exported to Excel.
# ------------------------------------------------------------------------------

ZERO_SELECTION <- "closest_to_median"


# ------------------------------------------------------------------------------
# SHAP representation for original missing values
#
# "MI_mean":
# Missing continuous predictor values are represented by the mean
# of the five FINAL saved MICE values.
#
# Observed values are NEVER changed.
# ------------------------------------------------------------------------------

MISSING_REPRESENTATION <- "MI_mean"


# ------------------------------------------------------------------------------
# Optional MI sensitivity analysis
#
# FALSE = recommended for first/main run
#
# TRUE:
# calculate an additional SHAP importance analysis separately
# in all five completed MI datasets.
#
# This is computationally expensive.
# ------------------------------------------------------------------------------

RUN_MI_SENSITIVITY <- FALSE

MI_SENS_NSIM_QUICK <- 5

MI_SENS_NSIM_FULL <- 100


MI_SENS_NSIM <- ifelse(
  RUN_MODE == "quick",
  MI_SENS_NSIM_QUICK,
  MI_SENS_NSIM_FULL
)


# ------------------------------------------------------------------------------
# PDF sizes
# ------------------------------------------------------------------------------

BAR_WIDTH <- 10
BAR_HEIGHT <- 8.5

BEE_WIDTH <- 10
BEE_HEIGHT <- 8.5

DEP_WIDTH <- 18
DEP_HEIGHT <- 13


# ==============================================================================
# ZERO-CROSSING LABEL POSITION SETTINGS
#
# X:
#   positive  -> move right
#   negative  -> move left
#
# Y:
#   positive  -> move upward
#   negative  -> move downward
#
# Values are proportions of each panel's own x/y range.
# ==============================================================================

# ------------------------------------------------------------------------------
# Standard right-lower label placement:
# Age, NPR, GGT (displayed as γ-GT), Fbg, TyG, HB
# ------------------------------------------------------------------------------

RIGHT_DOWN_X <- 0.035
RIGHT_DOWN_Y <- -0.10


# ------------------------------------------------------------------------------
# Standard left-lower label placement:
# APTT
# ------------------------------------------------------------------------------

LEFT_DOWN_X <- -0.035
LEFT_DOWN_Y <- -0.10


# ------------------------------------------------------------------------------
# PAR: farther left and lower than the standard placement
#
# Adjust these two values if label placement requires refinement:
#   more negative PAR_X -> farther left
#   more negative PAR_Y -> farther down
# ------------------------------------------------------------------------------

PAR_X <- -0.10
PAR_Y <- -0.16


# ------------------------------------------------------------------------------
# TC: farther right and lower than the standard placement
#
# Adjust these two values if label placement requires refinement:
#   more positive TC_X -> farther right
#   more negative TC_Y -> farther down
# ------------------------------------------------------------------------------

TC_X <- 0.10
TC_Y <- -0.16


# ------------------------------------------------------------------------------
# Other continuous predictors:
# Place labels above the zero-crossing with automatic horizontal offset.
# ------------------------------------------------------------------------------

DEFAULT_X <- 0.025
DEFAULT_Y <- 0.10





# ==============================================================================
# 1. PACKAGE CHECK
# ==============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "fastshap",
  "shapviz",
  "openxlsx",
  "tidymodels",
  "themis",
  "kernlab",
  "xgboost",
  "glmnet"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Please install missing packages first:\n",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
  
}


suppressPackageStartupMessages({
  
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  library(fastshap)
  library(shapviz)
  
  library(openxlsx)
  
  library(tidymodels)
  library(themis)
  library(kernlab)
  library(xgboost)
  library(glmnet)
  
})


tidymodels_prefer()

set.seed(SEED)



# ==============================================================================
# 2. READ FINAL FROZEN MODEL BUNDLE
# ==============================================================================

cat(
  "\n============================================================\n",
  "READING FINAL MODEL BUNDLE\n",
  "============================================================\n"
)


bundle <- readRDS(
  BUNDLE_PATH
)



# ------------------------------------------------------------------------------
# Check essential objects
# ------------------------------------------------------------------------------

required_bundle_fields <- c(
  "predictors",
  "binary_predictors",
  "continuous_predictors",
  "mice_settings",
  "final_svm_fits",
  "final_xgb_fits",
  "meta_scaler",
  "meta_fit",
  "platt_calibrator"
)


missing_bundle_fields <- setdiff(
  required_bundle_fields,
  names(bundle)
)


if (length(missing_bundle_fields) > 0) {
  
  stop(
    paste0(
      "The bundle is missing required objects:\n",
      paste(
        missing_bundle_fields,
        collapse = ", "
      )
    )
  )
}



predictor_vars <- bundle$predictors

binary_vars <- bundle$binary_predictors

continuous_vars <- bundle$continuous_predictors

MICE_M <- bundle$mice_settings$m



cat(
  "Number of predictors:",
  length(predictor_vars),
  "\n"
)


cat(
  "Predictors:\n",
  paste(
    predictor_vars,
    collapse = ", "
  ),
  "\n"
)


cat(
  "MICE m =",
  MICE_M,
  "\n"
)



if (length(predictor_vars) != 16) {
  
  stop(
    "The frozen model does not contain exactly 16 predictors."
  )
}


if (
  length(bundle$final_svm_fits) != MICE_M ||
  length(bundle$final_xgb_fits) != MICE_M
) {
  
  stop(
    "Number of final SVM/XGBoost models does not match MICE m."
  )
}



# ==============================================================================
# 3. FIND THE FIVE FINAL MI SHEETS AUTOMATICALLY
# ==============================================================================

all_sheets <- readxl::excel_sheets(
  RESULTS_PATH
)


cat(
  "\nExcel sheets found:\n"
)

print(
  all_sheets
)


mi_sheets <- grep(
  "^Final_MI_[0-9]{2}$",
  all_sheets,
  value = TRUE
)


mi_sheets <- sort(
  mi_sheets
)


expected_mi_sheets <- sprintf(
  "Final_MI_%02d",
  seq_len(MICE_M)
)


if (!identical(
  mi_sheets,
  expected_mi_sheets
)) {
  
  stop(
    paste0(
      "Final MI sheets do not match expected names.\n",
      "Expected: ",
      paste(
        expected_mi_sheets,
        collapse = ", "
      ),
      "\nFound: ",
      paste(
        mi_sheets,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nFinal MI sheets detected successfully:\n"
)

print(
  mi_sheets
)



# ==============================================================================
# 4. READ ALL FIVE FINAL MICE DATASETS
# ==============================================================================

read_final_mi <- function(sheet_name) {
  
  dat <- readxl::read_excel(
    RESULTS_PATH,
    sheet = sheet_name
  ) |>
    as.data.frame()
  
  
  required_cols <- c(
    "Row",
    "PDstage",
    predictor_vars
  )
  
  
  missing_cols <- setdiff(
    required_cols,
    names(dat)
  )
  
  
  if (length(missing_cols) > 0) {
    
    stop(
      paste0(
        "Sheet ",
        sheet_name,
        " is missing columns: ",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    )
  }
  
  
  # Ensure all predictors are numeric
  dat[predictor_vars] <- lapply(
    dat[predictor_vars],
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  )
  
  
  if (
    anyNA(
      dat[predictor_vars]
    )
  ) {
    
    stop(
      paste0(
        "Missing values remain in ",
        sheet_name
      )
    )
  }
  
  
  return(
    dat
  )
  
}



mi_list <- lapply(
  mi_sheets,
  read_final_mi
)


names(mi_list) <- mi_sheets



# ------------------------------------------------------------------------------
# Validate patient ordering across all MI datasets
# ------------------------------------------------------------------------------

reference_row <- mi_list[[1]]$Row

reference_stage <- as.character(
  mi_list[[1]]$PDstage
)


for (j in seq_along(mi_list)) {
  
  if (!identical(
    mi_list[[j]]$Row,
    reference_row
  )) {
    
    stop(
      paste0(
        "Patient Row order differs in ",
        mi_sheets[j]
      )
    )
  }
  
  
  if (!identical(
    as.character(
      mi_list[[j]]$PDstage
    ),
    reference_stage
  )) {
    
    stop(
      paste0(
        "PDstage order differs in ",
        mi_sheets[j]
      )
    )
  }
  
}


N <- nrow(
  mi_list[[1]]
)


cat(
  "\nNumber of patients in each MI dataset =",
  N,
  "\n"
)



# ==============================================================================
# 5. READ ORIGINAL 1321-PATIENT DATASET
# ==============================================================================

raw_data <- readxl::read_excel(
  
  RAW_DATA_PATH,
  
  na = c(
    "",
    "NA",
    "N/A",
    "na",
    "NaN",
    "NULL"
  )
  
) |>
  as.data.frame()



missing_raw_cols <- setdiff(
  predictor_vars,
  names(raw_data)
)


if (length(missing_raw_cols) > 0) {
  
  stop(
    paste0(
      "Original dataset is missing predictors:\n",
      paste(
        missing_raw_cols,
        collapse = ", "
      )
    )
  )
}



if (nrow(raw_data) != N) {
  
  stop(
    paste0(
      "Original data N = ",
      nrow(raw_data),
      ", but MI data N = ",
      N
    )
  )
}



X_raw <- raw_data[
  ,
  predictor_vars,
  drop = FALSE
]


X_raw[] <- lapply(
  X_raw,
  function(x) {
    suppressWarnings(
      as.numeric(x)
    )
  }
)



# ==============================================================================
# 6. CHECK ORIGINAL MISSING VALUES
# ==============================================================================

missing_summary <- data.frame(
  
  Feature = predictor_vars,
  
  Missing_N = vapply(
    X_raw,
    function(x) {
      sum(
        is.na(x)
      )
    },
    integer(1)
  )
  
)


missing_summary$Missing_Rate <-
  missing_summary$Missing_N / N


cat(
  "\nOriginal predictor missingness:\n"
)


print(
  missing_summary |>
    filter(
      Missing_N > 0
    )
)



# ==============================================================================
# 7. CONFIRM OBSERVED VALUES AGREE WITH FINAL MI DATA
# ==============================================================================

for (f in predictor_vars) {
  
  observed_idx <- which(
    !is.na(
      X_raw[[f]]
    )
  )
  
  
  if (length(observed_idx) > 0) {
    
    diff_value <- abs(
      X_raw[[f]][observed_idx] -
        mi_list[[1]][[f]][observed_idx]
    )
    
    
    if (
      any(
        diff_value > 1e-8,
        na.rm = TRUE
      )
    ) {
      
      warning(
        paste0(
          "Observed values of ",
          f,
          " do not completely agree between raw data and Final_MI_01."
        )
      )
    }
    
  }
  
}



# ==============================================================================
# 8. CREATE THE COMPLETE RAW-SCALE MATRIX USED FOR MAIN SHAP
#
# Observed values:
#     retain the ORIGINAL value
#
# Missing values:
#     use mean of the five FINAL saved MICE values
#
# This affects only originally missing cells.
# ==============================================================================

X_shap <- X_raw


imputation_representation <- list()



for (f in predictor_vars) {
  
  mi_matrix_f <- do.call(
    
    cbind,
    
    lapply(
      mi_list,
      function(dat) {
        dat[[f]]
      }
    )
    
  )
  
  
  missing_idx <- which(
    is.na(
      X_raw[[f]]
    )
  )
  
  
  if (length(missing_idx) == 0) {
    
    next
    
  }
  
  
  if (f %in% binary_vars) {
    
    # Sex/Smoke are expected to be complete in the analysis dataset.
    # This is included only as a safety check.
    
    mi_mean <- rowMeans(
      mi_matrix_f[
        missing_idx,
        ,
        drop = FALSE
      ]
    )
    
    
    replacement <- ifelse(
      mi_mean >= 0.5,
      1,
      0
    )
    
    
  } else {
    
    replacement <- rowMeans(
      mi_matrix_f[
        missing_idx,
        ,
        drop = FALSE
      ]
    )
    
  }
  
  
  X_shap[[f]][missing_idx] <- replacement
  
  
  for (
    k in seq_along(
      missing_idx
    )
  ) {
    
    idx <- missing_idx[k]
    
    
    imputation_representation[[length(imputation_representation) + 1]] <- data.frame(
      
      Row = idx,
      
      Feature = f,
      
      MI_01 = mi_matrix_f[idx, 1],
      
      MI_02 = mi_matrix_f[idx, 2],
      
      MI_03 = mi_matrix_f[idx, 3],
      
      MI_04 = mi_matrix_f[idx, 4],
      
      MI_05 = mi_matrix_f[idx, 5],
      
      SHAP_Representation = replacement[k]
      
    )
    
  }
  
}



if (
  length(
    imputation_representation
  ) > 0
) {
  
  imputation_representation <- bind_rows(
    imputation_representation
  )
  
} else {
  
  imputation_representation <- data.frame(
    Status = "No missing predictors"
  )
  
}



if (
  anyNA(
    X_shap
  )
) {
  
  stop(
    "Missing values remain in X_shap."
  )
}



cat(
  "\nValues used to represent originally missing predictors:\n"
)

print(
  imputation_representation
)



# ==============================================================================
# 9. FINAL MODEL PREDICTION FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# 9.1 Meta learner + Platt calibration
# ------------------------------------------------------------------------------

finish_stacking_prediction <- function(
    object,
    svm_probability,
    xgb_probability
) {
  
  
  meta_input <- data.frame(
    
    SVM_z =
      (
        svm_probability -
          object$meta_scaler$SVM_mean
      ) /
      object$meta_scaler$SVM_sd,
    
    
    XGB_z =
      (
        xgb_probability -
          object$meta_scaler$XGB_mean
      ) /
      object$meta_scaler$XGB_sd
    
  )
  
  
  stack_raw <- predict(
    
    object$meta_fit,
    
    new_data = meta_input,
    
    type = "prob"
    
  )$.pred_Stage1
  
  
  
  # IMPORTANT:
  # The final calibration model fits:
  #
  #     labels ~ probs
  #
  # not labels ~ logit(probs)
  #
  calibrated_probability <- predict(
    
    object$platt_calibrator,
    
    newdata = data.frame(
      probs = stack_raw
    ),
    
    type = "response"
    
  )
  
  
  return(
    as.numeric(
      calibrated_probability
    )
  )
  
}



# ------------------------------------------------------------------------------
# 9.2 Pooled final model for SHAP
#
# The SAME complete feature vector x is evaluated through all:
#
#     5 SVM models
#     5 XGBoost models
#
# then probabilities are averaged BEFORE the meta learner.
# ------------------------------------------------------------------------------

predict_pooled_final <- function(
    object,
    newdata
) {
  
  
  x <- as.data.frame(
    newdata
  )
  
  
  m <- object$mice_settings$m
  
  
  svm_prob_mat <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = m
  )
  
  
  xgb_prob_mat <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = m
  )
  
  
  
  for (j in seq_len(m)) {
    
    
    svm_prob_mat[, j] <- predict(
      
      object$final_svm_fits[[j]],
      
      new_data = x,
      
      type = "prob"
      
    )$.pred_Stage1
    
    
    
    xgb_prob_mat[, j] <- predict(
      
      object$final_xgb_fits[[j]],
      
      new_data = x,
      
      type = "prob"
      
    )$.pred_Stage1
    
  }
  
  
  
  svm_probability <- rowMeans(
    svm_prob_mat
  )
  
  
  xgb_probability <- rowMeans(
    xgb_prob_mat
  )
  
  
  
  finish_stacking_prediction(
    
    object = object,
    
    svm_probability = svm_probability,
    
    xgb_probability = xgb_probability
    
  )
  
}



# ------------------------------------------------------------------------------
# 9.3 Exact development-cohort prediction using the five saved MI datasets
#
# Imp01 data -> model 1
# Imp02 data -> model 2
# ...
# Imp05 data -> model 5
#
# This replicates the final TRUE-MI prediction structure.
# ------------------------------------------------------------------------------

predict_exact_saved_mi <- function(
    object,
    mi_data_list
) {
  
  
  m <- object$mice_settings$m
  
  n <- nrow(
    mi_data_list[[1]]
  )
  
  
  svm_prob_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = m
  )
  
  
  xgb_prob_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = m
  )
  
  
  
  for (j in seq_len(m)) {
    
    
    x_j <- mi_data_list[[j]][
      ,
      object$predictors,
      drop = FALSE
    ]
    
    
    svm_prob_mat[, j] <- predict(
      
      object$final_svm_fits[[j]],
      
      new_data = x_j,
      
      type = "prob"
      
    )$.pred_Stage1
    
    
    
    xgb_prob_mat[, j] <- predict(
      
      object$final_xgb_fits[[j]],
      
      new_data = x_j,
      
      type = "prob"
      
    )$.pred_Stage1
    
  }
  
  
  
  svm_probability <- rowMeans(
    svm_prob_mat
  )
  
  
  xgb_probability <- rowMeans(
    xgb_prob_mat
  )
  
  
  
  finish_stacking_prediction(
    
    object = object,
    
    svm_probability = svm_probability,
    
    xgb_probability = xgb_probability
    
  )
  
}



# ==============================================================================
# 10. PREDICTION DIAGNOSTIC
#
# Compare:
#
# A. Exact five-MI prediction
# B. SHAP single complete representation
#
# For patients without original missing values,
# these should essentially be identical.
# ==============================================================================

cat(
  "\n============================================================\n",
  "PREDICTION DIAGNOSTIC\n",
  "============================================================\n"
)



pred_exact_mi <- predict_exact_saved_mi(
  
  object = bundle,
  
  mi_data_list = mi_list
  
)



pred_shap_representation <- predict_pooled_final(
  
  object = bundle,
  
  newdata = X_shap
  
)



has_original_missing <- apply(
  is.na(
    X_raw
  ),
  1,
  any
)



prediction_diagnostic <- data.frame(
  
  Row = reference_row,
  
  Original_Missing = has_original_missing,
  
  Exact_MI_Prediction = pred_exact_mi,
  
  SHAP_Representation_Prediction =
    pred_shap_representation,
  
  Absolute_Difference =
    abs(
      pred_exact_mi -
        pred_shap_representation
    )
  
)



cat(
  "\nMean absolute prediction difference =",
  mean(
    prediction_diagnostic$Absolute_Difference
  ),
  "\n"
)


cat(
  "Maximum absolute prediction difference =",
  max(
    prediction_diagnostic$Absolute_Difference
  ),
  "\n"
)



if (
  any(
    !has_original_missing
  )
) {
  
  cat(
    "Maximum difference among originally complete patients =",
    max(
      prediction_diagnostic$Absolute_Difference[
        !has_original_missing
      ]
    ),
    "\n"
  )
  
}



cat(
  "\nPatients with original missing predictors:\n"
)


print(
  
  prediction_diagnostic |>
    filter(
      Original_Missing
    )
  
)



# ==============================================================================
# 11. QUICK MODEL PREDICTION TEST
# ==============================================================================

test_n <- min(
  10,
  N
)


test_prediction <- predict_pooled_final(
  
  bundle,
  
  X_shap[
    seq_len(test_n),
    ,
    drop = FALSE
  ]
  
)


cat(
  "\nFirst ",
  test_n,
  " model predictions:\n",
  sep = ""
)


print(
  test_prediction
)



if (
  length(
    test_prediction
  ) != test_n ||
  any(
    !is.finite(
      test_prediction
    )
  )
) {
  
  stop(
    "Prediction wrapper test failed."
  )
}



# ==============================================================================
# 12. FASTSHAP
# ==============================================================================

cat(
  "\n============================================================\n",
  "STARTING FASTSHAP\n",
  "RUN MODE = ",
  RUN_MODE,
  "\nNSIM = ",
  NSIM,
  "\n============================================================\n",
  sep = ""
)



set.seed(SEED)


shap_time <- system.time(
  
  
  shap_result <- fastshap::explain(
    
    object = bundle,
    
    X = X_shap,
    
    pred_wrapper = predict_pooled_final,
    
    nsim = NSIM,
    
    adjust = TRUE,
    
    shap_only = TRUE
    
  )
  
  
)



cat(
  "\nSHAP computation completed.\n"
)


print(
  shap_time
)



shap_baseline <- attr(
  shap_result,
  "baseline"
)


shap_matrix <- as.matrix(
  shap_result
)


colnames(
  shap_matrix
) <- predictor_vars



# The SHAP RDS contains patient-level feature and SHAP values.
# Keep it under private/ and do not commit it to the public repository.
saveRDS(

  list(
    
    SHAP = shap_matrix,
    
    X = X_shap,
    
    baseline = shap_baseline,
    
    NSIM = NSIM,
    
    run_mode = RUN_MODE,
    
    prediction_diagnostic =
      prediction_diagnostic
    
  ),
  
  file = file.path(
    OUTPUT_DIR,
    "shap_analysis_results.rds"
  )
  
)



# ==============================================================================
# 13. SHAP ADDITIVITY CHECK
# ==============================================================================

if (
  !is.null(
    shap_baseline
  ) &&
  length(
    shap_baseline
  ) == 1
) {
  
  
  shap_reconstructed_prediction <-
    shap_baseline +
    rowSums(
      shap_matrix
    )
  
  
  shap_additivity_difference <- abs(
    
    shap_reconstructed_prediction -
      pred_shap_representation
    
  )
  
  
  cat(
    "\nSHAP additivity check:\n"
  )
  
  
  cat(
    "Mean absolute difference =",
    mean(
      shap_additivity_difference
    ),
    "\n"
  )
  
  
  cat(
    "Maximum absolute difference =",
    max(
      shap_additivity_difference
    ),
    "\n"
  )
  
}



# ==============================================================================
# 14. DISPLAY LABELS
#
# Model variable:
#     AG_ratio
#     GGT
#
# Figure label:
#     A/G
#     γ-GT
#
# gamma is rendered using plotmath,
# avoiding Unicode PDF/font problems on macOS.
# ==============================================================================

feature_label_text <- function(x) {
  
  
  if (
    x == "AG_ratio"
  ) {
    
    return(
      "'A/G'"
    )
    
  }
  
  
  if (
    x == "GGT"
  ) {
    
    return(
      "gamma*'-GT'"
    )
    
  }
  
  
  return(
    paste0(
      "'",
      x,
      "'"
    )
  )
  
}



feature_axis_labels <- function(x) {
  
  parse(
    text = vapply(
      x,
      feature_label_text,
      character(1)
    )
  )
  
}



facet_label_map <- setNames(
  
  vapply(
    predictor_vars,
    feature_label_text,
    character(1)
  ),
  
  predictor_vars
  
)



feature_facet_labeller <- as_labeller(
  
  facet_label_map,
  
  default = label_parsed
  
)



# ==============================================================================
# 15. COMMON PUBLICATION THEME
# ==============================================================================

theme_shap <- theme_classic(
  
  base_family = "Helvetica",
  
  base_size = 12
  
) +
  
  theme(
    
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    )
    
  )



# ==============================================================================
# 16. PDF SAVE FUNCTION
#
# Uses base pdf()
# No cairo_pdf()
# No XQuartz dependency
# ==============================================================================

save_vector_pdf <- function(
    filename,
    plot_object,
    width,
    height
) {
  
  
  pdf(
    
    file = file.path(
      OUTPUT_DIR,
      filename
    ),
    
    width = width,
    
    height = height,
    
    family = "Helvetica",
    
    useDingbats = FALSE
    
  )
  
  
  print(
    plot_object
  )
  
  
  dev.off()
  
}



# ==============================================================================
# 17. MEAN ABSOLUTE SHAP
# ==============================================================================

mean_abs_shap <- colMeans(
  abs(
    shap_matrix
  )
)



importance_df <- data.frame(
  
  Feature = names(
    mean_abs_shap
  ),
  
  Mean_Absolute_SHAP =
    as.numeric(
      mean_abs_shap
    )
  
) |>
  
  arrange(
    desc(
      Mean_Absolute_SHAP
    )
  ) |>
  
  mutate(
    
    Rank = row_number()
    
  )



importance_order <- importance_df$Feature



importance_plot_df <- importance_df |>
  
  mutate(
    
    Feature = factor(
      
      Feature,
      
      levels = rev(
        importance_order
      )
      
    )
    
  )



# ==============================================================================
# 18. FIGURE 6-1
# Mean |SHAP| bar plot
# ==============================================================================

p_bar <- ggplot(
  
  importance_plot_df,
  
  aes(
    x = Mean_Absolute_SHAP,
    y = Feature
  )
  
) +
  
  geom_col(
    
    fill = "#4C78A8",
    
    width = 0.68
    
  ) +
  
  scale_y_discrete(
    
    labels = feature_axis_labels
    
  ) +
  
  theme_shap +
  
  theme(
    
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.line.x = element_line(color = "black", linewidth = 0.5)
    
  ) +
  
  labs(
    
    x = "Mean |SHAP value|",
    
    y = NULL
    
  )



print(
  p_bar
)



save_vector_pdf(
  
  filename = "shap_mean_absolute.pdf",
  
  plot_object = p_bar,
  
  width = BAR_WIDTH,
  
  height = BAR_HEIGHT
  
)



# ==============================================================================
# 19. CREATE SHAPVIZ OBJECT
# ==============================================================================

sv_final <- shapviz::shapviz(
  
  shap_matrix,
  
  X = X_shap
  
)



# ==============================================================================
# 20. FIGURE 6-2
# SHAP beeswarm summary plot
# ==============================================================================

p_bee <- shapviz::sv_importance(
  
  sv_final,
  
  kind = "beeswarm",
  
  max_display = 16,
  
  show_numbers = FALSE,
  
  size = 1.25,
  
  alpha = 0.80
  
) +
  
  scale_y_discrete(
    
    labels = feature_axis_labels
    
  ) +
  
  scale_color_gradientn(
    
    colors = c(
      "#4575B4",
      "#F7F7F7",
      "#D73027"
    ),
    
    values = c(
      0,
      0.5,
      1
    ),
    
    breaks = c(
      0,
      1
    ),
    
    labels = c(
      "Low",
      "High"
    ),
    
    name = "Feature value"
    
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    color = "grey40",
    
    linewidth = 0.45
    
  ) +
  
  theme_shap +
  
  theme(
    
    panel.grid.major.y = element_line(
      
      color = "grey93",
      
      linetype = "dotted",
      
      linewidth = 0.4
      
    ),
    
    legend.position = "right"
    
  ) +
  
  labs(
    
    x = "SHAP value (impact on model output)",
    
    y = NULL
    
  )



print(
  p_bee
)



save_vector_pdf(
  
  filename = "shap_beeswarm.pdf",
  
  plot_object = p_bee,
  
  width = BEE_WIDTH,
  
  height = BEE_HEIGHT
  
)



# ==============================================================================
# 21. PREPARE LONG DATA FOR DEPENDENCE PLOTS
# ==============================================================================

S_df <- as.data.frame(
  shap_matrix
)


X_df <- as.data.frame(
  X_shap
)



S_long <- S_df |>
  
  mutate(
    RowID = row_number()
  ) |>
  
  pivot_longer(
    
    cols = -RowID,
    
    names_to = "Feature",
    
    values_to = "SHAP_Value"
    
  )



X_long <- X_df |>
  
  mutate(
    RowID = row_number()
  ) |>
  
  pivot_longer(
    
    cols = -RowID,
    
    names_to = "Feature",
    
    values_to = "Original_Value"
    
  )



color_information <- X_df |>
  
  mutate(
    RowID = row_number()
  ) |>
  
  transmute(
    
    RowID,
    
    Age_Color = Age,
    
    Sex_Color = factor(
      
      Sex,
      
      levels = c(
        0,
        1
      ),
      
      labels = c(
        "Female",
        "Male"
      )
      
    )
    
  )



plot_df <- S_long |>
  
  left_join(
    
    X_long,
    
    by = c(
      "RowID",
      "Feature"
    )
    
  ) |>
  
  left_join(
    
    color_information,
    
    by = "RowID"
    
  ) |>
  
  mutate(
    
    Feature = factor(
      
      Feature,
      
      levels = importance_order
      
    )
    
  )



# ==============================================================================
# 22. ZERO-CROSSING FUNCTION
# ==============================================================================

find_zero_crossings <- function(
    x,
    y,
    span = LOESS_SPAN,
    quantiles = ZERO_QUANTILES
) {
  
  
  valid <- is.finite(x) &
    is.finite(y)
  
  
  x <- x[valid]
  
  y <- y[valid]
  
  
  
  if (
    length(x) < 20 ||
    length(
      unique(x)
    ) < 5
  ) {
    
    return(
      numeric(0)
    )
  }
  
  
  
  fit <- try(
    
    loess(
      
      y ~ x,
      
      span = span,
      
      degree = 2,
      
      control = loess.control(
        surface = "direct"
      )
      
    ),
    
    silent = TRUE
    
  )
  
  
  
  if (
    inherits(
      fit,
      "try-error"
    )
  ) {
    
    return(
      numeric(0)
    )
    
  }
  
  
  
  x_low <- as.numeric(
    
    quantile(
      
      x,
      
      probs = quantiles[1],
      
      na.rm = TRUE
      
    )
    
  )
  
  
  
  x_high <- as.numeric(
    
    quantile(
      
      x,
      
      probs = quantiles[2],
      
      na.rm = TRUE
      
    )
    
  )
  
  
  
  x_grid <- seq(
    
    x_low,
    
    x_high,
    
    length.out = 5000
    
  )
  
  
  
  y_grid <- predict(
    
    fit,
    
    newdata = data.frame(
      x = x_grid
    )
    
  )
  
  
  
  roots <- numeric(0)
  
  
  
  for (
    i in seq_len(
      length(
        x_grid
      ) - 1
    )
  ) {
    
    
    y1 <- y_grid[i]
    
    y2 <- y_grid[i + 1]
    
    
    
    if (
      !is.finite(y1) ||
      !is.finite(y2)
    ) {
      
      next
      
    }
    
    
    
    if (
      y1 == 0
    ) {
      
      roots <- c(
        roots,
        x_grid[i]
      )
      
      
    } else if (
      y1 * y2 < 0
    ) {
      
      
      # Linear interpolation between the two grid points
      
      root <- x_grid[i] -
        
        y1 *
        (
          x_grid[i + 1] -
            x_grid[i]
        ) /
        (
          y2 -
            y1
        )
      
      
      roots <- c(
        roots,
        root
      )
      
    }
    
  }
  
  
  
  roots <- sort(
    
    unique(
      
      round(
        roots,
        4
      )
      
    )
    
  )
  
  
  
  return(
    roots
  )
  
}



# ==============================================================================
# 23. CALCULATE ALL ZERO-CROSSINGS
# ==============================================================================

zero_list <- lapply(
  
  continuous_vars,
  
  function(f) {
    
    
    sub_data <- plot_df |>
      
      filter(
        as.character(
          Feature
        ) == f
      )
    
    
    
    roots <- find_zero_crossings(
      
      x = sub_data$Original_Value,
      
      y = sub_data$SHAP_Value
      
    )
    
    
    
    if (
      length(
        roots
      ) == 0
    ) {
      
      return(
        NULL
      )
      
    }
    
    
    
    data.frame(
      
      Feature = f,
      
      Zero_X = roots,
      
      Feature_Median = median(
        
        sub_data$Original_Value,
        
        na.rm = TRUE
        
      )
      
    )
    
  }
  
)



zero_all <- bind_rows(
  zero_list
)



# ==============================================================================
# 24. SELECT ONE ZERO-CROSSING FOR FIGURE ANNOTATION
#
# ALL crossings remain in zero_all and are exported.
# ==============================================================================

if (
  nrow(
    zero_all
  ) > 0
) {
  
  
  if (
    ZERO_SELECTION == "closest_to_median"
  ) {
    
    
    zero_primary <- zero_all |>
      
      group_by(
        Feature
      ) |>
      
      slice_min(
        
        order_by = abs(
          Zero_X -
            Feature_Median
        ),
        
        n = 1,
        
        with_ties = FALSE
        
      ) |>
      
      ungroup()
    
    
  } else {
    
    
    stop(
      "Unknown ZERO_SELECTION setting."
    )
    
  }
  
  
  
  zero_primary <- zero_primary |>
    
    mutate(
      
      Feature = factor(
        
        Feature,
        
        levels = importance_order
        
      ),
      
      Label = paste0(
        
        "Zero-crossing = ",
        
        format(
          
          round(
            Zero_X,
            2
          ),
          
          nsmall = 2
          
        )
        
      ),
      
      hjust_value = ifelse(
        
        Zero_X >= Feature_Median,
        
        1.05,
        
        -0.05
        
      )
      
    )
  
  
  
} else {
  
  
  zero_primary <- data.frame(
    
    Feature = factor(
      character(0),
      levels = importance_order
    ),
    
    Zero_X = numeric(0),
    
    Feature_Median = numeric(0),
    
    Label = character(0),
    
    hjust_value = numeric(0)
    
  )
  
}



# ==============================================================================
# 24B. CUSTOM ZERO-CROSSING LABEL COORDINATES
#
# These offsets change label placement only and do not change:
#   - the estimated zero-crossing
#   - the vertical reference line
#   - the zero-crossing marker
#
# With facet_wrap(scales = "free"), offsets are calculated from each panel's
# own x/y range rather than fixed hjust/vjust values.
# ==============================================================================

panel_ranges <- plot_df |>

  mutate(
    Feature_chr = as.character(Feature)
  ) |>

  group_by(
    Feature_chr
  ) |>

  summarise(

    X_Min = min(
      Original_Value,
      na.rm = TRUE
    ),

    X_Max = max(
      Original_Value,
      na.rm = TRUE
    ),

    Y_Min = min(
      SHAP_Value,
      na.rm = TRUE
    ),

    Y_Max = max(
      SHAP_Value,
      na.rm = TRUE
    ),

    .groups = "drop"

  ) |>

  mutate(

    X_Range = X_Max - X_Min,

    Y_Range = Y_Max - Y_Min,

    X_Range = ifelse(
      !is.finite(X_Range) | X_Range == 0,
      1,
      X_Range
    ),

    Y_Range = ifelse(
      !is.finite(Y_Range) | Y_Range == 0,
      1,
      Y_Range
    )

  )


# The model variable remains GGT; figure labels display it as γ-GT.
right_down_features <- c(
  "Age",
  "NPR",
  "GGT",
  "Fbg",
  "TyG",
  "HB"
)

left_down_features <- c(
  "APTT"
)


zero_primary_plot <- zero_primary |>

  mutate(
    Feature_chr = as.character(Feature)
  ) |>

  left_join(
    panel_ranges,
    by = "Feature_chr"
  ) |>

  mutate(

    # X coordinate
    Label_X = case_when(

      Feature_chr %in% right_down_features ~
        Zero_X + RIGHT_DOWN_X * X_Range,

      Feature_chr %in% left_down_features ~
        Zero_X + LEFT_DOWN_X * X_Range,

      Feature_chr == "PAR" ~
        Zero_X + PAR_X * X_Range,

      Feature_chr == "TC" ~
        Zero_X + TC_X * X_Range,

      Zero_X >= Feature_Median ~
        Zero_X - DEFAULT_X * X_Range,

      TRUE ~
        Zero_X + DEFAULT_X * X_Range

    ),


    # Y coordinate
    Label_Y = case_when(

      Feature_chr %in% right_down_features ~
        RIGHT_DOWN_Y * Y_Range,

      Feature_chr %in% left_down_features ~
        LEFT_DOWN_Y * Y_Range,

      Feature_chr == "PAR" ~
        PAR_Y * Y_Range,

      Feature_chr == "TC" ~
        TC_Y * Y_Range,

      TRUE ~
        DEFAULT_Y * Y_Range

    ),


    # Text alignment
    Label_Hjust = case_when(

      Feature_chr %in% right_down_features ~ 0,

      Feature_chr %in% left_down_features ~ 1,

      Feature_chr == "PAR" ~ 1,

      Feature_chr == "TC" ~ 0,

      Zero_X >= Feature_Median ~ 1,

      TRUE ~ 0

    )

  )


# ==============================================================================
# 25. DEPENDENCE PLOT FUNCTION
# ==============================================================================

make_dependence_plot <- function(
    color_by = c(
      "Sex",
      "Age"
    )
) {
  
  
  color_by <- match.arg(
    color_by
  )
  
  
  df_plot <- plot_df
  
  
  
  if (
    color_by == "Sex"
  ) {
    
    df_plot$ColorVariable <-
      df_plot$Sex_Color
    
  } else {
    
    df_plot$ColorVariable <-
      df_plot$Age_Color
    
  }
  
  
  
  p <- ggplot(
    
    df_plot,
    
    aes(
      x = Original_Value,
      y = SHAP_Value
    )
    
  ) +
    
    
    # --------------------------------------------------------------------------
  # Continuous feature points
  # --------------------------------------------------------------------------
  
  geom_point(
    
    data = df_plot |>
      
      filter(
        !as.character(
          Feature
        ) %in%
          binary_vars
      ),
    
    aes(
      color = ColorVariable
    ),
    
    size = 0.85,
    
    alpha = 0.50
    
  ) +
    
    
    # --------------------------------------------------------------------------
  # Binary variables: Sex / Smoke
  # --------------------------------------------------------------------------
  
  geom_jitter(
    
    data = df_plot |>
      
      filter(
        as.character(
          Feature
        ) %in%
          binary_vars
      ),
    
    aes(
      color = ColorVariable
    ),
    
    width = 0.10,
    
    height = 0,
    
    size = 0.85,
    
    alpha = 0.50
    
  ) +
    
    
    # --------------------------------------------------------------------------
  # Overall LOESS
  #
  # Coloring is descriptive only.
  # LOESS describes the overall complete development cohort.
  # --------------------------------------------------------------------------
  
  geom_smooth(
    
    data = df_plot |>
      
      filter(
        !as.character(
          Feature
        ) %in%
          binary_vars
      ),
    
    method = "loess",
    
    formula = y ~ x,
    
    span = LOESS_SPAN,
    
    se = FALSE,
    
    color = "black",
    
    linewidth = 0.85
    
  ) +
    
    
    # SHAP = 0
    
    geom_hline(
      
      yintercept = 0,
      
      color = "grey40",
      
      linewidth = 0.40
      
    ) +
    
    
    # Zero crossing vertical line
    
    geom_vline(
      
      data = zero_primary_plot,
      
      aes(
        xintercept = Zero_X
      ),
      
      inherit.aes = FALSE,
      
      linetype = "dotted",
      
      color = "black",
      
      linewidth = 0.40
      
    ) +
    
    
    # Zero crossing point
    
    geom_point(
      
      data = zero_primary_plot,
      
      aes(
        x = Zero_X,
        y = 0
      ),
      
      inherit.aes = FALSE,
      
      shape = 21,
      
      fill = "white",
      
      color = "black",
      
      stroke = 0.9,
      
      size = 2.4
      
    ) +
    
    
    # Zero crossing label
    
    geom_text(
      
      data = zero_primary_plot,
      
      aes(
        
        x = Label_X,
        
        y = Label_Y,
        
        label = Label,
        
        hjust = Label_Hjust
        
      ),
      
      inherit.aes = FALSE,
      
      vjust = 0.5,
      
      size = 2.6,
      
      family = "Helvetica",
      
      color = "black",
      
      check_overlap = FALSE
      
    ) +
    
    
    facet_wrap(
      
      ~ Feature,
      
      scales = "free",
      
      ncol = 4,
      
      labeller = feature_facet_labeller
      
    ) +
    
    
    theme_bw(
      
      base_family = "Helvetica",
      
      base_size = 10
      
    ) +
    
    
    theme(
      
      panel.grid = element_blank(),
      
      panel.border = element_rect(
        
        color = "black",
        
        fill = NA,
        
        linewidth = 0.55
        
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        
        size = 10.5,
        
        color = "black"
        
      ),
      
      axis.text = element_text(
        
        size = 8.5,
        
        color = "black"
        
      ),
      
      axis.title = element_text(
        
        size = 13,
        
        face = "bold",
        
        color = "black"
        
      ),
      
      panel.spacing = grid::unit(
        
        1.25,
        
        "lines"
        
      ),
      
      plot.margin = margin(
        
        10,
        25,
        10,
        10
        
      )
      
    ) +
    
    
    coord_cartesian(
      clip = "off"
    ) +
    
    
    labs(
      
      x = "Actual feature value",
      
      y = "SHAP value"
      
    )
  
  
  
  # ---------------------------------------------------------------------------
  # Sex coloring
  # ---------------------------------------------------------------------------
  
  if (
    color_by == "Sex"
  ) {
    
    
    p <- p +
      
      scale_color_viridis_d(
        
        option = "plasma",
        
        name = "Sex",
        
        end = 0.82,
        
        drop = FALSE
        
      ) +
      
      theme(
        
        legend.position = "bottom"
        
      )
    
  }
  
  
  
  # ---------------------------------------------------------------------------
  # Age coloring
  # ---------------------------------------------------------------------------
  
  if (
    color_by == "Age"
  ) {
    
    
    p <- p +
      
      scale_color_viridis_c(
        
        option = "plasma",
        
        name = "Age"
        
      ) +
      
      theme(
        
        legend.position = "right"
        
      )
    
  }
  
  
  
  return(
    p
  )
  
}



# ==============================================================================
# 26. FIGURE 6-3
# Dependence plots colored by Sex
# ==============================================================================

p_dep_sex <- make_dependence_plot(
  color_by = "Sex"
)


print(
  p_dep_sex
)



save_vector_pdf(
  
  filename = "shap_dependence_sex.pdf",
  
  plot_object = p_dep_sex,
  
  width = DEP_WIDTH,
  
  height = DEP_HEIGHT
  
)



# ==============================================================================
# 27. FIGURE 6-3
# Dependence plots colored by Age
# ==============================================================================

p_dep_age <- make_dependence_plot(
  color_by = "Age"
)


print(
  p_dep_age
)



save_vector_pdf(
  
  filename = "shap_dependence_age.pdf",
  
  plot_object = p_dep_age,
  
  width = DEP_WIDTH,
  
  height = DEP_HEIGHT
  
)



# ==============================================================================
# 28. EXPORT SHAP VALUES
# ==============================================================================

shap_export <- as.data.frame(
  shap_matrix
)


shap_export <- cbind(
  
  data.frame(
    Row = reference_row
  ),
  
  shap_export
  
)



# ==============================================================================
# 29. ZERO-CROSSING EXPORT
# ==============================================================================

if (
  nrow(
    zero_all
  ) > 0
) {
  
  
  zero_all_export <- zero_all |>
    
    arrange(
      Feature,
      Zero_X
    )
  
  
  zero_primary_export <- zero_primary |>
    
    transmute(
      
      Feature = as.character(
        Feature
      ),
      
      Primary_Zero_Crossing = Zero_X,
      
      Selection_Rule =
        ZERO_SELECTION
      
    )
  
  
} else {
  
  
  zero_all_export <- data.frame(
    Status = "No zero-crossings detected"
  )
  
  
  zero_primary_export <- data.frame(
    Status = "No zero-crossings detected"
  )
  
}



# ==============================================================================
# 30. SETTINGS EXPORT
# ==============================================================================

settings_export <- data.frame(
  
  Parameter = c(
    
    "RUN_MODE",
    "NSIM",
    "SEED",
    "MICE_M",
    "LOESS_SPAN",
    "ZERO_LOWER_QUANTILE",
    "ZERO_UPPER_QUANTILE",
    "ZERO_SELECTION",
    "MISSING_REPRESENTATION"
    
  ),
  
  Value = c(
    
    RUN_MODE,
    NSIM,
    SEED,
    MICE_M,
    LOESS_SPAN,
    ZERO_QUANTILES[1],
    ZERO_QUANTILES[2],
    ZERO_SELECTION,
    MISSING_REPRESENTATION
    
  )
  
)



# ==============================================================================
# 31. EXPORT MAIN RESULTS TO EXCEL
# ==============================================================================

openxlsx::write.xlsx(
  
  list(
    
    Settings =
      settings_export,
    
    Missing_Original =
      missing_summary,
    
    Missing_MI_Representation =
      imputation_representation,
    
    Prediction_Diagnostic =
      prediction_diagnostic,
    
    Mean_Absolute_SHAP =
      importance_df,
    
    SHAP_Values =
      shap_export,
    
    Zero_Crossings_All =
      zero_all_export,
    
    Zero_Primary =
      zero_primary_export
    
  ),
  
  file = file.path(
    
    OUTPUT_DIR,
    
    "shap_analysis_results.xlsx"
    
  ),
  
  rowNames = FALSE,
  
  overwrite = TRUE
  
)



# ==============================================================================
# 32. OPTIONAL FIVE-MI SHAP SENSITIVITY ANALYSIS
#
# OFF by default.
#
# This checks whether feature importance is stable if SHAP is evaluated
# separately on each of the five completed datasets.
# ==============================================================================

if (
  RUN_MI_SENSITIVITY
) {
  
  
  cat(
    "\n============================================================\n",
    "STARTING OPTIONAL 5-MI SHAP SENSITIVITY ANALYSIS\n",
    "NSIM = ",
    MI_SENS_NSIM,
    "\n============================================================\n",
    sep = ""
  )
  
  
  mi_importance_list <- list()
  
  
  
  for (
    j in seq_len(
      MICE_M
    )
  ) {
    
    
    cat(
      "\nMI sensitivity ",
      j,
      "/",
      MICE_M,
      "\n",
      sep = ""
    )
    
    
    X_j <- mi_list[[j]][
      
      ,
      predictor_vars,
      drop = FALSE
      
    ]
    
    
    set.seed(
      SEED + j
    )
    
    
    shap_j <- fastshap::explain(
      
      object = bundle,
      
      X = X_j,
      
      pred_wrapper =
        predict_pooled_final,
      
      nsim = MI_SENS_NSIM,
      
      adjust = TRUE,
      
      shap_only = TRUE
      
    )
    
    
    shap_j <- as.matrix(
      shap_j
    )
    
    
    imp_j <- colMeans(
      abs(
        shap_j
      )
    )
    
    
    mi_importance_list[[j]] <- data.frame(
      
      Imputation =
        paste0(
          "MI_",
          sprintf(
            "%02d",
            j
          )
        ),
      
      Feature =
        names(
          imp_j
        ),
      
      Mean_Absolute_SHAP =
        as.numeric(
          imp_j
        )
      
    )
    
  }
  
  
  
  mi_importance <- bind_rows(
    mi_importance_list
  )
  
  
  
  mi_importance_summary <- mi_importance |>
    
    group_by(
      Feature
    ) |>
    
    summarise(
      
      Mean =
        mean(
          Mean_Absolute_SHAP
        ),
      
      SD =
        sd(
          Mean_Absolute_SHAP
        ),
      
      Min =
        min(
          Mean_Absolute_SHAP
        ),
      
      Max =
        max(
          Mean_Absolute_SHAP
        ),
      
      .groups = "drop"
      
    ) |>
    
    arrange(
      desc(
        Mean
      )
    )
  
  
  
  openxlsx::write.xlsx(
    
    list(
      
      MI_Importance =
        mi_importance,
      
      MI_Importance_Summary =
        mi_importance_summary
      
    ),
    
    file = file.path(
      
      OUTPUT_DIR,
      
      "shap_mi_sensitivity.xlsx"
      
    ),
    
    rowNames = FALSE,
    
    overwrite = TRUE
    
  )
  
}



# ==============================================================================
# 33. SAVE SESSION INFO
# ==============================================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    
    OUTPUT_DIR,
    
    "shap_session_info.txt"
    
  )
  
)



# ==============================================================================
# FINISHED
# ==============================================================================

cat(
  "\n============================================================\n",
  "SHAP ANALYSIS COMPLETED\n",
  "============================================================\n",
  "\nOutput folder:\n",
  OUTPUT_DIR,
  "\n\nMain outputs:\n",
  "shap_mean_absolute.pdf\n",
  "shap_beeswarm.pdf\n",
  "shap_dependence_sex.pdf\n",
  "shap_dependence_age.pdf\n",
  "shap_analysis_results.xlsx\n",
  "SHAP_Main_Result.rds\n",
  "SHAP_R_SessionInfo.txt\n",
  "\n"
)



