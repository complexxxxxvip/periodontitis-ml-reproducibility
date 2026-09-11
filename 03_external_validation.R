# ==============================================================================
# External Validation: Frozen Stacking Model + Four Logistic Comparison Models
# Includes 95% CIs for sensitivity, specificity, Brier score,
# calibration intercept, and calibration slope
#
# External dataset:
#   private/data/external_validation_data.xlsx
#
# Models:
#   M1: Age alone
#   M2: Age + Sex + Smoke
#   M3: 13 biomarkers only
#   M4: Age + Sex + Smoke + 13 biomarkers
#   M5: Frozen SVM + XGBoost stacking model
#
# IMPORTANT PRINCIPLES
#   1) External outcome is NEVER used for imputation, model fitting, calibration,
#      threshold selection, or any other model-development step.
#   2) The frozen stacking model is NOT retrained.
#   3) M1-M4 are synchronized to the final baseline-model analysis.
#      The script first verifies that the baseline RDS was generated from the
#      exact same frozen stacking bundle (OOF Stack probabilities, fold map and
#      M5 threshold must agree). If the updated baseline RDS contains saved final
#      M1-M4 fits, those frozen fits are used directly. Otherwise the same final
#      M1-M4 fits are deterministically reconstructed from the frozen development
#      imputations; external data are used only for prediction/evaluation.
#   4) MICE operates ONLY on the same 18 BASIC/source variables used in development:
#        Sex, Smoke, Age, TP, ALB, CR, Fbg, PLT, TC, GGT, ALP,
#        NE_C, PT, HB, APTT, LDH, GLU, TG
#   5) The six derived predictors are NEVER imputed and are recalculated AFTER
#      each completed imputation:
#        AG_ratio = ALB/(TP-ALB)
#        eGFR     = CKD-EPI equation from CR, Age, Sex
#        PAR      = PLT/ALB
#        NPR      = NE_C/PLT
#        INR      = (PT/13.1)^1.31
#        TyG      = log(1594.26*TG*GLU/2)
#   6) Sex coding and eGFR race coefficient are read from the frozen bundle.
#      Expected: Sex=0 female, Sex=1 male, race coefficient=1.0.
#   7) All classification thresholds are frozen from the development cohort.
#      NO external Youden threshold is calculated for the primary analysis.
#   8) Platt calibration for the stacking model is frozen and is NOT refitted.
#
# Main outputs:
#   external_validation_results.xlsx
#   external_roc_comparison.png / .pdf
#   external_calibration_comparison.png / .pdf
#   external_validation_results.rds
#   external_validation_session_info.txt
# ==============================================================================

rm(list = ls())
gc()

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------

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
EXTERNAL_DATA_PATH <- file.path(DATA_DIR, "external_validation_data.xlsx")
STACK_BUNDLE_PATH <- file.path(MODEL_DIR, "final_stacking_model_bundle.rds")
BASELINE_RESULTS_PATH <- file.path(MODEL_DIR, "baseline_model_results.rds")

OUTPUT_DIR <- file.path(OUTPUT_ROOT, "external_validation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ECE_BINS <- 10

# Performance uncertainty settings.
PERFORMANCE_CI_LEVEL <- 0.95
PERFORMANCE_BOOT_N <- 2000L
PERFORMANCE_BOOT_SEED_OFFSET <- 90000L

# ------------------------------------------------------------------------------
# 1. PACKAGE CHECK AND LOADING
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "tidymodels", "themis", "mice", "readxl", "kernlab", "xgboost",
  "glmnet", "pROC", "openxlsx", "dplyr", "purrr", "tibble", "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them first, e.g. install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))"
    )
  )
}

suppressPackageStartupMessages({
  library(tidymodels)
  library(themis)
  library(mice)
  library(readxl)
  library(kernlab)
  library(xgboost)
  library(glmnet)
  library(pROC)
  library(openxlsx)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

tidymodels_prefer()

# ------------------------------------------------------------------------------
# 2. FILE CHECKS
# ------------------------------------------------------------------------------

required_files <- c(
  external_data = EXTERNAL_DATA_PATH,
  stacking_model_bundle = STACK_BUNDLE_PATH,
  baseline_model_results = BASELINE_RESULTS_PATH
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required private file(s) not found:\n",
      paste0("  - ", names(missing_files), ": ", unname(missing_files), collapse = "\n"),
      "\nPlace these files under private/ according to the repository README."
    )
  )
}

cat("\n==================== FILES ====================\n")
cat("External data:   ", EXTERNAL_DATA_PATH, "\n")
cat("Stack bundle:    ", STACK_BUNDLE_PATH, "\n")
cat("Baseline results:", BASELINE_RESULTS_PATH, "\n")

# ------------------------------------------------------------------------------
# 3. LOAD FROZEN DEVELOPMENT OBJECTS
# ------------------------------------------------------------------------------

bundle <- readRDS(STACK_BUNDLE_PATH)
baseline_results <- readRDS(BASELINE_RESULTS_PATH)

required_bundle_items <- c(
  "seed",
  "predictors",
  "binary_predictors",
  "continuous_predictors",
  "derived_predictors",
  "imputation_variables",
  "egfr_settings",
  "mice_settings",
  "raw_development_imputation_data",
  "svm_oof_predictions_by_imputation",
  "final_completed_basic_datasets",
  "final_svm_fits",
  "final_xgb_fits",
  "meta_scaler",
  "meta_fit",
  "platt_calibrator",
  "final_threshold",
  "meta_oof_predictions"
)

missing_bundle_items <- setdiff(required_bundle_items, names(bundle))
if (length(missing_bundle_items) > 0) {
  stop(
    paste0(
      "Frozen stack bundle is missing required elements: ",
      paste(missing_bundle_items, collapse = ", "),
      "\nPlease use the final derived-after-MICE frozen bundle."
    )
  )
}

if (is.null(baseline_results$performance_summary)) {
  stop(
    "baseline-model_Simple_Baseline_Results.rds does not contain performance_summary."
  )
}

SEED <- as.integer(bundle$seed)
set.seed(SEED)

MICE_M <- as.integer(bundle$mice_settings$m)
MICE_MAXIT <- as.integer(bundle$mice_settings$maxit)
MICE_EPS <- as.numeric(bundle$mice_settings$eps)
MICE_METHOD <- as.character(bundle$mice_settings$method)
MICE_PRINT_FLAG <- FALSE

predictor_vars <- as.character(bundle$predictors)
binary_vars <- as.character(bundle$binary_predictors)
continuous_vars <- as.character(bundle$continuous_predictors)
derived_vars <- as.character(bundle$derived_predictors)
imputation_vars <- as.character(bundle$imputation_variables)

expected_predictor_vars <- c(
  "Sex", "Smoke", "Age", "AG_ratio", "eGFR", "Fbg", "PAR", "TC",
  "GGT", "ALP", "NPR", "INR", "HB", "APTT", "LDH", "TyG"
)

expected_imputation_vars <- c(
  "Sex", "Smoke", "Age", "TP", "ALB", "CR", "Fbg", "PLT", "TC",
  "GGT", "ALP", "NE_C", "PT", "HB", "APTT", "LDH", "GLU", "TG"
)

expected_derived_vars <- c(
  "AG_ratio", "eGFR", "PAR", "NPR", "INR", "TyG"
)

if (!identical(predictor_vars, expected_predictor_vars)) {
  stop(
    paste0(
      "Frozen predictor set does not match the locked 16-variable specification.\n",
      "Expected: ", paste(expected_predictor_vars, collapse = ", "), "\n",
      "Found: ", paste(predictor_vars, collapse = ", ")
    )
  )
}

if (!identical(imputation_vars, expected_imputation_vars)) {
  stop(
    paste0(
      "Frozen MICE input does not match the locked 18 BASIC/source variables.\n",
      "Expected: ", paste(expected_imputation_vars, collapse = ", "), "\n",
      "Found: ", paste(imputation_vars, collapse = ", ")
    )
  )
}

if (!identical(derived_vars, expected_derived_vars)) {
  stop(
    paste0(
      "Frozen derived-predictor set is inconsistent.\n",
      "Expected: ", paste(expected_derived_vars, collapse = ", "), "\n",
      "Found: ", paste(derived_vars, collapse = ", ")
    )
  )
}

SEX_FEMALE_CODE <- as.numeric(bundle$egfr_settings$female_code)
SEX_MALE_CODE <- as.numeric(bundle$egfr_settings$male_code)
EGFR_RACE_COEFFICIENT <- as.numeric(bundle$egfr_settings$race_coefficient)

if (!identical(SEX_FEMALE_CODE, 0) ||
    !identical(SEX_MALE_CODE, 1) ||
    abs(EGFR_RACE_COEFFICIENT - 1.0) > 1e-12) {
  stop(
    paste0(
      "Frozen eGFR coding does not match the locked specification: ",
      "Sex=0 female, Sex=1 male, race coefficient=1.0."
    )
  )
}

if (length(bundle$final_svm_fits) != MICE_M ||
    length(bundle$final_xgb_fits) != MICE_M ||
    length(bundle$final_completed_basic_datasets) != MICE_M) {
  stop(
    "Frozen bundle does not contain exactly MICE_M final development datasets/models."
  )
}

# ------------------------------------------------------------------------------
# 3B. STRICT SYNCHRONIZATION AUDIT:
#     baseline-model object must belong to THIS exact frozen stack bundle.
# ------------------------------------------------------------------------------

cat("\n==================== BASELINE/STACK SYNCHRONIZATION AUDIT ====================\n")

# Baseline preprocessing specification must match the frozen stack specification.
if (is.null(baseline_results$imputation_variables) ||
    !identical(as.character(baseline_results$imputation_variables), imputation_vars)) {
  stop(
    paste0(
      "baseline-model object does not use the same 18 BASIC MICE variables as ",
      "the frozen stacking-model bundle. Rerun the baseline-model script using ",
      "the final stacking-model bundle before external validation."
    )
  )
}

if (is.null(baseline_results$derived_predictors) ||
    !setequal(as.character(baseline_results$derived_predictors), derived_vars)) {
  stop(
    "baseline-model object does not use the same six derived-after-MICE predictors."
  )
}

if (is.null(baseline_results$egfr_settings)) {
  stop("baseline-model object is missing egfr_settings.")
}

if (abs(as.numeric(baseline_results$egfr_settings$female_code) - SEX_FEMALE_CODE) > 1e-12 ||
    abs(as.numeric(baseline_results$egfr_settings$male_code) - SEX_MALE_CODE) > 1e-12 ||
    abs(as.numeric(baseline_results$egfr_settings$race_coefficient) -
        EGFR_RACE_COEFFICIENT) > 1e-12) {
  stop("baseline-model object eGFR coding differs from the frozen stacking bundle.")
}

# Predictor-set audit.
expected_simple_predictor_sets <- list(
  age_only = c("Age"),
  demographics = c("Age", "Sex", "Smoke"),
  biomarkers_only = setdiff(predictor_vars, c("Age", "Sex", "Smoke")),
  demographics_plus_biomarkers = predictor_vars
)

if (is.null(baseline_results$predictor_sets)) {
  stop("baseline-model object is missing predictor_sets.")
}

for (nm in names(expected_simple_predictor_sets)) {
  if (is.null(baseline_results$predictor_sets[[nm]]) ||
      !identical(
        as.character(baseline_results$predictor_sets[[nm]]),
        as.character(expected_simple_predictor_sets[[nm]])
      )) {
    stop(
      paste0(
        "baseline-model predictor set mismatch for '", nm,
        "'. Rerun the baseline-model script with the final stack bundle."
      )
    )
  }
}

# Strongest version check: the Stack OOF probabilities stored in the baseline-model
# result must be exactly the Stack OOF probabilities stored in the frozen stacking-model bundle.
if (is.null(baseline_results$oof_pooled)) {
  stop("baseline-model object is missing oof_pooled.")
}

baseline_oof_sync <- as_tibble(baseline_results$oof_pooled) |>
  arrange(.row)

stack_oof_sync <- as_tibble(bundle$meta_oof_predictions) |>
  arrange(.row)

required_baseline_oof_cols <- c(".row", "PDstage", "Stack")
required_stack_oof_cols <- c(".row", "PDstage", "Stack_calibrated")

if (!all(required_baseline_oof_cols %in% names(baseline_oof_sync))) {
  stop(
    "baseline-model oof_pooled lacks .row, PDstage, or Stack."
  )
}

if (!all(required_stack_oof_cols %in% names(stack_oof_sync))) {
  stop(
    "Frozen stack meta_oof_predictions lacks .row, PDstage, or Stack_calibrated."
  )
}

if (nrow(baseline_oof_sync) != nrow(stack_oof_sync) ||
    !identical(
      as.integer(baseline_oof_sync$.row),
      as.integer(stack_oof_sync$.row)
    ) ||
    !identical(
      as.character(baseline_oof_sync$PDstage),
      as.character(stack_oof_sync$PDstage)
    )) {
  stop(
    "baseline-model OOF rows/outcomes do not align with the current frozen stack bundle."
  )
}

stack_oof_max_abs_diff <- max(
  abs(
    as.numeric(baseline_oof_sync$Stack) -
      as.numeric(stack_oof_sync$Stack_calibrated)
  ),
  na.rm = TRUE
)

if (!is.finite(stack_oof_max_abs_diff) || stack_oof_max_abs_diff > 1e-12) {
  stop(
    paste0(
      "baseline-model object was not generated from the current frozen stack bundle: ",
      "Stack OOF probabilities differ (max absolute difference = ",
      signif(stack_oof_max_abs_diff, 6), "). ",
      "Rerun the baseline-model script before external validation."
    )
  )
}

# Fold-map audit. This guards against an inconsistent baseline-model run that happened to use
# the same final Stack predictions but a different simple-model OOF partition.
parse_stack_fold_number <- function(x) {
  suppressWarnings(
    as.integer(
      sub("^Fold([0-9]+)_Imp[0-9]+$", "\\1", as.character(x))
    )
  )
}

stack_fold_source <- as_tibble(bundle$svm_oof_predictions_by_imputation)

if (!all(c(".row", "id") %in% names(stack_fold_source))) {
  stop("Frozen stack SVM OOF object lacks .row/id needed for fold synchronization.")
}

current_stack_fold_map <- stack_fold_source |>
  transmute(
    .row = as.integer(.row),
    Fold = parse_stack_fold_number(id)
  ) |>
  group_by(.row) |>
  summarise(
    Fold = first(Fold),
    N_unique_folds = n_distinct(Fold),
    MI_N = n(),
    .groups = "drop"
  ) |>
  arrange(.row)

if (anyNA(current_stack_fold_map$Fold) ||
    any(current_stack_fold_map$N_unique_folds != 1) ||
    any(current_stack_fold_map$MI_N != MICE_M)) {
  stop("Could not reconstruct a valid frozen stack fold map.")
}

if (is.null(baseline_results$fold_map)) {
  stop("baseline-model object is missing fold_map.")
}

baseline_fold_map_sync <- as_tibble(baseline_results$fold_map) |>
  transmute(
    .row = as.integer(.row),
    Fold = as.integer(Fold)
  ) |>
  arrange(.row)

if (nrow(baseline_fold_map_sync) != nrow(current_stack_fold_map) ||
    !identical(
      baseline_fold_map_sync$.row,
      current_stack_fold_map$.row
    ) ||
    !identical(
      baseline_fold_map_sync$Fold,
      current_stack_fold_map$Fold
    )) {
  stop(
    "baseline-model fold map differs from the current frozen stacking outer folds."
  )
}

cat(
  "Synchronization audit PASSED:",
  "same MICE specification, predictor sets, Stack OOF probabilities and outer folds.\n"
)

# ------------------------------------------------------------------------------
# 4. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

to_stage_factor <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))

  out <- dplyr::case_when(
    x_chr %in% c("stage1", "1") ~ "Stage1",
    x_chr %in% c("stage0", "0") ~ "Stage0",
    TRUE ~ NA_character_
  )

  if (anyNA(out)) {
    stop("PDstage contains missing or unrecognized values.")
  }

  factor(out, levels = c("Stage1", "Stage0"))
}


derive_model_predictors <- function(
    df,
    strict = TRUE,
    female_code = SEX_FEMALE_CODE,
    male_code = SEX_MALE_CODE,
    race_coefficient = EGFR_RACE_COEFFICIENT
) {
  missing_basic <- setdiff(imputation_vars, names(df))

  if (length(missing_basic) > 0) {
    stop(
      paste0(
        "Cannot derive model predictors. Missing BASIC columns: ",
        paste(missing_basic, collapse = ", ")
      )
    )
  }

  out <- as_tibble(df) |>
    mutate(
      across(
        all_of(imputation_vars),
        ~ suppressWarnings(as.numeric(.x))
      )
    )

  sex_nonmissing <- out$Sex[!is.na(out$Sex)]
  if (length(sex_nonmissing) > 0 &&
      !all(sex_nonmissing %in% c(female_code, male_code))) {
    stop(
      paste0(
        "Sex contains values outside ",
        female_code, " (female) and ",
        male_code, " (male)."
      )
    )
  }

  n <- nrow(out)

  scr_mg_dl <- out$CR / 88.4
  is_female <- out$Sex == female_code
  is_male <- out$Sex == male_code

  kappa <- ifelse(is_female, 0.7, ifelse(is_male, 0.9, NA_real_))
  alpha <- ifelse(is_female, -0.329, ifelse(is_male, -0.411, NA_real_))
  sex_coefficient <- ifelse(
    is_female,
    1.018,
    ifelse(is_male, 1.0, NA_real_)
  )

  AG_ratio_calc <- rep(NA_real_, n)
  eGFR_calc <- rep(NA_real_, n)
  PAR_calc <- rep(NA_real_, n)
  NPR_calc <- rep(NA_real_, n)
  INR_calc <- rep(NA_real_, n)
  TyG_calc <- rep(NA_real_, n)

  ag_denom <- out$TP - out$ALB
  idx_ag <- !is.na(out$ALB) & !is.na(ag_denom) & ag_denom > 0
  AG_ratio_calc[idx_ag] <- out$ALB[idx_ag] / ag_denom[idx_ag]

  idx_egfr <- (
    !is.na(scr_mg_dl) & scr_mg_dl > 0 &
      !is.na(out$Age) & out$Age >= 0 &
      !is.na(kappa) & !is.na(alpha) & !is.na(sex_coefficient)
  )

  if (any(idx_egfr)) {
    scr_ratio_valid <- scr_mg_dl[idx_egfr] / kappa[idx_egfr]

    eGFR_calc[idx_egfr] <-
      141 *
      pmin(scr_ratio_valid, 1)^alpha[idx_egfr] *
      pmax(scr_ratio_valid, 1)^(-1.209) *
      0.993^out$Age[idx_egfr] *
      sex_coefficient[idx_egfr] *
      race_coefficient
  }

  idx_par <- !is.na(out$PLT) & !is.na(out$ALB) & out$ALB > 0
  PAR_calc[idx_par] <- out$PLT[idx_par] / out$ALB[idx_par]

  idx_npr <- !is.na(out$NE_C) & !is.na(out$PLT) & out$PLT > 0
  NPR_calc[idx_npr] <- out$NE_C[idx_npr] / out$PLT[idx_npr]

  idx_inr <- !is.na(out$PT) & out$PT > 0
  INR_calc[idx_inr] <- (out$PT[idx_inr] / 13.1)^1.31

  tyg_arg <- 1594.26 * out$TG * out$GLU / 2
  idx_tyg <- !is.na(tyg_arg) & tyg_arg > 0
  TyG_calc[idx_tyg] <- log(tyg_arg[idx_tyg])

  derived <- tibble(
    Sex = out$Sex,
    Smoke = out$Smoke,
    Age = out$Age,
    AG_ratio = AG_ratio_calc,
    eGFR = eGFR_calc,
    Fbg = out$Fbg,
    PAR = PAR_calc,
    TC = out$TC,
    GGT = out$GGT,
    ALP = out$ALP,
    NPR = NPR_calc,
    INR = INR_calc,
    HB = out$HB,
    APTT = out$APTT,
    LDH = out$LDH,
    TyG = TyG_calc
  ) |>
    select(all_of(predictor_vars))

  if (strict && anyNA(derived)) {
    bad_cols <- names(derived)[vapply(derived, anyNA, logical(1))]
    stop(
      paste0(
        "Derived model predictors contain missing/invalid values after MICE: ",
        paste(bad_cols, collapse = ", "),
        ". Check TP>ALB, CR>0, ALB>0, PLT>0, PT>0, TG>0, GLU>0 and Sex coding."
      )
    )
  }

  if (strict) {
    nonfinite_cols <- names(derived)[
      vapply(
        derived,
        function(x) any(!is.finite(as.numeric(x))),
        logical(1)
      )
    ]

    if (length(nonfinite_cols) > 0) {
      stop(
        paste0(
          "Non-finite derived/model values found after MICE: ",
          paste(nonfinite_cols, collapse = ", ")
        )
      )
    }
  }

  derived
}


mice_impute_predictors <- function(
    predictor_df,
    ignore_vec = NULL,
    seed = SEED,
    m = MICE_M,
    maxit = MICE_MAXIT,
    method_name = MICE_METHOD,
    print_flag = MICE_PRINT_FLAG,
    eps = MICE_EPS
) {
  dat <- as.data.frame(predictor_df)

  for (nm in names(dat)) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  }

  if (is.null(ignore_vec)) {
    ignore_vec <- rep(FALSE, nrow(dat))
  }

  if (length(ignore_vec) != nrow(dat)) {
    stop("ignore_vec length does not match number of rows.")
  }

  if (!anyNA(dat)) {
    completed_list <- lapply(
      seq_len(m),
      function(j) as_tibble(dat)
    )
    names(completed_list) <- paste0("Imp", seq_len(m))

    return(
      list(
        completed_list = completed_list,
        mids = NULL,
        method = setNames(rep("", ncol(dat)), names(dat)),
        predictorMatrix = NULL,
        loggedEvents = NULL,
        m = m
      )
    )
  }

  method_vec <- mice::make.method(dat)
  method_vec[] <- ""

  has_missing <- vapply(dat, anyNA, logical(1))
  method_vec[has_missing] <- method_name

  pred_mat <- mice::make.predictorMatrix(dat)
  diag(pred_mat) <- 0

  set.seed(seed)

  imp <- mice::mice(
    dat,
    method = method_vec,
    predictorMatrix = pred_mat,
    m = m,
    maxit = maxit,
    ignore = ignore_vec,
    seed = seed,
    printFlag = print_flag,
    eps = eps
  )

  completed_list <- lapply(
    seq_len(m),
    function(j) {
      out <- mice::complete(imp, action = j) |>
        as_tibble()

      if (anyNA(out)) {
        stop(
          sprintf(
            "MICE completed dataset %d still contains missing values.",
            j
          )
        )
      }
      out
    }
  )

  names(completed_list) <- paste0("Imp", seq_len(m))

  list(
    completed_list = completed_list,
    mids = imp,
    method = method_vec,
    predictorMatrix = pred_mat,
    loggedEvents = imp$loggedEvents,
    m = m
  )
}


make_roc <- function(truth, prob) {
  y <- ifelse(as.character(truth) == "Stage1", 1, 0)

  pROC::roc(
    response = y,
    predictor = prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}


safe_div <- function(num, den) {
  ifelse(den == 0, NA_real_, num / den)
}


exact_binomial_ci <- function(successes, total, conf_level = PERFORMANCE_CI_LEVEL) {
  if (!is.finite(total) || total <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  ci <- stats::binom.test(
    x = successes,
    n = total,
    conf.level = conf_level
  )$conf.int

  c(lower = as.numeric(ci[1]), upper = as.numeric(ci[2]))
}


fit_calibration_stats <- function(y, prob) {
  p_clip <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  lp <- qlogis(p_clip)

  cal_intercept_fit <- suppressWarnings(
    glm(
      y ~ 1,
      offset = lp,
      family = binomial(link = "logit")
    )
  )

  cal_slope_fit <- suppressWarnings(
    glm(
      y ~ lp,
      family = binomial(link = "logit")
    )
  )

  intercept_value <- unname(coef(cal_intercept_fit)[1])
  slope_value <- unname(coef(cal_slope_fit)["lp"])

  c(
    Calibration_Intercept = intercept_value,
    Calibration_Slope = slope_value
  )
}


bootstrap_brier_calibration_ci <- function(
    y,
    prob,
    n_boot = PERFORMANCE_BOOT_N,
    conf_level = PERFORMANCE_CI_LEVEL,
    seed = SEED + PERFORMANCE_BOOT_SEED_OFFSET
) {
  if (length(y) != length(prob)) {
    stop("bootstrap_brier_calibration_ci: y and prob lengths differ.")
  }

  n <- length(y)

  if (n < 2) {
    stop("bootstrap_brier_calibration_ci: at least 2 observations are required.")
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 200) {
    stop("PERFORMANCE_BOOT_N should be an integer >= 200; 2000 is recommended.")
  }

  set.seed(seed)

  boot_mat <- matrix(
    NA_real_,
    nrow = as.integer(n_boot),
    ncol = 3,
    dimnames = list(
      NULL,
      c("Brier", "Calibration_Intercept", "Calibration_Slope")
    )
  )

  for (b in seq_len(as.integer(n_boot))) {
    idx <- sample.int(n, size = n, replace = TRUE)

    y_b <- y[idx]
    p_b <- prob[idx]

    boot_mat[b, "Brier"] <- mean((p_b - y_b)^2)

    # With these sample sizes, a one-class bootstrap sample is essentially
    # impossible, but keep this guard so the script remains robust.
    if (length(unique(y_b)) < 2) {
      next
    }

    cal_b <- tryCatch(
      fit_calibration_stats(y_b, p_b),
      error = function(e) c(
        Calibration_Intercept = NA_real_,
        Calibration_Slope = NA_real_
      )
    )

    boot_mat[b, "Calibration_Intercept"] <- cal_b[["Calibration_Intercept"]]
    boot_mat[b, "Calibration_Slope"] <- cal_b[["Calibration_Slope"]]
  }

  alpha <- (1 - conf_level) / 2
  probs_ci <- c(alpha, 1 - alpha)

  percentile_ci <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) < max(100, 0.5 * n_boot)) {
      warning(
        paste0(
          "Too few valid bootstrap replicates for a stable percentile CI: ",
          length(x), " / ", n_boot
        )
      )
    }

    if (length(x) == 0) {
      return(c(lower = NA_real_, upper = NA_real_))
    }

    q <- stats::quantile(
      x,
      probs = probs_ci,
      na.rm = TRUE,
      names = FALSE,
      type = 6
    )

    c(lower = q[1], upper = q[2])
  }

  brier_ci <- percentile_ci(boot_mat[, "Brier"])
  intercept_ci <- percentile_ci(boot_mat[, "Calibration_Intercept"])
  slope_ci <- percentile_ci(boot_mat[, "Calibration_Slope"])

  list(
    Brier_CI = brier_ci,
    Calibration_Intercept_CI = intercept_ci,
    Calibration_Slope_CI = slope_ci,
    Valid_Brier_Boot = sum(is.finite(boot_mat[, "Brier"])),
    Valid_Intercept_Boot = sum(is.finite(boot_mat[, "Calibration_Intercept"])),
    Valid_Slope_Boot = sum(is.finite(boot_mat[, "Calibration_Slope"]))
  )
}


calculate_performance <- function(
    truth,
    prob,
    model_name,
    threshold,
    threshold_source,
    ece_bins = ECE_BINS,
    ci_level = PERFORMANCE_CI_LEVEL,
    n_boot = PERFORMANCE_BOOT_N
) {
  truth_chr <- as.character(truth)
  y <- ifelse(truth_chr == "Stage1", 1, 0)

  if (!is.finite(threshold)) {
    stop(paste0(model_name, ": frozen threshold is non-finite."))
  }

  roc_obj <- make_roc(truth, prob)
  auc_value <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(
    pROC::ci.auc(
      roc_obj,
      method = "delong",
      conf.level = ci_level
    )
  )

  pred <- ifelse(prob >= threshold, 1, 0)

  TP <- sum(pred == 1 & y == 1)
  TN <- sum(pred == 0 & y == 0)
  FP <- sum(pred == 1 & y == 0)
  FN <- sum(pred == 0 & y == 1)

  sensitivity <- safe_div(TP, TP + FN)
  specificity <- safe_div(TN, TN + FP)
  accuracy <- safe_div(TP + TN, TP + TN + FP + FN)
  ppv <- safe_div(TP, TP + FP)
  npv <- safe_div(TN, TN + FN)
  f1 <- safe_div(2 * ppv * sensitivity, ppv + sensitivity)
  gmean <- sqrt(sensitivity * specificity)

  # Exact binomial confidence intervals for threshold-based sensitivity/specificity.
  sensitivity_ci <- exact_binomial_ci(
    successes = TP,
    total = TP + FN,
    conf_level = ci_level
  )

  specificity_ci <- exact_binomial_ci(
    successes = TN,
    total = TN + FP,
    conf_level = ci_level
  )

  brier <- mean((prob - y)^2)

  cal_stats <- fit_calibration_stats(y, prob)
  cal_intercept <- unname(cal_stats[["Calibration_Intercept"]])
  cal_slope <- unname(cal_stats[["Calibration_Slope"]])

  # Bootstrap uncertainty for Brier and calibration parameters.
  # The same deterministic bootstrap seed is deliberately reused across models,
  # so the row resamples are aligned across M1-M5.
  boot_uncertainty <- bootstrap_brier_calibration_ci(
    y = y,
    prob = prob,
    n_boot = n_boot,
    conf_level = ci_level,
    seed = SEED + PERFORMANCE_BOOT_SEED_OFFSET
  )

  bins <- cut(
    prob,
    breaks = seq(0, 1, length.out = ece_bins + 1),
    include.lowest = TRUE,
    right = TRUE
  )

  ece_df <- tibble(
    prob = prob,
    y = y,
    bin = bins
  ) |>
    group_by(bin) |>
    summarise(
      n = n(),
      mean_prob = mean(prob),
      mean_obs = mean(y),
      .groups = "drop"
    )

  ece <- sum(
    (ece_df$n / length(prob)) *
      abs(ece_df$mean_prob - ece_df$mean_obs)
  )

  tibble(
    Model = model_name,
    External_N = length(y),
    Stage1_N = sum(y == 1),
    Stage0_N = sum(y == 0),
    Stage1_Prevalence = mean(y),

    AUC = auc_value,
    AUC_95CI_Lower = auc_ci[1],
    AUC_95CI_Upper = auc_ci[3],

    Threshold = threshold,
    Threshold_Source = threshold_source,

    Sensitivity = sensitivity,
    Sensitivity_95CI_Lower = sensitivity_ci[["lower"]],
    Sensitivity_95CI_Upper = sensitivity_ci[["upper"]],

    Specificity = specificity,
    Specificity_95CI_Lower = specificity_ci[["lower"]],
    Specificity_95CI_Upper = specificity_ci[["upper"]],

    Accuracy = accuracy,
    PPV = ppv,
    NPV = npv,
    F1 = f1,
    G_mean = gmean,

    Brier = brier,
    Brier_95CI_Lower = boot_uncertainty$Brier_CI[["lower"]],
    Brier_95CI_Upper = boot_uncertainty$Brier_CI[["upper"]],

    ECE = ece,

    Calibration_Intercept = cal_intercept,
    Calibration_Intercept_95CI_Lower =
      boot_uncertainty$Calibration_Intercept_CI[["lower"]],
    Calibration_Intercept_95CI_Upper =
      boot_uncertainty$Calibration_Intercept_CI[["upper"]],

    Calibration_Slope = cal_slope,
    Calibration_Slope_95CI_Lower =
      boot_uncertainty$Calibration_Slope_CI[["lower"]],
    Calibration_Slope_95CI_Upper =
      boot_uncertainty$Calibration_Slope_CI[["upper"]],

    Bootstrap_Replicates = n_boot,
    Bootstrap_Valid_Brier = boot_uncertainty$Valid_Brier_Boot,
    Bootstrap_Valid_Calibration_Intercept =
      boot_uncertainty$Valid_Intercept_Boot,
    Bootstrap_Valid_Calibration_Slope =
      boot_uncertainty$Valid_Slope_Boot,

    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  )
}


paired_delong_compare <- function(
    truth,
    prob_new,
    prob_old,
    new_name,
    old_name
) {
  roc_new <- make_roc(truth, prob_new)
  roc_old <- make_roc(truth, prob_old)

  test <- pROC::roc.test(
    roc_new,
    roc_old,
    method = "delong",
    paired = TRUE,
    conf.int = TRUE
  )

  auc_new <- as.numeric(pROC::auc(roc_new))
  auc_old <- as.numeric(pROC::auc(roc_old))

  ci_diff <- if (!is.null(test$conf.int)) {
    as.numeric(test$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    Comparison = paste0(new_name, " vs ", old_name),
    Old_Model = old_name,
    New_Model = new_name,
    AUC_Old = auc_old,
    AUC_New = auc_new,
    Delta_AUC_New_minus_Old = auc_new - auc_old,
    Delta_AUC_95CI_Lower = ci_diff[1],
    Delta_AUC_95CI_Upper = ci_diff[2],
    DeLong_P = as.numeric(test$p.value)
  )
}


fit_simple_logistic <- function(
    train_df,
    predictor_names,
    continuous_names
) {
  dat <- train_df

  vars_to_scale <- intersect(predictor_names, continuous_names)

  scaler <- list()

  if (length(vars_to_scale) > 0) {
    for (nm in vars_to_scale) {
      mu <- mean(dat[[nm]])
      sigma <- stats::sd(dat[[nm]])

      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
        stop(
          paste0(
            "Cannot standardize ", nm,
            " in full-development simple model."
          )
        )
      }

      scaler[[nm]] <- c(mean = mu, sd = sigma)
      dat[[nm]] <- (dat[[nm]] - mu) / sigma
    }
  }

  dat <- dat |>
    mutate(y = ifelse(PDstage == "Stage1", 1, 0))

  formula_obj <- reformulate(
    termlabels = predictor_names,
    response = "y"
  )

  fit <- suppressWarnings(
    glm(
      formula = formula_obj,
      data = dat,
      family = binomial(link = "logit")
    )
  )

  if (!isTRUE(fit$converged)) {
    warning(
      paste0(
        "Simple logistic model did not converge: ",
        paste(predictor_names, collapse = " + ")
      )
    )
  }

  list(
    fit = fit,
    scaler = scaler,
    predictors = predictor_names
  )
}


predict_simple_logistic <- function(fit_obj, new_df) {
  dat <- new_df

  for (nm in names(fit_obj$scaler)) {
    mu <- fit_obj$scaler[[nm]][["mean"]]
    sigma <- fit_obj$scaler[[nm]][["sd"]]
    dat[[nm]] <- (dat[[nm]] - mu) / sigma
  }

  prob <- suppressWarnings(
    predict(
      fit_obj$fit,
      newdata = dat,
      type = "response"
    )
  )

  if (any(!is.finite(prob))) {
    stop("Non-finite simple-model prediction detected.")
  }

  pmin(pmax(as.numeric(prob), 1e-12), 1 - 1e-12)
}


find_saved_simple_fit_container <- function(x) {
  candidates <- c(
    "final_simple_model_fits",
    "frozen_simple_model_fits",
    "simple_model_final_fits",
    "final_model_fits",
    "simple_model_fits"
  )

  for (nm in candidates) {
    obj <- x[[nm]]
    if (is.list(obj) && length(obj) > 0) {
      return(list(field = nm, fits = obj))
    }
  }

  NULL
}


validate_saved_simple_fit_container <- function(
    fit_container,
    simple_model_specs,
    m
) {
  if (is.null(fit_container) || !is.list(fit_container)) return(FALSE)

  if (!all(names(simple_model_specs) %in% names(fit_container))) {
    return(FALSE)
  }

  for (model_nm in names(simple_model_specs)) {
    model_fits <- fit_container[[model_nm]]

    if (!is.list(model_fits) || length(model_fits) != m) {
      return(FALSE)
    }

    for (j in seq_len(m)) {
      fit_obj <- model_fits[[j]]

      if (!is.list(fit_obj) ||
          !all(c("fit", "scaler", "predictors") %in% names(fit_obj))) {
        return(FALSE)
      }

      if (!identical(
        as.character(fit_obj$predictors),
        as.character(simple_model_specs[[model_nm]])
      )) {
        return(FALSE)
      }
    }
  }

  TRUE
}

# ------------------------------------------------------------------------------
# 5. READ AND VALIDATE EXTERNAL DATA
# ------------------------------------------------------------------------------

cat("\n==================== READ EXTERNAL DATA ====================\n")

external_raw <- readxl::read_excel(
  EXTERNAL_DATA_PATH,
  na = c("", "NA", "N/A", "na", "NaN", "NULL")
) |>
  as_tibble()

required_external_cols <- c("PDstage", imputation_vars)
missing_external_cols <- setdiff(required_external_cols, names(external_raw))

if (length(missing_external_cols) > 0) {
  stop(
    paste0(
      "External dataset is missing required columns: ",
      paste(missing_external_cols, collapse = ", ")
    )
  )
}

external_truth <- to_stage_factor(external_raw$PDstage)
n_external <- nrow(external_raw)

if (n_external < 2) {
  stop("External validation dataset must contain at least 2 rows.")
}

if (length(unique(external_truth)) < 2) {
  stop("External validation dataset must contain both Stage1 and Stage0.")
}

external_basic_raw <- external_raw |>
  select(all_of(imputation_vars)) |>
  mutate(
    across(
      everything(),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

for (nm in c("Sex", "Smoke")) {
  vals <- sort(unique(external_basic_raw[[nm]][!is.na(external_basic_raw[[nm]])]))
  if (!all(vals %in% c(0, 1))) {
    stop(
      paste0(
        "External ", nm,
        " must be coded 0/1. Observed nonmissing values: ",
        paste(vals, collapse = ", ")
      )
    )
  }
}

external_missingness <- tibble(
  Variable = imputation_vars,
  Missing_N = vapply(
    external_basic_raw,
    function(x) sum(is.na(x)),
    integer(1)
  )
) |>
  mutate(Missing_Rate = Missing_N / n_external)

external_class_distribution <- tibble(
  Class = c("Stage1", "Stage0"),
  N = c(
    sum(external_truth == "Stage1"),
    sum(external_truth == "Stage0")
  ),
  Proportion = N / n_external
)

cat("External N:", n_external, "\n")
print(external_class_distribution)
cat("\nExternal BASIC-variable missingness:\n")
print(external_missingness |> filter(Missing_N > 0))

# ------------------------------------------------------------------------------
# 6. CREATE EXTERNAL TRUE-MI BASIC DATASETS
#
# Development rows estimate the imputation model.
# External rows are ignore=TRUE and cannot alter imputation parameters.
# External PDstage is not supplied to mice().
# ------------------------------------------------------------------------------

cat("\n==================== EXTERNAL TRUE-MI ====================\n")

dev_imp_raw <- as_tibble(bundle$raw_development_imputation_data) |>
  select(all_of(imputation_vars)) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

needs_external_imputation <- anyNA(external_basic_raw)

external_mice_logged_events <- NULL
external_mice_method <- NULL
external_mice_predictor_matrix <- NULL

if (!needs_external_imputation) {

  cat("No missing BASIC variables in external cohort; no external MICE call is needed.\n")

  external_basic_complete_list <- lapply(
    seq_len(MICE_M),
    function(j) external_basic_raw
  )
  names(external_basic_complete_list) <- paste0("Imp", seq_len(MICE_M))

} else {

  cat(
    "Missing BASIC variables detected. External rows will be imputed using ",
    "development-estimated MICE relations only.\n",
    sep = ""
  )

  combined_imp_data <- bind_rows(
    dev_imp_raw,
    external_basic_raw
  )

  ignore_vec <- c(
    rep(FALSE, nrow(dev_imp_raw)),
    rep(TRUE, nrow(external_basic_raw))
  )

  ext_imp <- mice_impute_predictors(
    predictor_df = combined_imp_data,
    ignore_vec = ignore_vec,
    seed = SEED + 5000,
    m = MICE_M,
    maxit = MICE_MAXIT,
    method_name = MICE_METHOD,
    print_flag = MICE_PRINT_FLAG,
    eps = MICE_EPS
  )

  external_mice_logged_events <- ext_imp$loggedEvents
  external_mice_method <- ext_imp$method
  external_mice_predictor_matrix <- ext_imp$predictorMatrix

  external_basic_complete_list <- vector("list", MICE_M)

  for (j in seq_len(MICE_M)) {
    complete_j <- ext_imp$completed_list[[j]]

    external_basic_complete_list[[j]] <- as_tibble(
      complete_j[
        (nrow(dev_imp_raw) + 1):nrow(complete_j),
        ,
        drop = FALSE
      ]
    )
  }

  names(external_basic_complete_list) <- paste0("Imp", seq_len(MICE_M))
}

# Derive the six non-imputed predictors separately within every completed dataset.
external_model_complete_list <- lapply(
  external_basic_complete_list,
  function(x) {
    derive_model_predictors(
      x,
      strict = TRUE,
      female_code = SEX_FEMALE_CODE,
      male_code = SEX_MALE_CODE,
      race_coefficient = EGFR_RACE_COEFFICIENT
    )
  }
)

if (any(vapply(external_model_complete_list, anyNA, logical(1)))) {
  stop("At least one external completed model dataset still contains missing values.")
}

# ------------------------------------------------------------------------------
# 7. FROZEN STACKING MODEL PREDICTION
# ------------------------------------------------------------------------------

cat("\n==================== FROZEN STACK PREDICTION ====================\n")

svm_prob_mat <- matrix(
  NA_real_,
  nrow = n_external,
  ncol = MICE_M
)

xgb_prob_mat <- matrix(
  NA_real_,
  nrow = n_external,
  ncol = MICE_M
)

for (j in seq_len(MICE_M)) {

  svm_prob_mat[, j] <- predict(
    bundle$final_svm_fits[[j]],
    new_data = external_model_complete_list[[j]],
    type = "prob"
  )$.pred_Stage1

  xgb_prob_mat[, j] <- predict(
    bundle$final_xgb_fits[[j]],
    new_data = external_model_complete_list[[j]],
    type = "prob"
  )$.pred_Stage1
}

colnames(svm_prob_mat) <- paste0("Imp", seq_len(MICE_M))
colnames(xgb_prob_mat) <- paste0("Imp", seq_len(MICE_M))

svm_prob <- rowMeans(svm_prob_mat)
xgb_prob <- rowMeans(xgb_prob_mat)

svm_mi_sd <- if (MICE_M > 1) {
  apply(svm_prob_mat, 1, sd)
} else {
  rep(0, n_external)
}

xgb_mi_sd <- if (MICE_M > 1) {
  apply(xgb_prob_mat, 1, sd)
} else {
  rep(0, n_external)
}

meta_new <- tibble(
  SVM_z = (svm_prob - bundle$meta_scaler$SVM_mean) /
    bundle$meta_scaler$SVM_sd,
  XGB_z = (xgb_prob - bundle$meta_scaler$XGB_mean) /
    bundle$meta_scaler$XGB_sd
)

stack_raw <- predict(
  bundle$meta_fit,
  new_data = meta_new,
  type = "prob"
)$.pred_Stage1

stack_calibrated <- predict(
  bundle$platt_calibrator,
  newdata = data.frame(probs = stack_raw),
  type = "response"
)

stack_threshold <- as.numeric(bundle$final_threshold)

if (length(stack_threshold) != 1 || !is.finite(stack_threshold)) {
  stop("Frozen stacking threshold is missing or invalid.")
}

# ------------------------------------------------------------------------------
# 8. FROZEN M1-M4 SIMPLE LOGISTIC MODELS
#
# Preferred path:
#   use final M1-M4 fits saved by the UPDATED baseline-model object.
#
# Backward-compatible path:
#   if the RDS predates saved final fits, deterministically reconstruct the same
#   final models using only the frozen development MI datasets. This does NOT use
#   external outcome/predictors for model fitting and is therefore still a valid
#   external validation. The source is explicitly recorded in the audit output.
# ------------------------------------------------------------------------------

cat("\n==================== FROZEN SIMPLE MODELS ====================\n")

demographic_vars <- c("Age", "Sex", "Smoke")
biomarker_vars <- setdiff(predictor_vars, demographic_vars)

if (length(biomarker_vars) != 13) {
  stop(
    paste0(
      "Expected 13 biomarkers, found ",
      length(biomarker_vars), "."
    )
  )
}

simple_model_specs <- list(
  "M1 Age alone" = c("Age"),
  "M2 Age + Sex + Smoking" = demographic_vars,
  "M3 13 biomarkers only (logistic)" = biomarker_vars,
  "M4 Demographics + 13 biomarkers (logistic)" = predictor_vars
)

saved_simple <- find_saved_simple_fit_container(baseline_results)

use_saved_simple_fits <- FALSE
simple_fit_source <- NA_character_

if (!is.null(saved_simple) &&
    validate_saved_simple_fit_container(
      saved_simple$fits,
      simple_model_specs,
      MICE_M
    )) {

  simple_fit_lists <- saved_simple$fits
  use_saved_simple_fits <- TRUE
  simple_fit_source <- paste0(
    "Frozen M1-M4 fits loaded directly from baseline_results$",
    saved_simple$field
  )

  cat(simple_fit_source, "\n")

} else {

  simple_fit_source <- paste0(
    "Frozen-specification M1-M4 fits deterministically reconstructed from ",
    "the final development MI datasets because the baseline RDS does not contain ",
    "a recognized saved-fit container."
  )

  warning(
    paste0(
      simple_fit_source,
      " For a consistent reproducibility workflow, rerun the baseline-model script that ",
      "saves final simple-model fits, then rerun external validation."
    )
  )

  # Recover development outcomes in exact row order from the frozen stack object.
  dev_meta_oof <- as_tibble(bundle$meta_oof_predictions) |>
    arrange(.row)

  if (!all(c(".row", "PDstage") %in% names(dev_meta_oof))) {
    stop("bundle$meta_oof_predictions lacks .row or PDstage.")
  }

  n_dev <- nrow(dev_imp_raw)

  if (nrow(dev_meta_oof) != n_dev ||
      !identical(as.integer(dev_meta_oof$.row), seq_len(n_dev))) {
    stop("Development outcome rows cannot be aligned to final development imputations.")
  }

  development_truth <- factor(
    as.character(dev_meta_oof$PDstage),
    levels = c("Stage1", "Stage0")
  )

  simple_fit_lists <- setNames(
    lapply(simple_model_specs, function(x) vector("list", MICE_M)),
    names(simple_model_specs)
  )

  for (j in seq_len(MICE_M)) {

    dev_model_j <- derive_model_predictors(
      bundle$final_completed_basic_datasets[[j]],
      strict = TRUE,
      female_code = SEX_FEMALE_CODE,
      male_code = SEX_MALE_CODE,
      race_coefficient = EGFR_RACE_COEFFICIENT
    )

    dev_train_j <- bind_cols(
      tibble(PDstage = development_truth),
      dev_model_j
    )

    for (model_nm in names(simple_model_specs)) {

      simple_fit_lists[[model_nm]][[j]] <- fit_simple_logistic(
        train_df = dev_train_j,
        predictor_names = simple_model_specs[[model_nm]],
        continuous_names = continuous_vars
      )
    }
  }
}

# External prediction only. No fitting/calibration/threshold selection uses
# external outcome information.
simple_prob_mats <- setNames(
  lapply(
    simple_model_specs,
    function(x) matrix(
      NA_real_,
      nrow = n_external,
      ncol = MICE_M
    )
  ),
  names(simple_model_specs)
)

simple_coef_rows <- list()
coef_counter <- 0L

for (j in seq_len(MICE_M)) {

  ext_test_j <- external_model_complete_list[[j]]

  for (model_nm in names(simple_model_specs)) {

    fit_j <- simple_fit_lists[[model_nm]][[j]]

    simple_prob_mats[[model_nm]][, j] <- predict_simple_logistic(
      fit_obj = fit_j,
      new_df = ext_test_j
    )

    coef_counter <- coef_counter + 1L

    simple_coef_rows[[coef_counter]] <- tibble(
      Model = model_nm,
      Imputation = j,
      Term = names(coef(fit_j$fit)),
      Estimate = as.numeric(coef(fit_j$fit))
    )
  }
}

simple_coefficients <- bind_rows(simple_coef_rows)

for (model_nm in names(simple_prob_mats)) {
  colnames(simple_prob_mats[[model_nm]]) <- paste0(
    "Imp",
    seq_len(MICE_M)
  )
}

simple_pooled <- lapply(
  simple_prob_mats,
  rowMeans
)

simple_mi_sd <- lapply(
  simple_prob_mats,
  function(x) {
    if (MICE_M > 1) {
      apply(x, 1, sd)
    } else {
      rep(0, nrow(x))
    }
  }
)

# ------------------------------------------------------------------------------
# 9. LOAD FROZEN DEVELOPMENT THRESHOLDS FOR M1-M4
# ------------------------------------------------------------------------------

cat("\n==================== FROZEN THRESHOLDS ====================\n")

dev_perf <- as_tibble(baseline_results$performance_summary)

if (!all(c("Model", "Threshold") %in% names(dev_perf))) {
  stop(
    "Baseline performance_summary must contain Model and Threshold columns."
  )
}

get_frozen_simple_threshold <- function(model_name) {
  x <- dev_perf |>
    filter(Model == model_name) |>
    pull(Threshold)

  if (length(x) != 1 || !is.finite(as.numeric(x))) {
    stop(
      paste0(
        "Cannot recover one frozen development threshold for: ",
        model_name
      )
    )
  }

  as.numeric(x)
}

frozen_thresholds <- tibble(
  Model = c(
    names(simple_model_specs),
    "M5 SVM + XGBoost stacking"
  ),
  Threshold = c(
    vapply(
      names(simple_model_specs),
      get_frozen_simple_threshold,
      numeric(1)
    ),
    stack_threshold
  ),
  Threshold_Source = c(
    rep("Frozen development OOF Youden threshold", 4),
    "Frozen stacking development threshold"
  )
)

# Optional consistency check if M5 is present in the baseline result bundle.
baseline_stack_thr <- dev_perf |>
  filter(Model == "M5 SVM + XGBoost stacking") |>
  pull(Threshold)

if (length(baseline_stack_thr) != 1 ||
    !is.finite(as.numeric(baseline_stack_thr))) {
  stop(
    "baseline-model results do not contain one valid frozen M5 threshold."
  )
}

if (abs(as.numeric(baseline_stack_thr) - stack_threshold) > 1e-12) {
  stop(
    paste0(
      "baseline-model M5 threshold differs from bundle$final_threshold. ",
      "This indicates that the baseline RDS and frozen stack bundle are from ",
      "different fitted-model objects. Rerun the baseline-model script using the final ",
      "stack bundle before external validation."
    )
  )
}

print(frozen_thresholds)

# ------------------------------------------------------------------------------
# 10. BUILD EXTERNAL PREDICTION TABLE
# ------------------------------------------------------------------------------

prediction_table <- tibble(
  .row = seq_len(n_external),
  PDstage = external_truth,

  M1_Age_Only = simple_pooled[["M1 Age alone"]],
  M1_MI_SD = simple_mi_sd[["M1 Age alone"]],

  M2_Demographics = simple_pooled[["M2 Age + Sex + Smoking"]],
  M2_MI_SD = simple_mi_sd[["M2 Age + Sex + Smoking"]],

  M3_Biomarkers_Only = simple_pooled[["M3 13 biomarkers only (logistic)"]],
  M3_MI_SD = simple_mi_sd[["M3 13 biomarkers only (logistic)"]],

  M4_Demo_Biomarkers = simple_pooled[["M4 Demographics + 13 biomarkers (logistic)"]],
  M4_MI_SD = simple_mi_sd[["M4 Demographics + 13 biomarkers (logistic)"]],

  Stack_SVM_Base = svm_prob,
  Stack_SVM_MI_SD = svm_mi_sd,
  Stack_XGB_Base = xgb_prob,
  Stack_XGB_MI_SD = xgb_mi_sd,
  Stack_Raw = stack_raw,
  Stack_Calibrated = stack_calibrated
)

threshold_lookup <- setNames(
  frozen_thresholds$Threshold,
  frozen_thresholds$Model
)

prediction_table <- prediction_table |>
  mutate(
    M1_Class = ifelse(
      M1_Age_Only >= threshold_lookup[["M1 Age alone"]],
      "Stage1", "Stage0"
    ),
    M2_Class = ifelse(
      M2_Demographics >= threshold_lookup[["M2 Age + Sex + Smoking"]],
      "Stage1", "Stage0"
    ),
    M3_Class = ifelse(
      M3_Biomarkers_Only >= threshold_lookup[["M3 13 biomarkers only (logistic)"]],
      "Stage1", "Stage0"
    ),
    M4_Class = ifelse(
      M4_Demo_Biomarkers >= threshold_lookup[["M4 Demographics + 13 biomarkers (logistic)"]],
      "Stage1", "Stage0"
    ),
    M5_Stack_Class = ifelse(
      Stack_Calibrated >= stack_threshold,
      "Stage1", "Stage0"
    )
  )

# ------------------------------------------------------------------------------
# 11. EXTERNAL PERFORMANCE
# ------------------------------------------------------------------------------

cat("\n==================== EXTERNAL PERFORMANCE ====================\n")

performance_summary <- bind_rows(
  calculate_performance(
    truth = external_truth,
    prob = prediction_table$M1_Age_Only,
    model_name = "M1 Age alone",
    threshold = threshold_lookup[["M1 Age alone"]],
    threshold_source = "Frozen development OOF Youden"
  ),

  calculate_performance(
    truth = external_truth,
    prob = prediction_table$M2_Demographics,
    model_name = "M2 Age + Sex + Smoking",
    threshold = threshold_lookup[["M2 Age + Sex + Smoking"]],
    threshold_source = "Frozen development OOF Youden"
  ),

  calculate_performance(
    truth = external_truth,
    prob = prediction_table$M3_Biomarkers_Only,
    model_name = "M3 13 biomarkers only (logistic)",
    threshold = threshold_lookup[["M3 13 biomarkers only (logistic)"]],
    threshold_source = "Frozen development OOF Youden"
  ),

  calculate_performance(
    truth = external_truth,
    prob = prediction_table$M4_Demo_Biomarkers,
    model_name = "M4 Demographics + 13 biomarkers (logistic)",
    threshold = threshold_lookup[["M4 Demographics + 13 biomarkers (logistic)"]],
    threshold_source = "Frozen development OOF Youden"
  ),

  calculate_performance(
    truth = external_truth,
    prob = prediction_table$Stack_Calibrated,
    model_name = "M5 SVM + XGBoost stacking",
    threshold = stack_threshold,
    threshold_source = "Frozen stack development threshold"
  )
)

print(
  performance_summary |>
    select(
      Model,
      AUC,
      AUC_95CI_Lower,
      AUC_95CI_Upper,
      Sensitivity,
      Sensitivity_95CI_Lower,
      Sensitivity_95CI_Upper,
      Specificity,
      Specificity_95CI_Lower,
      Specificity_95CI_Upper,
      PPV,
      NPV,
      Brier,
      Brier_95CI_Lower,
      Brier_95CI_Upper,
      ECE,
      Calibration_Intercept,
      Calibration_Intercept_95CI_Lower,
      Calibration_Intercept_95CI_Upper,
      Calibration_Slope,
      Calibration_Slope_95CI_Lower,
      Calibration_Slope_95CI_Upper
    )
)

# ------------------------------------------------------------------------------
# 12. EXTERNAL PAIRED DELONG COMPARISONS
# ------------------------------------------------------------------------------

external_delong <- bind_rows(
  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$M2_Demographics,
    prob_old = prediction_table$M1_Age_Only,
    new_name = "M2 Age + Sex + Smoking",
    old_name = "M1 Age alone"
  ),

  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$M3_Biomarkers_Only,
    prob_old = prediction_table$M1_Age_Only,
    new_name = "M3 13 biomarkers only (logistic)",
    old_name = "M1 Age alone"
  ),

  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$M4_Demo_Biomarkers,
    prob_old = prediction_table$M2_Demographics,
    new_name = "M4 Demographics + 13 biomarkers (logistic)",
    old_name = "M2 Age + Sex + Smoking"
  ),

  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$M4_Demo_Biomarkers,
    prob_old = prediction_table$M3_Biomarkers_Only,
    new_name = "M4 Demographics + 13 biomarkers (logistic)",
    old_name = "M3 13 biomarkers only (logistic)"
  ),

  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$Stack_Calibrated,
    prob_old = prediction_table$M4_Demo_Biomarkers,
    new_name = "M5 SVM + XGBoost stacking",
    old_name = "M4 Demographics + 13 biomarkers (logistic)"
  ),

  paired_delong_compare(
    truth = external_truth,
    prob_new = prediction_table$Stack_Calibrated,
    prob_old = prediction_table$M3_Biomarkers_Only,
    new_name = "M5 SVM + XGBoost stacking",
    old_name = "M3 13 biomarkers only (logistic)"
  )
)

cat("\nExternal paired DeLong comparisons:\n")
print(external_delong)

# ------------------------------------------------------------------------------
# 13. CALIBRATION BIN TABLE
# ------------------------------------------------------------------------------

model_prob_map <- list(
  "M1 Age alone" = prediction_table$M1_Age_Only,
  "M2 Age + Sex + Smoking" = prediction_table$M2_Demographics,
  "M3 13 biomarkers only (logistic)" = prediction_table$M3_Biomarkers_Only,
  "M4 Demographics + 13 biomarkers (logistic)" = prediction_table$M4_Demo_Biomarkers,
  "M5 SVM + XGBoost stacking" = prediction_table$Stack_Calibrated
)

external_y <- ifelse(external_truth == "Stage1", 1, 0)

calibration_bins <- bind_rows(
  lapply(
    names(model_prob_map),
    function(model_nm) {
      prob <- model_prob_map[[model_nm]]

      # Quantile-based groups are used for plotting only.
      # Calibration intercept/slope and Brier are computed directly from all subjects.
      rank_group <- dplyr::ntile(prob, min(10, length(prob)))

      tibble(
        Model = model_nm,
        Probability = prob,
        Outcome = external_y,
        Group = rank_group
      ) |>
        group_by(Model, Group) |>
        summarise(
          N = n(),
          Mean_Predicted = mean(Probability),
          Observed_Rate = mean(Outcome),
          .groups = "drop"
        )
    }
  )
)

# ------------------------------------------------------------------------------
# 14. ROC PLOT
# ------------------------------------------------------------------------------

roc_plot_data <- bind_rows(
  lapply(
    names(model_prob_map),
    function(model_nm) {
      roc_obj <- make_roc(
        external_truth,
        model_prob_map[[model_nm]]
      )

      tibble(
        Model = model_nm,
        False_Positive_Rate = 1 - roc_obj$specificities,
        True_Positive_Rate = roc_obj$sensitivities
      )
    }
  )
)

auc_label_table <- performance_summary |>
  transmute(
    Model,
    Model_Label = sprintf(
      "%s (AUC %.3f, 95%% CI %.3f-%.3f)",
      Model,
      AUC,
      AUC_95CI_Lower,
      AUC_95CI_Upper
    )
  )

label_lookup <- setNames(
  auc_label_table$Model_Label,
  auc_label_table$Model
)

roc_plot_data <- roc_plot_data |>
  mutate(
    Model_Label = factor(
      label_lookup[Model],
      levels = label_lookup[names(model_prob_map)]
    )
  )

roc_plot <- ggplot(
  roc_plot_data,
  aes(
    x = False_Positive_Rate,
    y = True_Positive_Rate,
    group = Model_Label,
    linetype = Model_Label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "External Validation ROC Comparison",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    linetype = "Model"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  ) +
  guides(
    linetype = guide_legend(ncol = 1, byrow = TRUE)
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "external_roc_comparison.png"
  ),
  plot = roc_plot,
  width = 10,
  height = 8,
  dpi = 320
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "external_roc_comparison.pdf"
  ),
  plot = roc_plot,
  width = 10,
  height = 8
)

# ------------------------------------------------------------------------------
# 15. CALIBRATION PLOT
# ------------------------------------------------------------------------------

calibration_plot <- ggplot(
  calibration_bins,
  aes(
    x = Mean_Predicted,
    y = Observed_Rate,
    group = Model,
    linetype = Model,
    shape = Model
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    title = "External Validation Calibration",
    x = "Mean Predicted Probability",
    y = "Observed Stage1 Proportion",
    linetype = "Model",
    shape = "Model"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "external_calibration_comparison.png"
  ),
  plot = calibration_plot,
  width = 10,
  height = 8,
  dpi = 320
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "external_calibration_comparison.pdf"
  ),
  plot = calibration_plot,
  width = 10,
  height = 8
)

# ------------------------------------------------------------------------------
# 16. IMPUTATION-SPECIFIC PREDICTIONS FOR AUDIT
# ------------------------------------------------------------------------------

external_predictions_by_mi <- bind_rows(
  lapply(
    seq_len(MICE_M),
    function(j) {
      tibble(
        .row = seq_len(n_external),
        PDstage = external_truth,
        Imputation = j,

        M1_Age_Only = simple_prob_mats[["M1 Age alone"]][, j],
        M2_Demographics = simple_prob_mats[["M2 Age + Sex + Smoking"]][, j],
        M3_Biomarkers_Only =
          simple_prob_mats[["M3 13 biomarkers only (logistic)"]][, j],
        M4_Demo_Biomarkers =
          simple_prob_mats[["M4 Demographics + 13 biomarkers (logistic)"]][, j],

        SVM_Base = svm_prob_mat[, j],
        XGB_Base = xgb_prob_mat[, j]
      )
    }
  )
)

# ------------------------------------------------------------------------------
# 17. SETTINGS / AUDIT TABLES
# ------------------------------------------------------------------------------

settings_table <- tibble(
  Parameter = c(
    "External_data_path",
    "Stack_bundle_path",
    "Baseline_results_path",
    "External_N",
    "SEED",
    "MICE_M",
    "MICE_MAXIT",
    "MICE_EPS",
    "MICE_METHOD",
    "External_outcome_used_for_imputation",
    "External_rows_estimate_MICE_parameters",
    "Derived_predictors_imputed",
    "Sex_female_code",
    "Sex_male_code",
    "eGFR_race_coefficient",
    "External_threshold_reoptimization",
    "External_stack_recalibration",
    "Baseline_stack_sync_check",
    "Baseline_stack_OOF_max_abs_diff",
    "Simple_model_fit_source",
    "Saved_simple_fits_used",
    "MI_pooling_rule"
  ),
  Value = as.character(c(
    EXTERNAL_DATA_PATH,
    STACK_BUNDLE_PATH,
    BASELINE_RESULTS_PATH,
    n_external,
    SEED,
    MICE_M,
    MICE_MAXIT,
    MICE_EPS,
    MICE_METHOD,
    FALSE,
    FALSE,
    FALSE,
    SEX_FEMALE_CODE,
    SEX_MALE_CODE,
    EGFR_RACE_COEFFICIENT,
    FALSE,
    FALSE,
    "PASSED",
    stack_oof_max_abs_diff,
    simple_fit_source,
    use_saved_simple_fits,
    "Mean predicted probability across imputations per patient"
  ))
)


performance_ci_settings <- tibble(
  Parameter = c(
    "CI_level",
    "Sensitivity_specificity_CI_method",
    "Brier_CI_method",
    "Calibration_intercept_CI_method",
    "Calibration_slope_CI_method",
    "Bootstrap_replicates",
    "Bootstrap_seed",
    "Bootstrap_sampling_unit"
  ),
  Value = c(
    as.character(PERFORMANCE_CI_LEVEL),
    "Exact Clopper-Pearson binomial interval",
    "Nonparametric patient-level bootstrap percentile interval",
    "Nonparametric patient-level bootstrap percentile interval",
    "Nonparametric patient-level bootstrap percentile interval",
    as.character(PERFORMANCE_BOOT_N),
    as.character(SEED + PERFORMANCE_BOOT_SEED_OFFSET),
    "Patients sampled with replacement; same resample seed reused across M1-M5"
  )
)

derivation_table <- tibble(
  Predictor = c(
    "AG_ratio",
    "eGFR",
    "PAR",
    "NPR",
    "INR",
    "TyG"
  ),
  Rule = c(
    "ALB/(TP-ALB)",
    paste0(
      "141*min((CR/88.4)/kappa,1)^alpha*max((CR/88.4)/kappa,1)^(-1.209)",
      "*0.993^Age*sex_coefficient*1.0"
    ),
    "PLT/ALB",
    "NE_C/PLT",
    "(PT/13.1)^1.31",
    "log(1594.26*TG*GLU/2)"
  )
)

mice_log_table <- if (
  is.null(external_mice_logged_events) ||
  nrow(as.data.frame(external_mice_logged_events)) == 0
) {
  tibble(Status = "No external MICE logged events")
} else {
  as_tibble(external_mice_logged_events)
}

# ------------------------------------------------------------------------------
# 18. EXPORT EXCEL
# ------------------------------------------------------------------------------

cat("\n==================== EXPORT ====================\n")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "External_Performance")
openxlsx::writeData(
  wb,
  "External_Performance",
  performance_summary
)

openxlsx::addWorksheet(wb, "Frozen_Thresholds")
openxlsx::writeData(
  wb,
  "Frozen_Thresholds",
  frozen_thresholds
)

openxlsx::addWorksheet(wb, "External_DeLong")
openxlsx::writeData(
  wb,
  "External_DeLong",
  external_delong
)

openxlsx::addWorksheet(wb, "External_Predictions")
openxlsx::writeData(
  wb,
  "External_Predictions",
  prediction_table
)

openxlsx::addWorksheet(wb, "Predictions_By_MI")
openxlsx::writeData(
  wb,
  "Predictions_By_MI",
  external_predictions_by_mi
)

openxlsx::addWorksheet(wb, "Calibration_Bins")
openxlsx::writeData(
  wb,
  "Calibration_Bins",
  calibration_bins
)

openxlsx::addWorksheet(wb, "Simple_Coefficients_By_MI")
openxlsx::writeData(
  wb,
  "Simple_Coefficients_By_MI",
  simple_coefficients
)

openxlsx::addWorksheet(wb, "External_Class_Distribution")
openxlsx::writeData(
  wb,
  "External_Class_Distribution",
  external_class_distribution
)

openxlsx::addWorksheet(wb, "External_Missingness")
openxlsx::writeData(
  wb,
  "External_Missingness",
  external_missingness
)

openxlsx::addWorksheet(wb, "External_MICE_Log")
openxlsx::writeData(
  wb,
  "External_MICE_Log",
  mice_log_table
)

openxlsx::addWorksheet(wb, "Derivation_Rules")
openxlsx::writeData(
  wb,
  "Derivation_Rules",
  derivation_table
)

openxlsx::addWorksheet(wb, "Performance_CI_Settings")
openxlsx::writeData(
  wb,
  "Performance_CI_Settings",
  performance_ci_settings
)

openxlsx::addWorksheet(wb, "Settings")
openxlsx::writeData(
  wb,
  "Settings",
  settings_table
)

openxlsx::saveWorkbook(
  wb,
  file = file.path(
    OUTPUT_DIR,
    "external_validation_results.xlsx"
  ),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 19. SAVE PRIVATE R OBJECT + SESSION INFO + SUMMARY TEXT
# ------------------------------------------------------------------------------
# The RDS object contains patient-level completed data and predictions.
# Keep all generated files under private/ and do not commit them to the public repository.

external_results <- list(
  created_at = Sys.time(),
  external_data_path = EXTERNAL_DATA_PATH,
  stack_bundle_path = STACK_BUNDLE_PATH,
  baseline_results_path = BASELINE_RESULTS_PATH,
  settings = settings_table,
  performance_ci_settings = performance_ci_settings,
  synchronization = list(
    status = "PASSED",
    stack_oof_max_abs_diff = stack_oof_max_abs_diff,
    simple_model_fit_source = simple_fit_source,
    saved_simple_fits_used = use_saved_simple_fits
  ),
  derivation_rules = derivation_table,
  frozen_thresholds = frozen_thresholds,
  external_class_distribution = external_class_distribution,
  external_missingness = external_missingness,
  external_mice_logged_events = external_mice_logged_events,
  external_basic_completed_datasets = external_basic_complete_list,
  external_model_completed_datasets = external_model_complete_list,
  simple_model_fits = simple_fit_lists,
  simple_model_coefficients = simple_coefficients,
  predictions_by_imputation = external_predictions_by_mi,
  pooled_predictions = prediction_table,
  performance_summary = performance_summary,
  paired_delong = external_delong,
  calibration_bins = calibration_bins
)

saveRDS(
  external_results,
  file = file.path(
    OUTPUT_DIR,
    "external_validation_results.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "external_validation_session_info.txt"
  )
)

summary_lines <- c(
  "EXTERNAL VALIDATION - FROZEN MODELS",
  "===================================",
  paste0("Created: ", Sys.time()),
  paste0("External N: ", n_external),
  paste0(
    "Stage1 / Stage0: ",
    sum(external_truth == "Stage1"),
    " / ",
    sum(external_truth == "Stage0")
  ),
  "",
  "Validation rule:",
  "No external tuning, no threshold optimization, no stack recalibration.",
  "",
  paste0("Baseline/stack synchronization: PASSED; max Stack OOF abs diff = ",
         signif(stack_oof_max_abs_diff, 6)),
  paste0("M1-M4 final-fit source: ", simple_fit_source),
  "",
  "Frozen thresholds:",
  paste(
    capture.output(print(frozen_thresholds)),
    collapse = "\n"
  ),
  "",
  "Performance uncertainty:",
  paste0("Sensitivity/specificity CI: exact Clopper-Pearson; level = ", PERFORMANCE_CI_LEVEL),
  paste0("Brier/calibration CI: ", PERFORMANCE_BOOT_N, "-replicate nonparametric percentile bootstrap"),
  "",
  "External performance:",
  paste(
    capture.output(
      print(
        performance_summary |>
          select(
            Model,
            AUC,
            AUC_95CI_Lower,
            AUC_95CI_Upper,
            Sensitivity,
            Sensitivity_95CI_Lower,
            Sensitivity_95CI_Upper,
            Specificity,
            Specificity_95CI_Lower,
            Specificity_95CI_Upper,
            PPV,
            NPV,
            Brier,
            Brier_95CI_Lower,
            Brier_95CI_Upper,
            ECE,
            Calibration_Intercept,
            Calibration_Intercept_95CI_Lower,
            Calibration_Intercept_95CI_Upper,
            Calibration_Slope,
            Calibration_Slope_95CI_Lower,
            Calibration_Slope_95CI_Upper
          )
      )
    ),
    collapse = "\n"
  )
)

writeLines(
  summary_lines,
  con = file.path(
    OUTPUT_DIR,
    "external_validation_summary.txt"
  )
)

cat("\n============================================================\n")
cat("EXTERNAL VALIDATION COMPLETED.\n\n")
cat("Outputs:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "external_validation_results.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "external_roc_comparison.png"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "external_calibration_comparison.png"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "external_validation_results.rds"
  ),
  "\n"
)
cat("============================================================\n")
