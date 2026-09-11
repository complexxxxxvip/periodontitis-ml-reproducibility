# ==============================================================================
# Baseline Model Comparison (M1-M5)
# Full development cohort, TRUE multiple imputation, OOF evaluation
#
# Models
#   M1: Age alone                                  -> logistic regression
#   M2: Age + Sex + Smoke                          -> logistic regression
#   M3: 13 biomarkers only                         -> logistic regression
#   M4: Age + Sex + Smoke + 13 biomarkers          -> logistic regression
#   M5: Existing locked SVM + XGBoost stacking     -> read directly from bundle
#
# Key design choices
#   1) Do NOT retrain the stacking model.
#   2) Read the frozen stacking bundle and reuse its exact outer-fold assignment.
#   3) Reproduce the same fold-wise TRUE-MI procedure used by the stacking script.
#   4) Outcome is excluded from imputation.
#   5) MICE imputes ONLY 18 BASIC/source variables. The six derived predictors
#      (AG_ratio, eGFR, PAR, NPR, INR, TyG) are NEVER stochastically imputed.
#      They are deterministically recalculated after EACH completed imputation.
#   6) Logistic baselines are intentionally simple:
#        - fold-specific Z-standardization of continuous predictors only
#        - no SMOTE/Tomek
#        - no hyperparameter tuning
#        - no Platt recalibration
#        - no nonlinear terms in this first-pass baseline analysis
#   7) OOF probabilities are averaged across imputations for each patient.
#   8) Each M1-M4 threshold is explicitly derived ONCE from pooled development OOF
#      probabilities using Youden's index and then frozen for external validation.
#   9) After OOF evaluation, one final M1-M4 logistic model is fitted in EACH of the
#      final full-development completed MI datasets stored in the frozen stack bundle.
#      These fitted models and their development-only scaling parameters are saved
#      in baseline_model_results.rds for external validation.
#  10) Main comparisons:
#        M2 vs M1  : added sex + smoking beyond age
#        M3 vs M1  : biomarkers-only model vs age alone
#        M4 vs M2  : added biomarkers beyond demographics
#        M4 vs M3  : added age + sex + smoking beyond biomarkers only
#        M5 vs M4  : added ML architecture beyond same 16 predictors
#
# NOTE
#   This script assumes the frozen bundle was produced by:
#   01_development_model.R
# ==============================================================================

rm(list = ls())
gc()

# Project-relative paths. Run this script from the repository root.
PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PRIVATE_DIR <- file.path(PROJECT_DIR, "private")
MODEL_DIR <- file.path(PRIVATE_DIR, "model_objects")
OUTPUT_ROOT <- file.path(PRIVATE_DIR, "analysis_outputs")
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# The private/ directory is intended for local use only and should not be
# committed to the public repository.

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------

STACK_BUNDLE_PATH <- file.path(
  MODEL_DIR,
  "final_stacking_model_bundle.rds"
)

BASELINE_RESULTS_PATH <- file.path(
  MODEL_DIR,
  "baseline_model_results.rds"
)

OUTPUT_DIR <- file.path(
  OUTPUT_ROOT,
  "baseline_models"
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ECE bins.
ECE_BINS <- 10

# ------------------------------------------------------------------------------
# 1. PACKAGE CHECK AND LOADING
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "mice", "pROC", "openxlsx", "dplyr", "tibble", "purrr", "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them first, e.g. install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))"
    )
  )
}

suppressPackageStartupMessages({
  library(mice)
  library(pROC)
  library(openxlsx)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 2. READ THE EXISTING LOCKED STACKING BUNDLE
# ------------------------------------------------------------------------------

cat("\n==================== 1. READ STACK BUNDLE ====================\n")

if (!file.exists(STACK_BUNDLE_PATH)) {
  stop(
    paste0(
      "Cannot find STACK_BUNDLE_PATH:\n",
      STACK_BUNDLE_PATH,
      "\nPlease update STACK_BUNDLE_PATH at the top of this script."
    )
  )
}

bundle <- readRDS(STACK_BUNDLE_PATH)

# File fingerprint used to guarantee that later external validation uses the
# exact same frozen stack bundle that generated the baseline thresholds/models.
STACK_BUNDLE_MD5 <- unname(tools::md5sum(STACK_BUNDLE_PATH))

required_bundle_items <- c(
  "seed",
  "predictors",
  "binary_predictors",
  "continuous_predictors",
  "derived_predictors",
  "imputation_variables",
  "egfr_settings",
  "mice_settings",
  "raw_development_predictors",
  "raw_development_imputation_data",
  "final_completed_basic_datasets",
  "svm_oof_predictions_by_imputation",
  "meta_oof_predictions",
  "final_threshold"
)

missing_bundle_items <- setdiff(required_bundle_items, names(bundle))

if (length(missing_bundle_items) > 0) {
  stop(
    paste0(
      "The bundle is missing required elements: ",
      paste(missing_bundle_items, collapse = ", ")
    )
  )
}

SEED <- bundle$seed
MICE_M <- bundle$mice_settings$m
MICE_MAXIT <- bundle$mice_settings$maxit

# Require the eps-fixed stack bundle so M1-M5 use the same MICE tolerance.
if (is.null(bundle$mice_settings$eps)) {
  stop(
    paste0(
      "The stacking bundle does not contain mice_settings$eps. ",
      "Please rerun the final development-model script first, ",
      "then rerun this baseline comparison."
    )
  )
}
MICE_EPS <- as.numeric(bundle$mice_settings$eps)

if (length(MICE_EPS) != 1 || !is.finite(MICE_EPS) || MICE_EPS < 0) {
  stop("Invalid bundle$mice_settings$eps.")
}

MICE_METHOD <- bundle$mice_settings$method
MICE_PRINT_FLAG <- bundle$mice_settings$printFlag
CV_FOLDS <- bundle$mice_settings$cv_folds

if (length(bundle$final_completed_basic_datasets) != MICE_M) {
  stop(
    paste0(
      "The frozen stack bundle must contain exactly ", MICE_M,
      " final_completed_basic_datasets for final M1-M4 deployment fitting."
    )
  )
}

predictor_vars <- bundle$predictors
binary_vars <- bundle$binary_predictors
continuous_vars <- bundle$continuous_predictors
derived_vars <- as.character(bundle$derived_predictors)
imputation_vars <- as.character(bundle$imputation_variables)

# Locked derived-after-MICE specification.
expected_imputation_vars <- c(
  "Sex", "Smoke", "Age", "TP", "ALB", "CR", "Fbg", "PLT", "TC",
  "GGT", "ALP", "NE_C", "PT", "HB", "APTT", "LDH", "GLU", "TG"
)
expected_derived_vars <- c("AG_ratio", "eGFR", "PAR", "NPR", "INR", "TyG")

if (!identical(imputation_vars, expected_imputation_vars)) {
  stop(
    paste0(
      "The stack bundle does not use the locked 18-variable BASIC MICE input.\n",
      "Expected: ", paste(expected_imputation_vars, collapse = ", "), "\n",
      "Found: ", paste(imputation_vars, collapse = ", ")
    )
  )
}

if (!setequal(derived_vars, expected_derived_vars)) {
  stop(
    paste0(
      "The stack bundle does not contain the expected six derived predictors: ",
      paste(expected_derived_vars, collapse = ", ")
    )
  )
}

if (isTRUE(bundle$mice_settings$derived_predictors_imputed)) {
  stop("The stack bundle indicates that derived predictors were imputed. Use the new derived-after-MICE bundle.")
}

SEX_FEMALE_CODE <- as.numeric(bundle$egfr_settings$female_code)
SEX_MALE_CODE <- as.numeric(bundle$egfr_settings$male_code)
EGFR_RACE_COEFFICIENT <- as.numeric(bundle$egfr_settings$race_coefficient)

if (!identical(SEX_FEMALE_CODE, 0) || !identical(SEX_MALE_CODE, 1)) {
  stop("This analysis requires Sex=0 for female and Sex=1 for male.")
}

if (length(EGFR_RACE_COEFFICIENT) != 1 ||
    !is.finite(EGFR_RACE_COEFFICIENT) ||
    abs(EGFR_RACE_COEFFICIENT - 1.0) > 1e-12) {
  stop("This analysis requires the eGFR race coefficient to be 1.0.")
}

demographic_vars <- c("Age", "Sex", "Smoke")
biomarker_vars <- setdiff(predictor_vars, demographic_vars)

if (!all(demographic_vars %in% predictor_vars)) {
  stop("The frozen predictor set does not contain Age, Sex, and Smoke.")
}

if (length(predictor_vars) != 16) {
  warning(
    paste0(
      "Expected 16 frozen predictors, but bundle contains ",
      length(predictor_vars), "."
    )
  )
}

if (length(biomarker_vars) != 13) {
  warning(
    paste0(
      "Expected 13 biomarkers after removing Age/Sex/Smoke, but found ",
      length(biomarker_vars), "."
    )
  )
}

cat("Seed:", SEED, "\n")
cat("CV folds:", CV_FOLDS, "\n")
cat("MICE m:", MICE_M, "\n")
cat("MICE maxit:", MICE_MAXIT, "\n")
cat("MICE eps:", format(MICE_EPS, scientific = TRUE), "\n")
cat("MICE method:", MICE_METHOD, "\n")
cat("Frozen predictors:", paste(predictor_vars, collapse = ", "), "\n")
cat("Biomarkers:", paste(biomarker_vars, collapse = ", "), "\n")
cat("MICE BASIC/source variables (18):", paste(imputation_vars, collapse = ", "), "\n")
cat("Derived AFTER MICE (never imputed):", paste(derived_vars, collapse = ", "), "\n")
cat("eGFR coding: Sex 0=female, 1=male; race coefficient=", EGFR_RACE_COEFFICIENT, "\n")

# ------------------------------------------------------------------------------
# 3. RECONSTRUCT THE DEVELOPMENT DATA FROM THE BUNDLE
# ------------------------------------------------------------------------------

cat("\n==================== 2. RECONSTRUCT DEVELOPMENT DATA ====================\n")

stack_oof <- bundle$meta_oof_predictions |>
  arrange(.row)

if (!all(c(".row", "PDstage", "Stack_calibrated") %in% names(stack_oof))) {
  stop(
    "bundle$meta_oof_predictions must contain .row, PDstage, and Stack_calibrated."
  )
}

raw_predictors <- as_tibble(bundle$raw_development_predictors)
raw_imputation_df <- as_tibble(bundle$raw_development_imputation_data)

if (ncol(raw_imputation_df) != length(imputation_vars) ||
    !identical(names(raw_imputation_df), imputation_vars)) {
  stop(
    "raw_development_imputation_data is not the locked 18-column BASIC MICE input."
  )
}

if (any(derived_vars %in% names(raw_imputation_df))) {
  stop("A derived predictor is present in raw_development_imputation_data; derived variables must never enter mice().")
}

n_dev <- nrow(raw_imputation_df)

if (nrow(raw_predictors) != n_dev) {
  stop("raw_development_predictors and raw_development_imputation_data have different row counts.")
}

if (nrow(stack_oof) != n_dev) {
  stop(
    paste0(
      "Row-count mismatch: MICE source data N=", n_dev,
      ", stack OOF N=", nrow(stack_oof), "."
    )
  )
}

if (!identical(as.integer(stack_oof$.row), seq_len(n_dev))) {
  stop(
    "stack_oof$.row is not exactly 1:N. ",
    "Please inspect row alignment before proceeding."
  )
}

# At this point model_raw deliberately contains only row ID and outcome.
# Final 16 predictors are rebuilt from each completed 18-variable MICE dataset
# inside every fold/imputation below.
model_raw <- tibble(
  .row = seq_len(n_dev),
  PDstage = factor(
    as.character(stack_oof$PDstage),
    levels = c("Stage1", "Stage0")
  )
)

if (anyNA(model_raw$PDstage)) {
  stop("PDstage contains missing or unrecognized values.")
}

raw_imputation_df <- raw_imputation_df |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

# Binary coding check on the actual MICE source variables.
for (nm in c("Sex", "Smoke")) {
  vals <- sort(unique(raw_imputation_df[[nm]][!is.na(raw_imputation_df[[nm]])]))
  if (!all(vals %in% c(0, 1))) {
    stop(paste0(nm, " must be coded as 0/1."))
  }
}

class_distribution <- tibble(
  Class = c("Stage1", "Stage0"),
  N = c(
    sum(model_raw$PDstage == "Stage1"),
    sum(model_raw$PDstage == "Stage0")
  )
)

# Missingness audit now refers to the 18 BASIC/source variables entering mice().
missingness <- tibble(
  Variable = imputation_vars,
  Missing_N = vapply(
    raw_imputation_df |> select(all_of(imputation_vars)),
    function(x) sum(is.na(x)),
    integer(1)
  )
) |>
  mutate(Missing_Rate = Missing_N / n_dev)

cat("Development N:", n_dev, "\n")
print(class_distribution)
cat("Missingness in the 18 BASIC/source MICE variables:\n")
print(missingness |> filter(Missing_N > 0))

# ------------------------------------------------------------------------------
# 4. RECOVER THE EXACT OUTER-FOLD ASSIGNMENT USED BY THE STACKING MODEL
# ------------------------------------------------------------------------------

cat("\n==================== 3. RECOVER EXACT STACK FOLDS ====================\n")

svm_oof_by_mi <- bundle$svm_oof_predictions_by_imputation

if (!all(c(".row", "id") %in% names(svm_oof_by_mi))) {
  stop(
    "bundle$svm_oof_predictions_by_imputation must contain .row and id."
  )
}

parse_fold_number <- function(x) {
  out <- suppressWarnings(
    as.integer(
      sub("^Fold([0-9]+)_Imp[0-9]+$", "\\1", as.character(x))
    )
  )
  out
}

fold_long <- svm_oof_by_mi |>
  transmute(
    .row = as.integer(.row),
    id = as.character(id),
    Fold = parse_fold_number(id)
  )

if (anyNA(fold_long$Fold)) {
  bad_ids <- unique(fold_long$id[is.na(fold_long$Fold)])
  stop(
    paste0(
      "Could not parse fold number from these IDs: ",
      paste(bad_ids, collapse = ", ")
    )
  )
}

fold_map <- fold_long |>
  group_by(.row) |>
  summarise(
    Fold = first(Fold),
    N_unique_folds = n_distinct(Fold),
    MI_N = n(),
    .groups = "drop"
  ) |>
  arrange(.row)

if (nrow(fold_map) != n_dev) {
  stop("Could not recover one fold assignment for every development patient.")
}

if (any(fold_map$N_unique_folds != 1)) {
  stop("At least one patient appears in more than one outer assessment fold.")
}

if (any(fold_map$MI_N != MICE_M)) {
  stop(
    paste0(
      "At least one patient does not have exactly ",
      MICE_M, " imputation-specific stack OOF predictions."
    )
  )
}

if (length(unique(fold_map$Fold)) != CV_FOLDS) {
  stop(
    paste0(
      "Recovered ", length(unique(fold_map$Fold)),
      " folds, but bundle reports ", CV_FOLDS, "."
    )
  )
}

fold_sizes <- fold_map |>
  count(Fold, name = "Assessment_N") |>
  arrange(Fold)

print(fold_sizes)

# ------------------------------------------------------------------------------
# 5A. HELPER: DERIVE THE SIX NON-IMPUTED MODEL PREDICTORS AFTER MICE
#
# These six variables NEVER enter mice():
#   AG_ratio = ALB / (TP - ALB)
#   eGFR     = 141 * min(Scr/kappa,1)^alpha * max(Scr/kappa,1)^(-1.209) *
#              0.993^Age * sex_coefficient * race_coefficient,
#              where Scr = CR/88.4 mg/dL;
#              female Sex=0: kappa=0.7, alpha=-0.329, coefficient=1.018;
#              male   Sex=1: kappa=0.9, alpha=-0.411, coefficient=1.0;
#              race coefficient = 1.0.
#   PAR      = PLT / ALB
#   NPR      = NE_C / PLT
#   INR      = (PT / 13.1)^1.31
#   TyG      = ln(1594.26 * TG * GLU / 2)
# ------------------------------------------------------------------------------

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
    stop("Sex contains values outside 0=female and 1=male.")
  }

  n <- nrow(out)
  scr_mg_dl <- out$CR / 88.4
  is_female <- out$Sex == female_code
  is_male <- out$Sex == male_code

  kappa <- ifelse(is_female, 0.7, ifelse(is_male, 0.9, NA_real_))
  alpha <- ifelse(is_female, -0.329, ifelse(is_male, -0.411, NA_real_))
  sex_coefficient <- ifelse(is_female, 1.018, ifelse(is_male, 1.0, NA_real_))

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
        ". Check source values (especially TP>ALB, CR>0, ALB>0, PLT>0, ",
        "PT>0, TG>0, GLU>0) and Sex coding."
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

# ------------------------------------------------------------------------------
# 5. HELPER: TRUE MICE-PMM IMPUTATION
#
# This reproduces the same basic MICE logic as the stacking script:
#   - ONLY the 18 BASIC/source variables are supplied to mice()
#   - the six derived predictors are absent from mice() and are calculated afterward
#   - outcome is NOT supplied to mice()
#   - assessment rows have ignore=TRUE
#   - assessment rows are still imputed, but cannot estimate imputation parameters
#   - every completed imputation is retained
# ------------------------------------------------------------------------------

mice_impute_predictors <- function(
    predictor_df,
    ignore_vec,
    seed,
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

# ------------------------------------------------------------------------------
# 6. HELPER: FIT A SIMPLE LOGISTIC MODEL AND PREDICT SEVERE PERIODONTITIS
#
# IMPORTANT:
# PDstage is stored with Stage1 as the first factor level for tidymodels.
# Base R glm() treats factor levels differently, so we explicitly construct:
#   y = 1 for Stage1 (severe)
#   y = 0 for Stage0 (mild)
# ------------------------------------------------------------------------------

fit_predict_logistic <- function(
    train_df,
    test_df,
    predictor_names,
    model_label,
    continuous_names = continuous_vars
) {

  # Fold-specific Z-standardization:
  #   - estimate mean/SD from analysis rows only;
  #   - apply those parameters unchanged to assessment rows;
  #   - keep binary Sex/Smoke as 0/1.
  #
  # For ordinary unpenalized logistic regression, linear rescaling should not
  # materially change fitted probabilities, but this keeps preprocessing
  # consistent with the stacking pipeline and improves numerical stability.
  dat_train <- train_df
  dat_test <- test_df

  vars_to_scale <- intersect(predictor_names, continuous_names)

  if (length(vars_to_scale) > 0) {
    for (nm in vars_to_scale) {

      mu <- mean(dat_train[[nm]])
      sigma <- stats::sd(dat_train[[nm]])

      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
        stop(
          paste0(
            model_label,
            ": cannot standardize predictor ", nm,
            " because the analysis-fold SD is non-positive/non-finite."
          )
        )
      }

      dat_train[[nm]] <- (dat_train[[nm]] - mu) / sigma
      dat_test[[nm]] <- (dat_test[[nm]] - mu) / sigma
    }
  }

  dat_train <- dat_train |>
    mutate(y = ifelse(PDstage == "Stage1", 1, 0))

  formula_obj <- reformulate(
    termlabels = predictor_names,
    response = "y"
  )

  fit <- suppressWarnings(
    glm(
      formula = formula_obj,
      data = dat_train,
      family = binomial(link = "logit")
    )
  )

  if (!isTRUE(fit$converged)) {
    warning(paste0(model_label, ": glm did not converge."))
  }

  prob <- suppressWarnings(
    predict(
      fit,
      newdata = dat_test,
      type = "response"
    )
  )

  if (length(prob) != nrow(dat_test)) {
    stop(paste0(model_label, ": prediction length mismatch."))
  }

  if (any(!is.finite(prob))) {
    stop(paste0(model_label, ": non-finite predicted probabilities detected."))
  }

  prob <- pmin(pmax(as.numeric(prob), 1e-12), 1 - 1e-12)

  prob
}


# ------------------------------------------------------------------------------
# 6B. HELPER: FIT / PREDICT A FROZEN FULL-DEVELOPMENT SIMPLE MODEL
#
# Unlike fit_predict_logistic(), this helper returns the fitted glm plus the
# development-only scaling parameters so the exact model can be deployed later.
# No external data are used here.
# ------------------------------------------------------------------------------

fit_frozen_simple_logistic <- function(
    train_df,
    predictor_names,
    model_label,
    continuous_names = continuous_vars
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
            model_label,
            ": cannot estimate final development scaler for ", nm,
            " because mean/SD is non-finite or SD <= 0."
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
    warning(paste0(model_label, ": final full-development glm did not converge."))
  }

  list(
    fit = fit,
    scaler = scaler,
    predictors = predictor_names,
    model_label = model_label
  )
}


predict_frozen_simple_logistic <- function(fit_obj, new_df) {

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

  if (length(prob) != nrow(dat) || any(!is.finite(prob))) {
    stop(
      paste0(
        fit_obj$model_label,
        ": invalid prediction from frozen simple model."
      )
    )
  }

  pmin(pmax(as.numeric(prob), 1e-12), 1 - 1e-12)
}


# ------------------------------------------------------------------------------
# 7. FOLD-WISE TRUE-MI -> DERIVE PREDICTORS -> OOF PREDICTIONS FOR FOUR LOGISTIC BASELINES
# ------------------------------------------------------------------------------

cat("\n==================== 4. FIT SIMPLE BASELINE MODELS ====================\n")

prediction_list <- list()
imputation_log_list <- list()
counter <- 0L

all_rows <- seq_len(n_dev)
fold_values <- sort(unique(fold_map$Fold))

for (i in fold_values) {

  cat(
    sprintf(
      "\nFold %d/%d: TRUE-MI on 18 BASIC variables + derive 6 predictors + 4 logistic baselines ...\n",
      i, CV_FOLDS
    )
  )

  assessment_idx <- fold_map$.row[fold_map$Fold == i]
  analysis_idx <- setdiff(all_rows, assessment_idx)

  # Exact same leakage-control logic:
  # analysis rows build the imputation model;
  # assessment rows are imputed but cannot affect imputation-model parameters.
  ignore_vec <- rep(TRUE, n_dev)
  ignore_vec[analysis_idx] <- FALSE

  imp_i <- mice_impute_predictors(
    predictor_df = raw_imputation_df,
    ignore_vec = ignore_vec,
    seed = SEED + i,
    m = MICE_M,
    maxit = MICE_MAXIT,
    method_name = MICE_METHOD,
    print_flag = MICE_PRINT_FLAG,
    eps = MICE_EPS
  )

  imputation_log_list[[i]] <- list(
    Fold = i,
    Analysis_N = length(analysis_idx),
    Assessment_N = length(assessment_idx),
    MICE_method = imp_i$method,
    MICE_predictorMatrix = imp_i$predictorMatrix,
    MICE_eps = MICE_EPS,
    MICE_loggedEvents = imp_i$loggedEvents
  )

  for (j in seq_len(MICE_M)) {

    completed_basic <- imp_i$completed_list[[j]]

    # IMPORTANT: all six derived predictors are rebuilt from the completed BASIC
    # variables. No observed/precalculated AG_ratio/eGFR/PAR/NPR/INR/TyG value is
    # carried forward into modeling.
    completed_predictors <- derive_model_predictors(
      completed_basic,
      strict = TRUE
    )

    fold_data <- bind_cols(
      model_raw |>
        select(.row, PDstage),
      completed_predictors
    )

    if (anyNA(fold_data)) {
      stop(
        sprintf(
          "Missing values remain in Fold %d, Imputation %d.",
          i, j
        )
      )
    }

    train_df <- fold_data[analysis_idx, , drop = FALSE]
    test_df <- fold_data[assessment_idx, , drop = FALSE]

    p_age <- fit_predict_logistic(
      train_df = train_df,
      test_df = test_df,
      predictor_names = c("Age"),
      model_label = "M1 Age alone"
    )

    p_demo <- fit_predict_logistic(
      train_df = train_df,
      test_df = test_df,
      predictor_names = demographic_vars,
      model_label = "M2 Age + Sex + Smoke"
    )

    p_biomarkers_only <- fit_predict_logistic(
      train_df = train_df,
      test_df = test_df,
      predictor_names = biomarker_vars,
      model_label = "M3 13 biomarkers only"
    )

    p_full_logit <- fit_predict_logistic(
      train_df = train_df,
      test_df = test_df,
      predictor_names = predictor_vars,
      model_label = "M4 Demographics + biomarkers"
    )

    counter <- counter + 1L

    prediction_list[[counter]] <- tibble(
      .row = test_df$.row,
      PDstage = test_df$PDstage,
      Fold = i,
      Imputation = j,
      Age_Only = p_age,
      Demographics = p_demo,
      Biomarkers_Only = p_biomarkers_only,
      Demo_Biomarkers_Logistic = p_full_logit
    )
  }
}

baseline_oof_by_mi <- bind_rows(prediction_list) |>
  arrange(.row, Imputation)

required_baseline_prob_cols <- c(
  "Age_Only",
  "Demographics",
  "Biomarkers_Only",
  "Demo_Biomarkers_Logistic"
)

if (!all(required_baseline_prob_cols %in% names(baseline_oof_by_mi))) {
  stop("One or more baseline OOF probability columns are missing.")
}

expected_rows <- n_dev * MICE_M

if (nrow(baseline_oof_by_mi) != expected_rows) {
  stop(
    paste0(
      "Expected ", expected_rows,
      " imputation-specific OOF rows, but obtained ",
      nrow(baseline_oof_by_mi), "."
    )
  )
}

# ------------------------------------------------------------------------------
# 7B. MICE DIAGNOSTIC AUDIT
# ------------------------------------------------------------------------------

mice_log_summary <- purrr::map_dfr(
  seq_along(imputation_log_list),
  function(i) {
    x <- imputation_log_list[[i]]$MICE_loggedEvents
    tibble(
      Fold = i,
      Analysis_N = imputation_log_list[[i]]$Analysis_N,
      Assessment_N = imputation_log_list[[i]]$Assessment_N,
      MICE_eps = MICE_EPS,
      Logged_events = if (is.null(x)) 0L else nrow(x)
    )
  }
)

mice_logged_events_all <- purrr::map_dfr(
  seq_along(imputation_log_list),
  function(i) {
    x <- imputation_log_list[[i]]$MICE_loggedEvents
    if (is.null(x) || nrow(x) == 0) return(tibble())
    as_tibble(x) |> mutate(Fold = i, .before = 1)
  }
)

print(mice_log_summary)

if (nrow(mice_logged_events_all) > 0) {
  warning(
    paste0(
      "MICE logged events remain after MICE_EPS=", MICE_EPS,
      ". Inspect MICE_LoggedEvents in the output workbook."
    )
  )
} else {
  cat("No MICE logged events after eps fix.\n")
}

# ------------------------------------------------------------------------------
# 8. POOL OOF PROBABILITIES ACROSS IMPUTATIONS
# ------------------------------------------------------------------------------

cat("\n==================== 5. POOL OOF PROBABILITIES ====================\n")

baseline_oof <- baseline_oof_by_mi |>
  group_by(.row, PDstage) |>
  summarise(
    MI_N = n(),

    Age_Only = mean(Age_Only),
    Age_Only_MI_SD = ifelse(n() > 1, sd(Age_Only), 0),

    Demographics = mean(Demographics),
    Demographics_MI_SD = ifelse(n() > 1, sd(Demographics), 0),

    Biomarkers_Only = mean(Biomarkers_Only),
    Biomarkers_Only_MI_SD =
      ifelse(n() > 1, sd(Biomarkers_Only), 0),

    Demo_Biomarkers_Logistic = mean(Demo_Biomarkers_Logistic),
    Demo_Biomarkers_Logistic_MI_SD =
      ifelse(n() > 1, sd(Demo_Biomarkers_Logistic), 0),

    .groups = "drop"
  ) |>
  arrange(.row)

if (any(baseline_oof$MI_N != MICE_M)) {
  stop("At least one patient does not have the expected number of MI predictions.")
}

if (!identical(as.integer(baseline_oof$.row), seq_len(n_dev))) {
  stop("Baseline pooled OOF rows are not aligned to 1:N.")
}

# ------------------------------------------------------------------------------
# 9. ADD THE EXISTING STACKING OOF PROBABILITY
# ------------------------------------------------------------------------------

stack_for_join <- stack_oof |>
  transmute(
    .row = as.integer(.row),
    PDstage_stack = factor(
      as.character(PDstage),
      levels = c("Stage1", "Stage0")
    ),
    Stack = as.numeric(Stack_calibrated)
  )

comparison_oof <- baseline_oof |>
  left_join(
    stack_for_join,
    by = ".row"
  )

if (anyNA(comparison_oof$Stack)) {
  stop("Stack OOF probabilities could not be aligned to the baseline OOF table.")
}

if (!all(
  as.character(comparison_oof$PDstage) ==
    as.character(comparison_oof$PDstage_stack)
)) {
  stop("Outcome mismatch between baseline OOF and stack OOF after row alignment.")
}

comparison_oof <- comparison_oof |>
  select(-PDstage_stack)

# ------------------------------------------------------------------------------
# 10. PERFORMANCE FUNCTIONS
# ------------------------------------------------------------------------------

make_roc <- function(truth, prob) {

  truth_num <- ifelse(as.character(truth) == "Stage1", 1, 0)

  pROC::roc(
    response = truth_num,
    predictor = prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

safe_div <- function(num, den) {
  ifelse(den == 0, NA_real_, num / den)
}

get_youden_threshold <- function(truth, prob) {

  roc_obj <- make_roc(truth, prob)

  coord <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c(
      "threshold",
      "sensitivity",
      "specificity"
    ),
    transpose = FALSE
  )

  as.numeric(coord$threshold[1])
}

calculate_performance <- function(
    truth,
    prob,
    model_name,
    threshold = NULL,
    threshold_source = "OOF Youden",
    ece_bins = ECE_BINS
) {

  truth_chr <- as.character(truth)
  y <- ifelse(truth_chr == "Stage1", 1, 0)

  roc_obj <- make_roc(truth, prob)

  auc_value <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(
    pROC::ci.auc(
      roc_obj,
      method = "delong"
    )
  )

  if (is.null(threshold) || !is.finite(threshold)) {
    threshold <- get_youden_threshold(truth, prob)
    threshold_source <- "OOF Youden"
  }

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

  brier <- mean((prob - y)^2)

  # Calibration intercept and slope.
  # Protect against logit(0) / logit(1).
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

  cal_intercept <- unname(coef(cal_intercept_fit)[1])
  cal_slope <- unname(coef(cal_slope_fit)["lp"])

  # 10-bin ECE by default.
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

    AUC = auc_value,
    AUC_95CI_Lower = auc_ci[1],
    AUC_95CI_Upper = auc_ci[3],

    Threshold = threshold,
    Threshold_Source = threshold_source,

    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    PPV = ppv,
    NPV = npv,
    F1 = f1,
    G_mean = gmean,

    Brier = brier,
    ECE = ece,

    Calibration_Intercept = cal_intercept,
    Calibration_Slope = cal_slope,

    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  )
}

# ------------------------------------------------------------------------------
# 11. EXPLICITLY FREEZE DEVELOPMENT OOF THRESHOLDS + PERFORMANCE
# ------------------------------------------------------------------------------

cat("\n==================== 6. FREEZE OOF THRESHOLDS ====================\n")

# M1-M4 thresholds are derived ONLY from each model's pooled development OOF
# probabilities. They are calculated once here and then passed explicitly to all
# downstream performance calculations and saved for external validation.
simple_frozen_thresholds <- tibble(
  Model = c(
    "M1 Age alone",
    "M2 Age + Sex + Smoking",
    "M3 13 biomarkers only (logistic)",
    "M4 Demographics + 13 biomarkers (logistic)"
  ),
  Probability_Column = c(
    "Age_Only",
    "Demographics",
    "Biomarkers_Only",
    "Demo_Biomarkers_Logistic"
  ),
  Threshold = c(
    get_youden_threshold(comparison_oof$PDstage, comparison_oof$Age_Only),
    get_youden_threshold(comparison_oof$PDstage, comparison_oof$Demographics),
    get_youden_threshold(comparison_oof$PDstage, comparison_oof$Biomarkers_Only),
    get_youden_threshold(
      comparison_oof$PDstage,
      comparison_oof$Demo_Biomarkers_Logistic
    )
  ),
  Threshold_Source = "Frozen development pooled OOF Youden"
)

if (any(!is.finite(simple_frozen_thresholds$Threshold))) {
  stop("At least one M1-M4 frozen OOF threshold is non-finite.")
}

# The stacking threshold MUST already be frozen in the stack bundle.
# Do not silently recalculate it in this baseline script.
stack_threshold <- suppressWarnings(as.numeric(bundle$final_threshold))

if (length(stack_threshold) != 1 || !is.finite(stack_threshold)) {
  stop(
    paste0(
      "The frozen stack bundle does not contain one valid final_threshold. ",
      "Rerun/fix the final stacking model before running this baseline script."
    )
  )
}

frozen_thresholds <- bind_rows(
  simple_frozen_thresholds |>
    select(Model, Threshold, Threshold_Source),
  tibble(
    Model = "M5 SVM + XGBoost stacking",
    Threshold = stack_threshold,
    Threshold_Source = "Frozen stacking development OOF threshold"
  )
)

print(frozen_thresholds)

threshold_lookup <- setNames(
  frozen_thresholds$Threshold,
  frozen_thresholds$Model
)

performance_summary <- bind_rows(

  calculate_performance(
    truth = comparison_oof$PDstage,
    prob = comparison_oof$Age_Only,
    model_name = "M1 Age alone",
    threshold = threshold_lookup[["M1 Age alone"]],
    threshold_source = "Frozen development pooled OOF Youden"
  ),

  calculate_performance(
    truth = comparison_oof$PDstage,
    prob = comparison_oof$Demographics,
    model_name = "M2 Age + Sex + Smoking",
    threshold = threshold_lookup[["M2 Age + Sex + Smoking"]],
    threshold_source = "Frozen development pooled OOF Youden"
  ),

  calculate_performance(
    truth = comparison_oof$PDstage,
    prob = comparison_oof$Biomarkers_Only,
    model_name = "M3 13 biomarkers only (logistic)",
    threshold = threshold_lookup[["M3 13 biomarkers only (logistic)"]],
    threshold_source = "Frozen development pooled OOF Youden"
  ),

  calculate_performance(
    truth = comparison_oof$PDstage,
    prob = comparison_oof$Demo_Biomarkers_Logistic,
    model_name = "M4 Demographics + 13 biomarkers (logistic)",
    threshold = threshold_lookup[["M4 Demographics + 13 biomarkers (logistic)"]],
    threshold_source = "Frozen development pooled OOF Youden"
  ),

  calculate_performance(
    truth = comparison_oof$PDstage,
    prob = comparison_oof$Stack,
    model_name = "M5 SVM + XGBoost stacking",
    threshold = stack_threshold,
    threshold_source = "Frozen stacking development OOF threshold"
  )
)

print(
  performance_summary |>
    select(
      Model,
      AUC,
      AUC_95CI_Lower,
      AUC_95CI_Upper,
      Threshold,
      Threshold_Source,
      Brier,
      ECE,
      Calibration_Intercept,
      Calibration_Slope,
      Sensitivity,
      Specificity
    )
)

# ------------------------------------------------------------------------------
# 11B. FIT AND FREEZE FINAL M1-M4 MODELS ON FULL DEVELOPMENT MI DATASETS
#
# These models are deployment objects. The external cohort must NEVER be used to
# refit them. One model is fitted per final completed development imputation.
# ------------------------------------------------------------------------------

cat("\n==================== 6B. FREEZE FINAL M1-M4 MODELS ====================\n")

simple_model_specs <- list(
  "M1 Age alone" = c("Age"),
  "M2 Age + Sex + Smoking" = demographic_vars,
  "M3 13 biomarkers only (logistic)" = biomarker_vars,
  "M4 Demographics + 13 biomarkers (logistic)" = predictor_vars
)

final_simple_model_fits <- setNames(
  lapply(simple_model_specs, function(x) vector("list", MICE_M)),
  names(simple_model_specs)
)

final_simple_coef_rows <- list()
final_simple_scaler_rows <- list()
coef_counter <- 0L
scaler_counter <- 0L

for (j in seq_len(MICE_M)) {

  completed_basic_j <- as_tibble(
    bundle$final_completed_basic_datasets[[j]]
  )

  if (nrow(completed_basic_j) != n_dev) {
    stop(
      sprintf(
        "Final completed BASIC dataset %d has %d rows; expected %d.",
        j, nrow(completed_basic_j), n_dev
      )
    )
  }

  completed_predictors_j <- derive_model_predictors(
    completed_basic_j,
    strict = TRUE
  )

  final_dev_j <- bind_cols(
    model_raw |> select(PDstage),
    completed_predictors_j
  )

  if (anyNA(final_dev_j)) {
    stop(
      sprintf(
        "Missing values remain in final simple-model development dataset %d.",
        j
      )
    )
  }

  for (model_nm in names(simple_model_specs)) {

    fit_obj <- fit_frozen_simple_logistic(
      train_df = final_dev_j,
      predictor_names = simple_model_specs[[model_nm]],
      model_label = model_nm,
      continuous_names = continuous_vars
    )

    final_simple_model_fits[[model_nm]][[j]] <- fit_obj

    coef_counter <- coef_counter + 1L
    final_simple_coef_rows[[coef_counter]] <- tibble(
      Model = model_nm,
      Imputation = j,
      Term = names(coef(fit_obj$fit)),
      Estimate = as.numeric(coef(fit_obj$fit))
    )

    if (length(fit_obj$scaler) > 0) {
      for (nm in names(fit_obj$scaler)) {
        scaler_counter <- scaler_counter + 1L
        final_simple_scaler_rows[[scaler_counter]] <- tibble(
          Model = model_nm,
          Imputation = j,
          Variable = nm,
          Mean = unname(fit_obj$scaler[[nm]][["mean"]]),
          SD = unname(fit_obj$scaler[[nm]][["sd"]])
        )
      }
    }
  }
}

names_by_imp <- paste0("Imp", seq_len(MICE_M))
for (model_nm in names(final_simple_model_fits)) {
  names(final_simple_model_fits[[model_nm]]) <- names_by_imp
}

final_simple_coefficients <- bind_rows(final_simple_coef_rows)
final_simple_scalers <- bind_rows(final_simple_scaler_rows)

cat(
  "Frozen simple models:",
  length(simple_model_specs), "models x", MICE_M,
  "development imputations =", length(simple_model_specs) * MICE_M,
  "fitted glm objects.\n"
)

# ------------------------------------------------------------------------------
# 12. PRESPECIFIED PAIRED DELONG COMPARISONS
#
# pROC paired DeLong CI is for AUC(new) - AUC(old) because roc_new is passed first.
# ------------------------------------------------------------------------------

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

auc_incremental_comparisons <- bind_rows(

  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Demographics,
    prob_old = comparison_oof$Age_Only,
    new_name = "M2 Age + Sex + Smoking",
    old_name = "M1 Age alone"
  ),

  # Biomarkers-only versus age alone.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Biomarkers_Only,
    prob_old = comparison_oof$Age_Only,
    new_name = "M3 13 biomarkers only (logistic)",
    old_name = "M1 Age alone"
  ),

  # Primary comparison: biomarkers beyond demographics.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Demo_Biomarkers_Logistic,
    prob_old = comparison_oof$Demographics,
    new_name = "M4 Demographics + 13 biomarkers (logistic)",
    old_name = "M2 Age + Sex + Smoking"
  ),

  # Added demographics beyond biomarkers-only.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Demo_Biomarkers_Logistic,
    prob_old = comparison_oof$Biomarkers_Only,
    new_name = "M4 Demographics + 13 biomarkers (logistic)",
    old_name = "M3 13 biomarkers only (logistic)"
  ),

  # ML architecture beyond the same 16 predictors.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Stack,
    prob_old = comparison_oof$Demo_Biomarkers_Logistic,
    new_name = "M5 SVM + XGBoost stacking",
    old_name = "M4 Demographics + 13 biomarkers (logistic)"
  ),

  # Direct stack versus biomarkers-only comparison.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Stack,
    prob_old = comparison_oof$Biomarkers_Only,
    new_name = "M5 SVM + XGBoost stacking",
    old_name = "M3 13 biomarkers only (logistic)"
  ),

  # Useful global comparisons.
  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Demo_Biomarkers_Logistic,
    prob_old = comparison_oof$Age_Only,
    new_name = "M4 Demographics + 13 biomarkers (logistic)",
    old_name = "M1 Age alone"
  ),

  paired_delong_compare(
    truth = comparison_oof$PDstage,
    prob_new = comparison_oof$Stack,
    prob_old = comparison_oof$Age_Only,
    new_name = "M5 SVM + XGBoost stacking",
    old_name = "M1 Age alone"
  )
)

cat("\nPaired DeLong comparisons:\n")
print(auc_incremental_comparisons)

# ------------------------------------------------------------------------------
# 13. CREATE ROC PLOT
# ------------------------------------------------------------------------------

cat("\n==================== 7. ROC FIGURE ====================\n")

model_prob_map <- list(
  "Age alone" = comparison_oof$Age_Only,
  "Age + Sex + Smoking" = comparison_oof$Demographics,
  "13 biomarkers only (logistic)" =
    comparison_oof$Biomarkers_Only,
  "Demographics + 13 biomarkers (logistic)" =
    comparison_oof$Demo_Biomarkers_Logistic,
  "SVM + XGBoost stacking" = comparison_oof$Stack
)

roc_plot_data <- bind_rows(
  lapply(
    names(model_prob_map),
    function(model_nm) {

      roc_obj <- make_roc(
        comparison_oof$PDstage,
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
    Model = case_when(
      Model == "M1 Age alone" ~ "Age alone",
      Model == "M2 Age + Sex + Smoking" ~ "Age + Sex + Smoking",
      Model == "M3 13 biomarkers only (logistic)" ~
        "13 biomarkers only (logistic)",
      Model == "M4 Demographics + 13 biomarkers (logistic)" ~
        "Demographics + 13 biomarkers (logistic)",
      Model == "M5 SVM + XGBoost stacking" ~
        "SVM + XGBoost stacking",
      TRUE ~ Model
    ),
    AUC_Label = sprintf(
      "%s (AUC %.3f, 95%% CI %.3f-%.3f)",
      Model,
      AUC,
      AUC_95CI_Lower,
      AUC_95CI_Upper
    )
  )

label_lookup <- setNames(
  auc_label_table$AUC_Label,
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
    title = "ROC Comparison of Baseline and Stacking Models",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    linetype = "Model"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  ) +
  guides(
    linetype = guide_legend(ncol = 1, byrow = TRUE)
  )

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "ROC_Baseline_vs_Stacking.png"
  ),
  plot = roc_plot,
  width = 10,
  height = 8,
  dpi = 320
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "ROC_Baseline_vs_Stacking.pdf"
  ),
  plot = roc_plot,
  width = 10,
  height = 8
)

# ------------------------------------------------------------------------------
# 14. KEY COMPARISON TABLE
# ------------------------------------------------------------------------------

key_comparison_table <- performance_summary |>
  select(
    Model,
    AUC,
    AUC_95CI_Lower,
    AUC_95CI_Upper,
    Brier,
    ECE,
    Calibration_Intercept,
    Calibration_Slope,
    Sensitivity,
    Specificity,
    Accuracy,
    PPV,
    NPV,
    F1
  ) |>
  mutate(
    AUC_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      AUC,
      AUC_95CI_Lower,
      AUC_95CI_Upper
    )
  ) |>
  select(
    Model,
    AUC_95CI,
    everything(),
    -AUC_95CI_Lower,
    -AUC_95CI_Upper
  )

# Add prespecified AUC differences using clinically meaningful reference models.
# M2 and M3 are parallel simple-model branches, so "vs previous row" would be misleading.
key_comparison_table <- key_comparison_table |>
  mutate(
    Reference_Model = c(
      NA_character_,
      "M1 Age alone",
      "M1 Age alone",
      "M2 Age + Sex + Smoking",
      "M4 Demographics + 13 biomarkers (logistic)"
    ),
    Delta_AUC_vs_Reference = c(
      NA_real_,
      AUC[2] - AUC[1],
      AUC[3] - AUC[1],
      AUC[4] - AUC[2],
      AUC[5] - AUC[4]
    )
  )

# ------------------------------------------------------------------------------
# 15. EXPORT EXCEL WORKBOOK
# ------------------------------------------------------------------------------

cat("\n==================== 8. EXPORT RESULTS ====================\n")

settings_table <- tibble(
  Parameter = c(
    "Stack_bundle",
    "Development_N",
    "SEED",
    "CV_FOLDS",
    "MICE_M",
    "MICE_MAXIT",
    "MICE_EPS",
    "MICE_METHOD",
    "Outcome_used_for_imputation",
    "MICE_input_variables",
    "Derived_predictors_imputed",
    "Derived_predictors",
    "Derived_formula_rules",
    "eGFR_sex_race_rule",
    "Baseline_algorithm",
    "Baseline_standardization",
    "Baseline_standardization_rule",
    "Baseline_SMOTE_Tomek",
    "Baseline_recalibration",
    "MI_pooling_rule",
    "Threshold_rule_M1_M4",
    "Final_simple_model_fit_rule",
    "Stack_bundle_MD5"
  ),
  Value = c(
    STACK_BUNDLE_PATH,
    n_dev,
    SEED,
    CV_FOLDS,
    MICE_M,
    MICE_MAXIT,
    MICE_EPS,
    MICE_METHOD,
    FALSE,
    paste(imputation_vars, collapse = ", "),
    FALSE,
    paste(derived_vars, collapse = ", "),
    paste(
      c(
        "AG_ratio=ALB/(TP-ALB)",
        "INR=(PT/13.1)^1.31",
        "NPR=NE_C/PLT",
        "PAR=PLT/ALB",
        "TyG=ln(1594.26*TG*GLU/2)"
      ),
      collapse = "; "
    ),
    "eGFR: CR/88.4 to mg/dL; Sex=0 female (kappa=0.7, alpha=-0.329, x1.018); Sex=1 male (kappa=0.9, alpha=-0.411); race coefficient=1.0",
    "Standard unpenalized logistic regression",
    TRUE,
    "Z-score continuous predictors within each analysis fold/imputation; apply analysis mean/SD to assessment rows only",
    FALSE,
    FALSE,
    "Mean predicted probability across imputations per patient",
    "One Youden threshold per M1-M4 model from pooled development OOF probabilities; frozen before external validation",
    "One unpenalized logistic model per M1-M4 per final full-development completed imputation; scaler estimated only from that development imputation",
    STACK_BUNDLE_MD5
  )
)

predictor_set_table <- bind_rows(
  tibble(
    Model = "M1 Age alone",
    Predictor = "Age"
  ),
  tibble(
    Model = "M2 Age + Sex + Smoking",
    Predictor = demographic_vars
  ),
  tibble(
    Model = "M3 13 biomarkers only (logistic)",
    Predictor = biomarker_vars
  ),
  tibble(
    Model = "M4 Demographics + 13 biomarkers (logistic)",
    Predictor = predictor_vars
  ),
  tibble(
    Model = "M5 SVM + XGBoost stacking",
    Predictor = predictor_vars
  )
)

mice_input_table <- tibble(
  Order = seq_along(imputation_vars),
  Variable = imputation_vars,
  Role = "BASIC/source variable entered into mice()"
)

derivation_rules_table <- tibble(
  Derived_Predictor = c("AG_ratio", "eGFR", "PAR", "NPR", "INR", "TyG"),
  Imputed_by_MICE = FALSE,
  Rule = c(
    "ALB / (TP - ALB)",
    "141 * min((CR/88.4)/kappa,1)^alpha * max((CR/88.4)/kappa,1)^(-1.209) * 0.993^Age * sex_coefficient * 1.0",
    "PLT / ALB",
    "NE_C / PLT",
    "(PT / 13.1)^1.31",
    "ln(1594.26 * TG * GLU / 2)"
  )
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Key_Comparisons")
openxlsx::writeData(
  wb,
  "Key_Comparisons",
  key_comparison_table
)

openxlsx::addWorksheet(wb, "Performance_Full")
openxlsx::writeData(
  wb,
  "Performance_Full",
  performance_summary
)

openxlsx::addWorksheet(wb, "Incremental_AUC_DeLong")
openxlsx::writeData(
  wb,
  "Incremental_AUC_DeLong",
  auc_incremental_comparisons
)

openxlsx::addWorksheet(wb, "OOF_Pooled")
openxlsx::writeData(
  wb,
  "OOF_Pooled",
  comparison_oof
)

openxlsx::addWorksheet(wb, "OOF_By_Imputation")
openxlsx::writeData(
  wb,
  "OOF_By_Imputation",
  baseline_oof_by_mi
)

openxlsx::addWorksheet(wb, "Fold_Map")
openxlsx::writeData(
  wb,
  "Fold_Map",
  fold_map
)

openxlsx::addWorksheet(wb, "Fold_Sizes")
openxlsx::writeData(
  wb,
  "Fold_Sizes",
  fold_sizes
)

openxlsx::addWorksheet(wb, "Predictor_Sets")
openxlsx::writeData(
  wb,
  "Predictor_Sets",
  predictor_set_table
)

openxlsx::addWorksheet(wb, "MICE_Input_Variables")
openxlsx::writeData(
  wb,
  "MICE_Input_Variables",
  mice_input_table
)

openxlsx::addWorksheet(wb, "Derivation_Rules")
openxlsx::writeData(
  wb,
  "Derivation_Rules",
  derivation_rules_table
)

openxlsx::addWorksheet(wb, "MICE_Input_Missingness")
openxlsx::writeData(
  wb,
  "MICE_Input_Missingness",
  missingness
)

openxlsx::addWorksheet(wb, "Class_Distribution")
openxlsx::writeData(
  wb,
  "Class_Distribution",
  class_distribution
)

openxlsx::addWorksheet(wb, "MICE_Log_Summary")
openxlsx::writeData(
  wb,
  "MICE_Log_Summary",
  mice_log_summary
)

openxlsx::addWorksheet(wb, "MICE_LoggedEvents")
if (nrow(mice_logged_events_all) > 0) {
  openxlsx::writeData(
    wb,
    "MICE_LoggedEvents",
    mice_logged_events_all
  )
} else {
  openxlsx::writeData(
    wb,
    "MICE_LoggedEvents",
    tibble(Status = "No logged events")
  )
}

openxlsx::addWorksheet(wb, "Frozen_Thresholds")
openxlsx::writeData(
  wb,
  "Frozen_Thresholds",
  frozen_thresholds
)

openxlsx::addWorksheet(wb, "Final_Simple_Coefficients")
openxlsx::writeData(
  wb,
  "Final_Simple_Coefficients",
  final_simple_coefficients
)

openxlsx::addWorksheet(wb, "Final_Simple_Scalers")
openxlsx::writeData(
  wb,
  "Final_Simple_Scalers",
  final_simple_scalers
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
    "baseline_model_comparison.xlsx"
  ),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 16. SAVE PRIVATE R OBJECT FOR EXTERNAL VALIDATION
# ------------------------------------------------------------------------------
# This object contains patient-level OOF predictions and fitted models.
# Keep it under private/ and do not commit it to the public repository.

baseline_results_bundle <- list(
  created_at = Sys.time(),
  stack_bundle_path = STACK_BUNDLE_PATH,
  stack_bundle_md5 = STACK_BUNDLE_MD5,
  settings = settings_table,
  imputation_variables = imputation_vars,
  derived_predictors = derived_vars,
  egfr_settings = list(
    female_code = SEX_FEMALE_CODE,
    male_code = SEX_MALE_CODE,
    race_coefficient = EGFR_RACE_COEFFICIENT
  ),
  predictor_sets = list(
    age_only = c("Age"),
    demographics = demographic_vars,
    biomarkers_only = biomarker_vars,
    demographics_plus_biomarkers = predictor_vars
  ),
  fold_map = fold_map,
  imputation_logs = imputation_log_list,
  mice_log_summary = mice_log_summary,
  mice_logged_events = mice_logged_events_all,
  oof_by_imputation = baseline_oof_by_mi,
  oof_pooled = comparison_oof,

  # Explicitly frozen classification thresholds.
  frozen_thresholds = frozen_thresholds,

  # Frozen final M1-M4 deployment objects: one fit per development imputation.
  final_simple_model_fits = final_simple_model_fits,
  final_simple_coefficients = final_simple_coefficients,
  final_simple_scalers = final_simple_scalers,

  performance_summary = performance_summary,
  incremental_auc_delong = auc_incremental_comparisons
)

saveRDS(
  baseline_results_bundle,
  file = BASELINE_RESULTS_PATH
)

# ------------------------------------------------------------------------------
# 17. CONSOLE SUMMARY
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SIMPLE BASELINE COMPARISON COMPLETED.\n\n")

cat("Primary primary comparison:\n")
print(
  key_comparison_table |>
    select(
      Model,
      AUC_95CI,
      Reference_Model,
      Delta_AUC_vs_Reference,
      Brier,
      Calibration_Intercept,
      Calibration_Slope
    )
)

cat("\nPrimary incremental comparison:\n")
print(
  auc_incremental_comparisons |>
    filter(
      Comparison ==
        "M4 Demographics + 13 biomarkers (logistic) vs M2 Age + Sex + Smoking"
    )
)

cat("\nAdditional biomarkers-only comparisons:\n")
print(
  auc_incremental_comparisons |>
    filter(
      grepl("M3 13 biomarkers only", Comparison, fixed = TRUE)
    )
)

cat("\nOutputs:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "baseline_model_comparison.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "ROC_Baseline_vs_Stacking.png"
  ),
  "\n"
)
cat(
  BASELINE_RESULTS_PATH,
  "\n"
)

cat("============================================================\n")




