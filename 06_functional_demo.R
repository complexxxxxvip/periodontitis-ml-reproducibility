# ==============================================================================
# Functional demonstration using synthetic data only
#
# This script demonstrates:
#   1) the expected source-variable structure;
#   2) multiple imputation of source variables;
#   3) deterministic calculation of derived predictors;
#   4) the SVM + XGBoost stacking architecture;
#   5) Elastic Net meta-learning and Platt calibration.
#
# IMPORTANT:
# - This script does NOT contain the fitted study model.
# - It does NOT use real participant data.
# - It does NOT reproduce the performance estimates reported in the manuscript.
# ==============================================================================

rm(list = ls())
gc()

SEED <- 2025
set.seed(SEED)

# ------------------------------------------------------------------------------
# 0. Resolve paths relative to this script
# ------------------------------------------------------------------------------
get_script_dir <- function() {
  script_file <- tryCatch(
    sys.frame(1)$ofile,
    error = function(e) NULL
  )

  if (!is.null(script_file) && length(script_file) == 1 && nzchar(script_file)) {
    return(dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- get_script_dir()
DATA_PATH <- file.path(SCRIPT_DIR, "synthetic_example_data.csv")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "demo_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DATA_PATH)) {
  stop(
    paste0(
      "Cannot find synthetic_example_data.csv.\n",
      "Expected location: ", DATA_PATH, "\n",
      "Place the CSV in the same folder as 06_functional_demo.R."
    )
  )
}

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "mice", "dplyr", "tibble", "rsample", "recipes", "themis",
  "parsnip", "workflows", "kernlab", "xgboost", "glmnet"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall them before running this demonstration."
    )
  )
}

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
  library(tibble)
  library(rsample)
  library(recipes)
  library(themis)
  library(parsnip)
  library(workflows)
})

# ------------------------------------------------------------------------------
# 2. Variable definitions
# ------------------------------------------------------------------------------
source_vars <- c(
  "Sex", "Smoke", "Age", "TP", "ALB", "CR", "Fbg", "PLT", "TC",
  "GGT", "ALP", "NE_C", "PT", "HB", "APTT", "LDH", "GLU", "TG"
)

final_predictors <- c(
  "Sex", "Smoke", "Age", "AG_ratio", "eGFR", "Fbg", "PAR", "TC",
  "GGT", "ALP", "NPR", "INR", "HB", "APTT", "LDH", "TyG"
)

raw <- read.csv(DATA_PATH, check.names = FALSE)

required_cols <- c("PDstage", source_vars)
missing_cols <- setdiff(required_cols, names(raw))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (!all(raw$PDstage %in% c("Stage1", "Stage0"))) {
  stop("PDstage must contain only Stage1 or Stage0.")
}

raw$PDstage <- factor(raw$PDstage, levels = c("Stage1", "Stage0"))

for (nm in source_vars) {
  raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))
}

if (!all(na.omit(raw$Sex) %in% c(0, 1))) {
  stop("Sex must be coded 0/1.")
}

if (!all(na.omit(raw$Smoke) %in% c(0, 1))) {
  stop("Smoke must be coded 0/1.")
}

if (any(!is.na(raw$Age) & raw$Age != round(raw$Age))) {
  stop("Age must be recorded in whole years in the synthetic demonstration.")
}

# ------------------------------------------------------------------------------
# 3. Multiple imputation of source variables only
# ------------------------------------------------------------------------------
imp_data <- raw[, source_vars, drop = FALSE]

if (anyNA(imp_data)) {

  method_vec <- mice::make.method(imp_data)
  method_vec[] <- ""

  has_missing <- vapply(imp_data, anyNA, logical(1))
  method_vec[has_missing] <- "pmm"

  pred_mat <- mice::make.predictorMatrix(imp_data)
  diag(pred_mat) <- 0

  imp <- mice::mice(
    imp_data,
    method = method_vec,
    predictorMatrix = pred_mat,
    m = 5,
    maxit = 50,
    seed = SEED,
    printFlag = FALSE,
    eps = 1e-8
  )

  completed_source <- lapply(
    seq_len(5),
    function(j) mice::complete(imp, action = j)
  )

} else {

  # Keep the same five-dataset structure if the example file is complete.
  completed_source <- replicate(
    5,
    as.data.frame(imp_data),
    simplify = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Deterministic derived predictors
# ------------------------------------------------------------------------------
derive_predictors <- function(df) {

  df$Sex <- round(df$Sex)
  df$Smoke <- round(df$Smoke)

  scr_mg_dl <- df$CR / 88.4
  female <- df$Sex == 0
  male <- df$Sex == 1

  kappa <- ifelse(female, 0.7, ifelse(male, 0.9, NA_real_))
  alpha <- ifelse(female, -0.329, ifelse(male, -0.411, NA_real_))
  sex_coef <- ifelse(female, 1.018, ifelse(male, 1.0, NA_real_))

  out <- data.frame(
    Sex = df$Sex,
    Smoke = df$Smoke,
    Age = df$Age,
    AG_ratio = df$ALB / (df$TP - df$ALB),
    eGFR = 141 *
      pmin(scr_mg_dl / kappa, 1)^alpha *
      pmax(scr_mg_dl / kappa, 1)^(-1.209) *
      0.993^df$Age * sex_coef,
    Fbg = df$Fbg,
    PAR = df$PLT / df$ALB,
    TC = df$TC,
    GGT = df$GGT,
    ALP = df$ALP,
    NPR = df$NE_C / df$PLT,
    INR = (df$PT / 13.1)^1.31,
    HB = df$HB,
    APTT = df$APTT,
    LDH = df$LDH,
    TyG = log(1594.26 * df$TG * df$GLU / 2)
  )

  out <- out[, final_predictors, drop = FALSE]

  bad <- !vapply(
    out,
    function(x) all(is.finite(as.numeric(x))),
    logical(1)
  )

  if (any(bad)) {
    stop(
      paste0(
        "Non-finite values remain after imputation/derivation in: ",
        paste(names(out)[bad], collapse = ", ")
      )
    )
  }

  out
}

completed <- lapply(completed_source, derive_predictors)

# ------------------------------------------------------------------------------
# 5. Functional architecture demonstration
#
# The meta-learner below is fitted to in-sample base predictions to keep this
# synthetic demonstration concise. These predictions are NOT validation results.
# The full study analysis uses leakage-controlled out-of-fold predictions.
# ------------------------------------------------------------------------------
split_source <- tibble(
  .row = seq_len(nrow(raw)),
  PDstage = raw$PDstage
)

set.seed(SEED)
sp <- initial_split(split_source, prop = 0.80, strata = PDstage)

train_idx <- sort(analysis(sp)$.row)
test_idx <- sort(assessment(sp)$.row)

svm_train_mat <- matrix(NA_real_, nrow = length(train_idx), ncol = 5)
svm_test_mat  <- matrix(NA_real_, nrow = length(test_idx),  ncol = 5)
xgb_train_mat <- matrix(NA_real_, nrow = length(train_idx), ncol = 5)
xgb_test_mat  <- matrix(NA_real_, nrow = length(test_idx),  ncol = 5)

for (j in seq_len(5)) {

  dat_j <- completed[[j]]
  dat_j$PDstage <- raw$PDstage

  train_j <- dat_j[train_idx, c("PDstage", final_predictors)]
  test_j  <- dat_j[test_idx,  c("PDstage", final_predictors)]

  rec <- recipe(PDstage ~ ., data = train_j) |>
    step_normalize(all_numeric_predictors()) |>
    step_smote(PDstage) |>
    step_tomek(PDstage)

  svm_spec <- svm_rbf(
    cost = 4.10,
    rbf_sigma = 0.00145
  ) |>
    set_engine("kernlab") |>
    set_mode("classification")

  xgb_spec <- boost_tree(
    mtry = 0.774,
    trees = 635,
    min_n = 31,
    tree_depth = 4,
    learn_rate = 0.0149,
    loss_reduction = 0.00435,
    sample_size = 0.513
  ) |>
    set_engine(
      "xgboost",
      nthread = 1,
      verbosity = 0,
      counts = FALSE
    ) |>
    set_mode("classification")

  svm_fit <- workflow() |>
    add_recipe(rec) |>
    add_model(svm_spec) |>
    fit(train_j)

  xgb_fit <- workflow() |>
    add_recipe(rec) |>
    add_model(xgb_spec) |>
    fit(train_j)

  svm_train_mat[, j] <- predict(
    svm_fit, train_j, type = "prob"
  )$.pred_Stage1

  svm_test_mat[, j] <- predict(
    svm_fit, test_j, type = "prob"
  )$.pred_Stage1

  xgb_train_mat[, j] <- predict(
    xgb_fit, train_j, type = "prob"
  )$.pred_Stage1

  xgb_test_mat[, j] <- predict(
    xgb_fit, test_j, type = "prob"
  )$.pred_Stage1
}

prediction_objects <- list(
  svm_train_mat = svm_train_mat,
  svm_test_mat = svm_test_mat,
  xgb_train_mat = xgb_train_mat,
  xgb_test_mat = xgb_test_mat
)

for (nm in names(prediction_objects)) {
  if (any(!is.finite(prediction_objects[[nm]]))) {
    stop(
      paste0(
        "Non-finite base-learner predictions detected in ", nm, "."
      )
    )
  }
}

svm_train <- rowMeans(svm_train_mat)
svm_test  <- rowMeans(svm_test_mat)
xgb_train <- rowMeans(xgb_train_mat)
xgb_test  <- rowMeans(xgb_test_mat)

# Small synthetic datasets can occasionally produce a nearly constant base-model
# probability vector. Guard against division by zero when forming meta features.
safe_standardize <- function(train_values, test_values, label) {

  mu <- mean(train_values)
  sigma <- sd(train_values)

  if (!is.finite(mu)) {
    stop(paste0(label, " mean is non-finite."))
  }

  if (!is.finite(sigma) || sigma < 1e-12) {
    warning(
      paste0(
        label,
        " predictions have near-zero variance in this synthetic run; ",
        "the standardized demonstration feature is set to zero."
      )
    )

    return(
      list(
        train = rep(0, length(train_values)),
        test = rep(0, length(test_values)),
        mean = mu,
        sd = NA_real_
      )
    )
  }

  list(
    train = (train_values - mu) / sigma,
    test = (test_values - mu) / sigma,
    mean = mu,
    sd = sigma
  )
}

svm_scaled <- safe_standardize(
  svm_train, svm_test, "SVM"
)

xgb_scaled <- safe_standardize(
  xgb_train, xgb_test, "XGBoost"
)

meta_train <- data.frame(
  PDstage = raw$PDstage[train_idx],
  SVM_z = svm_scaled$train,
  XGB_z = xgb_scaled$train
)

meta_test <- data.frame(
  SVM_z = svm_scaled$test,
  XGB_z = xgb_scaled$test
)

if (any(!is.finite(as.matrix(meta_train[, c("SVM_z", "XGB_z")])))) {
  stop("Non-finite values remain in meta-learner training inputs.")
}

if (any(!is.finite(as.matrix(meta_test[, c("SVM_z", "XGB_z")])))) {
  stop("Non-finite values remain in meta-learner test inputs.")
}

meta_spec <- logistic_reg(
  penalty = 0.0117,
  mixture = 0.0429
) |>
  set_engine("glmnet") |>
  set_mode("classification")

meta_fit <- workflow() |>
  add_formula(PDstage ~ SVM_z + XGB_z) |>
  add_model(meta_spec) |>
  fit(meta_train)

stack_train <- predict(
  meta_fit, meta_train, type = "prob"
)$.pred_Stage1

stack_test <- predict(
  meta_fit, meta_test, type = "prob"
)$.pred_Stage1

platt_df <- data.frame(
  labels = ifelse(
    raw$PDstage[train_idx] == "Stage1",
    1,
    0
  ),
  probs = stack_train
)

platt <- glm(
  labels ~ probs,
  data = platt_df,
  family = binomial()
)

calibrated_test <- predict(
  platt,
  newdata = data.frame(probs = stack_test),
  type = "response"
)

# Frozen study threshold is shown only to demonstrate the classification step.
# Synthetic outputs are not study predictions.
threshold <- 0.6519343132

demo_predictions <- data.frame(
  Synthetic_row = test_idx,
  Calibrated_probability = calibrated_test,
  Demonstration_class = ifelse(
    calibrated_test >= threshold,
    "Stage III-IV",
    "Stage I-II"
  )
)

write.csv(
  demo_predictions,
  file.path(
    OUTPUT_DIR,
    "synthetic_demo_predictions.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFunctional demonstration completed successfully.\n",
  "Input data are fully synthetic and contain no real participant records.\n",
  "Synthetic predictions are for code demonstration only and are not study results.\n",
  "Output: ",
  file.path(OUTPUT_DIR, "synthetic_demo_predictions.csv"),
  "\n",
  sep = ""
)
