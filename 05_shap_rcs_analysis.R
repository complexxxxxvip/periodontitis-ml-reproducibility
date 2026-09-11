# ==============================================================================
# Exploratory SHAP and Restricted Cubic Spline Subgroup Analysis
# Final development cohort: n = 1321
#
# Analyses:
#   1. Sex-stratified SHAP |value| Mann-Whitney U test + BH-FDR
#   2. Age-stratified SHAP |value| Mann-Whitney U test + BH-FDR
#   3. Sex-stratified RCS logistic regression
#   4. Age-stratified RCS logistic regression
#   5. Family-specific BH-FDR correction for RCS interaction P values:
#        Sex family = 14 prespecified RCS features
#        Age family = 13 prespecified RCS features
#      Outputs both Pinteraction and Qinteraction.
#
# Model variable names and display labels:
#   AG_ratio -> displayed as A/G
#   GGT      -> displayed as γ-GT
#
# Age subgroup:
#   Younger <= 39.46 years
#   Older   > 39.46 years
# ==============================================================================


rm(list = ls())
gc()


# ==============================================================================
# 0. USER-ADJUSTABLE PARAMETERS
# ==============================================================================

SEED <- 2025


# ------------------------------------------------------------------------------
# Working directory
# ------------------------------------------------------------------------------

# Project-relative paths. Run this script from the repository root.
PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PRIVATE_DIR <- file.path(PROJECT_DIR, "private")
DATA_DIR <- file.path(PRIVATE_DIR, "data")
OUTPUT_ROOT <- file.path(PRIVATE_DIR, "analysis_outputs")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# The private/ directory is intended for local use only and should not be
# committed to the public repository.


# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------

RAW_DATA_PATH <- file.path(DATA_DIR, "development_data.xlsx")


# SHAP_new.R produced this workbook
SHAP_RESULTS_PATH <- file.path(
  OUTPUT_ROOT,
  "shap_analysis",
  "shap_analysis_results.xlsx"
)


SHAP_SHEET <- "SHAP_Values"



# ------------------------------------------------------------------------------
# Output folders
# ------------------------------------------------------------------------------

OUTPUT_DIR <- file.path(OUTPUT_ROOT, "shap_rcs_analysis")

SHAP_U_DIR <- file.path(
  OUTPUT_DIR,
  "SHAP_U"
)

RCS_DIR <- file.path(
  OUTPUT_DIR,
  "RCS"
)


dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


dir.create(
  SHAP_U_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


dir.create(
  RCS_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)



# ------------------------------------------------------------------------------
# Age subgroup cutoff
#
# IMPORTANT:
# raw age in years
#
# This is the raw-age cutoff, not a standardized value.
# ------------------------------------------------------------------------------

AGE_CUTOFF <- 39.46


AGE_GROUP_LABELS <- c(
  
  paste0(
    "Younger (age ≤ 39)"
  ),
  
  paste0(
    "Older (age > 39)"
  )
  
)



# ------------------------------------------------------------------------------
# SHAP U test
#
# TRUE is recommended:
#
# Sex subgroup analysis:
#     remove Sex SHAP itself
#
# Age subgroup analysis:
#     remove Age SHAP itself
# ------------------------------------------------------------------------------

EXCLUDE_GROUPING_FEATURE <- TRUE


# Multiple-testing correction
FDR_METHOD <- "BH"


# Minimum N required in each group
MIN_GROUP_N <- 3



# ------------------------------------------------------------------------------
# RCS settings
# ------------------------------------------------------------------------------

RCS_KNOT_PROBS <- c(
  0.05,
  0.275,
  0.50,
  0.725,
  0.95
)


RCS_BOOT_B <- 1000


RCS_PRED_POINTS <- 500



# ------------------------------------------------------------------------------
# Winsorization
#
# The reported analysis uses 1%-99% winsorization.
#
# TRUE  = preserve that analysis strategy
# FALSE = do not winsorize
#
# Keep this setting aligned with the Methods.
# ------------------------------------------------------------------------------

RCS_WINSORIZE <- TRUE


RCS_WINSOR_PROBS <- c(
  0.01,
  0.99
)



# ------------------------------------------------------------------------------
# RCS missing values
#
# complete_case:
# use original raw data and exclude only patients missing the
# target variable for that specific RCS.
#
# The raw dataset has minimal missingness,
# this may mean N=1320 rather than 1321 for a few variables.
# ------------------------------------------------------------------------------

RCS_MISSING_MODE <- "complete_case"


if (
  RCS_MISSING_MODE != "complete_case"
) {
  
  stop(
    "Currently RCS_MISSING_MODE must be 'complete_case'."
  )
  
}



# ==============================================================================
# RCS TARGET VARIABLES
# ==============================================================================
#
# Batch mode:
#   - Sex-stratified RCS: all 14 continuous predictors
#   - Age-stratified RCS: all continuous predictors except Age itself
#
# Sex and Smoke are binary variables and therefore are not modeled with RCS.
# Age is not used as the RCS target in the age-stratified analysis because
# Age itself defines the subgroup.
# ==============================================================================


RCS_SEX_TARGETS <- c(
  "Age",
  "AG_ratio",
  "eGFR",
  "Fbg",
  "PAR",
  "TC",
  "GGT",
  "ALP",
  "NPR",
  "INR",
  "HB",
  "APTT",
  "LDH",
  "TyG"
)


RCS_AGE_TARGETS <- c(
  "AG_ratio",
  "eGFR",
  "Fbg",
  "PAR",
  "TC",
  "GGT",
  "ALP",
  "NPR",
  "INR",
  "HB",
  "APTT",
  "LDH",
  "TyG"
)


# If TRUE, failure of one individual RCS model will be logged and
# the script will continue with the remaining variables.
RCS_CONTINUE_ON_ERROR <- TRUE



# ------------------------------------------------------------------------------
# Figure size
# ------------------------------------------------------------------------------

SHAP_U_WIDTH <- 8

SHAP_U_HEIGHT <- 5


RCS_WIDTH <- 12

RCS_HEIGHT <- 9

# ------------------------------------------------------------------------------
# Combined RCS figure settings
# ------------------------------------------------------------------------------

RCS_COMPOSITE_NCOL <- 4
RCS_COMPOSITE_WIDTH <- 18
RCS_COMPOSITE_HEIGHT <- 15
RCS_COMPOSITE_PNG_DPI <- 600

RCS_INTERACTION_LABEL_SIZE <- 2.8
RCS_INTERACTION_LABEL_FILL <- "white"
RCS_INTERACTION_LABEL_ALPHA <- 0.90
RCS_INTERACTION_LABEL_BORDER <- "grey60"




# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "rms",
  "scales",
  "patchwork"
  
)



missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]



if (
  length(
    missing_packages
  ) > 0
) {
  
  stop(
    paste0(
      "Please install missing packages:\n",
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
  
  library(openxlsx)
  
  library(rms)
  
  library(scales)
  
  library(patchwork)
  
})


set.seed(
  SEED
)



# ==============================================================================
# 2. VARIABLE DEFINITIONS
# ==============================================================================

predictor_vars <- c(
  
  "Sex",
  "Smoke",
  "Age",
  "AG_ratio",
  "eGFR",
  "Fbg",
  "PAR",
  "TC",
  "GGT",
  "ALP",
  "NPR",
  "INR",
  "HB",
  "APTT",
  "LDH",
  "TyG"
  
)



binary_vars <- c(
  
  "Sex",
  "Smoke"
  
)



continuous_vars <- setdiff(
  
  predictor_vars,
  
  binary_vars
  
)



# ==============================================================================
# 3. DISPLAY LABELS
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



feature_x_label <- function(x) {
  
  
  if (
    x == "AG_ratio"
  ) {
    
    return(
      "A/G"
    )
    
  }
  
  
  if (
    x == "GGT"
  ) {
    
    return(
      expression(
        gamma*"-GT"
      )
    )
    
  }
  
  
  return(
    x
  )
  
}



# ==============================================================================
# 4. PDF SAVE FUNCTION
#
# Mac-safe:
# Helvetica
# no cairo_pdf()
# no Arial.ttf
# no showtext
# ==============================================================================

save_pdf <- function(
    filename,
    plot_object,
    width,
    height
) {
  
  
  pdf(
    
    filename,
    
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
# 5. OUTCOME CONVERSION
# ==============================================================================

to_stage_binary <- function(x) {
  
  
  x_chr <- trimws(
    
    tolower(
      as.character(x)
    )
    
  )
  
  
  out <- case_when(
    
    x_chr %in%
      c(
        "stage1",
        "1"
      ) ~ 1,
    
    
    x_chr %in%
      c(
        "stage0",
        "0"
      ) ~ 0,
    
    
    TRUE ~ NA_real_
    
  )
  
  
  if (
    anyNA(out)
  ) {
    
    stop(
      "PDstage contains missing or unrecognized values."
    )
    
  }
  
  
  return(
    out
  )
  
}



# ==============================================================================
# 6. READ RAW 1321-PATIENT DATA
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



required_raw_cols <- c(
  
  "PDstage",
  
  predictor_vars
  
)



missing_raw_cols <- setdiff(
  
  required_raw_cols,
  
  names(
    raw_data
  )
  
)



if (
  length(
    missing_raw_cols
  ) > 0
) {
  
  stop(
    paste0(
      "Missing raw-data columns:\n",
      paste(
        missing_raw_cols,
        collapse = ", "
      )
    )
  )
  
}



raw_data <- raw_data |>
  
  mutate(
    
    .Row = seq_len(
      n()
    ),
    
    Y = to_stage_binary(
      PDstage
    ),
    
    across(
      
      all_of(
        predictor_vars
      ),
      
      ~ suppressWarnings(
        as.numeric(.x)
      )
      
    )
    
  )



if (
  nrow(
    raw_data
  ) != 1321
) {
  
  warning(
    paste0(
      "Expected N=1321, but current N=",
      nrow(
        raw_data
      )
    )
  )
  
}



# ------------------------------------------------------------------------------
# Check Sex / Smoke
# ------------------------------------------------------------------------------

for (
  nm in binary_vars
) {
  
  
  vals <- sort(
    
    unique(
      
      raw_data[[nm]][
        !is.na(
          raw_data[[nm]]
        )
      ]
      
    )
    
  )
  
  
  if (
    !all(
      vals %in%
      c(
        0,
        1
      )
    )
  ) {
    
    stop(
      paste0(
        nm,
        " must be coded 0/1."
      )
    )
    
  }
  
}



# ------------------------------------------------------------------------------
# Sex / Age groups
# ------------------------------------------------------------------------------

raw_data <- raw_data |>
  
  mutate(
    
    
    Sex_Group = factor(
      
      Sex,
      
      levels = c(
        0,
        1
      ),
      
      labels = c(
        "Female",
        "Male"
      )
      
    ),
    
    
    Age_Group = factor(
      
      ifelse(
        
        Age <= AGE_CUTOFF,
        
        AGE_GROUP_LABELS[1],
        
        AGE_GROUP_LABELS[2]
        
      ),
      
      levels =
        AGE_GROUP_LABELS
      
    )
    
  )



cat(
  
  "\n============================================================\n",
  
  "DEVELOPMENT COHORT\n",
  
  "============================================================\n",
  
  "Total N = ",
  nrow(raw_data),
  "\n",
  
  "Female N = ",
  sum(
    raw_data$Sex_Group == "Female",
    na.rm = TRUE
  ),
  "\n",
  
  "Male N = ",
  sum(
    raw_data$Sex_Group == "Male",
    na.rm = TRUE
  ),
  "\n",
  
  AGE_GROUP_LABELS[1],
  " N = ",
  sum(
    raw_data$Age_Group ==
      AGE_GROUP_LABELS[1],
    na.rm = TRUE
  ),
  "\n",
  
  AGE_GROUP_LABELS[2],
  " N = ",
  sum(
    raw_data$Age_Group ==
      AGE_GROUP_LABELS[2],
    na.rm = TRUE
  ),
  "\n",
  
  sep = ""
  
)



# ==============================================================================
# 7. READ SHAP VALUES
# ==============================================================================

if (
  !file.exists(
    SHAP_RESULTS_PATH
  )
) {
  
  stop(
    paste0(
      "Cannot find:\n",
      SHAP_RESULTS_PATH
    )
  )
  
}



shap_sheets <- readxl::excel_sheets(
  SHAP_RESULTS_PATH
)



if (
  !SHAP_SHEET %in%
  shap_sheets
) {
  
  stop(
    paste0(
      "Cannot find sheet: ",
      SHAP_SHEET
    )
  )
  
}



shap_data <- readxl::read_excel(
  
  SHAP_RESULTS_PATH,
  
  sheet = SHAP_SHEET
  
) |>
  
  as.data.frame()



if (
  !("Row" %in%
    names(
      shap_data
    ))
) {
  
  stop(
    "SHAP_Values must contain the Row column."
  )
  
}



missing_shap_features <- setdiff(
  
  predictor_vars,
  
  names(
    shap_data
  )
  
)



if (
  length(
    missing_shap_features
  ) > 0
) {
  
  stop(
    paste0(
      "Missing SHAP features:\n",
      paste(
        missing_shap_features,
        collapse = ", "
      )
    )
  )
  
}



if (
  nrow(
    shap_data
  ) !=
  nrow(
    raw_data
  )
) {
  
  stop(
    "SHAP rows do not match the raw dataset."
  )
  
}



shap_data$Row <- as.integer(
  shap_data$Row
)



if (
  !identical(
    shap_data$Row,
    raw_data$.Row
  )
) {
  
  stop(
    paste0(
      "SHAP Row index does not match the ",
      "original 1321-patient order."
    )
  )
  
}



shap_values <- shap_data |>
  
  select(
    all_of(
      predictor_vars
    )
  ) |>
  
  mutate(
    
    across(
      
      everything(),
      
      ~ suppressWarnings(
        as.numeric(.x)
      )
      
    )
    
  )



if (
  anyNA(
    shap_values
  )
) {
  
  stop(
    "NA values exist in SHAP matrix."
  )
  
}



# ==============================================================================
# 8. MANN-WHITNEY U FUNCTION
# ==============================================================================

safe_mann_whitney <- function(
    x,
    group,
    level1,
    level2
) {
  
  
  keep <- is.finite(x) &
    !is.na(group)
  
  
  x <- x[keep]
  
  group <- group[keep]
  
  
  x1 <- x[
    group == level1
  ]
  
  
  x2 <- x[
    group == level2
  ]
  
  
  x1 <- x1[
    is.finite(x1)
  ]
  
  
  x2 <- x2[
    is.finite(x2)
  ]
  
  
  n1 <- length(
    x1
  )
  
  
  n2 <- length(
    x2
  )
  
  
  
  if (
    n1 < MIN_GROUP_N ||
    n2 < MIN_GROUP_N
  ) {
    
    return(
      
      data.frame(
        
        N1 = n1,
        
        N2 = n2,
        
        Mean1 =
          ifelse(
            n1 > 0,
            mean(x1),
            NA
          ),
        
        Mean2 =
          ifelse(
            n2 > 0,
            mean(x2),
            NA
          ),
        
        Median1 =
          ifelse(
            n1 > 0,
            median(x1),
            NA
          ),
        
        Median2 =
          ifelse(
            n2 > 0,
            median(x2),
            NA
          ),
        
        U = NA,
        
        Rank_Biserial = NA,
        
        P_Raw = NA
        
      )
      
    )
    
  }
  
  
  
  test_result <- suppressWarnings(
    
    wilcox.test(
      
      x1,
      
      x2,
      
      exact = FALSE
      
    )
    
  )
  
  
  
  U <- as.numeric(
    test_result$statistic
  )
  
  
  
  # Rank-biserial effect size
  #
  # Positive:
  # level1 tends to have larger values
  #
  # Negative:
  # level2 tends to have larger values
  
  rank_biserial <-
    
    2 * U /
    (
      n1 * n2
    ) - 1
  
  
  
  return(
    
    data.frame(
      
      N1 = n1,
      
      N2 = n2,
      
      Mean1 = mean(x1),
      
      Mean2 = mean(x2),
      
      Median1 = median(x1),
      
      Median2 = median(x2),
      
      U = U,
      
      Rank_Biserial =
        rank_biserial,
      
      P_Raw =
        test_result$p.value
      
    )
    
  )
  
}



add_fdr <- function(df) {
  
  
  df |>
    
    mutate(
      
      P_Adj = p.adjust(
        
        P_Raw,
        
        method = FDR_METHOD
        
      ),
      
      stars = case_when(
        
        is.na(P_Adj) ~ "N/A",
        
        P_Adj < 0.001 ~ "***",
        
        P_Adj < 0.01 ~ "**",
        
        P_Adj < 0.05 ~ "*",
        
        TRUE ~ "ns"
        
      )
      
    )
  
}



# ==============================================================================
# 9. SEX-STRATIFIED SHAP U TEST
# ==============================================================================

sex_features <- predictor_vars



if (
  EXCLUDE_GROUPING_FEATURE
) {
  
  sex_features <- setdiff(
    
    sex_features,
    
    "Sex"
    
  )
  
}



sex_results <- lapply(
  
  sex_features,
  
  function(feat) {
    
    
    tmp <- safe_mann_whitney(
      
      x = abs(
        shap_values[[feat]]
      ),
      
      group =
        raw_data$Sex_Group,
      
      level1 =
        "Female",
      
      level2 =
        "Male"
      
    )
    
    
    
    data.frame(
      
      Feature = feat,
      
      N_Female =
        tmp$N1,
      
      N_Male =
        tmp$N2,
      
      Mean_Female =
        tmp$Mean1,
      
      Mean_Male =
        tmp$Mean2,
      
      Median_Female =
        tmp$Median1,
      
      Median_Male =
        tmp$Median2,
      
      U_Statistic =
        tmp$U,
      
      Rank_Biserial_Female_vs_Male =
        tmp$Rank_Biserial,
      
      P_Raw =
        tmp$P_Raw
      
    )
    
  }
  
) |>
  
  bind_rows() |>
  
  add_fdr() |>
  
  mutate(
    
    Total_Imp =
      
      (
        Mean_Female +
          Mean_Male
      ) / 2
    
  ) |>
  
  arrange(
    Total_Imp
  ) |>
  
  mutate(
    
    Feature = factor(
      
      Feature,
      
      levels =
        unique(
          Feature
        )
      
    )
    
  )



sex_plot_data <- sex_results |>
  
  pivot_longer(
    
    cols = c(
      
      Mean_Female,
      
      Mean_Male
      
    ),
    
    names_to =
      "Group",
    
    values_to =
      "Mean_Abs_SHAP"
    
  ) |>
  
  mutate(
    
    Group = factor(
      
      Group,
      
      levels = c(
        
        "Mean_Female",
        
        "Mean_Male"
        
      ),
      
      labels = c(
        
        "Female",
        
        "Male"
        
      )
      
    ),
    
    Feature = factor(
      
      Feature,
      
      levels =
        levels(
          sex_results$Feature
        )
      
    )
    
  )



sex_colors <- c(
  
  "Female" =
    "#D73027",
  
  "Male" =
    "#4575B4"
  
)



sex_offset <-
  
  max(
    sex_plot_data$Mean_Abs_SHAP,
    na.rm = TRUE
  ) * 0.08



p_shap_sex <- ggplot(
  
  sex_plot_data,
  
  aes(
    
    x =
      Mean_Abs_SHAP,
    
    y =
      Feature
    
  )
  
) +
  
  
  geom_segment(
    
    data =
      sex_results,
    
    aes(
      
      x =
        Mean_Female,
      
      xend =
        Mean_Male,
      
      y =
        Feature,
      
      yend =
        Feature
      
    ),
    
    inherit.aes = FALSE,
    
    color = "#EEEEEE",
    
    linewidth = 2.5,
    
    alpha = 0.8
    
  ) +
  
  
  geom_point(
    
    aes(
      color = Group
    ),
    
    size = 4.2
    
  ) +
  
  
  geom_text(
    
    data =
      sex_results,
    
    aes(
      
      x =
        pmax(
          Mean_Female,
          Mean_Male
        ) +
        sex_offset,
      
      y =
        Feature,
      
      label =
        stars
      
    ),
    
    inherit.aes = FALSE,
    
    family = "Helvetica",
    
    fontface = "bold",
    
    size = 4.3
    
  ) +
  
  
  scale_y_discrete(
    
    labels =
      feature_axis_labels
    
  ) +
  
  
  scale_color_manual(
    
    values =
      sex_colors
    
  ) +
  
  
  scale_x_continuous(
    
    expand =
      expansion(
        
        mult = c(
          0.02,
          0.16
        )
        
      )
    
  ) +
  
  
  theme_minimal(
    
    base_family =
      "Helvetica",
    
    base_size = 11
    
  ) +
  
  
  theme(
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    axis.title.y =
      element_text(
        
        size = 12,
        
        face = "bold",
        
        margin =
          margin(
            r = 15
          )
        
      ),
    
    axis.text.y =
      element_text(
        
        size = 11,
        
        face = "bold",
        
        color =
          "grey20"
        
      ),
    
    axis.title.x =
      element_text(
        
        size = 11,
        
        margin =
          margin(
            t = 10
          )
        
      ),
    
    legend.position =
      "top",
    
    legend.title =
      element_blank(),
    
    plot.margin =
      margin(
        20,
        25,
        20,
        20
      )
    
  ) +
  
  
  labs(
    
    x =
      "Average impact on model output (Mean |SHAP|)",
    
    y =
      "Feature"
    
  )



print(
  p_shap_sex
)



save_pdf(
  
  filename =
    file.path(
      
      SHAP_U_DIR,
      
      "SHAP_U_Sex.pdf"
      
    ),
  
  plot_object =
    p_shap_sex,
  
  width =
    SHAP_U_WIDTH,
  
  height =
    SHAP_U_HEIGHT
  
)



# ==============================================================================
# 10. AGE-STRATIFIED SHAP U TEST
# ==============================================================================

age_features <- predictor_vars



if (
  EXCLUDE_GROUPING_FEATURE
) {
  
  age_features <- setdiff(
    
    age_features,
    
    "Age"
    
  )
  
}



age_results <- lapply(
  
  age_features,
  
  function(feat) {
    
    
    tmp <- safe_mann_whitney(
      
      x = abs(
        shap_values[[feat]]
      ),
      
      group =
        raw_data$Age_Group,
      
      level1 =
        AGE_GROUP_LABELS[1],
      
      level2 =
        AGE_GROUP_LABELS[2]
      
    )
    
    
    
    data.frame(
      
      Feature = feat,
      
      N_Younger =
        tmp$N1,
      
      N_Older =
        tmp$N2,
      
      Mean_Younger =
        tmp$Mean1,
      
      Mean_Older =
        tmp$Mean2,
      
      Median_Younger =
        tmp$Median1,
      
      Median_Older =
        tmp$Median2,
      
      U_Statistic =
        tmp$U,
      
      Rank_Biserial_Younger_vs_Older =
        tmp$Rank_Biserial,
      
      P_Raw =
        tmp$P_Raw
      
    )
    
  }
  
) |>
  
  bind_rows() |>
  
  add_fdr() |>
  
  mutate(
    
    Total_Imp =
      
      (
        Mean_Younger +
          Mean_Older
      ) / 2
    
  ) |>
  
  arrange(
    Total_Imp
  ) |>
  
  mutate(
    
    Feature = factor(
      
      Feature,
      
      levels =
        unique(
          Feature
        )
      
    )
    
  )



age_plot_data <- age_results |>
  
  pivot_longer(
    
    cols = c(
      
      Mean_Younger,
      
      Mean_Older
      
    ),
    
    names_to =
      "Group",
    
    values_to =
      "Mean_Abs_SHAP"
    
  ) |>
  
  mutate(
    
    Group = factor(
      
      Group,
      
      levels = c(
        
        "Mean_Younger",
        
        "Mean_Older"
        
      ),
      
      labels =
        AGE_GROUP_LABELS
      
    ),
    
    Feature = factor(
      
      Feature,
      
      levels =
        levels(
          age_results$Feature
        )
      
    )
    
  )



age_colors <- setNames(
  
  c(
    "#00AFBB",
    "#E18727"
  ),
  
  AGE_GROUP_LABELS
  
)



age_offset <-
  
  max(
    age_plot_data$Mean_Abs_SHAP,
    na.rm = TRUE
  ) * 0.08



p_shap_age <- ggplot(
  
  age_plot_data,
  
  aes(
    
    x =
      Mean_Abs_SHAP,
    
    y =
      Feature
    
  )
  
) +
  
  
  geom_segment(
    
    data =
      age_results,
    
    aes(
      
      x =
        Mean_Younger,
      
      xend =
        Mean_Older,
      
      y =
        Feature,
      
      yend =
        Feature
      
    ),
    
    inherit.aes = FALSE,
    
    color =
      "#EEEEEE",
    
    linewidth = 2.5,
    
    alpha = 0.8
    
  ) +
  
  
  geom_point(
    
    aes(
      color = Group
    ),
    
    size = 4.2
    
  ) +
  
  
  geom_text(
    
    data =
      age_results,
    
    aes(
      
      x =
        pmax(
          Mean_Younger,
          Mean_Older
        ) +
        age_offset,
      
      y =
        Feature,
      
      label =
        stars
      
    ),
    
    inherit.aes = FALSE,
    
    family =
      "Helvetica",
    
    fontface =
      "bold",
    
    size = 4.3
    
  ) +
  
  
  scale_y_discrete(
    
    labels =
      feature_axis_labels
    
  ) +
  
  
  scale_color_manual(
    
    values =
      age_colors
    
  ) +
  
  
  scale_x_continuous(
    
    expand =
      expansion(
        
        mult = c(
          0.02,
          0.16
        )
        
      )
    
  ) +
  
  
  theme_minimal(
    
    base_family =
      "Helvetica",
    
    base_size = 11
    
  ) +
  
  
  theme(
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    axis.title.y =
      element_text(
        
        size = 12,
        
        face = "bold",
        
        margin =
          margin(
            r = 15
          )
        
      ),
    
    axis.text.y =
      element_text(
        
        size = 11,
        
        face = "bold",
        
        color =
          "grey20"
        
      ),
    
    axis.title.x =
      element_text(
        
        size = 11,
        
        margin =
          margin(
            t = 10
          )
        
      ),
    
    legend.position =
      "top",
    
    legend.title =
      element_blank(),
    
    plot.margin =
      margin(
        20,
        25,
        20,
        20
      )
    
  ) +
  
  
  labs(
    
    x =
      "Average impact on model output (Mean |SHAP|)",
    
    y =
      "Feature"
    
  )



print(
  p_shap_age
)



save_pdf(
  
  filename =
    file.path(
      
      SHAP_U_DIR,
      
      "SHAP_U_Age.pdf"
      
    ),
  
  plot_object =
    p_shap_age,
  
  width =
    SHAP_U_WIDTH,
  
  height =
    SHAP_U_HEIGHT
  
)



# ==============================================================================
# 11. EXPORT SHAP U RESULTS
# ==============================================================================

SHAP_U_settings <- data.frame(
  
  Parameter = c(
    
    "Development_N",
    
    "Age_cutoff_years",
    
    "Exclude_grouping_feature",
    
    "FDR_method",
    
    "Minimum_group_N"
    
  ),
  
  Value = c(
    
    nrow(
      raw_data
    ),
    
    AGE_CUTOFF,
    
    EXCLUDE_GROUPING_FEATURE,
    
    FDR_METHOD,
    
    MIN_GROUP_N
    
  )
  
)



openxlsx::write.xlsx(
  
  list(
    
    
    Settings =
      SHAP_U_settings,
    
    
    Sex_SHAP_U =
      sex_results |>
      
      mutate(
        Feature =
          as.character(
            Feature
          )
      ),
    
    
    Age_SHAP_U =
      age_results |>
      
      mutate(
        Feature =
          as.character(
            Feature
          )
      )
    
    
  ),
  
  file =
    file.path(
      
      SHAP_U_DIR,
      
      "SHAP_U_Subgroup_Results.xlsx"
      
    ),
  
  rowNames = FALSE,
  
  overwrite = TRUE
  
)



# ==============================================================================
# 12. RCS DATA PREPARATION
# ==============================================================================

prepare_rcs_data <- function(
    target_var,
    subgroup_type = c(
      "sex",
      "age"
    )
) {
  
  
  subgroup_type <-
    match.arg(
      subgroup_type
    )
  
  
  
  if (
    !(target_var %in%
      continuous_vars)
  ) {
    
    stop(
      paste0(
        target_var,
        " is not an eligible continuous RCS variable."
      )
    )
    
  }
  
  
  
  if (
    subgroup_type == "age" &&
    target_var == "Age"
  ) {
    
    stop(
      paste0(
        "Age cannot be the target ",
        "when Age defines the subgroup."
      )
    )
    
  }
  
  
  
  if (
    subgroup_type == "sex"
  ) {
    
    
    group_var <-
      "Sex_Group"
    
    
    group_levels <- c(
      
      "Female",
      
      "Male"
      
    )
    
    
    reference_group <-
      "Female"
    
    
  } else {
    
    
    group_var <-
      "Age_Group"
    
    
    group_levels <-
      AGE_GROUP_LABELS
    
    
    reference_group <-
      AGE_GROUP_LABELS[1]
    
  }
  
  
  
  dat <- raw_data |>
    
    transmute(
      
      .Row =
        .Row,
      
      Y =
        Y,
      
      .Target =
        .data[[target_var]],
      
      .Group =
        .data[[group_var]]
      
    ) |>
    
    filter(
      
      !is.na(Y),
      
      !is.na(.Target),
      
      !is.na(.Group)
      
    )
  
  
  
  dat$.Group <- factor(
    
    dat$.Group,
    
    levels =
      group_levels
    
  )
  
  
  
  analysis_n <-
    nrow(
      dat
    )
  
  
  
  # ============================================================================
  # Winsorization
  # ============================================================================
  
  if (
    RCS_WINSORIZE
  ) {
    
    
    lower_bound <-
      as.numeric(
        
        quantile(
          
          dat$.Target,
          
          probs =
            RCS_WINSOR_PROBS[1],
          
          na.rm = TRUE
          
        )
        
      )
    
    
    
    upper_bound <-
      as.numeric(
        
        quantile(
          
          dat$.Target,
          
          probs =
            RCS_WINSOR_PROBS[2],
          
          na.rm = TRUE
          
        )
        
      )
    
    
    
    dat$.Target <- pmax(
      
      pmin(
        
        dat$.Target,
        
        upper_bound
        
      ),
      
      lower_bound
      
    )
    
    
  } else {
    
    
    lower_bound <-
      min(
        dat$.Target,
        na.rm = TRUE
      )
    
    
    upper_bound <-
      max(
        dat$.Target,
        na.rm = TRUE
      )
    
  }
  
  
  
  # ============================================================================
  # Prespecified knots
  # ============================================================================
  
  reference_values <-
    
    dat$.Target[
      
      dat$.Group ==
        reference_group
      
    ]
  
  
  
  knots_val <-
    
    as.numeric(
      
      quantile(
        
        reference_values,
        
        probs =
          RCS_KNOT_PROBS,
        
        na.rm = TRUE,
        
        names = FALSE
        
      )
      
    )
  
  
  
  if (
    length(
      unique(
        knots_val
      )
    ) !=
    length(
      knots_val
    )
  ) {
    
    stop(
      paste0(
        "Duplicated RCS knot values for ",
        target_var,
        "."
      )
    )
    
  }
  
  
  
  dat$bin <- cut(
    
    dat$.Target,
    
    breaks = c(
      
      -Inf,
      
      knots_val,
      
      Inf
      
    ),
    
    include.lowest = TRUE
    
  )
  
  
  
  return(
    
    list(
      
      data = dat,
      
      group_levels =
        group_levels,
      
      reference_group =
        reference_group,
      
      knots =
        knots_val,
      
      lower_bound =
        lower_bound,
      
      upper_bound =
        upper_bound,
      
      analysis_n =
        analysis_n
      
    )
    
  )
  
}



# ==============================================================================
# 13. RCS BIN PREDICTED-PROBABILITY U TEST
#
# Previous version: within each knot-defined interval, the U test compared
# the RAW predictor (.Target) between subgroups.
#
# Within each interval, compare patient-level RCS model-
# predicted probabilities of Stage1 between subgroups.
#
# Five RCS knots define SIX intervals:
# (-Inf,k1], (k1,k2], (k2,k3], (k3,k4], (k4,k5], (k5,Inf].
#
# These interval tests are secondary/descriptive. Formal global inference
# remains the RCS ANOVA, especially the subgroup-by-RCS interaction.
# ==============================================================================

rcs_bin_u_test <- function(dat, group_levels) {

  if (!(".Predicted_Probability" %in% names(dat))) {
    stop("rcs_bin_u_test requires .Predicted_Probability.")
  }

  bin_levels <- levels(dat$bin)

  result <- lapply(
    bin_levels,
    function(bin_i) {

      dat_i <- dat |> filter(bin == bin_i)

      tmp <- safe_mann_whitney(
        x = dat_i$.Predicted_Probability,
        group = dat_i$.Group,
        level1 = group_levels[1],
        level2 = group_levels[2]
      )

      data.frame(
        Bin = bin_i,
        N_Group1 = tmp$N1,
        N_Group2 = tmp$N2,
        Mean_Probability_Group1 = tmp$Mean1,
        Mean_Probability_Group2 = tmp$Mean2,
        Median_Probability_Group1 = tmp$Median1,
        Median_Probability_Group2 = tmp$Median2,
        Median_Probability_Difference_Group1_minus_Group2 =
          tmp$Median1 - tmp$Median2,
        U_Statistic = tmp$U,
        Rank_Biserial = tmp$Rank_Biserial,
        P_Raw = tmp$P_Raw,
        Mid_Point = ifelse(
          nrow(dat_i) > 0,
          median(dat_i$.Target, na.rm = TRUE),
          NA
        ),
        Tested_Quantity =
          "Patient-level RCS predicted probability of Stage1",
        stringsAsFactors = FALSE
      )
    }
  ) |>
    bind_rows() |>
    mutate(
      P_Adj = p.adjust(P_Raw, method = FDR_METHOD),
      stars = case_when(
        is.na(P_Adj) ~ "N/A",
        P_Adj < 0.001 ~ "***",
        P_Adj < 0.01 ~ "**",
        P_Adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  return(result)
}




# ------------------------------------------------------------------------------
# Robust extractor for rms::lrm fit statistics
# ------------------------------------------------------------------------------

extract_lrm_named_stat <- function(fit_object, aliases, statistic_label) {
  s <- fit_object$stats
  if (is.null(s) || is.null(names(s))) {
    stop(paste0("Cannot extract ", statistic_label, ": fit$stats unavailable."))
  }

  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  idx <- which(norm(names(s)) %in% norm(aliases))

  if (length(idx) < 1) {
    stop(
      paste0(
        "Cannot identify ", statistic_label,
        ". Available fit$stats names: ",
        paste(names(s), collapse = ", ")
      )
    )
  }

  value <- suppressWarnings(as.numeric(s[idx[1]]))
  if (length(value) != 1 || !is.finite(value)) {
    stop(paste0(statistic_label, " is not a finite scalar."))
  }
  value
}



# ==============================================================================
# 14. SINGLE RCS MODEL + FIGURE FUNCTION
# ==============================================================================

run_one_rcs <- function(
    target_var,
    subgroup_type = c(
      "sex",
      "age"
    )
) {
  
  
  subgroup_type <-
    match.arg(
      subgroup_type
    )
  
  
  
  prep <- prepare_rcs_data(
    
    target_var =
      target_var,
    
    subgroup_type =
      subgroup_type
    
  )
  
  
  
  dat <-
    prep$data
  
  
  
  cat(
    
    "\n============================================================\n",
    
    "RCS: ",
    toupper(
      subgroup_type
    ),
    " | ",
    target_var,
    "\n",
    
    "============================================================\n",
    
    "Analysis N = ",
    prep$analysis_n,
    "\n",
    
    "Knots = ",
    paste(
      round(
        prep$knots,
        5
      ),
      collapse = ", "
    ),
    "\n",
    
    sep = ""
    
  )
  
  
  
  # ============================================================================
  # datadist
  #
  # rms stores only the NAME of the datadist object in options(datadist=...).
  # Because run_one_rcs() is a function, a local object such as "dd" is not
  # reliably visible to rms::Design(). Therefore the temporary datadist object
  # is deliberately placed in .GlobalEnv and removed automatically on exit.
  # ============================================================================

  old_datadist <-
    getOption(
      "datadist"
    )


  dd_object_name <-
    ".dd_rcs_temp"


  if (
    exists(
      dd_object_name,
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {

    rm(
      list = dd_object_name,
      envir = .GlobalEnv
    )

  }


  assign(

    dd_object_name,

    rms::datadist(
      dat
    ),

    envir = .GlobalEnv

  )


  options(
    datadist = dd_object_name
  )


  on.exit(

    {

      options(
        datadist = old_datadist
      )


      if (
        exists(
          dd_object_name,
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {

        rm(
          list = dd_object_name,
          envir = .GlobalEnv
        )

      }

    },

    add = TRUE

  )



  # ============================================================================
  # RCS logistic model
  #
  # IMPORTANT:
  # rms::rcs() takes the user-specified knot vector as its SECOND argument.
  # Use the rcs() argument structure below for the fixed knot specification.
  # and can cause rms to fall back to its default 5-knot placement.
  #
  # Here the exact numeric knot values are written into the formula so that
  # the prespecified 5th, 27.5th, 50th, 72.5th and 95th percentile knots
  # calculated above are the knots actually used by the model.
  # ============================================================================

  knots_text <- paste(

    format(
      prep$knots,
      digits = 15,
      scientific = FALSE,
      trim = TRUE
    ),

    collapse = ", "

  )


  formula_text <- paste0(

    "Y ~ rcs(.Target, c(",

    knots_text,

    ")) * .Group"

  )


  rcs_formula <- as.formula(
    formula_text
  )


  cat(

    "RCS formula: ",

    formula_text,

    "\n",

    sep = ""

  )


  fit_rcs <- rms::lrm(

    formula = rcs_formula,

    data = dat,

    x = TRUE,

    y = TRUE

  )


  # ============================================================================
  # GLOBAL SUBGROUP-BY-RCS INTERACTION TEST
  #
  # FIX: do not search rms::anova() for a literal "Full RCS*subgroup model versus reduced RCS+subgroup model" row.
  # Compare two nested models instead:
  #   Full    : RCS * subgroup
  #   Reduced : RCS + subgroup
  # ============================================================================

  reduced_formula_text <- paste0(
    "Y ~ rcs(.Target, c(",
    knots_text,
    ")) + .Group"
  )

  fit_rcs_reduced <- rms::lrm(
    formula = as.formula(reduced_formula_text),
    data = dat,
    x = TRUE,
    y = TRUE
  )

  full_lr_chisq <- extract_lrm_named_stat(
    fit_rcs,
    c("Model L.R.", "Model LR", "Model LR Chi-Square"),
    "full-model LR chi-square"
  )

  reduced_lr_chisq <- extract_lrm_named_stat(
    fit_rcs_reduced,
    c("Model L.R.", "Model LR", "Model LR Chi-Square"),
    "reduced-model LR chi-square"
  )

  full_model_df <- extract_lrm_named_stat(
    fit_rcs,
    c("d.f.", "df"),
    "full-model degrees of freedom"
  )

  reduced_model_df <- extract_lrm_named_stat(
    fit_rcs_reduced,
    c("d.f.", "df"),
    "reduced-model degrees of freedom"
  )

  interaction_chisq <- full_lr_chisq - reduced_lr_chisq
  interaction_df <- full_model_df - reduced_model_df

  if (interaction_chisq < 0 && abs(interaction_chisq) < 1e-8) {
    interaction_chisq <- 0
  }

  if (
    !is.finite(interaction_chisq) ||
    interaction_chisq < 0 ||
    !is.finite(interaction_df) ||
    interaction_df <= 0
  ) {
    stop(
      paste0(
        "Invalid global interaction LR test for ", target_var,
        ": chi-square=", interaction_chisq,
        ", df=", interaction_df
      )
    )
  }

  interaction_p <- stats::pchisq(
    interaction_chisq,
    df = interaction_df,
    lower.tail = FALSE
  )

  interaction_test_raw <- data.frame(
    Target = target_var,
    Stratification = ifelse(subgroup_type == "sex", "Sex", "Age"),
    Interaction_Term = "Global subgroup-by-RCS interaction",
    Test_Method = "Nested-model likelihood-ratio test",
    Interaction_ChiSquare = as.numeric(interaction_chisq),
    Interaction_df = as.numeric(interaction_df),
    Pinteraction = as.numeric(interaction_p),
    Full_Model = formula_text,
    Reduced_Model = reduced_formula_text,
    stringsAsFactors = FALSE
  )

  cat(
    "\nGLOBAL INTERACTION TEST: chi-square=",
    signif(interaction_chisq, 7),
    ", df=", interaction_df,
    ", Pinteraction=",
    format.pval(interaction_p, digits = 4, eps = 0.001),
    "\n",
    sep = ""
  )




  

  # ============================================================================
  # PATIENT-LEVEL RCS PREDICTED PROBABILITIES
  # Used for the knot-interval subgroup U tests.
  # ============================================================================

  dat$.Predicted_Probability <- plogis(
    as.numeric(
      predict(
        fit_rcs,
        newdata = dat,
        type = "lp"
      )
    )
  )

  if (
    any(!is.finite(dat$.Predicted_Probability)) ||
    any(dat$.Predicted_Probability < 0 | dat$.Predicted_Probability > 1)
  ) {
    stop("Invalid patient-level predicted probabilities from RCS model.")
  }


# ============================================================================
  # Bootstrap covariance
  # ============================================================================
  
  set.seed(
    SEED
  )
  
  
  fit_boot <- bootcov(
    
    fit_rcs,
    
    B =
      RCS_BOOT_B
    
  )
  
  
  
  # ============================================================================
  # Prediction
  # ============================================================================
  
  pred_data <- Predict(
    
    fit_boot,
    
    .Target,
    
    .Group,
    
    fun = plogis,
    
    np =
      RCS_PRED_POINTS
    
  ) |>
    
    as.data.frame() |>
    
    mutate(
      
      lower =
        pmax(
          lower,
          0
        ),
      
      upper =
        pmin(
          upper,
          1
        )
      
    )
  
  
  
  # ============================================================================
  # Knot-interval U tests on patient-level predicted probabilities
  # ============================================================================
  
  sig_data <- rcs_bin_u_test(
    
    dat =
      dat,
    
    group_levels =
      prep$group_levels
    
  )
  
  
  
  # ============================================================================
  # RCS ANOVA
  # ============================================================================
  
  anova_obj <-
    anova(
      fit_rcs
    )
  
  
  
  anova_df <- data.frame(
    
    Term =
      rownames(
        anova_obj
      ),
    
    as.data.frame(
      
      anova_obj,
      
      check.names =
        FALSE
      
    ),
    
    row.names =
      NULL,
    
    check.names =
      FALSE
    
  )
  
  
  
  # ============================================================================
  # Colors
  # ============================================================================
  
  if (
    subgroup_type ==
    "sex"
  ) {
    
    
    plot_colors <- c(
      
      "Female" =
        "#D73027",
      
      "Male" =
        "#4575B4"
      
    )
    
    
    legend_title <-
      "Sex"
    
    
  } else {
    
    
    plot_colors <- setNames(
      
      c(
        "#00AFBB",
        "#E18727"
      ),
      
      prep$group_levels
      
    )
    
    
    legend_title <-
      "Age group"
    
  }
  
  
  
  # ============================================================================
  # TOP RCS PLOT
  # ============================================================================
  
  p_top <- ggplot(
    
    pred_data,
    
    aes(
      
      x =
        .Target,
      
      y =
        yhat,
      
      group =
        .Group
      
    )
    
  ) +
    
    
    geom_vline(
      
      xintercept =
        prep$knots,
      
      linetype =
        "dashed",
      
      color =
        "grey88",
      
      linewidth =
        0.45
      
    ) +
    
    
    geom_ribbon(
      
      aes(
        
        ymin =
          lower,
        
        ymax =
          upper,
        
        fill =
          .Group
        
      ),
      
      alpha =
        0.15,
      
      color =
        NA
      
    ) +
    
    
    geom_line(
      
      aes(
        color =
          .Group
      ),
      
      linewidth =
        1.05
      
    ) +
    
    
    scale_color_manual(
      
      values =
        plot_colors,
      
      name =
        legend_title
      
    ) +
    
    
    scale_fill_manual(
      
      values =
        plot_colors,
      
      name =
        legend_title
      
    ) +
    
    
    scale_y_continuous(
      
      labels =
        scales::percent_format(
          accuracy = 1
        ),
      
      expand =
        expansion(
          
          mult = c(
            0.01,
            0.04
          )
          
        )
      
    ) +
    
    
    coord_cartesian(
      
      xlim = c(
        
        prep$lower_bound,
        
        prep$upper_bound
        
      )
      
    ) +
    
    
    theme_bw(
      
      base_family =
        "Helvetica",
      
      base_size =
        11
      
    ) +
    
    
    theme(
      
      legend.position =
        "top",
      
      axis.title.x =
        element_blank(),
      
      axis.text.x =
        element_blank(),
      
      axis.ticks.x =
        element_blank(),
      
      panel.grid =
        element_blank(),
      
      plot.margin =
        margin(
          b = 2,
          t = 10,
          l = 10,
          r = 10
        )
      
    ) +
    
    
    labs(
      
      y =
        "Estimated probability of Stage III–IV periodontitis"
      
    )
  
  
  
  # ============================================================================
  # BOTTOM BOXPLOT
  # ============================================================================
  
  p_bottom <- ggplot(
    
    dat,
    
    aes(
      
      x =
        .Target,
      
      y =
        .Group,
      
      fill =
        .Group
      
    )
    
  ) +
    
    
    geom_vline(
      
      xintercept =
        prep$knots,
      
      linetype =
        "dashed",
      
      color =
        "grey88",
      
      linewidth =
        0.45
      
    ) +
    
    
    geom_boxplot(
      
      aes(
        color =
          .Group
      ),
      
      width =
        0.5,
      
      alpha =
        0.30,
      
      outlier.shape =
        NA
      
    ) +
    
    
    geom_text(
      
      data =
        sig_data,
      
      aes(
        
        x =
          Mid_Point,
        
        label =
          stars
        
      ),
      
      y =
        1.5,
      
      inherit.aes =
        FALSE,
      
      family =
        "Helvetica",
      
      fontface =
        "bold",
      
      size =
        4.3,
      
      color =
        "black",
      
      na.rm =
        TRUE
      
    ) +
    
    
    scale_color_manual(
      
      values =
        plot_colors
      
    ) +
    
    
    scale_fill_manual(
      
      values =
        plot_colors
      
    ) +
    
    
    coord_cartesian(
      
      xlim = c(
        
        prep$lower_bound,
        
        prep$upper_bound
        
      )
      
    ) +
    
    
    theme_bw(
      
      base_family =
        "Helvetica",
      
      base_size =
        11
      
    ) +
    
    
    theme(
      
      legend.position =
        "none",
      
      panel.grid =
        element_blank(),
      
      axis.title.y =
        element_blank(),
      
      plot.margin =
        margin(
          t = 2,
          b = 10,
          l = 10,
          r = 10
        )
      
    ) +
    
    
    labs(
      
      x =
        feature_x_label(
          target_var
        )
      
    )
  
  
  
  # ============================================================================
  # Combine
  # ============================================================================
  
  final_plot <-
    
    p_top /
    
    p_bottom +
    
    plot_layout(
      
      heights =
        c(
          3.5,
          1
        )
      
    )
  
  
  
  # ============================================================================
  # File names
  # ============================================================================
  
  if (
    subgroup_type ==
    "sex"
  ) {
    
    prefix <-
      "RCS_Sex"
    
  } else {
    
    prefix <-
      "RCS_Age"
    
  }
  
  
  
  pdf_file <- file.path(
    
    RCS_DIR,
    
    paste0(
      
      prefix,
      
      "_",
      
      target_var,
      
      ".pdf"
      
    )
    
  )
  
  
  
  excel_file <- file.path(
    
    RCS_DIR,
    
    paste0(
      
      prefix,
      
      "_",
      
      target_var,
      
      "_Results.xlsx"
      
    )
    
  )
  
  
  
  summary_file <- file.path(
    
    RCS_DIR,
    
    paste0(
      
      prefix,
      
      "_",
      
      target_var,
      
      "_ModelSummary.txt"
      
    )
    
  )
  
  
  
  # ============================================================================
  # Save PDF
  # ============================================================================
  
  save_pdf(
    
    filename =
      pdf_file,
    
    plot_object =
      final_plot,
    
    width =
      RCS_WIDTH,
    
    height =
      RCS_HEIGHT
    
  )
  
  
  
  # ============================================================================
  # Settings table
  # ============================================================================
  
  # ============================================================================
  # Patient-level fitted probability audit
  # ============================================================================

  patient_probability_df <- dat |>
    transmute(
      Row = .Row,
      Target_Value = .Target,
      Subgroup = as.character(.Group),
      Interval = as.character(bin),
      Predicted_Probability_Stage1 = .Predicted_Probability
    )


  settings_df <- data.frame(
    
    Parameter = c(
      
      "Target",
      
      "Subgroup",
      
      "Development_cohort_N",
      
      "Analysis_N",
      
      "Age_cutoff",
      
      "Winsorization",
      
      "Lower_plot_bound",
      
      "Upper_plot_bound",
      
      "Knot_percentiles",
      
      "Knot_values",
      
      "Knot_reference_group",
      
      "Bootstrap_B",
      
      "FDR_method",

      "Interval_test_quantity",

      "Number_of_RCS_knots",

      "Number_of_knot_defined_intervals"
      
    ),
    
    
    Value = c(
      
      target_var,
      
      subgroup_type,
      
      nrow(
        raw_data
      ),
      
      prep$analysis_n,
      
      ifelse(
        subgroup_type ==
          "age",
        AGE_CUTOFF,
        NA
      ),
      
      RCS_WINSORIZE,
      
      prep$lower_bound,
      
      prep$upper_bound,
      
      paste(
        RCS_KNOT_PROBS,
        collapse = ", "
      ),
      
      paste(
        prep$knots,
        collapse = ", "
      ),
      
      prep$reference_group,
      
      RCS_BOOT_B,
      
      FDR_METHOD,

      "Patient-level model-predicted probability of Stage1",

      length(prep$knots),

      length(levels(dat$bin))
      
    )
    
  )
  
  
  
  # ============================================================================
  # Save Excel
  # ============================================================================
  
  openxlsx::write.xlsx(
    
    list(
      
      
      Settings =
        settings_df,
      
      
      Bin_Probability_U_Test =
        sig_data,

      Patient_Fitted_Probabilities =
        patient_probability_df,
      
      
      RCS_ANOVA =
        anova_df,


      RCS_Interaction_Raw =
        interaction_test_raw,

      
      
      RCS_Predictions =
        pred_data
      
      
    ),
    
    file =
      excel_file,
    
    rowNames =
      FALSE,
    
    overwrite =
      TRUE
    
  )
  
  
  
  # ============================================================================
  # Save model summary
  # ============================================================================
  
  capture.output(
    
    
    cat(
      "============================================================\n"
    ),
    
    
    cat(
      "RCS MODEL SUMMARY\n"
    ),
    
    
    cat(
      "============================================================\n"
    ),
    
    
    cat(
      "Target: ",
      target_var,
      "\n"
    ),
    
    
    cat(
      "Subgroup: ",
      subgroup_type,
      "\n"
    ),
    
    
    cat(
      "Analysis N: ",
      prep$analysis_n,
      "\n"
    ),
    
    
    cat(
      "Knots: ",
      paste(
        prep$knots,
        collapse = ", "
      ),
      "\n\n"
    ),
    
    
    print(
      fit_rcs
    ),
    
    
    cat(
      "\n\nKNOT-INTERVAL SUBGROUP TEST:\n"
    ),

    cat(
      "Tested quantity: patient-level model-predicted probability of Stage1\n"
    ),

    cat(
      "Number of RCS knots: ", length(prep$knots), "\n"
    ),

    cat(
      "Number of knot-defined intervals: ", length(levels(dat$bin)), "\n"
    ),

    cat(
      "Mann-Whitney U with BH-FDR across intervals; formal global inference is from RCS ANOVA.\n"
    ),

    cat(
      "\nRCS ANOVA:\n"
    ),
    
    
    print(
      anova_obj
    ),
    
    
    file =
      summary_file
    
  )
  
  
  
  print(
    final_plot
  )
  
  
  
  invisible(
    
    list(
      
      fit =
        fit_rcs,


      reduced_fit =
        fit_rcs_reduced,

      interaction_test_raw =
        interaction_test_raw,

      
      bootstrap_fit =
        fit_boot,
      
      predictions =
        pred_data,
      
      bin_probability_U =
        sig_data,
      
      anova =
        anova_df,
      
      knots =
        prep$knots,
      

      top_plot =
        p_top,

      bottom_plot =
        p_bottom,
      
      plot =
        final_plot
      
    )
    
  )
  
}



# ==============================================================================
# 15. BATCH RCS RUNNER
# ==============================================================================

run_rcs_batch <- function(
    targets,
    subgroup_type = c(
      "sex",
      "age"
    )
) {


  subgroup_type <- match.arg(
    subgroup_type
  )


  results <- list()


  error_log <- data.frame(

    Subgroup = character(0),

    Target = character(0),

    Error = character(0),

    stringsAsFactors = FALSE

  )


  total_targets <- length(
    targets
  )


  for (
    i in seq_along(
      targets
    )
  ) {


    target_var <- targets[i]


    cat(

      "\n\n############################################################\n",

      "BATCH RCS ",

      toupper(
        subgroup_type
      ),

      ": ",

      i,

      "/",

      total_targets,

      " | ",

      target_var,

      "\n",

      "############################################################\n",

      sep = ""

    )


    one_result <- tryCatch(

      {

        run_one_rcs(

          target_var = target_var,

          subgroup_type = subgroup_type

        )

      },

      error = function(e) {


        error_message <- conditionMessage(
          e
        )


        cat(

          "\n*** RCS ERROR ***\n",

          "Subgroup: ",

          subgroup_type,

          "\nTarget: ",

          target_var,

          "\nMessage: ",

          error_message,

          "\n",

          sep = ""

        )


        error_log <<- bind_rows(

          error_log,

          data.frame(

            Subgroup = subgroup_type,

            Target = target_var,

            Error = error_message,

            stringsAsFactors = FALSE

          )

        )


        if (
          !RCS_CONTINUE_ON_ERROR
        ) {

          stop(
            e
          )

        }


        return(
          NULL
        )

      }

    )


    results[[target_var]] <- one_result

  }


  list(

    results = results,

    errors = error_log

  )

}



# ==============================================================================
# 16. RUN ALL SEX- AND AGE-STRATIFIED RCS MODELS
#
# One source() call automatically generates:
#
#   Sex-stratified:
#       14 continuous-variable RCS figures
#
#   Age-stratified:
#       13 continuous-variable RCS figures
#
# No manual target-variable editing is required.
# ==============================================================================

sex_rcs_batch <- run_rcs_batch(

  targets = RCS_SEX_TARGETS,

  subgroup_type = "sex"

)


sex_rcs_results <- sex_rcs_batch$results

sex_rcs_errors <- sex_rcs_batch$errors



age_rcs_batch <- run_rcs_batch(

  targets = RCS_AGE_TARGETS,

  subgroup_type = "age"

)


age_rcs_results <- age_rcs_batch$results

age_rcs_errors <- age_rcs_batch$errors



rcs_error_log <- bind_rows(

  sex_rcs_errors,

  age_rcs_errors

)


if (
  nrow(
    rcs_error_log
  ) > 0
) {


  openxlsx::write.xlsx(

    rcs_error_log,

    file = file.path(

      RCS_DIR,

      "RCS_Batch_Error_Log.xlsx"

    ),

    rowNames = FALSE,

    overwrite = TRUE

  )


  cat(

    "\n============================================================\n",

    "RCS BATCH COMPLETED WITH ",

    nrow(
      rcs_error_log
    ),

    " ERROR(S).\n",

    "See: ",

    file.path(
      RCS_DIR,
      "RCS_Batch_Error_Log.xlsx"
    ),

    "\n============================================================\n",

    sep = ""

  )


} else {


  cat(

    "\n============================================================\n",

    "ALL BATCH RCS MODELS COMPLETED SUCCESSFULLY\n",

    "Sex RCS figures: ",

    length(
      RCS_SEX_TARGETS
    ),

    "\nAge RCS figures: ",

    length(
      RCS_AGE_TARGETS
    ),

    "\n============================================================\n",

    sep = ""

  )

}




# ==============================================================================
# 16A. RCS INTERACTION P/Q VALUES WITH FAMILY-SPECIFIC BH-FDR CORRECTION
#
# Pinteraction:
#   Raw GLOBAL interaction P value from rms::anova(fit_rcs), using the
#   "TOTAL INTERACTION" row.
#
# Qinteraction:
#   Benjamini-Hochberg FDR-adjusted P value.
#
# Prespecified multiple-testing families:
#   1) Sex RCS family: 14 features
#   2) Age RCS family: 13 features
#
# The two families are adjusted SEPARATELY. They are NOT pooled into one
# 27-test family.
# ==============================================================================


format_interaction_p <- function(p) {

  out <- rep(
    NA_character_,
    length(p)
  )

  ok <- is.finite(
    p
  )

  if (
    any(
      ok
    )
  ) {

    p2 <- p[
      ok
    ]

    out[
      ok
    ] <- ifelse(

      p2 < 0.001,

      "<0.001",

      ifelse(

        p2 < 0.01,

        sprintf(
          "%.4f",
          p2
        ),

        sprintf(
          "%.3f",
          p2
        )

      )

    )

  }

  out

}


extract_global_interaction <- function(
    one_result,
    target_var,
    family_name
) {

  empty_row <- function(status_text) {
    data.frame(
      Family = family_name,
      Target = target_var,
      Interaction_Term = NA_character_,
      Interaction_ChiSquare = NA_real_,
      Interaction_df = NA_real_,
      Pinteraction = NA_real_,
      Extraction_Status = status_text,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(one_result)) {
    return(empty_row("RCS model unavailable / failed"))
  }

  if (is.null(one_result$interaction_test_raw)) {
    return(
      empty_row(
        "interaction_test_raw unavailable; rerun this script from the beginning"
      )
    )
  }

  z <- as.data.frame(one_result$interaction_test_raw)

  needed <- c(
    "Interaction_Term",
    "Interaction_ChiSquare",
    "Interaction_df",
    "Pinteraction"
  )

  missing_cols <- setdiff(needed, names(z))
  if (length(missing_cols) > 0 || nrow(z) != 1) {
    return(empty_row("Invalid raw interaction-test table"))
  }

  chi <- suppressWarnings(as.numeric(z$Interaction_ChiSquare[1]))
  dfv <- suppressWarnings(as.numeric(z$Interaction_df[1]))
  pv <- suppressWarnings(as.numeric(z$Pinteraction[1]))

  data.frame(
    Family = family_name,
    Target = target_var,
    Interaction_Term = as.character(z$Interaction_Term[1]),
    Interaction_ChiSquare = chi,
    Interaction_df = dfv,
    Pinteraction = pv,
    Extraction_Status = ifelse(
      is.finite(chi) && is.finite(dfv) && is.finite(pv),
      "OK - nested-model LR interaction test",
      "Non-finite interaction statistic"
    ),
    stringsAsFactors = FALSE
  )
}


build_interaction_family <- function(
    result_list,
    targets,
    family_name
) {

  family_size <- length(
    targets
  )


  out <- lapply(

    targets,

    function(
      target_var
    ) {

      extract_global_interaction(

        one_result =
          result_list[[target_var]],

        target_var =
          target_var,

        family_name =
          family_name

      )

    }

  ) |>

    bind_rows()


  out$Family_Size_Prespecified <-
    family_size


  out$Qinteraction <-
    NA_real_


  valid_idx <- which(
    is.finite(
      out$Pinteraction
    )
  )


  if (
    length(
      valid_idx
    ) > 0
  ) {

    # n = family_size keeps the full prespecified family denominator even
    # if one model failed and its Pinteraction is unavailable.
    out$Qinteraction[
      valid_idx
    ] <- p.adjust(

      out$Pinteraction[
        valid_idx
      ],

      method =
        FDR_METHOD,

      n =
        family_size

    )

  }


  out$Pinteraction_Display <-
    format_interaction_p(
      out$Pinteraction
    )


  out$Qinteraction_Display <-
    format_interaction_p(
      out$Qinteraction
    )


  out$Qinteraction_Significance <- case_when(

    is.na(
      out$Qinteraction
    ) ~
      "N/A",

    out$Qinteraction < 0.001 ~
      "***",

    out$Qinteraction < 0.01 ~
      "**",

    out$Qinteraction < 0.05 ~
      "*",

    TRUE ~
      "ns"

  )


  out$FDR_Method <-
    FDR_METHOD


  out$FDR_Family <- paste0(

    family_name,

    " RCS interaction family (",

    family_size,

    " prespecified tests)"

  )


  out

}


# ------------------------------------------------------------------------------
# Sex family: 14 prespecified features
# ------------------------------------------------------------------------------

if (
  length(
    RCS_SEX_TARGETS
  ) != 14
) {

  stop(
    paste0(
      "Expected exactly 14 Sex-family RCS targets, found ",
      length(
        RCS_SEX_TARGETS
      ),
      "."
    )
  )

}


sex_interaction_fdr <- build_interaction_family(

  result_list =
    sex_rcs_results,

  targets =
    RCS_SEX_TARGETS,

  family_name =
    "Sex"

)


# ------------------------------------------------------------------------------
# Age family: 13 prespecified features
# ------------------------------------------------------------------------------

if (
  length(
    RCS_AGE_TARGETS
  ) != 13
) {

  stop(
    paste0(
      "Expected exactly 13 Age-family RCS targets, found ",
      length(
        RCS_AGE_TARGETS
      ),
      "."
    )
  )

}


age_interaction_fdr <- build_interaction_family(

  result_list =
    age_rcs_results,

  targets =
    RCS_AGE_TARGETS,

  family_name =
    "Age"

)


interaction_fdr_combined <- bind_rows(

  sex_interaction_fdr,

  age_interaction_fdr

)



if (any(!is.finite(interaction_fdr_combined$Pinteraction))) {

  bad_rows <- interaction_fdr_combined[
    !is.finite(interaction_fdr_combined$Pinteraction),
    c("Family", "Target", "Extraction_Status"),
    drop = FALSE
  ]

  print(bad_rows)

  stop(
    "At least one Pinteraction is unavailable; stopping instead of exporting NA qinteraction."
  )
}



# ------------------------------------------------------------------------------
# Console audit
# ------------------------------------------------------------------------------

cat(

  "\n============================================================\n",

  "RCS Pinteraction / Qinteraction\n",

  "============================================================\n",

  "Adjustment: ",

  FDR_METHOD,

  "\nSex family: 14 tests\n",

  "Age family: 13 tests\n",

  "Families adjusted separately.\n",

  "============================================================\n",

  sep = ""

)


cat(
  "\nSEX FAMILY\n"
)

print(

  sex_interaction_fdr |>

    select(

      Target,

      Pinteraction,

      Qinteraction,

      Pinteraction_Display,

      Qinteraction_Display,

      Qinteraction_Significance,

      Extraction_Status

    )

)


cat(
  "\nAGE FAMILY\n"
)

print(

  age_interaction_fdr |>

    select(

      Target,

      Pinteraction,

      Qinteraction,

      Pinteraction_Display,

      Qinteraction_Display,

      Qinteraction_Significance,

      Extraction_Status

    )

)


if (
  any(
    !is.finite(
      sex_interaction_fdr$Pinteraction
    )
  )
) {

  warning(
    paste0(
      "At least one Sex-family interaction P value was unavailable. ",
      "BH correction still uses n=14, preserving the prespecified family."
    )
  )

}


if (
  any(
    !is.finite(
      age_interaction_fdr$Pinteraction
    )
  )
) {

  warning(
    paste0(
      "At least one Age-family interaction P value was unavailable. ",
      "BH correction still uses n=13, preserving the prespecified family."
    )
  )

}


# ------------------------------------------------------------------------------
# Export family-specific Pinteraction/Qinteraction summary
# ------------------------------------------------------------------------------

interaction_settings <- data.frame(

  Parameter = c(

    "Interaction_source",

    "Interaction_ANOVA_row",

    "FDR_method",

    "Sex_family_size",

    "Sex_family_targets",

    "Age_family_size",

    "Age_family_targets",

    "Adjustment_rule"

  ),

  Value = c(

    "Nested-model likelihood-ratio test",

    "TOTAL INTERACTION",

    FDR_METHOD,

    length(
      RCS_SEX_TARGETS
    ),

    paste(
      RCS_SEX_TARGETS,
      collapse = ", "
    ),

    length(
      RCS_AGE_TARGETS
    ),

    paste(
      RCS_AGE_TARGETS,
      collapse = ", "
    ),

    paste0(
      "BH-FDR independently within Sex family (14) and Age family (13); ",
      "not pooled across all 27 tests."
    )

  ),

  stringsAsFactors =
    FALSE

)


interaction_summary_file <- file.path(

  RCS_DIR,

  "RCS_Interaction_Pinteraction_Qinteraction_FDR.xlsx"

)


openxlsx::write.xlsx(

  list(

    Settings =
      interaction_settings,

    Sex_14_Features =
      sex_interaction_fdr,

    Age_13_Features =
      age_interaction_fdr,

    Combined =
      interaction_fdr_combined

  ),

  file =
    interaction_summary_file,

  rowNames =
    FALSE,

  overwrite =
    TRUE

)


utils::write.csv(

  interaction_fdr_combined,

  file =
    file.path(
      RCS_DIR,
      "RCS_Interaction_Pinteraction_Qinteraction_FDR.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ------------------------------------------------------------------------------
# Append each feature's corrected P/Q values to its own RCS result workbook.
# ------------------------------------------------------------------------------

append_interaction_fdr_sheet <- function(
    family_df,
    prefix
) {

  for (
    i in seq_len(
      nrow(
        family_df
      )
    )
  ) {

    target_var <-
      family_df$Target[
        i
      ]


    result_file <- file.path(

      RCS_DIR,

      paste0(
        prefix,
        "_",
        target_var,
        "_Results.xlsx"
      )

    )


    if (
      !file.exists(
        result_file
      )
    ) {

      next

    }


    one_row <- family_df[
      i,
      ,
      drop = FALSE
    ]


    tryCatch(

      {

        wb <- openxlsx::loadWorkbook(
          result_file
        )


        current_sheets <- names(
          wb
        )


        if (
          "Interaction_FDR" %in%
            current_sheets
        ) {

          openxlsx::removeWorksheet(
            wb,
            "Interaction_FDR"
          )

        }


        openxlsx::addWorksheet(
          wb,
          "Interaction_FDR"
        )


        openxlsx::writeData(
          wb,
          "Interaction_FDR",
          one_row
        )


        openxlsx::saveWorkbook(
          wb,
          result_file,
          overwrite = TRUE
        )

      },

      error = function(e) {

        warning(
          paste0(
            "Could not append Interaction_FDR to ",
            result_file,
            ": ",
            conditionMessage(
              e
            )
          )
        )

      }

    )

  }

}


append_interaction_fdr_sheet(

  family_df =
    sex_interaction_fdr,

  prefix =
    "RCS_Sex"

)


append_interaction_fdr_sheet(

  family_df =
    age_interaction_fdr,

  prefix =
    "RCS_Age"

)


cat(

  "\nInteraction P/Q summary saved to:\n",

  interaction_summary_file,

  "\n",

  sep = ""

)




# ==============================================================================
# 16B. TWO COMBINED RCS FIGURES WITH Pinteraction + qinteraction
# ==============================================================================

feature_display_plain <- function(x) {
  out <- as.character(x)
  out[out == "AG_ratio"] <- "A/G"
  out[out == "GGT"] <- "γ-GT"
  out
}


make_annotated_rcs_panel <- function(
    one_result,
    interaction_row,
    target_var,
    panel_index,
    ncol = RCS_COMPOSITE_NCOL
) {

  if (is.null(one_result) || is.null(one_result$top_plot)) {
    return(NULL)
  }

  if (is.null(interaction_row) || nrow(interaction_row) != 1) {
    p_text <- "NA"
    q_text <- "NA"
  } else {
    p_text <- as.character(interaction_row$Pinteraction_Display[1])
    q_text <- as.character(interaction_row$Qinteraction_Display[1])
  }

  annotation_text <- paste0(
    "Pinteraction = ", p_text,
    "\nqinteraction = ", q_text
  )

  p <- one_result$top_plot +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.04,
      vjust = 1.10,
      size = RCS_INTERACTION_LABEL_SIZE,
      family = "Helvetica",
      color = "black",
      fill = RCS_INTERACTION_LABEL_FILL,
      alpha = RCS_INTERACTION_LABEL_ALPHA,
      linewidth = 0.28,
      label.padding = grid::unit(0.16, "lines")
    ) +
    labs(
      x = feature_x_label(target_var),
      y = if (((panel_index - 1) %% ncol) == 0) {
        "Predicted probability"
      } else {
        NULL
      }
    ) +
    theme(
      axis.title.x = element_text(
        size = 9.5,
        face = "bold",
        margin = margin(t = 6)
      ),
      axis.text.x = element_text(
        size = 8.0,
        color = "black"
      ),
      axis.ticks.x = element_line(
        color = "black",
        linewidth = 0.35
      ),
      axis.title.y = element_text(
        size = 9.0,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8.0,
        color = "black"
      ),
      legend.position = "top",
      legend.title = element_text(
        size = 9.0,
        face = "bold"
      ),
      legend.text = element_text(
        size = 8.5
      ),
      plot.margin = margin(6, 6, 6, 6)
    )

  p
}


build_rcs_composite <- function(
    result_list,
    interaction_df,
    targets,
    figure_title
) {

  plot_list <- vector("list", length(targets))

  for (i in seq_along(targets)) {

    target_var <- targets[i]

    one_result <- result_list[[target_var]]

    interaction_row <- interaction_df[
      interaction_df$Target == target_var,
      ,
      drop = FALSE
    ]

    plot_list[[i]] <- make_annotated_rcs_panel(
      one_result = one_result,
      interaction_row = interaction_row,
      target_var = target_var,
      panel_index = i,
      ncol = RCS_COMPOSITE_NCOL
    )
  }

  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    stop(
      paste0(
        "No valid RCS panels available for: ",
        figure_title
      )
    )
  }

  combined <- patchwork::wrap_plots(
    plotlist = plot_list,
    ncol = RCS_COMPOSITE_NCOL,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = figure_title,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(
          family = "Helvetica",
          size = 15,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 8)
        ),
        plot.tag = element_text(
          family = "Helvetica",
          size = 10,
          face = "bold"
        )
      )
    )

  combined <- combined &
    theme(
      legend.position = "top",
      legend.box = "horizontal"
    )

  combined
}


sex_rcs_composite <- build_rcs_composite(
  result_list = sex_rcs_results,
  interaction_df = sex_interaction_fdr,
  targets = RCS_SEX_TARGETS,
  figure_title = "Sex-stratified restricted cubic spline analyses"
)

age_rcs_composite <- build_rcs_composite(
  result_list = age_rcs_results,
  interaction_df = age_interaction_fdr,
  targets = RCS_AGE_TARGETS,
  figure_title = "Age-stratified restricted cubic spline analyses"
)

print(sex_rcs_composite)
print(age_rcs_composite)


sex_composite_pdf <- file.path(
  RCS_DIR,
  "RCS_Sex_Combined_Pinteraction_Qinteraction.pdf"
)

save_pdf(
  filename = sex_composite_pdf,
  plot_object = sex_rcs_composite,
  width = RCS_COMPOSITE_WIDTH,
  height = RCS_COMPOSITE_HEIGHT
)

sex_composite_png <- file.path(
  RCS_DIR,
  "RCS_Sex_Combined_Pinteraction_Qinteraction.png"
)

ggsave(
  filename = sex_composite_png,
  plot = sex_rcs_composite,
  width = RCS_COMPOSITE_WIDTH,
  height = RCS_COMPOSITE_HEIGHT,
  dpi = RCS_COMPOSITE_PNG_DPI,
  bg = "white"
)


age_composite_pdf <- file.path(
  RCS_DIR,
  "RCS_Age_Combined_Pinteraction_Qinteraction.pdf"
)

save_pdf(
  filename = age_composite_pdf,
  plot_object = age_rcs_composite,
  width = RCS_COMPOSITE_WIDTH,
  height = RCS_COMPOSITE_HEIGHT
)

age_composite_png <- file.path(
  RCS_DIR,
  "RCS_Age_Combined_Pinteraction_Qinteraction.png"
)

ggsave(
  filename = age_composite_png,
  plot = age_rcs_composite,
  width = RCS_COMPOSITE_WIDTH,
  height = RCS_COMPOSITE_HEIGHT,
  dpi = RCS_COMPOSITE_PNG_DPI,
  bg = "white"
)


# ==============================================================================
# 16C. EXPORT REQUESTED RCS STATISTICAL TABLES
# ==============================================================================

interaction_export <- bind_rows(

  sex_interaction_fdr |>
    transmute(
      Stratification = "Sex",
      Feature = feature_display_plain(Target),
      Chi_square = Interaction_ChiSquare,
      df = Interaction_df,
      Raw_P_for_interaction = Pinteraction,
      BH_adjusted_P_for_interaction = Qinteraction
    ),

  age_interaction_fdr |>
    transmute(
      Stratification = "Age",
      Feature = feature_display_plain(Target),
      Chi_square = Interaction_ChiSquare,
      df = Interaction_df,
      Raw_P_for_interaction = Pinteraction,
      BH_adjusted_P_for_interaction = Qinteraction
    )
)

names(interaction_export) <- c(
  "Stratification",
  "Feature",
  "χ²",
  "df",
  "Raw P for interaction",
  "BH-adjusted P for interaction"
)


collect_interval_results <- function(
    result_list,
    targets
) {

  rows <- lapply(
    targets,
    function(target_var) {

      one_result <- result_list[[target_var]]

      if (is.null(one_result) || is.null(one_result$bin_probability_U)) {
        return(NULL)
      }

      one_result$bin_probability_U |>
        transmute(
          Feature = feature_display_plain(target_var),
          Interval = as.character(Bin),
          N_Group1 = N_Group1,
          N_Group2 = N_Group2,
          Median_probability_Group1 = Median_Probability_Group1,
          Median_probability_Group2 = Median_Probability_Group2,
          difference = Median_Probability_Difference_Group1_minus_Group2,
          raw_P = P_Raw,
          BH_adjusted_P = P_Adj
        )
    }
  ) |>
    bind_rows()

  names(rows) <- c(
    "Feature",
    "Interval",
    "N Group1",
    "N Group2",
    "Median probability Group1",
    "Median probability Group2",
    "difference",
    "raw P",
    "BH-adjusted P"
  )

  rows
}


sex_interval_export <- collect_interval_results(
  result_list = sex_rcs_results,
  targets = RCS_SEX_TARGETS
)

age_interval_export <- collect_interval_results(
  result_list = age_rcs_results,
  targets = RCS_AGE_TARGETS
)

combined_interval_export <- bind_rows(
  data.frame(
    Stratification = "Sex",
    sex_interval_export,
    check.names = FALSE
  ),
  data.frame(
    Stratification = "Age",
    age_interval_export,
    check.names = FALSE
  )
)


rcs_export_settings <- data.frame(
  Item = c(
    "Sex interaction family",
    "Age interaction family",
    "Interaction multiple-testing correction",
    "Sex interval Group1",
    "Sex interval Group2",
    "Age interval Group1",
    "Age interval Group2",
    "Interval tested quantity",
    "Interval multiple-testing correction",
    "Difference definition",
    "Number of RCS knots",
    "Number of knot-defined intervals"
  ),
  Value = c(
    "14 prespecified RCS features",
    "13 prespecified RCS features",
    paste0(
      FDR_METHOD,
      " separately within Sex and Age interaction families"
    ),
    "Female",
    "Male",
    AGE_GROUP_LABELS[1],
    AGE_GROUP_LABELS[2],
    "Patient-level RCS predicted probability of Stage III-IV periodontitis",
    paste0(
      FDR_METHOD,
      " separately across the 6 intervals within each feature"
    ),
    "Median probability Group1 minus Median probability Group2",
    length(RCS_KNOT_PROBS),
    length(RCS_KNOT_PROBS) + 1
  ),
  stringsAsFactors = FALSE
)


rcs_statistics_file <- file.path(
  RCS_DIR,
  "RCS_Interaction_and_Interval_Statistics.xlsx"
)

wb_rcs <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb_rcs, "Interaction_Summary")
openxlsx::writeData(
  wb_rcs,
  "Interaction_Summary",
  interaction_export
)

openxlsx::addWorksheet(wb_rcs, "Sex_Intervals")
openxlsx::writeData(
  wb_rcs,
  "Sex_Intervals",
  sex_interval_export
)

openxlsx::addWorksheet(wb_rcs, "Age_Intervals")
openxlsx::writeData(
  wb_rcs,
  "Age_Intervals",
  age_interval_export
)

openxlsx::addWorksheet(wb_rcs, "Combined_Intervals")
openxlsx::writeData(
  wb_rcs,
  "Combined_Intervals",
  combined_interval_export
)

openxlsx::addWorksheet(wb_rcs, "Settings")
openxlsx::writeData(
  wb_rcs,
  "Settings",
  rcs_export_settings
)


rcs_header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fgFill = "#EAF2F8",
  border = "Bottom",
  borderStyle = "thin"
)

sheet_col_counts <- c(
  Interaction_Summary = ncol(interaction_export),
  Sex_Intervals = ncol(sex_interval_export),
  Age_Intervals = ncol(age_interval_export),
  Combined_Intervals = ncol(combined_interval_export),
  Settings = ncol(rcs_export_settings)
)

for (sheet_nm in names(sheet_col_counts)) {

  n_cols <- sheet_col_counts[[sheet_nm]]

  openxlsx::addStyle(
    wb_rcs,
    sheet = sheet_nm,
    style = rcs_header_style,
    rows = 1,
    cols = seq_len(n_cols),
    gridExpand = TRUE
  )

  openxlsx::freezePane(
    wb_rcs,
    sheet = sheet_nm,
    firstActiveRow = 2
  )

  openxlsx::setColWidths(
    wb_rcs,
    sheet = sheet_nm,
    cols = seq_len(n_cols),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  wb_rcs,
  file = rcs_statistics_file,
  overwrite = TRUE
)


cat(
  "\n============================================================\n",
  "COMBINED RCS FIGURES + REQUESTED TABLES COMPLETED\n",
  "============================================================\n",
  "Sex combined RCS PDF:\n", sex_composite_pdf, "\n\n",
  "Sex combined RCS PNG:\n", sex_composite_png, "\n\n",
  "Age combined RCS PDF:\n", age_composite_pdf, "\n\n",
  "Age combined RCS PNG:\n", age_composite_png, "\n\n",
  "Requested Excel tables:\n", rcs_statistics_file, "\n",
  "============================================================\n",
  sep = ""
)



# ==============================================================================
# 17. OVERALL SETTINGS
# ==============================================================================

overall_settings <- data.frame(
  
  Parameter = c(
    
    "Development_N",
    
    "Age_cutoff",
    
    "SHAP_U_exclude_group_variable",
    
    "FDR_method",
    
    "RCS_knot_percentiles",
    
    "RCS_bootstrap_B",
    
    "RCS_winsorization",
    
    "RCS_winsor_percentiles",
    
    "RCS_sex_targets",
    
    "RCS_age_targets"
    
  ),
  
  
  Value = c(
    
    nrow(
      raw_data
    ),
    
    AGE_CUTOFF,
    
    EXCLUDE_GROUPING_FEATURE,
    
    FDR_METHOD,
    
    paste(
      RCS_KNOT_PROBS,
      collapse = ", "
    ),
    
    RCS_BOOT_B,
    
    RCS_WINSORIZE,
    
    paste(
      RCS_WINSOR_PROBS,
      collapse = ", "
    ),
    
    paste(
      RCS_SEX_TARGETS,
      collapse = ", "
    ),
    
    paste(
      RCS_AGE_TARGETS,
      collapse = ", "
    )
    
  )
  
)



raw_missingness <- data.frame(
  
  Feature =
    predictor_vars,
  
  Missing_N = vapply(
    
    raw_data[
      predictor_vars
    ],
    
    function(x) {
      
      sum(
        is.na(x)
      )
      
    },
    
    integer(1)
    
  )
  
)



openxlsx::write.xlsx(
  
  list(
    
    Settings =
      overall_settings,
    
    Raw_Missingness =
      raw_missingness
    
  ),
  
  file =
    file.path(
      
      OUTPUT_DIR,
      
      "shap_rcs_settings.xlsx"
      
    ),
  
  rowNames =
    FALSE,
  
  overwrite =
    TRUE
  
)



# ==============================================================================
# 18. SESSION INFO
# ==============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      
      OUTPUT_DIR,
      
      "R_SessionInfo.txt"
      
    )
  
)



# ==============================================================================
# FINISHED
# ==============================================================================

cat(
  
  "\n============================================================\n",
  
  "RCS + SHAP SUBGROUP ANALYSES COMPLETED\n",
  
  "============================================================\n",
  
  "\nOutput folder:\n",
  
  OUTPUT_DIR,
  
  "\n\nSHAP U outputs:\n",
  
  file.path(
    SHAP_U_DIR,
    "SHAP_U_Sex.pdf"
  ),
  
  "\n",
  
  file.path(
    SHAP_U_DIR,
    "SHAP_U_Age.pdf"
  ),
  
  "\n",
  
  file.path(
    SHAP_U_DIR,
    "SHAP_U_Subgroup_Results.xlsx"
  ),
  
  "\n\nRCS outputs:\n",
  
  RCS_DIR,
  
  "\n"
  
)
