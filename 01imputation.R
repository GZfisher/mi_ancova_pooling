# Load the R workspace file (assumed to contain 'wide_data' and related objects)
# load("dummy.RData")

# Load necessary libraries for imputation and data manipulation
library(norm)         # For normal-based multiple imputation models
library(dplyr)        # For data manipulation (pipes, mutate, etc.)
library(tidyverse)    # Includes data wrangling and visualization tools

# Source custom R functions (user-defined; step-wise imputation and Fortran routines)
source("function/two_step_imputation.R")
source("function/Fortran_function.R")

# Add an indicator variable for imputation number (0 = original data) to the original data
ori_data <- wide_data %>% mutate(impno = 0)

# Function to convert imputed and original (with missing values) data from wide to long format,
# and compute analysis variables and a unique subject-visit ID
convert_long <- function(imputed_data, original_data, visit_prefix = "VISIT", subj_col = "SUBJID",
                        base_col = "BASE", target_var = "CHG", imp_col = "impno") {
  # Combine original data (with missing values) and imputed data
  combined <- rbind(original_data, imputed_data)
  
  # Pivot to long format: each row is a subject-visit combination
  long <- combined %>%
    pivot_longer(
      cols = starts_with(visit_prefix), # Visit columns (e.g., VISIT1, VISIT2, ...)
      names_to = "AVISIT",              # Name for visit variable
      values_to = target_var            # Name for change-from-baseline (or other measure)
    ) %>%
    mutate(
      AVAL = !!sym(base_col) + !!sym(target_var),                 # Calculate analysis value
      id = paste0(!!sym(subj_col), "_", AVISIT),           # Create unique ID per subject-visit
      AVISITN = gsub(visit_prefix, "", AVISIT)
    )
  return(long)
}

# Select required columns and convert to numeric matrix for imputation
# Includes covariates (TRT01PN, REGIONN, BMBLIG1N, BASE) and visit variables
adqs_mono <- data.matrix(wide_data[,c("TRT01PN","REGIONN","BLBMIG1N","BASE",colnames(wide_data)[8:35])])

# First step imputation: user-defined function, generates initial imputed datasets
mono100 <- step1(adqs_mono, 100, 200, 100, seed = 13141)

# Merge back key identifying variables to the imputed data
mono100_added <- cbind(wide_data[,c("SUBJID","TRT01P","DCTFL")], mono100)  %>% 
  mutate(TRT01PN = as.factor(TRT01PN),
         REGIONN = as.factor(REGIONN),
         BLBMIG1N = as.factor(BLBMIG1N))

# Define covariate columns for modeling
cov_cols <- c('TRT01PN','REGIONN','BLBMIG1N','BASE')
# Identify visit variables that need imputation or modeling
visit_cols = colnames(mono100)[5:ncol(mono100)]

# Initialize model formula list for each step
formula_list <- c()

# Dynamically construct the formula for each visit variable, 
# each time adding the current visit as a predictor in subsequent formulas
for (i in 1:length(visit_cols)) {
  now_visit_col <- visit_cols[i]
  
  # If the current visit column has missing values:
  if (any(is.na(mono100[[now_visit_col]]))) {
    response <- now_visit_col
    formula_list <- c(
      formula_list, 
      as.formula(paste0(response,"~",paste(cov_cols,collapse = '+')))
    )
    # Add this column to covariates for the next steps
    cov_cols <- c(cov_cols, now_visit_col)
  }
  # If the column is complete, just add it to covariates
  else {
    cov_cols <- c(cov_cols, now_visit_col)
  }
}

# Second step imputation: conditionally impute using formulas constructed above
complete100 <- step2(mono100_added, 100, 'norm', formula_list, seed = 23424)

# Example usage:
# original_data should be the initial dataset containing missing values
# imputed_data is the completed data after imputation steps
complete_long <- convert_long(imputed_data = complete100, original_data = ori_data,
                         visit_prefix = "WEEK", subj_col = "SUBJID", base_col = "BASE", target_var = "CHG")

# Save final imputed/long data for further use (commented out here)
# save(complete_long, complete100, mono100_added, file = "imputation.RData")