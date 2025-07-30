# rm(list = ls())
load("imputation.RData")

library(dplyr)
# Set reference levels for treatment variables as required for analysis
complete_long$TRT01PN <- as.character(complete_long$TRT01PN)
complete_long$TRT01P <- as.character(complete_long$TRT01P)

# Convert to 'mids' object for multiply imputed analysis
complete_imp <- mice::as.mids(complete_long, .imp = "impno", .id = "id")

covariates <- c("TRT01PN", "BASE", "BLBMIG1N","REGIONN")
covariates_part <- paste(covariates, collapse = " + ")
formula <- paste0("CHG", " ~ ", covariates_part)

ancova_res <- function(imp, formula, trt_var, ref) {
  lm_fit <- with(data=imp,exp=stats::lm(
    formula = stats::as.formula(formula)
  )) # Fit model based on each imputation
  emmeans_fit <- emmeans::emmeans(
    lm_fit,
    # Specify here the group variable over which EMM are desired.
    specs = trt_var,
    weights = "proportional"
  )
  emmeans_contrasts <- emmeans::contrast(
    emmeans_fit,
    # Compare dummy arms versus the control arm.
    method = "trt.vs.ctrl",
    # Take the arm factor from param "ref" as the control arm.
    ref = ref,
    level = 0.95
  )
  sum_contrasts <- summary(
    emmeans_contrasts,
    # Derive confidence intervals, t-tests and p-values.
    infer = TRUE,
    # Do not adjust the p-values for multiplicity.
    adjust = "none"
  )
  return(list(contrasts=sum_contrasts, lm_fit=lm_fit))
}


for (vv in seq(2,56,2)) {
  res_r100 <- ancova_res(complete_imp %>% filter(AVISIT==paste0("WEEK",vv)), formula = formula, 
                         trt_var = "TRT01PN", ref = 2)$contrasts %>% mutate(n=100, source='R', visit=vv)
  res_temp <- res_r100
  if (vv==2) {
    res_com100 <- res_temp
  } else {
    res_com100 <- rbind(res_com100, res_temp)
  }
  print(vv)
}
diff_r <- res_com100
save(diff_r, file = "diff_r.RData")
