# Models to find association between various characteristics and number of substitutions as well as (non)synonymous changes

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(pscl); library(MASS); library(pscl)

#read in data
all_without_ct_changes_first_last_placements <- readRDS("")
all_without_ct_changes_first_last_placements$num_vacc_at_infection <- as.numeric(as.character(all_without_ct_changes_first_last_placements$num_vacc_at_infection))
first_last_dnds_all_without_ct <- readRDS(file= "")

# Code for number of substitutions -------------

#Zero inflated Poisson model
zinb_model_count <- zeroinfl(total_change_count ~ total_days_between_sequences + KOEN + num_vacc_at_infection + age_group | 1, 
                       data = all_without_ct_changes_first_last_placements,
                       dist = "poisson")
summary(zinb_model_count)


#  Code for number of (non)synonymous changes -------------
first_last_dnds_all_without_ct <- readRDS(file = "")
first_last_dnds_all_without_ct$num_vacc_at_infection <- as.numeric(first_last_dnds_all_without_ct$num_vacc_at_infection)

# Stepwise selection of variables
data <- first_last_dnds_all_without_ct
outcomes <- c("nonsynonymous_count", "synonymous_count")
fixed_var <- "days_between_samples"
optional_vars <- c("KOEN", "age_group", "major_scorpio_call","charlson.index.5yrs_prepandemic", "num_vacc_at_infection")
optional_vars_contemporary <- c("KOEN", "age_group", "major_scorpio_call","charlson.index.5yrs")
diagnoses <- c("immunosuppression_diagnosis", 
               "autoimmune_diagnosis", 
               "diabetes_diagnosis", 
               "IMID_diagnosis")

#Baseline models ------

# Using contemporary diagnoses -----

# Baseline Model for Non-Synonymous Changes with Stepwise Selection
formula_nonsyn_baseline <- as.formula(paste("nonsynonymous_count ~", paste(c(fixed_var, optional_vars), collapse = " + "), "| 1"))
model_full_nonsyn_baseline <- zeroinfl(formula_nonsyn_baseline, data = data, dist = "negbin")
stepwise_model_nonsyn_baseline <- stepAIC(model_full_nonsyn_baseline, direction = "both", trace = TRUE)
summary(stepwise_model_nonsyn_baseline)

model_final_nonsyn_baseline_contemporary <- zeroinfl(nonsynonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs | 1,
                                                     data = data, dist = "negbin")
summary(model_final_nonsyn_baseline_contemporary)

# Baseline Model for Synonymous Changes with Stepwise Selection
model_full_syn_baseline_contemporary <- zeroinfl(synonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs | 1,
                                                 data = data, dist = "negbin")
summary(model_full_syn_baseline_contemporary)


# Supplementary: Using prepandemic diagnoses -----
model_final_nonsyn_baseline <- zeroinfl(nonsynonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs_prepandemic + num_vacc_at_infection | 1,
                                        data = data, dist = "negbin")
summary(model_final_nonsyn_baseline)

# Baseline Model for Synonymous Changes
model_full_syn_baseline <- zeroinfl(synonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs_prepandemic | 1,
                                    data = data, dist = "negbin")
summary(model_full_syn_baseline)






# Association between diagnostic groups and number of (non)synonymous changes ------
baseline_variables <- c("age_group", "major_scorpio_call")

# Using contemporary diagnoses ------
baseline_variables <- c("age_group", "major_scorpio_call")

# SIG! Non-Synonymous Changes with Immunosuppression (fixed)
formula_nonsyn_imm_contemporary <- as.formula(paste("nonsynonymous_count ~", 
                                       paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                       "+ immunosuppression_diagnosis | 1"))
model_full_nonsyn_imm_contemporary <- zeroinfl(formula_nonsyn_imm_contemporary, data = data, dist = "negbin")
summary(model_full_nonsyn_imm_contemporary)

# Non-Synonymous Changes with Autoimmunity (fixed)
formula_nonsyn_auto_contemporary <- as.formula(paste("nonsynonymous_count ~", 
                                        paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                        "+ autoimmune_diagnosis | 1"))
model_full_nonsyn_auto_contemporary <- zeroinfl(formula_nonsyn_auto_contemporary, data = data, dist = "negbin")
summary(model_full_nonsyn_auto_contemporary)

# Non-Synonymous Changes with Diabetes (fixed)
formula_nonsyn_diabetes_contemporary <- as.formula(paste("nonsynonymous_count ~", 
                                            paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                            "+ diabetes_diagnosis | 1"))
model_full_nonsyn_diabetes_contemporary <- zeroinfl(formula_nonsyn_diabetes_contemporary, data = data, dist = "negbin")
summary(model_full_nonsyn_diabetes_contemporary)

# Non-Synonymous Changes with IMID (fixed)
formula_nonsyn_imid_contemporary <- as.formula(paste("nonsynonymous_count ~", 
                                        paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                        "+ IMID_diagnosis | 1"))
model_full_nonsyn_imid_contemporary <- zeroinfl(formula_nonsyn_imid_contemporary, data = data, dist = "negbin")
summary(model_full_nonsyn_imid_contemporary)

# Synonymous Changes with Immunosuppression (fixed)
formula_syn_imm_contemporary <- as.formula(paste("synonymous_count ~", 
                                    paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                    "+ immunosuppression_diagnosis | 1"))
model_full_syn_imm_contemporary <- zeroinfl(formula_syn_imm_contemporary, data = data, dist = "negbin")
summary(model_full_syn_imm_contemporary)

# Synonymous Changes with Autoimmunity (fixed)
formula_syn_auto_contemporary <- as.formula(paste("synonymous_count ~", 
                                     paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                     "+ autoimmune_diagnosis | 1"))
model_full_syn_auto_contemporary <- zeroinfl(formula_syn_auto_contemporary, data = data, dist = "negbin")
summary(model_full_syn_auto_contemporary)

# Synonymous Changes with Diabetes (fixed)
formula_syn_diabetes_contemporary <- as.formula(paste("synonymous_count ~", 
                                         paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                         "+ diabetes_diagnosis | 1"))
model_full_syn_diabetes_contemporary <- zeroinfl(formula_syn_diabetes_contemporary, data = data, dist = "negbin")
summary(model_full_syn_diabetes_contemporary)

# Synonymous Changes with IMID (fixed)
formula_syn_imid_contemporary <- as.formula(paste("synonymous_count ~", 
                                     paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                     "+ IMID_diagnosis | 1"))
model_full_syn_imid_contemporary <- zeroinfl(formula_syn_imid_contemporary, data = data, dist = "negbin")
summary(model_full_syn_imid_contemporary)





# Supplementary: Using prepandemic diagnoses ----------
# Non-Synonymous Changes with Immunosuppression (fixed)
formula_nonsyn_imm <- as.formula(paste("nonsynonymous_count ~", 
                                       paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                       "+ immunosuppression_diagnosis_prepandemic | 1"))
model_full_nonsyn_imm <- zeroinfl(formula_nonsyn_imm, data = data, dist = "negbin")
summary(model_full_nonsyn_imm)

# Non-Synonymous Changes with Autoimmunity (fixed)
formula_nonsyn_auto <- as.formula(paste("nonsynonymous_count ~", 
                                        paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                        "+ autoimmune_diagnosis_prepandemic | 1"))
model_full_nonsyn_auto <- zeroinfl(formula_nonsyn_auto, data = data, dist = "negbin")
summary(model_full_nonsyn_auto)

# Non-Synonymous Changes with Diabetes (fixed)
formula_nonsyn_diabetes <- as.formula(paste("nonsynonymous_count ~", 
                                            paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                            "+ diabetes_diagnosis_prepandemic | 1"))
model_full_nonsyn_diabetes <- zeroinfl(formula_nonsyn_diabetes, data = data, dist = "negbin")
summary(model_full_nonsyn_diabetes)

# Non-Synonymous Changes with IMID (fixed)
formula_nonsyn_imid <- as.formula(paste("nonsynonymous_count ~", 
                                        paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                        "+ IMID_diagnosis_prepandemic | 1"))
model_full_nonsyn_imid <- zeroinfl(formula_nonsyn_imid, data = data, dist = "negbin")
summary(model_full_nonsyn_imid)

# Synonymous Changes with Immunosuppression (fixed)
formula_syn_imm <- as.formula(paste("synonymous_count ~", 
                                    paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                    "+ immunosuppression_diagnosis_prepandemic | 1"))
model_full_syn_imm <- zeroinfl(formula_syn_imm, data = data, dist = "negbin")
summary(model_full_syn_imm)

# Synonymous Changes with Autoimmunity (fixed)
formula_syn_auto <- as.formula(paste("synonymous_count ~", 
                                     paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                     "+ autoimmune_diagnosis_prepandemic | 1"))
model_full_syn_auto <- zeroinfl(formula_syn_auto, data = data, dist = "negbin")
summary(model_full_syn_auto)

# Synonymous Changes with Diabetes (fixed)
formula_syn_diabetes <- as.formula(paste("synonymous_count ~", 
                                         paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                         "+ diabetes_diagnosis_prepandemic | 1"))
model_full_syn_diabetes <- zeroinfl(formula_syn_diabetes, data = data, dist = "negbin")
summary(model_full_syn_diabetes)

# Synonymous Changes with IMID (fixed)
formula_syn_imid <- as.formula(paste("synonymous_count ~", 
                                     paste(c(fixed_var, baseline_variables), collapse = " + "), 
                                     "+ IMID_diagnosis_prepandemic | 1"))
model_full_syn_imid <- zeroinfl(formula_syn_imid, data = data, dist = "negbin")
summary(model_full_syn_imid)





# Supplementary analysis to explore the role of vaccination --------

#Vaccination and number of changes
quasi_model <- glm(
  rate_per_site ~ total_days_between_sequences + KOEN + num_vacc_at_infection + age_group, 
  family = quasipoisson(link = "log"), 
  data = all_without_ct_changes_first_last_placements
)

summary(quasi_model)

# Extract coefficients, standard errors, and p-values
coef_summary <- summary(quasi_model)$coefficients
coef_estimates <- coef_summary[, 1]  # Coefficients
std_errors <- coef_summary[, 2]      # Standard errors
p_values <- coef_summary[, 4]        # P-values
exp_coef <- exp(coef_estimates)
lower_ci <- exp(coef_estimates - 1.96 * std_errors)
upper_ci <- exp(coef_estimates + 1.96 * std_errors)

# Create a summary table
summary_table <- data.frame(
  Variable = rownames(coef_summary),
  Estimate = coef_estimates,
  `Exp(Estimate)` = exp_coef,
  `Lower 95% CI` = lower_ci,
  `Upper 95% CI` = upper_ci,
  `P-value` = p_values
)

# Print the summary table
print(summary_table)

# Summary statistics
non_zero_rate_per_site <- all_without_ct_changes_first_last_placements$rate_per_site[all_without_ct_changes_first_last_placements$rate_per_site > 0]
summary_stats <- data.frame(
  mean = mean(non_zero_rate_per_site),
  median = median(non_zero_rate_per_site),
  IQR = IQR(non_zero_rate_per_site),
  range_min = min(non_zero_rate_per_site),
  range_max = max(non_zero_rate_per_site),
  q95 = quantile(non_zero_rate_per_site, 0.95),
  q05 = quantile(non_zero_rate_per_site, 0.05)
)
print(summary_stats)


#Vaccination analysis for (non)synonymous substitutions
first_last_dnds_all_without_ct <- readRDS(file = "")
first_last_dnds_all_without_ct$num_vacc_at_infection <- as.numeric(first_last_dnds_all_without_ct$num_vacc_at_infection)

vaccination_model_final_nonsyn_baseline <- zeroinfl(nonsynonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs_prepandemic + num_vacc_at_infection | 1,
                                                    data = first_last_dnds_all_without_ct, dist = "negbin")
summary(vaccination_model_final_nonsyn_baseline)

vaccination_model_full_syn_baseline <- zeroinfl(synonymous_count ~ days_between_samples + age_group + major_scorpio_call + charlson.index.5yrs_prepandemic + num_vacc_at_infection | 1,
                                                data = first_last_dnds_all_without_ct, dist = "negbin")
summary(vaccination_model_full_syn_baseline)



