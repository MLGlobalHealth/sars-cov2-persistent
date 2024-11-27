# Scripts for running logistic regression for case-control study ---------
# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(broom); library(mediation); library(brms); library(RcppEigen); library(ggplot2); library(brms);
library(bayestestR); library(mediation)

# Reading dataframes from RDS files
subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson <- readRDS(file="")
subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson$num_vacc_at_infection <- as.numeric(as.character(subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson$num_vacc_at_infection))

# BRMS functions ---------
run_bayesian_logistic_regression_fixed_data <- function(predictors, 
                                                        data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                                        outcome = "outcome", 
                                                        chains = 4, iter = 2000, warmup = 1000, seed = 123) {
  
  # Construct the formula for the model using the specified predictors
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  
  # Define prior for the model
  prior <- set_prior("normal(0, 1)", class = "b")
  
  # Fit the Bayesian logistic regression model
  model <- brm(
    formula = formula,
    data = data,
    family = bernoulli(link = "logit"),
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed
  )
  
  return(model)
}

get_bayesian_odds_ratios <- function(brm_model, ci_level = 0.95) {
  # Calculate the credible interval probabilities
  ci_lower <- (1 - ci_level) / 2
  ci_upper <- 1 - ci_lower
  
  # Get the posterior summary
  summary_brms <- posterior_summary(brm_model, robust = TRUE, probs = c(ci_lower, ci_upper))
  
  # Exponentiate the estimate and credible intervals
  summary_brms_exp <- summary_brms
  summary_brms_exp[, "Estimate"] <- exp(summary_brms_exp[, "Estimate"])
  summary_brms_exp[, paste0("Q", round(ci_lower * 100, 2))] <- exp(summary_brms_exp[, paste0("Q", round(ci_lower * 100, 2))])
  summary_brms_exp[, paste0("Q", round(ci_upper * 100, 2))] <- exp(summary_brms_exp[, paste0("Q", round(ci_upper * 100, 2))])
  
  # Rename columns for clarity
  colnames(summary_brms_exp)[colnames(summary_brms_exp) == "Estimate"] <- "Adjusted Odds Ratio"
  colnames(summary_brms_exp)[colnames(summary_brms_exp) == paste0("Q", round(ci_lower * 100, 2))] <- paste0(ci_level * 100, "% CrI (Lower)")
  colnames(summary_brms_exp)[colnames(summary_brms_exp) == paste0("Q", round(ci_upper * 100, 2))] <- paste0(ci_level * 100, "% CrI (Upper)")
  
  # Return the result as a data frame
  return(as.data.frame(summary_brms_exp))
}



# Main logistic regression analyses -------------------

# Baseline: Age group, sex, major variant, charlson, number of vaccinations at infection
predictors_model_geo <- c("age_group", "KOEN", "charlson.index.5yrs_prepandemic", "num_vacc_at_infection")
brm_model_geo <- run_bayesian_logistic_regression_fixed_data(predictors)
brm_model_geo_output <- get_bayesian_odds_ratios(brm_model_geo)
brm_model_geo_output
brm_model_geo$data %>%
  group_by(outcome) %>%  # Group by the outcome
  summarize(
    median_charlson = median(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Median value
    lower_IQR = quantile(charlson.index.5yrs_prepandemic, 0.25, na.rm = TRUE),  # Lower quartile (Q1)
    upper_IQR = quantile(charlson.index.5yrs_prepandemic, 0.75, na.rm = TRUE),  # Upper quartile (Q3)
    mean_charlson = mean(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Mean value
    range_min = min(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Minimum value (start of range)
    range_max = max(charlson.index.5yrs_prepandemic, na.rm = TRUE)   # Maximum value (end of range)
  )

#Univariate
predictors_model_geo_uni <- c("charlson.index.5yrs_prepandemic")
brm_model_geo_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_uni)
brm_model_geo_output_uni <- get_bayesian_odds_ratios(brm_model_geo_uni)
brm_model_geo_output_uni
brm_model_geo_uni$data %>%
  group_by(outcome) %>%  # Group by the outcome
  summarize(
    median_charlson = median(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Median value
    lower_IQR = quantile(charlson.index.5yrs_prepandemic, 0.25, na.rm = TRUE),  # Lower quartile (Q1)
    upper_IQR = quantile(charlson.index.5yrs_prepandemic, 0.75, na.rm = TRUE),  # Upper quartile (Q3)
    mean_charlson = mean(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Mean value
    range_min = min(charlson.index.5yrs_prepandemic, na.rm = TRUE),  # Minimum value (start of range)
    range_max = max(charlson.index.5yrs_prepandemic, na.rm = TRUE)   # Maximum value (end of range)
  )

#Age
glm_model_geo_age_uni <- run_bayesian_logistic_regression_fixed_data("age_group")
glm_model_geo_age_uni_output <- get_bayesian_odds_ratios(glm_model_geo_age_uni)
glm_model_geo_age_uni_output
table(glm_model_geo_age_uni$data$age_group, glm_model_geo_age_uni$data$outcome)

#Sex
glm_model_geo_sex_uni <- run_bayesian_logistic_regression_fixed_data("KOEN")
glm_model_geo_sex_uni_output <- get_bayesian_odds_ratios(glm_model_geo_sex_uni)
glm_model_geo_sex_uni_output
table(glm_model_geo_sex_uni$data$KOEN, glm_model_geo_sex_uni$data$outcome)

#Vacc
glm_model_geo_vacc_uni <- run_bayesian_logistic_regression_fixed_data("num_vacc_at_infection")
glm_model_geo_vacc_uni_output <- get_bayesian_odds_ratios(glm_model_geo_vacc_uni)
glm_model_geo_vacc_uni_output
glm_model_geo_vacc_uni$data %>%
  group_by(outcome) %>%  # Group by the outcome
  summarize(
    median_vacc = median(num_vacc_at_infection, na.rm = TRUE),  # Median value
    lower_IQR = quantile(num_vacc_at_infection, 0.25, na.rm = TRUE),  # Lower quartile (Q1)
    upper_IQR = quantile(num_vacc_at_infection, 0.75, na.rm = TRUE),  # Upper quartile (Q3)
    mean_charlson = mean(num_vacc_at_infection, na.rm = TRUE),  # Mean value
    range_min = min(num_vacc_at_infection, na.rm = TRUE),  # Minimum value (start of range)
    range_max = max(num_vacc_at_infection, na.rm = TRUE)   # Maximum value (end of range)
  )




# DIABETES: Age group, sex, major variant, diabetes, number of vaccinations at infection
predictors_model_geo_diabetes <- c("age_group", "KOEN", "diabetes_diagnosis_prepandemic", "num_vacc_at_infection")
brm_model_geo_diabetes <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes)
brm_model_geo_diabetes_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes)
brm_model_geo_diabetes_output

predictors_model_geo_diabetes_uni <- c("diabetes_diagnosis_prepandemic")
brm_model_geo_diabetes_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes_uni)
brm_model_geo_diabetes_uni_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes_uni)
brm_model_geo_diabetes_uni_output
contingency_table <- table(brm_model_geo_diabetes$data$diabetes_diagnosis_prepandemic, 
                           brm_model_geo_diabetes$data$outcome)
contingency_table

# IMMUNOSUPPRESSION
predictors_model_geo_immunosuppression <- c("age_group", "immunosuppression_diagnosis_prepandemic", "KOEN", "num_vacc_at_infection")
brm_model_geo_immunosuppression <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression)
brm_model_geo_immunosuppression_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression)
print(brm_model_geo_immunosuppression_output)

predictors_model_geo_immunosuppression_uni <- c("immunosuppression_diagnosis_prepandemic")
brm_model_geo_immunosuppression_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression_uni)
brm_model_geo_immunosuppression_uni_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression_uni)
print(brm_model_geo_immunosuppression_uni_output)

contingency_table_immunosuppression <- table(brm_model_geo_immunosuppression$data$immunosuppression_diagnosis_prepandemic, 
                                             brm_model_geo_immunosuppression$data$outcome)
print(contingency_table_immunosuppression)

# AUTOIMMUNE
predictors_model_geo_autoimmune <- c("age_group", "KOEN", "autoimmune_diagnosis_prepandemic", "num_vacc_at_infection")
brm_model_geo_autoimmune <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune)
brm_model_geo_autoimmune_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune)
print(brm_model_geo_autoimmune_output)

predictors_model_geo_autoimmune_uni <- c("autoimmune_diagnosis_prepandemic")
brm_model_geo_autoimmune_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune_uni)
brm_model_geo_autoimmune_uni_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune_uni)
print(brm_model_geo_autoimmune_uni_output)

contingency_table_autoimmune <- table(brm_model_geo_autoimmune$data$autoimmune_diagnosis_prepandemic, 
                                      brm_model_geo_autoimmune$data$outcome)
print(contingency_table_autoimmune)

# IMID
predictors_model_geo_IMID <- c("age_group", "KOEN", "IMID_diagnosis_prepandemic", "num_vacc_at_infection")
brm_model_geo_IMID <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID)
brm_model_geo_IMID_output <- get_bayesian_odds_ratios(brm_model_geo_IMID)
print(brm_model_geo_IMID_output)

predictors_model_geo_IMID_uni <- c("IMID_diagnosis_prepandemic")
brm_model_geo_IMID_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID_uni)
brm_model_geo_IMID_uni_output <- get_bayesian_odds_ratios(brm_model_geo_IMID_uni)
print(brm_model_geo_IMID_uni_output)

contingency_table_IMID <- table(brm_model_geo_IMID$data$IMID_diagnosis_prepandemic, 
                                brm_model_geo_IMID$data$outcome)
print(contingency_table_IMID)








# Sensitivity analyses -------------------

# Calculating Charlson score for previous 10 years
predictors_model_geo_10yrs <- c("age_group", "KOEN", "charlson.index.10yrs_prepandemic", "num_vacc_at_infection")
brm_model_geo_10yrs <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_10yrs)
brm_model_geo_10yrs_output <- get_bayesian_odds_ratios(brm_model_geo_10yrs)
print(brm_model_geo_10yrs_output)

summary_10yrs <- brm_model_geo_10yrs$data %>%
  group_by(outcome) %>% 
  summarize(
    median_charlson = median(charlson.index.10yrs_prepandemic, na.rm = TRUE), 
    lower_IQR = quantile(charlson.index.10yrs_prepandemic, 0.25, na.rm = TRUE), 
    upper_IQR = quantile(charlson.index.10yrs_prepandemic, 0.75, na.rm = TRUE), 
    mean_charlson = mean(charlson.index.10yrs_prepandemic, na.rm = TRUE), 
    range_min = min(charlson.index.10yrs_prepandemic, na.rm = TRUE), 
    range_max = max(charlson.index.10yrs_prepandemic, na.rm = TRUE)
  )
print(summary_10yrs)

# Calculating Charlson score for contemporary diagnoses
predictors_model_geo_5yrs <- c("age_group", "KOEN", "charlson.index.5yrs", "num_vacc_at_infection")
brm_model_geo_5yrs <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_5yrs)
brm_model_geo_5yrs_output <- get_bayesian_odds_ratios(brm_model_geo_5yrs)
print(brm_model_geo_5yrs_output)

summary_5yrs <- brm_model_geo_5yrs$data %>%
  group_by(outcome) %>% 
  summarize(
    median_charlson = median(charlson.index.5yrs, na.rm = TRUE), 
    lower_IQR = quantile(charlson.index.5yrs, 0.25, na.rm = TRUE), 
    upper_IQR = quantile(charlson.index.5yrs, 0.75, na.rm = TRUE), 
    mean_charlson = mean(charlson.index.5yrs, na.rm = TRUE), 
    range_min = min(charlson.index.5yrs, na.rm = TRUE), 
    range_max = max(charlson.index.5yrs, na.rm = TRUE)
  )
print(summary_5yrs)









# Mediation analysis for role of vaccinations ------
library(effectsize)

predictors_model_geo_immunosuppression <- c("age_group", "immunosuppression_diagnosis_prepandemic", "KOEN", "num_vacc_at_infection")

f1 <- bf(num_vacc_at_infection ~ immunosuppression_diagnosis_prepandemic + age_group + KOEN)

# Define the Bayesian outcome model (outcome) - specify as logistic
f2 <- bf(outcome ~ immunosuppression_diagnosis_prepandemic + age_group + KOEN + num_vacc_at_infection,
         family = bernoulli(link = "logit"))

standardized_data <- subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson %>%
  mutate(across(c(immunosuppression_diagnosis_prepandemic, age_group, KOEN), 
                ~ standardize(.x)))

# Fit the model using the standardized data
brm_model_geo_immunosuppression <- brm(
  f1 + f2 + set_rescor(FALSE), 
  data = standardized_data, 
  refresh = 0
)
params <- insight::get_parameters(brm_model_geo_immunosuppression)

mediation_results <- bayestestR::mediation(
  brm_model_geo_immunosuppression,
  treatment = "immunosuppression_diagnosis_prepandemic.L", # Replace with your treatment variable
  mediator = "num_vacc_at_infection"
)
mediation_results







# Using both pre- and post-pandemic ICD-10 diagnoses --------

# DIABETES: Age group, sex, major variant, diabetes, number of vaccinations at infection
predictors_model_geo_diabetes <- c("age_group", "KOEN", "diabetes_diagnosis", "num_vacc_at_infection")
brm_model_geo_diabetes <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes)
brm_model_geo_diabetes_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes)
print(brm_model_geo_diabetes_output)

brm_model_geo_diabetes_uni <- run_bayesian_logistic_regression_fixed_data(c("diabetes_diagnosis"))
brm_model_geo_diabetes_uni_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes_uni)
print(brm_model_geo_diabetes_uni_output)

table(brm_model_geo_diabetes$data$diabetes_diagnosis, brm_model_geo_diabetes$data$outcome)

# IMMUNOSUPPRESSION: Age group, sex, immunosuppression, number of vaccinations at infection
predictors_model_geo_immunosuppression <- c("age_group", "KOEN", "immunosuppression_diagnosis", "num_vacc_at_infection")
brm_model_geo_immunosuppression <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression)
brm_model_geo_immunosuppression_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression)
print(brm_model_geo_immunosuppression_output)

brm_model_geo_immunosuppression_uni <- run_bayesian_logistic_regression_fixed_data(c("immunosuppression_diagnosis"))
brm_model_geo_immunosuppression_uni_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression_uni)
print(brm_model_geo_immunosuppression_uni_output)

table(brm_model_geo_immunosuppression$data$immunosuppression_diagnosis, brm_model_geo_immunosuppression$data$outcome)

# AUTOIMMUNE: Age group, sex, autoimmune, number of vaccinations at infection
predictors_model_geo_autoimmune <- c("age_group", "KOEN", "autoimmune_diagnosis", "num_vacc_at_infection")
brm_model_geo_autoimmune <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune)
brm_model_geo_autoimmune_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune)
print(brm_model_geo_autoimmune_output)

brm_model_geo_autoimmune_uni <- run_bayesian_logistic_regression_fixed_data(c("autoimmune_diagnosis"))
brm_model_geo_autoimmune_uni_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune_uni)
print(brm_model_geo_autoimmune_uni_output)

table(brm_model_geo_autoimmune$data$autoimmune_diagnosis, brm_model_geo_autoimmune$data$outcome)

# IMID: Age group, sex, IMID, number of vaccinations
predictors_model_geo_IMID <- c("age_group", "KOEN", "IMID_diagnosis", "num_vacc_at_infection")
brm_model_geo_IMID <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID)
brm_model_geo_IMID_output <- get_bayesian_odds_ratios(brm_model_geo_IMID)
print(brm_model_geo_IMID_output)

brm_model_geo_IMID_uni <- run_bayesian_logistic_regression_fixed_data(c("IMID_diagnosis"))
brm_model_geo_IMID_uni_output <- get_bayesian_odds_ratios(brm_model_geo_IMID_uni)
print(brm_model_geo_IMID_uni_output)

table(brm_model_geo_IMID$data$IMID_diagnosis, brm_model_geo_IMID$data$outcome)











# Testing for interactions between age and diagnoses ---------
# DIABETES
glm_model_geo_diabetes_age_interaction <- glm(outcome ~ KOEN + age_group*diabetes_diagnosis_prepandemic + num_vacc_at_infection, 
                                              data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                              family = binomial(link = "logit"))
tidy(glm_model_geo_diabetes_age_interaction, exponentiate = TRUE, conf.int=TRUE)

# IMMUNOSUPPRESSION: Age group, age, immunosuppression, number of vaccinations at infection
glm_model_geo_immunosuppression_age_interaction <- glm(outcome ~ age_group*immunosuppression_diagnosis_prepandemic + KOEN + num_vacc_at_infection, 
                                                       data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                                       family = binomial(link = "logit"))
tidy_model_geo_immunosuppression_age_interaction <- tidy(glm_model_geo_immunosuppression_age_interaction, exponentiate = TRUE, conf.int = TRUE)
print(tidy_model_geo_immunosuppression_age_interaction)

# AUTOIMMUNE: Age group, age, autoimmune, number of vaccinations at infection
glm_model_geo_autoimmune_age_interaction <- glm(outcome ~ KOEN + age_group*autoimmune_diagnosis + num_vacc_at_infection, 
                                                data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                                family = binomial(link = "logit"))
tidy_model_geo_autoimmune_age_interaction <- tidy(glm_model_geo_autoimmune_age_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_autoimmune_age_interaction)
# IMID: Age group, age, IMID, number of vaccinations
glm_model_geo_IMID_age_interaction <- glm(outcome ~ KOEN + age_group*IMID_diagnosis + num_vacc_at_infection, 
                                          data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                          family = binomial(link = "logit"))
tidy_model_geo_IMID_age_interaction <- tidy(glm_model_geo_IMID_age_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_IMID_age_interaction)


# Testing for interactions between sex and diagnoses -------

# DIABETES: Age group, sex*diabetes, major variant, number of vaccinations at infection
glm_model_geo_diabetes_sex_interaction <- glm(outcome ~ age_group + KOEN*diabetes_diagnosis + num_vacc_at_infection, 
                              data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                              family = binomial(link = "logit"))
tidy_model_geo_diabetes_sex_interaction <- tidy(glm_model_geo_diabetes_sex_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_diabetes_sex_interaction)
# IMMUNOSUPPRESSION: Age group, sex, immunosuppression, number of vaccinations at infection
glm_model_geo_immunosuppression_sex_interaction <- glm(outcome ~ age_group + immunosuppression_diagnosis*KOEN + num_vacc_at_infection, 
                                       data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                       family = binomial(link = "logit"))
tidy_model_geo_immunosuppression_sex_interaction <- tidy(glm_model_geo_immunosuppression_sex_interaction, exponentiate = TRUE, conf.int = TRUE)
print(tidy_model_geo_immunosuppression_sex_interaction)
# AUTOIMMUNE: Age group, sex, autoimmune, number of vaccinations at infection
glm_model_geo_autoimmune_sex_interaction <- glm(outcome ~ age_group + KOEN*autoimmune_diagnosis + num_vacc_at_infection, 
                                data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                family = binomial(link = "logit"))
tidy_model_geo_autoimmune_sex_interaction <- tidy(glm_model_geo_autoimmune_sex_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_autoimmune_sex_interaction)
# IMID: Age group, sex, IMID, number of vaccinations
glm_model_geo_IMID_sex_interaction <- glm(outcome ~ age_group + KOEN*IMID_diagnosis + num_vacc_at_infection, 
                          data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                          family = binomial(link = "logit"))
tidy_model_geo_IMID_sex_interaction <- tidy(glm_model_geo_IMID_sex_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_IMID_sex_interaction)


# Testing for interactions between number of vaccinations and diagnoses -------

# DIABETES: Age group, sex, diabetes, number of vaccinations at infection
glm_model_geo_diabetes_numvacc_interaction <- glm(outcome ~ age_group + KOEN + diabetes_diagnosis*num_vacc_at_infection, 
                                              data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                              family = binomial(link = "logit"))
tidy_model_geo_diabetes_numvacc_interaction <- tidy(glm_model_geo_diabetes_numvacc_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_diabetes_numvacc_interaction)
# IMMUNOSUPPRESSION: Age group, sex, immunosuppression, number of vaccinations at infection
glm_model_geo_immunosuppression_numvacc_interaction <- glm(outcome ~ age_group + immunosuppression_diagnosis*num_vacc_at_infection + KOEN, 
                                                       data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                                       family = binomial(link = "logit"))
tidy_model_geo_immunosuppression_numvacc_interaction <- tidy(glm_model_geo_immunosuppression_numvacc_interaction, exponentiate = TRUE, conf.int = TRUE)
print(tidy_model_geo_immunosuppression_numvacc_interaction)
# AUTOIMMUNE: Age group, sex, autoimmune, number of vaccinations at infection
glm_model_geo_autoimmune_numvacc_interaction <- glm(outcome ~ age_group + KOEN + autoimmune_diagnosis*num_vacc_at_infection, 
                                                data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                                family = binomial(link = "logit"))
tidy_model_geo_autoimmune_numvacc_interaction <- tidy(glm_model_geo_autoimmune_numvacc_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_autoimmune_numvacc_interaction)
# IMID: Age group, sex, IMID, number of vaccinations
glm_model_geo_IMID_numvacc_interaction <- glm(outcome ~ age_group + KOEN + IMID_diagnosis*num_vacc_at_infection, 
                                          data = subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson, 
                                          family = binomial(link = "logit"))
tidy_model_geo_IMID_numvacc_interaction <- tidy(glm_model_geo_IMID_numvacc_interaction, exponentiate = TRUE, conf.int=TRUE)
print(tidy_model_geo_IMID_numvacc_interaction)























