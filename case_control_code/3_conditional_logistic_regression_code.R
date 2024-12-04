# Scripts for running conditional logistic regression for case-control study ---------
# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(broom); library(mediation); library(brms); library(RcppEigen); library(ggplot2); library(brms);
library(bayestestR); library(mediation); library(rstanarm); library(survival); library(arm)



# Set working directory
# Reading dataframes from RDS files
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- readRDS(file="")
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$event_id <- as.factor(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$event_id)
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$vaccinated <- 
  ifelse(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$num_vacc_at_infection > 0, "yes", "no")
# Convert DATESAMPLING and vaccination dates to Date format
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$DATESAMPLING <- ymd(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$DATESAMPLING)
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$FIRST_VACCINEDATE <- ymd(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$FIRST_VACCINEDATE)
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$SECOND_VACCINEDATE <- ymd(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$SECOND_VACCINEDATE)
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$THIRD_VACCINEDATE <- ymd(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$THIRD_VACCINEDATE)
# Create a time_since_vaccination column for vaccinated individuals
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$time_since_vaccination <- NA
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson %>%
  mutate(
    time_since_vaccination = case_when(
      vaccinated == "yes" & 
        (
          DATESAMPLING >= (FIRST_VACCINEDATE + days(7)) |
            DATESAMPLING >= (SECOND_VACCINEDATE + days(7)) |
            DATESAMPLING >= (THIRD_VACCINEDATE + days(7))
        ) ~ pmin(
          # Calculate time since vaccination for each dose, considering the grace period
          ifelse(DATESAMPLING >= (FIRST_VACCINEDATE + days(7)), as.numeric(interval(FIRST_VACCINEDATE, DATESAMPLING) / days(1)), NA),
          ifelse(DATESAMPLING >= (SECOND_VACCINEDATE + days(7)), as.numeric(interval(SECOND_VACCINEDATE, DATESAMPLING) / days(1)), NA),
          ifelse(DATESAMPLING >= (THIRD_VACCINEDATE + days(7)), as.numeric(interval(THIRD_VACCINEDATE, DATESAMPLING) / days(1)), NA),
          na.rm = TRUE  # Ignore NA values when calculating the minimum
        ) - 7,  # Subtract 7 days grace period to account for the post-vaccination window
      
      TRUE ~ NA_real_  # If no vaccination or conditions aren't met, set as NA
    )
  )
#Factorize
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson %>%
  mutate(
    time_category = case_when(
      is.na(time_since_vaccination) ~ NA_character_,    # If no vaccination, set as NA
      time_since_vaccination < 84 ~ "< 12 weeks",        # Less than 12 weeks
      time_since_vaccination >= 84 & time_since_vaccination < 168 ~ "12-24 weeks", # 12-24 weeks
      time_since_vaccination >= 168 ~ "> 24 weeks",      # More than 24 weeks
      TRUE ~ NA_character_  # Any other case, should be NA
    ),
    time_category = factor(time_category, levels = c("< 12 weeks", "12-24 weeks", "> 24 weeks"))
  )
# Transform NA values in time_category to "Unknown"
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$time_category <-
  as.factor(ifelse(is.na(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$time_category), "No Vaccination", 
                   as.character(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson$time_category)))

# Create a combined category variable
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- 
  conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson %>%
  mutate(
    combined_category = case_when(
      # No vaccination
      num_vacc_at_infection == 0 ~ "No Vaccination",
      
      # One vaccination
      num_vacc_at_infection == 1 ~ "1",

      # Two or Three vaccinations
      (num_vacc_at_infection == 2 | num_vacc_at_infection == 3) & time_since_vaccination < 126 ~ "2_or_3_less_18",
      (num_vacc_at_infection == 2 | num_vacc_at_infection == 3) & time_since_vaccination >= 126 ~ "2_or_3_greater_18",
      
      # Any other case (if data is missing or invalid)
      TRUE ~ NA_character_
    ),
    # Convert to factor with levels in logical order
    combined_category = factor(
      combined_category, 
      levels = c(
        "No Vaccination", 
        "1",
        "2_or_3_less_18", "2_or_3_greater_18"
      )
    )
  )





# Conditional logistic regression functions ---------
options(mc.cores = 4)

run_bayesian_logistic_regression_fixed_data <- function(
    predictors, 
    data, 
    prior = normal(0, 1, autoscale = FALSE), 
    prior_covariance = decov(), 
    chains = 4, 
    iter = 2000, 
    warmup = 1000, 
    seed = 123,
    adapt_delta = 0.95,
    ...
) {
  # Create formula
  formula <- as.formula(paste("outcome", "~", paste(predictors, collapse = " + ")))
  
  # Fit model
  stan_clogit(
    formula = formula,
    data = conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson,
    strata = event_id,
    prior = prior,
    prior_covariance = prior_covariance,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    adapt_delta = adapt_delta,
    ...
  )
}

get_bayesian_odds_ratios <- function(clogit_model, ci_level = 0.95) {
  # Calculate the credible interval probabilities
  ci_lower <- (1 - ci_level) / 2
  ci_upper <- 1 - ci_lower
  
  # Extract posterior summary from the model
  summary_clogit <- posterior_summary(clogit_model, robust = TRUE, probs = c(ci_lower, ci_upper))
  
  # Exponentiate the estimate and credible intervals
  summary_clogit_exp <- summary_clogit
  summary_clogit_exp[, "Estimate"] <- exp(summary_clogit_exp[, "Estimate"])
  summary_clogit_exp[, paste0("Q", round(ci_lower * 100, 2))] <- exp(summary_clogit_exp[, paste0("Q", round(ci_lower * 100, 2))])
  summary_clogit_exp[, paste0("Q", round(ci_upper * 100, 2))] <- exp(summary_clogit_exp[, paste0("Q", round(ci_upper * 100, 2))])
  
  # Rename columns for clarity
  colnames(summary_clogit_exp)[colnames(summary_clogit_exp) == "Estimate"] <- "Adjusted Odds Ratio"
  colnames(summary_clogit_exp)[colnames(summary_clogit_exp) == paste0("Q", round(ci_lower * 100, 2))] <- paste0(ci_level * 100, "% CrI (Lower)")
  colnames(summary_clogit_exp)[colnames(summary_clogit_exp) == paste0("Q", round(ci_upper * 100, 2))] <- paste0(ci_level * 100, "% CrI (Upper)")
  
  # Return the result as a data frame
  return(as.data.frame(summary_clogit_exp))
}

# Main logistic regression analyses -------------------
predictors_baseline <- c("age_group", "KOEN", "charlson.index.5yrs_prepandemic", "combined_category")

# Baseline: Age group, sex, major variant, charlson, number of vaccinations at infection
brm_model_geo <- run_bayesian_logistic_regression_fixed_data(predictors_baseline)
saveRDS(brm_model_geo, file="brm_model_geo.RDS")
brm_model_geo <- readRDS(file="brm_model_geo.RDS")
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
saveRDS(brm_model_geo_uni, file="brm_model_geo_uni.RDS")
brm_model_geo_uni <- readRDS(file="brm_model_geo_uni.RDS")

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
saveRDS(glm_model_geo_age_uni, file="glm_model_geo_age_uni.RDS")
glm_model_geo_age_uni <- readRDS(file="glm_model_geo_age_uni.RDS")

glm_model_geo_age_uni_output <- get_bayesian_odds_ratios(glm_model_geo_age_uni)
glm_model_geo_age_uni_output
table(glm_model_geo_age_uni$data$age_group, glm_model_geo_age_uni$data$outcome)

#Sex
glm_model_geo_sex_uni <- run_bayesian_logistic_regression_fixed_data("KOEN")
saveRDS(glm_model_geo_sex_uni, file="glm_model_geo_sex_uni.RDS")
glm_model_geo_sex_uni <- readRDS(file="glm_model_geo_sex_uni.RDS")

glm_model_geo_sex_uni_output <- get_bayesian_odds_ratios(glm_model_geo_sex_uni)
glm_model_geo_sex_uni_output
table(glm_model_geo_sex_uni$data$KOEN, glm_model_geo_sex_uni$data$outcome)

#Vacc
glm_model_geo_vacc_uni <- run_bayesian_logistic_regression_fixed_data("combined_category")
saveRDS(glm_model_geo_vacc_uni, file="glm_model_geo_vacc_uni.RDS")
glm_model_geo_vacc_uni <- readRDS(file="glm_model_geo_vacc_uni.RDS")


glm_model_geo_vacc_uni_output <- get_bayesian_odds_ratios(glm_model_geo_vacc_uni)
glm_model_geo_vacc_uni_output
table(glm_model_geo_vacc_uni$data$combined_category, glm_model_geo_vacc_uni$data$outcome)





# DIABETES: Age group, sex, major variant, diabetes, number of vaccinations at infection
predictors_model_geo_diabetes <- c("age_group", "KOEN", "diabetes_diagnosis_prepandemic", "combined_category")
brm_model_geo_diabetes <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes)
saveRDS(brm_model_geo_diabetes, file="brm_model_geo_diabetes.RDS")

brm_model_geo_diabetes_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes)
brm_model_geo_diabetes_output

predictors_model_geo_diabetes_uni <- c("diabetes_diagnosis_prepandemic")
brm_model_geo_diabetes_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes_uni)
saveRDS(brm_model_geo_diabetes_uni, file="brm_model_geo_diabetes_uni.RDS")

brm_model_geo_diabetes_uni_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes_uni)
brm_model_geo_diabetes_uni_output
contingency_table <- table(brm_model_geo_diabetes$data$diabetes_diagnosis_prepandemic, 
                           brm_model_geo_diabetes$data$outcome)
contingency_table

# IMMUNOSUPPRESSION
predictors_model_geo_immunosuppression <- c("age_group", "immunosuppression_diagnosis_prepandemic", "KOEN", "combined_category")
brm_model_geo_immunosuppression <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression)
saveRDS(brm_model_geo_immunosuppression, file="brm_model_geo_immunosuppression.RDS")

brm_model_geo_immunosuppression_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression)
print(brm_model_geo_immunosuppression_output)

predictors_model_geo_immunosuppression_uni <- c("immunosuppression_diagnosis_prepandemic")
brm_model_geo_immunosuppression_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression_uni)
saveRDS(brm_model_geo_immunosuppression_uni, file="brm_model_geo_immunosuppression_uni.RDS")

brm_model_geo_immunosuppression_uni_output <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression_uni)
print(brm_model_geo_immunosuppression_uni_output)

contingency_table_immunosuppression <- table(brm_model_geo_immunosuppression$data$immunosuppression_diagnosis_prepandemic, 
                                             brm_model_geo_immunosuppression$data$outcome)
print(contingency_table_immunosuppression)

# AUTOIMMUNE
predictors_model_geo_autoimmune <- c("age_group", "KOEN", "autoimmune_diagnosis_prepandemic", "combined_category")
brm_model_geo_autoimmune <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune)
saveRDS(brm_model_geo_autoimmune, file="brm_model_geo_autoimmune.RDS")

brm_model_geo_autoimmune_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune)
print(brm_model_geo_autoimmune_output)

predictors_model_geo_autoimmune_uni <- c("autoimmune_diagnosis_prepandemic")
brm_model_geo_autoimmune_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune_uni)
saveRDS(brm_model_geo_autoimmune_uni, file="brm_model_geo_autoimmune_uni.RDS")

brm_model_geo_autoimmune_uni_output <- get_bayesian_odds_ratios(brm_model_geo_autoimmune_uni)
print(brm_model_geo_autoimmune_uni_output)

contingency_table_autoimmune <- table(brm_model_geo_autoimmune$data$autoimmune_diagnosis_prepandemic, 
                                      brm_model_geo_autoimmune$data$outcome)
print(contingency_table_autoimmune)

# IMID
predictors_model_geo_IMID <- c("age_group", "KOEN", "IMID_diagnosis_prepandemic", "combined_category")
brm_model_geo_IMID <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID)
saveRDS(brm_model_geo_IMID, file="brm_model_geo_IMID.RDS")

brm_model_geo_IMID_output <- get_bayesian_odds_ratios(brm_model_geo_IMID)
print(brm_model_geo_IMID_output)

predictors_model_geo_IMID_uni <- c("IMID_diagnosis_prepandemic")
brm_model_geo_IMID_uni <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID_uni)
saveRDS(brm_model_geo_IMID_uni, file="brm_model_geo_IMID_uni.RDS")

brm_model_geo_IMID_uni_output <- get_bayesian_odds_ratios(brm_model_geo_IMID_uni)
print(brm_model_geo_IMID_uni_output)

contingency_table_IMID <- table(brm_model_geo_IMID$data$IMID_diagnosis_prepandemic, 
                                brm_model_geo_IMID$data$outcome)
print(contingency_table_IMID)








# Sensitivity analyses -------------------

# Calculating Charlson score for previous 10 years
predictors_model_geo_10yrs <- c("age_group", "KOEN", "charlson.index.10yrs_prepandemic", "combined_category")
brm_model_geo_10yrs <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_10yrs)
saveRDS(brm_model_geo_10yrs, file="brm_model_geo_10yrs.RDS")

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
predictors_model_geo_5yrs <- c("age_group", "KOEN", "charlson.index.5yrs", "combined_category")
brm_model_geo_5yrs <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_5yrs)
saveRDS(brm_model_geo_5yrs, file="brm_model_geo_5yrs.RDS")

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











# Using both pre- and post-pandemic ICD-10 diagnoses --------

# DIABETES: Age group, sex, major variant, diabetes, number of vaccinations at infection
predictors_model_geo_diabetes_postpandemic <- c("age_group", "KOEN", "diabetes_diagnosis", "combined_category")
brm_model_geo_diabetes_postpandemic <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_diabetes_postpandemic)
saveRDS(brm_model_geo_diabetes_postpandemic, file="brm_model_geo_diabetes_postpandemic.RDS")

brm_model_geo_diabetes_postpandemic_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes_postpandemic)
print(brm_model_geo_diabetes_postpandemic_output)

brm_model_geo_diabetes_uni_postpandemic <- run_bayesian_logistic_regression_fixed_data(c("diabetes_diagnosis"))
saveRDS(brm_model_geo_diabetes_uni_postpandemic, file="brm_model_geo_diabetes_uni_postpandemic.RDS")

brm_model_geo_diabetes_uni_output <- get_bayesian_odds_ratios(brm_model_geo_diabetes_uni_postpandemic)
print(brm_model_geo_diabetes_uni_postpandemic)

table(brm_model_geo_diabetes_postpandemic$data$diabetes_diagnosis, brm_model_geo_diabetes_postpandemic$data$outcome)

# IMMUNOSUPPRESSION: Age group, sex, immunosuppression, number of vaccinations at infection
predictors_model_geo_immunosuppression_postpandemic <- c("age_group", "KOEN", "immunosuppression_diagnosis", "combined_category")
brm_model_geo_immunosuppression_postpandemic <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_immunosuppression_postpandemic)
saveRDS(brm_model_geo_immunosuppression_postpandemic, file="brm_model_geo_immunosuppression_postpandemic.RDS")

brm_model_geo_immunosuppression_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression_postpandemic)
print(brm_model_geo_immunosuppression_output_postpandemic)

brm_model_geo_immunosuppression_uni_postpandemic <- run_bayesian_logistic_regression_fixed_data(c("immunosuppression_diagnosis"))
saveRDS(brm_model_geo_immunosuppression_uni_postpandemic, file="brm_model_geo_immunosuppression_uni_postpandemic.RDS")

brm_model_geo_immunosuppression_uni_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_immunosuppression_uni_postpandemic)
print(brm_model_geo_immunosuppression_uni_output_postpandemic)

table(brm_model_geo_immunosuppression_postpandemic$data$immunosuppression_diagnosis, brm_model_geo_immunosuppression_postpandemic$data$outcome)

# AUTOIMMUNE: Age group, sex, autoimmune, number of vaccinations at infection
predictors_model_geo_autoimmune_postpandemic <- c("age_group", "KOEN", "autoimmune_diagnosis", "combined_category")
brm_model_geo_autoimmune_postpandemic <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_autoimmune_postpandemic)
saveRDS(brm_model_geo_autoimmune_postpandemic, file="brm_model_geo_autoimmune_postpandemic.RDS")

brm_model_geo_autoimmune_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_autoimmune_postpandemic)
print(brm_model_geo_autoimmune_output_postpandemic)

brm_model_geo_autoimmune_uni_postpandemic <- run_bayesian_logistic_regression_fixed_data(c("autoimmune_diagnosis"))
saveRDS(brm_model_geo_autoimmune_uni_postpandemic, file="brm_model_geo_autoimmune_uni_postpandemic.RDS")

brm_model_geo_autoimmune_uni_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_autoimmune_uni_postpandemic)
print(brm_model_geo_autoimmune_uni_output_postpandemic)

table(brm_model_geo_autoimmune_postpandemic$data$autoimmune_diagnosis, brm_model_geo_autoimmune_postpandemic$data$outcome)

# IMID: Age group, sex, IMID, number of vaccinations
predictors_model_geo_IMID_postpandemic <- c("age_group", "KOEN", "IMID_diagnosis", "combined_category")
brm_model_geo_IMID_postpandemic <- run_bayesian_logistic_regression_fixed_data(predictors_model_geo_IMID_postpandemic)
saveRDS(brm_model_geo_IMID_postpandemic, file="brm_model_geo_IMID_postpandemic.RDS")

brm_model_geo_IMID_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_IMID_postpandemic)
print(brm_model_geo_IMID_output_postpandemic)

brm_model_geo_IMID_uni_postpandemic <- run_bayesian_logistic_regression_fixed_data(c("IMID_diagnosis"))
saveRDS(brm_model_geo_IMID_uni_postpandemic, file="brm_model_geo_IMID_uni_postpandemic.RDS")

brm_model_geo_IMID_uni_output_postpandemic <- get_bayesian_odds_ratios(brm_model_geo_IMID_uni_postpandemic)
print(brm_model_geo_IMID_uni_output_postpandemic)

table(brm_model_geo_IMID_postpandemic$data$IMID_diagnosis, brm_model_geo_IMID_postpandemic$data$outcome)


# Regression analyses using all diagnoses  -------------------

predictors_all_diagnoses <- c("age_group", "KOEN", "combined_category", "autoimmune_diagnosis_prepandemic", "immunosuppression_diagnosis_prepandemic", 
                              "diabetes_diagnosis_prepandemic", "IMID_diagnosis_prepandemic")
brm_model_all_diagnoses <- run_bayesian_logistic_regression_fixed_data(predictors_all_diagnoses)
saveRDS(brm_model_all_diagnoses, file="brm_model_all_diagnoses.RDS")
brm_model_all_diagnoses_output <- get_bayesian_odds_ratios(brm_model_all_diagnoses)
print(brm_model_all_diagnoses_output)
