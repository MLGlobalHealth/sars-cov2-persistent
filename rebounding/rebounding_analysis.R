# Rebounding analysis

# read packages
library(ape);library(phytools);library(TreeTools);library(dplyr);library(tidyverse);library(data.table);library(dbplyr);library(lubridate)
library(rlang);library(foreach);library(doParallel);library(DSTora);library(ROracle);library(DSTcolectica);library(DSTdb);library(DBI);library(parallel)
library(doParallel);library(foreach);library(ggsignif);library(Rcpp);library(purrr);library(tidyr);library(furrr);library(future);library(future.apply);library(lubridate);library(seqinr);library(adegenet)
library(ggplot2);library(viridis);library(patchwork); library(knitr);library(brms)
library(pscl); library(MASS);

# Testing, vaccination and general metadata
COVID_TEST_filtered <- readRDS(file = "")
COVID_VACC_filtered <- readRDS(file = "")
lifelines_koen_filtered <- readRDS(file = "")
lifelines_filtered <- readRDS(file = "")

# Persistent infection metadata
all_without_ct_metadata <- readRDS(file = "")

# Reading in more data
all_without_ct_changes_first_last_placements <- readRDS("")
all_without_ct_changes_first_last_placements$num_vacc_at_infection <- as.numeric(as.character(all_without_ct_changes_first_last_placements$num_vacc_at_infection))
first_last_dnds_all_without_ct <- readRDS(file= "")


# Prepare negative tests
negative_tests <- COVID_TEST_filtered %>%
  filter(SVARRESULTAT == 0, CASEDEF == "SARS2") %>%
  dplyr::select(PERSON_ID, PRDATE_ADJUSTED, STATUS) %>%
  rename(DATESAMPLING = PRDATE_ADJUSTED)

# Finding rebounders:
# Ensure DATESAMPLING and infection_episode_start_date/infection_episode_end_date are Date objects
negative_tests$DATESAMPLING <- as.Date(negative_tests$DATESAMPLING)
all_without_ct_changes_first_last_placements$infection_episode_start_date <- as.Date(all_without_ct_changes_first_last_placements$infection_episode_start_date)
all_without_ct_changes_first_last_placements$infection_episode_end_date <- as.Date(all_without_ct_changes_first_last_placements$infection_episode_end_date)

# Initialize columns for the results
all_without_ct_changes_first_last_placements$negative_test_count <- 0
all_without_ct_changes_first_last_placements$rebound <- 0
all_without_ct_changes_first_last_placements$negative_test_dates <- ""  
all_without_ct_changes_first_last_placements$negative_test_times <- ""  

# Loop through each row of all_without_ct_changes_first_last_placements
for (i in 1:nrow(all_without_ct_changes_first_last_placements)) {
  # Extract current row information
  person_id <- all_without_ct_changes_first_last_placements$PERSON_ID[i]
  start_date <- all_without_ct_changes_first_last_placements$infection_episode_start_date[i]
  end_date <- all_without_ct_changes_first_last_placements$infection_episode_end_date[i]
  
  # Define the time window for negative tests
  valid_start_date <- start_date + 2  # >= 2 days after the start date
  valid_end_date <- end_date - 2      # <= 2 days before the end date
  
  # Filter negative tests for the current person within the time window
  negative_tests_in_window <- negative_tests %>%
    filter(
      PERSON_ID == person_id,                  # Match PERSON_ID
      DATESAMPLING >= valid_start_date,        # >= 2 days after start
      DATESAMPLING <= valid_end_date           # <= 2 days before end
    )
  
  # Count the number of negative tests
  num_negative_tests <- nrow(negative_tests_in_window)
  
  # Extract the dates of negative tests
  negative_test_dates <- negative_tests_in_window$DATESAMPLING
  
  # Calculate the time in days from the start of the infection to each negative test
  time_from_start <- as.numeric(negative_test_dates - start_date)
  
  # Update columns for the current row
  all_without_ct_changes_first_last_placements$negative_test_count[i] <- num_negative_tests
  all_without_ct_changes_first_last_placements$rebound[i] <- ifelse(num_negative_tests > 0, 1, 0)
  
  # Save the dates of the negative tests (as a comma-separated string)
  all_without_ct_changes_first_last_placements$negative_test_dates[i] <- paste(negative_test_dates, collapse = ", ")
  
  # Save the times from the infection start to the negative tests (as a comma-separated string)
  all_without_ct_changes_first_last_placements$negative_test_times[i] <- paste(time_from_start, collapse = ", ")
}

# Reorganize columns to place the new ones next to the end date column
all_without_ct_changes_first_last_placements <- all_without_ct_changes_first_last_placements %>%
  relocate(
    negative_test_count, rebound, negative_test_dates, negative_test_times,
    .after = infection_episode_end_date
  )

#Counting rebounders and those with immunosuppression -------
total_rebounders <- nrow(all_without_ct_changes_first_last_placements %>% filter(rebound == 1))

# Calculate counts for each combination of immunosuppression diagnosis and rebound
summary_table <- all_without_ct_changes_first_last_placements %>%
  group_by(immunosuppression_diagnosis, rebound) %>%
  summarise(
    count = n(),
    .groups = "drop" # Avoid grouped output
  ) %>%
  mutate(
    proportion = count / sum(count) # Calculate proportion for each group
  )

kable(summary_table, caption = "Counts and Proportions of Rebounders by Immunosuppression Status")



#Statistical testing -------------


# Immunosuppresion status and likelihood of rebounding
model_rebound <- brm(
  rebound ~ immunosuppression_diagnosis,  # Predictor is immunosuppression status
  data = all_without_ct_changes_first_last_placements,
  family = bernoulli(),                   # Binary outcome
  prior = c(
    prior(normal(0, 1), class = "b"),    # Prior for regression coefficients
    prior(cauchy(0, 2), class = "Intercept")  # Prior for intercept
  ),
  iter = 4000, chains = 4, cores = 4, seed = 42
)
summary(model_rebound)
posterior_samples <- as_draws(model_rebound)
posterior_summary <- posterior_summary(model_rebound)
beta_immuno <- posterior_summary["b_immunosuppression_diagnosis.L", ]
prob_less_than_zero <- mean(as.matrix(model_rebound) %>% 
                              as.data.frame() %>% 
                              select(b_immunosuppression_diagnosis.L) < 0)
print(beta_immuno)
cat("Probability that immunosuppression_diagnosis coefficient < 0:", prob_less_than_zero, "\n")


# Immunosuppresion status and likelihood of rebounding, adjusting for time
all_without_ct_changes_first_last_placements <- all_without_ct_changes_first_last_placements %>%
  mutate(infection_duration = as.numeric(infection_episode_end_date - infection_episode_start_date))
model_rebound_adjusted <- brm(
  rebound ~ immunosuppression_diagnosis + infection_duration,
  data = all_without_ct_changes_first_last_placements,
  family = bernoulli(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(cauchy(0, 2), class = "Intercept")
  ),
  iter = 4000, chains = 4, cores = 4, seed = 42
)
summary(model_rebound_adjusted)
posterior_summary_adjusted <- posterior_summary(model_rebound)
beta_immuno_adjusted <- posterior_summary_adjusted["b_immunosuppression_diagnosis.L", ]
prob_less_than_zero_adjusted <- mean(as.matrix(model_rebound_adjusted) %>% 
                              as.data.frame() %>% 
                              select(b_immunosuppression_diagnosis.L) < 0)
print(beta_immuno_adjusted)
cat("Probability that immunosuppression_diagnosis coefficient < 0:", prob_less_than_zero_adjusted, "\n")

proportions <- all_without_ct_changes_first_last_placements %>%
  group_by(immunosuppression_diagnosis) %>%
  summarize(
    total = n(),
    rebounders = sum(rebound == 1, na.rm = TRUE),
    proportion_rebounders = rebounders / total
  )
print(proportions)



# Association between rebounding and substitution rate
quasi_model <- glm(
  rate_per_site ~ total_days_between_sequences + age_group + rebound + charlson.index.5yrs, 
  family = quasipoisson(link = "log"), 
  data = all_without_ct_changes_first_last_placements
)

summary(quasi_model)

# Extract coefficients, standard errors, and p-values
coef_summary <- summary(quasi_model)$coefficients
coef_estimates <- coef_summary[, 1]  # Coefficients
std_errors <- coef_summary[, 2]      # Standard errors
p_values <- coef_summary[, 4]        # P-values

# Exponentiate coefficients and calculate confidence intervals
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




