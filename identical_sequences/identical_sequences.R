# Script to identify sequences with 0, 1 and 2 differences compared to persistent infections within 11 days

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(broom); library(mediation); library(brms); library(RcppEigen); library(ggplot2); library(brms);
library(bayestestR); library(mediation); library(rstanarm); library(survival); library(arm); library(dplyr)

# Read in relevant data
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- readRDS(file="")
conditional_subset <- conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson
first_last_dnds_all_without_ct <- readRDS(file= "")

#Read more data in
final_merged_new_id_QC_unique_metadata <- read.csv(file="")

#Getting relevant sequences
result <- first_last_dnds_all_without_ct %>%
  filter(nonsynonymous_count >= 1) %>%  # Only keep rows with nonsynonymous_count >= 1
  group_by(pair_id, PERSON_ID) %>%
  slice(1) %>% # Select the first row within each group
  ungroup()
#128 rows
result <- result[, c("pair_id", "PERSON_ID", "last_strain", "last_DATESAMPLING")]
last_strain_values <- result$last_strain
writeLines(as.character(last_strain_values), "last_strain_values.txt")

# Finding cases (persistent) and controls ----------

# Load case strains (you mentioned last_strain_values.txt)
case_strains <- readLines("last_strain_values.txt")

# Filter the 'cases' dataframe
cases <- final_merged_new_id_QC_unique_metadata[
  final_merged_new_id_QC_unique_metadata$strain %in% case_strains, 
]

# Ensure DATESAMPLING is in Date format (yyyy-mm-dd)
final_merged_new_id_QC_unique_metadata$DATESAMPLING <- ymd(final_merged_new_id_QC_unique_metadata$DATESAMPLING)
cases$DATESAMPLING <- ymd(cases$DATESAMPLING)

# Ensure that cases dataframe is not empty
if (nrow(cases) == 0) {
  stop("The cases dataframe is empty, check your filtering conditions!")
}

# Assign unique case codes to cases (this will work as long as `cases` has rows)
cases$case_code <- seq_len(nrow(cases))

# Store all PERSON_ID values from the cases to avoid using them as controls
case_person_ids <- unique(cases$PERSON_ID)

# Initialize an empty list for storing controls
controls_list <- list()

# Iterate over each case
for (i in seq_len(nrow(cases))) {
  # Get the current case row
  case_row <- cases[i, ]
  
  # Debugging output to track progress
  cat("Processing case:", i, "Strain:", case_row$strain, "Date:", case_row$DATESAMPLING, "\n")
  
  # Find potential controls:
  # - Match scorpio_call
  # - Match sampling date within 11 days
  # - Exclude the case itself and exclude any rows where PERSON_ID is the same as the case's PERSON_ID
  potential_controls <- final_merged_new_id_QC_unique_metadata[
    final_merged_new_id_QC_unique_metadata$scorpio_call == case_row$scorpio_call &
      final_merged_new_id_QC_unique_metadata$DATESAMPLING >= case_row$DATESAMPLING &
      final_merged_new_id_QC_unique_metadata$DATESAMPLING <= (case_row$DATESAMPLING + days(11)) &
      final_merged_new_id_QC_unique_metadata$strain != case_row$strain &  # Exclude the case itself
      !final_merged_new_id_QC_unique_metadata$PERSON_ID %in% case_person_ids,  # Exclude same PERSON_ID
  ]
  
  # Check if potential_controls is empty, if it is, skip this iteration
  if (nrow(potential_controls) == 0) {
    cat("No potential controls found for case", case_row$case_code, "\n")
    next  # Skip to the next iteration (case)
  }
  
  # Add the case_code to the potential controls
  potential_controls$case_code <- case_row$case_code
  
  # Add controls to the list
  controls_list[[i]] <- potential_controls
}

# Combine all controls into a single dataframe (only non-empty controls will be added)
controls <- do.call(rbind, controls_list)

# Combine cases and controls into a final dataframe
final_df <- rbind(
  cases[, c("strain", "case_code", "DATESAMPLING", "scorpio_call", "PERSON_ID")], 
  controls[, c("strain", "case_code", "DATESAMPLING", "scorpio_call", "PERSON_ID")]
)

# Add a column to distinguish between cases and controls
final_df$type <- ifelse(final_df$strain %in% case_strains, "case", "control")

# Debugging: Check final case count and ensure proper combination
cat("Final dataframe case count:", sum(final_df$type == "case"), "\n")

# View a summary of the final dataframe
summary(final_df)

saveRDS(final_df, file="case_control_df.RDS")

cases <- final_df[final_df$type == "case", ]
controls <- final_df[final_df$type == "control", ]

# Save case strains to a text file
writeLines(cases$strain, "case_strains.txt")

# Save control strains to a text file
writeLines(controls$strain, "control_strains.txt")

# Create an empty list to store the case-to-controls mapping
case_to_controls <- list()

# For each case, get its corresponding controls
for (i in seq_len(nrow(cases))) {
  case_row <- cases[i, ]
  
  # Get the control strains for this specific case by matching on case_code
  control_strains <- controls[controls$case_code == case_row$case_code, "strain"]
  
  # Store the case strain and its associated control strains
  case_to_controls[[case_row$strain]] <- control_strains
}

# Write the case-to-controls mapping to a text file
writeLines(
  sapply(names(case_to_controls), function(case) {
    paste(case, paste(case_to_controls[[case]], collapse = "\t"), sep = "\t")
  }),
  "case_to_controls.txt"
)

# Run find_identical_sequences.py

#Reading back output ------
zero_diff <- read.csv(file="zero_diff_matches.csv")
one_diff <- read.csv(file="one_diff_matches.csv")
two_diff <- read.csv(file="two_diff_matches.csv")

# Create a summary table for each dataframe
summaries <- data.frame(
  Dataset = c("Zero Differences", "One Difference", "Two Differences"),
  Total_Rows = c(nrow(zero_diff), nrow(one_diff), nrow(two_diff)),
  Unique_Controls = c(length(unique(zero_diff$case_strain)), 
                      length(unique(one_diff$case_strain)), 
                      length(unique(two_diff$case_strain)))
)

# Print the summaries
print(summaries)