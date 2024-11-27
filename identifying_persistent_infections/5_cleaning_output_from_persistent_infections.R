# Collecting and cleaning data for persistent infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate);
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle);
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(doParallel);
library(foreach); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(furrr); library(future); library(future.apply); library(lubridate); library(seqinr);
library(adegenet); library(ggplot2); library(viridis)

#Establish connection
drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej('', dbuser = '')

# To access data tables 
lifelines <- dbReadTable(conn = conn,
                         name = "",
                         schema = '')

lifelines_koen <- dbReadTable(conn = conn,
                              name = "",
                              schema = '')

COVID_VACC <- dbReadTable(conn = conn,
                          name = "",
                          schema = '')

COVID_TEST <- dbReadTable(conn = conn,
                          name = "",
                          schema = '')




# Code to process input data ------------------------------------------------

# Code to process fasta files
process_sequences <- function(sequences, chunk_size = 1000) {
  
  # Convert to data.table
  df <- data.table(
    strain = names(sequences),
    sequence = sapply(sequences, function(x) toupper(paste(x, collapse = "")))
  )
  
  # Find the maximum sequence length
  max_len <- max(nchar(df$sequence))
  
  # Initialize data.table with position columns
  for (start_pos in seq(1, max_len, by = chunk_size)) {
    end_pos <- min(start_pos + chunk_size - 1, max_len)
    df <- add_position_columns(df, start_pos, end_pos)
    cat("Processed columns:", start_pos, "to", end_pos, "\n")
  }
  
  # remove the sequence column
  df[, sequence := NULL]
  
  # Ensure 'strain' is the first column
  setcolorder(df, c("strain", setdiff(names(df), "strain")))
  return(df)
}

# Function to add position columns
add_position_columns <- function(df, start_pos, end_pos) {
  for (i in start_pos:end_pos) {
    df[[paste0("Position_", i)]] <- substring(df$sequence, i, i)
  }
  return(df)
}



# Importing data from each folder ---------------------------------------------

# Without CT
results_all_possible_persistent_without_CT <- readRDS(file="results_all_possible_persistent_without_CT.RDS")
results_all_possible_persistent_without_CT <- results_all_possible_persistent_without_CT[results_all_possible_persistent_without_CT$rare_snp_count > 0, ]
all_possible_persistent_infections_without_CT <- readRDS(file="final_possible_persistent_infections_with_infection_episodes_without_CT.RDS")

# Getting subset metadata for each category and adding birthday, sex, location --------------------------------
all_without_ct_metadata <- all_possible_persistent_infections_without_CT %>%
  filter((PERSON_ID %in% results_all_possible_persistent_without_CT$PERSON_ID) & 
           (interaction(PERSON_ID, infection_episode) %in% interaction(results_all_possible_persistent_without_CT$PERSON_ID, 
                                                                       results_all_possible_persistent_without_CT$infection_episode)))


relevant_person_ids_all_without_ct <- all_without_ct_metadata %>%
  pull(PERSON_ID) %>%
  unique()

# Filter lifelines to include only relevant PERSON_IDs
lifelines_filtered <- lifelines %>%
  filter(PERSON_ID %in% c(
    relevant_person_ids_all_without_ct
  )) %>%
  group_by(PERSON_ID) %>%
  slice(1) %>%  # Keep the first row for each PERSON_ID
  select(PERSON_ID, BIRTHDAY)

# Filter lifelines_koen to include only relevant PERSON_IDs
lifelines_koen_filtered <- lifelines_koen %>%
  filter(PERSON_ID %in% c(
    relevant_person_ids_all_without_ct
  )) %>%
  group_by(PERSON_ID) %>%
  slice(1) %>%  # Keep the first row for each PERSON_ID
  select(PERSON_ID, KOEN)

# Filter COVID_VACC to include only relevant PERSON_IDs
COVID_VACC_filtered <- COVID_VACC %>%
  filter(PERSON_ID %in% c(
    relevant_person_ids_all_without_ct
  )) %>%
  group_by(PERSON_ID) %>%
  slice(1)

# Filter lifelines_koen to include only relevant PERSON_IDs
COVID_TEST_filtered <- COVID_TEST %>%
  filter(PERSON_ID %in% c(
    relevant_person_ids_all_without_ct
  ))

#Save for later use
saveRDS(COVID_TEST_filtered, file = "")
saveRDS(COVID_VACC_filtered, file = "")
saveRDS(lifelines_koen_filtered, file = "")
saveRDS(lifelines_filtered, file = "")

all_without_ct_metadata <- all_without_ct_metadata %>%
  left_join(lifelines_filtered, by = "PERSON_ID") %>%
  left_join(lifelines_koen_filtered, by = "PERSON_ID") %>%
  left_join(COVID_VACC_filtered, by = "PERSON_ID") %>%
  arrange(PERSON_ID, DATESAMPLING)

# Ensure unique strain values by keeping the first occurrence
all_without_ct_metadata <- all_without_ct_metadata %>%
  arrange(strain) %>%  # Optionally arrange by strain or another column to determine which row to keep
  distinct(strain, .keep_all = TRUE) %>%
  arrange(PERSON_ID, DATESAMPLING)

# Calculating at age infection
# Calculate age_at_infection and age_at_first_test for each dataframe
compute_ages <- function(df) {
  df %>%
    mutate(
      DATESAMPLING = ymd(DATESAMPLING),  # Convert DATESAMPLING to Date format
      BIRTHDAY = ymd(BIRTHDAY)  # Convert BIRTHDAY to Date format
    ) %>%
    mutate(
      age_at_infection = as.numeric(difftime(DATESAMPLING, BIRTHDAY, units = "days")) / 365.25  # Calculate age at infection in years
    ) %>%
    mutate(
      age_at_infection = ifelse(is.infinite(age_at_infection), NA_real_, age_at_infection)  # Replace Inf with NA
    ) %>%
    group_by(PERSON_ID, infection_episode) %>%
    mutate(
      age_at_first_test = ifelse(all(is.na(age_at_infection)), NA_real_, min(age_at_infection, na.rm = TRUE))  # Calculate minimum age_at_infection; handle all NA case
    ) %>%
    ungroup()
}

# Apply the function to each dataframe
all_without_ct_metadata <- compute_ages(all_without_ct_metadata)
saveRDS(all_without_ct_metadata, file = "all_without_ct_metadata.rds")
