
# Preparing data for case-control study ----------------------
# Full metadata contains sequence data for all individuals; all_without_ct_metadata is metadata for persistent infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(data.table)

full_metadata <- read.csv(file = "")
all_without_ct_metadata <- readRDS(file = "")

#Slicing data to get first positive row
all_without_ct_metadata_unique <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  slice_min(DATESAMPLING) %>%
  ungroup() %>% 
  select(PERSON_ID, infection_episode, DATESAMPLING, strain, CT, major_scorpio_call, major_variant)

full_metadata_clean <- full_metadata %>%
  select(PERSON_ID, DATESAMPLING, strain, CT, scorpio_call, major_variant) %>%
  mutate(DATESAMPLING = as.Date(DATESAMPLING))

variant_mapping <- c(
  "Alpha (B.1.1.7-like)" = "Alpha",
  "Delta (AY.4-like)" = "Delta",
  "Delta (AY.4.2-like)" = "Delta",
  "Delta (B.1.617.2-like)" = "Delta",
  "Omicron (BA.1-like)" = "Omicron BA.1",
  "Omicron (BA.2-like)" = "Omicron BA.2",
  "Omicron (BA.5-like)" = "Omicron BA.5"
)

full_metadata_clean <- full_metadata_clean %>%
  mutate(major_scorpio_call = recode(scorpio_call, !!!variant_mapping))

# Creating dataframe to find potential matches for controls -------

# Match cases with potential controls
matched <- all_without_ct_metadata_unique %>%
  mutate(event_id = row_number()) %>%  # Create unique event_id for each case
  cross_join(full_metadata_clean, suffix = c("_case", "_control")) %>%  # Perform cross join
  filter(abs(difftime(DATESAMPLING_case, DATESAMPLING_control, units = "days")) <= 1,
         PERSON_ID_case != PERSON_ID_control,
         major_scorpio_call_case == major_scorpio_call_control) %>%  # Ensure same major_scorpio_call
  mutate(type = "control")  # Add column to indicate it's a control

# Prepare cases dataframe
cases <- all_without_ct_metadata_unique %>%
  mutate(event_id = row_number()) %>%  # Recreate event_id to align with matched data
  mutate(type = "case") %>%  # Add column to indicate it's a case
  select(event_id, PERSON_ID, infection_episode, DATESAMPLING, strain, CT, major_scorpio_call, major_variant, type)

# Prepare potential controls dataframe
potential_controls <- matched %>%
  select(event_id, 
         PERSON_ID = PERSON_ID_control, 
         DATESAMPLING = DATESAMPLING_control, 
         strain = strain_control, 
         CT = CT_control, 
         major_scorpio_call = major_scorpio_call_control, 
         major_variant = major_variant_control, 
         type)

# Combine Cases and Controls
potential_cases_and_controls_all <- bind_rows(cases, potential_controls) %>%
  arrange(event_id)

#Loading extra test metadata in
drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej()
COVID_TEST <- dbReadTable(conn = ,
                          name = "",
                          schema = '')
unique_person_ids <- potential_cases_and_controls_all %>%
  pull(PERSON_ID) %>%
  unique()
COVID_TEST <- COVID_TEST %>%
  filter(PERSON_ID %in% unique_person_ids) %>%
  filter(CASEDEF == "SARS2")

potential_cases_and_controls_all$DATESAMPLING <- as.Date(potential_cases_and_controls_all$DATESAMPLING)
COVID_test_subset_potential_cases_and_controls_all$PRDATE_ADJUSTED <- as.Date(COVID_test_subset_potential_cases_and_controls_all$PRDATE_ADJUSTED)


# Confirming that the cases are in fact not persistent infections ------------
# Ensuring that the following conditions are met:
# No other positive tests AFTER within >= 26 days AND <120 days 
# AND 
# no other positive tests within >= 26 days AND <120 days BEFORE the test
# AND
# negative test within >= 26 days AND <120 days AFTER


# Load your data into data.tables
potential_cases_and_controls_all <- as.data.table(potential_cases_and_controls_all)
controls <- potential_cases_and_controls_all[type == "control"]
cases <- potential_cases_and_controls_all[type != "control"]

COVID_test_subset_potential_cases_and_controls_all <- as.data.table(COVID_test_subset_potential_cases_and_controls_all)
positive_tests <- COVID_test_subset_potential_cases_and_controls_all[SVARRESULTAT == 1]
negative_tests <- COVID_test_subset_potential_cases_and_controls_all[SVARRESULTAT == 0]

# Check if there is a valid negative test
has_valid_negative_test <- function(person_id, sampling_date) {
  nrow(negative_tests[
    PERSON_ID == person_id &
      PRDATE_ADJUSTED > (sampling_date + 26) & 
      PRDATE_ADJUSTED < (sampling_date + 120)
  ]) > 0
}

# Define function to check for invalid positive tests
has_invalid_positive_tests <- function(person_id, sampling_date) {
  pos_before <- nrow(positive_tests[
    PERSON_ID == person_id &
      PRDATE_ADJUSTED >= (sampling_date - 120) & 
      PRDATE_ADJUSTED < (sampling_date - 26)
  ])
  
  pos_after <- nrow(positive_tests[
    PERSON_ID == person_id &
      PRDATE_ADJUSTED >= (sampling_date + 26) & 
      PRDATE_ADJUSTED < (sampling_date + 120)
  ])
  
  # Return TRUE if there are any positive tests in either range (invalid), FALSE otherwise
  pos_before > 0 | pos_after > 0
}

# Define function to validate controls
is_valid_control_with_neg <- function(control_row) {
  person_id <- control_row$PERSON_ID
  sampling_date <- control_row$DATESAMPLING
  
  # Check for valid negative test and invalid positive tests
  valid_neg_test <- has_valid_negative_test(person_id, sampling_date)
  invalid_pos_tests <- has_invalid_positive_tests(person_id, sampling_date)
  
  # Control is valid if there is a valid negative test and no invalid positive tests
  valid_neg_test & !invalid_pos_tests
}

# Prepare for parallel processing
plan(multisession, workers = 32)
options(future.globals.maxSize = 1000 * 1024^2)  # Set to 1 GB, or adjust as needed

# Split the controls data into chunks
n_chunks <- 200
chunk_size <- ceiling(nrow(controls) / n_chunks)
control_chunks <- split(controls, ceiling(seq_len(nrow(controls)) / chunk_size))

# Use future_lapply to process each chunk in parallel
valid_controls_list <- future_lapply(seq_along(control_chunks), function(i) {
  chunk <- control_chunks[[i]]
  
  # Apply validation function to each row of the chunk
  valid_indices <- chunk[, sapply(1:.N, function(idx) is_valid_control_with_neg(.SD[idx, ])), .SDcols = names(chunk)]
  
  # Print progress
  cat(sprintf("Processed chunk %d of %d\n", i, length(control_chunks)))
  
  # Return the valid controls in the current chunk
  chunk[valid_indices]
}, future.seed = TRUE)

# Combine valid controls from all chunks
valid_controls <- rbindlist(valid_controls_list)

# Combine valid controls with cases
subset_cases_and_controls_all <- rbindlist(list(cases, valid_controls))

# Arrange by event_id and save
setorder(subset_cases_and_controls_all, event_id)
saveRDS(subset_cases_and_controls_all, file="subset_cases_and_controls_all.RDS")


# Loading extra metadata ---------

drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej('STATPROD', dbuser = 'XP2')

COVID_VACC <- dbReadTable(conn = conn,
                          name = "COVID_VAC202301",
                          schema = 'D222301')


lifelines <- dbReadTable(conn = conn,
                         name = "LIFELINES",
                         schema = 'D221911')


lifelines_koen <- dbReadTable(conn = conn,
                              name = "LIFELINES_KOEN",
                              schema = 'D221911')


relevant_person_ids <- unique(subset_cases_and_controls_all$PERSON_ID)

COVID_VACC_filtered <- COVID_VACC %>%
  filter(PERSON_ID %in% relevant_person_ids) %>%     # Filter for relevant PERSON_IDs
  group_by(PERSON_ID) %>%
  slice(1) %>%                                       # Keep the first row for each PERSON_ID
  ungroup()

lifelines_filtered <- lifelines %>%
  filter(PERSON_ID %in% relevant_person_ids) %>%     # Filter for relevant PERSON_IDs
  group_by(PERSON_ID) %>%
  slice(1) %>%                                       # Keep the first row for each PERSON_ID
  ungroup()

lifelines_koen_filtered <- lifelines_koen %>%
  filter(PERSON_ID %in% relevant_person_ids) %>%     # Filter for relevant PERSON_IDs
  group_by(PERSON_ID) %>%
  slice(1) %>%                                       # Keep the first row for each PERSON_ID
  ungroup()

subset_cases_and_controls_all_with_metadata <- subset_cases_and_controls_all %>%
  left_join(COVID_VACC_filtered, by = "PERSON_ID") %>%    
  left_join(lifelines_koen_filtered, by = "PERSON_ID") %>%    
  left_join(lifelines_filtered, by = "PERSON_ID")     

saveRDS(subset_cases_and_controls_all_with_metadata, file="subset_cases_and_controls_all_with_metadata.RDS")



# Geography data ---------------------------
query <- ""
geo_address <- dbGetQuery(conn = con2, statement = query)

filter_geo_address_row_by_row_dt <- function(original_df, geo_address_df) {
  # Convert data frames to data.table
  original_dt <- as.data.table(original_df)
  geo_dt <- as.data.table(geo_address_df)
  
  # Keep only the relevant columns in geo_dt
  geo_dt <- geo_dt[, .(PERSON_ID, KOM, ID, BOP_VFRA, BOP_VTIL)]
  
  # Ensure column classes for geo_address_df match
  geo_dt[, PERSON_ID := as.numeric(PERSON_ID)]
  geo_dt[, KOM := as.character(KOM)]
  geo_dt[, ID := as.character(ID)]
  geo_dt[, BOP_VFRA := as.POSIXct(BOP_VFRA)]
  geo_dt[, BOP_VTIL := as.POSIXct(BOP_VTIL)]
  
  # Create an empty list to store results
  results_list <- vector("list", nrow(original_dt))
  
  # Loop through each row in the original data.table
  for (i in seq_len(nrow(original_dt))) {
    # Extract the current row
    row <- original_dt[i, ]
    
    # Print progress for every 10,000 rows
    if (i %% 10000 == 0) {
      print(i)
    }
    
    # Get the current PERSON_ID and DATESAMPLING
    person_id <- row$PERSON_ID
    datesampling <- row$DATESAMPLING
    
    # Find matches in geo_address for this PERSON_ID
    geo_match <- geo_dt[PERSON_ID == person_id]
    
    # Further filter by DATESAMPLING falling between BOP_VFRA and BOP_VTIL
    geo_match <- geo_match[datesampling >= BOP_VFRA & datesampling <= BOP_VTIL]
    
    if (nrow(geo_match) > 0) {
      # Use the first matching row and exclude duplicates
      result <- cbind(row, geo_match[1, ])
    } else {
      # Create an NA row with the correct column names and types
      geo_na <- data.table(
        PERSON_ID = NA_real_,  # Numeric
        KOM = NA_character_,   # Character
        ID = NA_character_,    # Character
        BOP_VFRA = as.POSIXct(NA),  # Date
        BOP_VTIL = as.POSIXct(NA)   # Date
      )
      
      # Combine the current row with the NA row
      result <- cbind(row, geo_na)
    }
    
    # Store the result
    results_list[[i]] <- result
    
  }
  
  # Combine the results list into a final data.table
  final_result <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  
  return(final_result)
}

# Apply the function to your data
subset_cases_and_controls_all_with_metadata_geo <- filter_geo_address_row_by_row_dt(subset_cases_and_controls_all_with_metadata, geo_address)

# Adding number of vaccinations as a categorical variable ---
# Define a function to process each dataframe
process_vaccination_data <- function(df) {
  # Ensure the date columns are in Date format
  df$DATESAMPLING <- as.Date(df$DATESAMPLING)
  df$FIRST_VACCINEDATE <- as.Date(df$FIRST_VACCINEDATE)
  df$SECOND_VACCINEDATE <- as.Date(df$SECOND_VACCINEDATE)
  df$THIRD_VACCINEDATE <- as.Date(df$THIRD_VACCINEDATE)
  
  # Define a function to count vaccinations with a 7-day lag
  vaccination_count <- function(sample_date, vacc_date) {
    ifelse(!is.na(vacc_date) & (sample_date - vacc_date) >= 7, 1, 0)
  }
  
  # Calculate the number of vaccinations at the time of infection
  df$num_vacc_at_infection <- 
    vaccination_count(df$DATESAMPLING, df$FIRST_VACCINEDATE) +
    vaccination_count(df$DATESAMPLING, df$SECOND_VACCINEDATE) +
    vaccination_count(df$DATESAMPLING, df$THIRD_VACCINEDATE)
  
  # Convert the number of vaccinations into a categorical variable
  df$num_vacc_at_infection_cat <- cut(
    df$num_vacc_at_infection, 
    breaks = c(-Inf, 0, 1, 2, 3), # breaks define intervals for 0, 1, 2, and 3+ vaccinations
    labels = c("0", "1", "2", "3"), # assign category labels
    right = TRUE # include the upper limit (i.e., the right boundary) in the interval
  )
  
  return(df)
}

# Apply the function to both dataframes
subset_cases_and_controls_all_with_metadata_geo <- process_vaccination_data(subset_cases_and_controls_all_with_metadata_geo)

process_age_data <- function(df) {
  # Ensure BIRTHDAY is in Date format
  df$BIRTHDAY <- as.Date(df$BIRTHDAY, format = "%Y-%m-%d")
  df$DATESAMPLING <- as.Date(df$DATESAMPLING)
  
  # Calculate age at infection
  df$age_at_infection <- as.numeric(difftime(df$DATESAMPLING, df$BIRTHDAY, units = "weeks")) / 52.25
  
  # Create age group as a categorical variable
  df$age_group <- cut(
    df$age_at_infection, 
    breaks = c(-Inf, 15, 30, 45, 60, 75, Inf), # breaks for the specified age groups
    labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"), # labels for the groups
    right = FALSE # age intervals are left-closed, so [0-15) means age >= 0 and < 15
  )
  
  return(df)
}

# Apply the function to both dataframes
subset_cases_and_controls_all_with_metadata_geo <- process_age_data(subset_cases_and_controls_all_with_metadata_geo)
colnames(subset_cases_and_controls_all_with_metadata_geo)[42] <- "PERSON_ID_2"



# Loading diagnosis database to add diagnosis metadata to all rows -------------------------

# IMID --------------------
IMID_codes <- data.frame(icd_10 = c("K50", "K50.0", "K50.1", "K50.8", "K50.9",
                                    "K51", "K51.0", "K51.2", "K51.3", "K51.4", "K51.5", "K51.8", "K51.9",
                                    "M45", "M46", "M46.0", "M46.2", "M46.3", "M46.4", "M46.5", "M46.8", "M46.9",
                                    "M05", "M05.0", "M05.1", "M05.2", "M05.3", "M05.8", "M05.9",
                                    "M06", "M06.0", "M06.1", "M06.2", "M06.3", "M06.4", "M06.8", "M06.9",
                                    "M07", "M07.0", "M07.1", "M07.2", "M07.3", "M07.4", "M07.5", "M07.6",
                                    "L40", "L40.0", "L40.1", "L40.2", "L40.3", "L40.4", "L40.5", "L40.8", "L40.9")
)
IMID_codes$diag_code <- gsub("\\.", "", paste0("D", IMID_codes$icd_10))

IMID_diag_codes_list <- paste0("'", IMID_codes$diag_code, "'", collapse = ", ")
IMID_query_LPR3 <- paste0("",
                     "", IMID_diag_codes_list, ") ")
IMID_diagnoses_LPR3 <- dbGetQuery(conn, IMID_query_LPR3)

IMID_query_LPR2 <- paste0("",
                     "")
IMID_diagnoses_LPR2 <- dbGetQuery(conn, IMID_query_LPR2)
IMID_diagnoses_LPR2 <- IMID_diagnoses_LPR2 %>%
  rename(AKTIONSDIAGNOSE = C_DIAG)

# Immunosuppression ---------------------

# Codes taken from here: https://datacompass.lshtm.ac.uk/id/eprint/1284/

immunosuppression_codes <- data.frame(
  icd_10 = c("B20", "B20.0", "B20.1", "B20.2", "B20.3", "B20.4", "B20.5", "B20.6", "B20.7", "B20.8", "B20.9",
             "B21", "B21.0", "B21.1", "B21.2", "B21.3", "B21.7", "B21.8", "B21.9", "B22", "B22.0", "B22.1",
             "B22.2", "B22.7", "B23", "B23.0", "B23.1", "B23.2", "B23.8", "B24", "C81", "C81.0", "C81.1",
             "C81.2", "C81.3", "C81.7", "C81.9", "C82", "C82.0", "C82.1", "C82.2", "C82.7", "C82.9", "C83",
             "C83.0", "C83.1", "C83.2", "C83.3", "C83.4", "C83.5", "C83.6", "C83.7", "C83.8", "C83.9", "C84",
             "C84.0", "C84.1", "C84.2", "C84.3", "C84.4", "C84.5", "C85", "C85.0", "C85.1", "C85.7", "C85.9",
             "C88", "C88.0", "C88.1", "C88.2", "C88.3", "C88.7", "C88.9", "C90", "C90.0", "C90.1", "C90.2",
             "C91", "C91.0", "C91.1", "C91.2", "C91.3", "C91.4", "C91.5", "C91.7", "C91.9", "C92", "C92.0",
             "C92.1", "C92.2", "C92.3", "C92.4", "C92.5", "C92.7", "C92.9", "C93", "C93.0", "C93.1", "C93.2",
             "C93.7", "C93.9", "C94", "C94.0", "C94.1", "C94.2", "C94.3", "C94.4", "C94.5", "C94.7", "C95",
             "C95.0", "C95.1", "C95.2", "C95.7", "C95.9", "C96", "C96.0", "C96.1", "C96.2", "C96.3", "C96.7",
             "C96.9", "D61.1", "D61.2", "D61.3", "D61.8", "D61.9", "D81", "D81.0", "D81.1", "D81.2", "D81.3",
             "D81.4", "D81.5", "D81.6", "D81.7", "D81.8", "D81.9", "D82.0", "D82.1", "D82.2", "D83", "D83.0",
             "D83.1", "D83.2", "D83.8", "D83.9", "D89.3", "F02.4", "T86.0", "T86.2", "T86.3", "T86.4", "Y83.0",
             "Z21", "Z94.0", "Z94.1", "Z94.2", "Z94.3", "Z94.4"))
immunosuppression_codes$diag_code <- gsub("\\.", "", paste0("D", immunosuppression_codes$icd_10))

immunosuppression_diag_codes_list <- paste0("'", immunosuppression_codes$diag_code, "'", collapse = ", ")
immunosuppression_query_LPR3 <- paste0("",
                     "")
immunosuppression_diagnoses_LPR3 <- dbGetQuery(conn, immunosuppression_query_LPR3)

immunosuppression_query_LPR2 <- paste0("",
                                       "")
immunosuppression_diagnoses_LPR2 <- dbGetQuery(conn, immunosuppression_query_LPR2)
immunosuppression_diagnoses_LPR2 <- immunosuppression_diagnoses_LPR2 %>%
  rename(AKTIONSDIAGNOSE = C_DIAG)

# Diabetes --------------

# Codes taken from here: https://datacompass.lshtm.ac.uk/id/eprint/1273/

diabetes_codes <- data.frame(
  icd_10 = c("E10", "E10.0", "E10.1", "E10.2", "E10.3", "E10.4", "E10.5", "E10.6", "E10.7", "E10.8", "E10.9",
          "E11", "E11.0", "E11.1", "E11.2", "E11.3", "E11.4", "E11.5", "E11.6", "E11.7", "E11.8", "E11.9",
          "E13", "E13.0", "E13.1", "E13.2", "E13.3", "E13.4", "E13.5", "E13.6", "E13.7", "E13.8", "E13.9",
          "E14", "E14.0", "E14.1", "E14.2", "E14.3", "E14.4", "E14.5", "E14.6", "E14.7", "E14.8", "E14.9",
          "G59.0", "G63.2", "H28.0", "H36.0", "M14.2", "N08.3", "O24.0", "O24.1", "O24.3", "Y42.3")
)
diabetes_codes$diag_code <- gsub("\\.", "", paste0("D", diabetes_codes$icd_10))

diabetes_diag_codes_list <- paste0("'", diabetes_codes$diag_code, "'", collapse = ", ")
diabetes_query_LPR3 <- paste0("")
diabetes_diagnoses_LPR3 <- dbGetQuery(conn, diabetes_query_LPR3)

diabetes_query_LPR2 <- paste0("")
diabetes_diagnoses_LPR2 <- dbGetQuery(conn, diabetes_query_LPR2)
diabetes_diagnoses_LPR2 <- diabetes_diagnoses_LPR2 %>%
  rename(AKTIONSDIAGNOSE = C_DIAG)

# Autoimmune diseases ---------------------

# Codes taken from here: https://datacompass.lshtm.ac.uk/id/eprint/2853/

autoimmune_codes <- data.frame(
  icd_10 = c("M08.1", "M08.20", "M08.211", "M08.212", "M08.219", "M08.221", "M08.222", "M08.229", "M08.231",
             "M08.232", "M08.239", "M08.241", "M08.242", "M08.249", "M08.251", "M08.252", "M08.259",
             "M08.261", "M08.262", "M08.269", "M08.271", "M08.272", "M08.279", "M08.28", "M08.29",
             "M08.2A", "M08.3", "M08.40", "M08.411", "M08.412", "M08.419", "M08.421", "M08.422",
             "M08.429", "M08.431", "M08.432", "M08.439", "M08.441", "M08.442", "M08.449", "M08.451",
             "M08.452", "M08.459", "M08.461", "M08.462", "M08.469", "M08.471", "M08.472", "M08.479",
             "M08.48", "M08.4A", "M08.80", "M08.811", "M08.812", "M08.819", "M08.821", "M08.822",
             "M08.829", "M08.831", "M08.832", "M08.839", "M08.841", "M08.842", "M08.849", "M08.851",
             "M08.852", "M08.859", "M08.861", "M08.862", "M08.869", "M08.871", "M08.872", "M08.879",
             "M08.88", "M08.89", "M08.90", "M08.911", "M08.912", "M08.919", "M08.921", "M08.922",
             "M08.929", "M08.931", "M08.932", "M08.939", "M08.941", "M08.942", "M08.949", "M08.951",
             "M08.952", "M08.959", "M08.961", "M08.962", "M08.969", "M08.971", "M08.972", "M08.979",
             "M08.98", "M08.99", "M08.9A", "M09.8", "M31.30", "M31.31", "M33.00", "M33.01", "M33.02",
             "M33.03", "M33.09", "M33.10", "M33.11", "M33.12", "M33.13", "M33.19", "M33.20", "M33.21",
             "M33.22", "M33.29", "M33.90", "M33.91", "M33.92", "M33.93", "M33.99", "M31.5", "M31.6",
             "M35.3", "G70.00", "G70.01", "M34.0", "M34.1", "M34.81", "M34.82", "M34.83", "M34.89",
             "M34.9", "M32.10", "M32.11", "M32.12", "M32.13", "M32.14", "M32.15", "M32.19", "M32.8",
             "M32.9", "M35.00", "M35.01", "M35.02", "M35.03", "M35.04", "M35.05", "M35.06", "M35.07",
             "M35.08", "M35.09", "M35.0A", "M35.0B", "M35.0C", "M45.0", "M45.1", "M45.2", "M45.3",
             "M45.4", "M45.5", "M45.6", "M45.7", "M45.8", "M45.9", "M45.A0", "M45.A1", "M45.A2",
             "M45.A3", "M45.A4", "M45.A5", "M45.A6", "M45.A7", "M45.A8", "M45.AB", "M35.2", "D86.0",
             "D86.1", "D86.2", "D86.3", "D86.81", "D86.82", "D86.83", "D86.84", "D86.85", "D86.86",
             "D86.87", "D86.89", "D86.9", "M63.3", "G53.2"))
autoimmune_codes$diag_code <- gsub("\\.", "", paste0("D", autoimmune_codes$icd_10))

autoimmune_diag_codes_list <- paste0("'", autoimmune_codes$diag_code, "'", collapse = ", ")
autoimmune_query_LPR3 <- paste0("")
autoimmune_diagnoses_LPR3 <- dbGetQuery(conn, autoimmune_query_LPR3)

autoimmune_query_LPR2 <- paste0("")
autoimmune_diagnoses_LPR2 <- dbGetQuery(conn, autoimmune_query_LPR2)
autoimmune_diagnoses_LPR2 <- autoimmune_diagnoses_LPR2 %>%
  rename(AKTIONSDIAGNOSE = C_DIAG)

# Subsetting diagnosis dataframes to only include relevant people
ids_with_metadata_geo <- subset_cases_and_controls_all_with_metadata_geo$PERSON_ID
common_ids <- ids_with_metadata_geo

IMID_diagnoses_LPR3_subset <- IMID_diagnoses_LPR3[IMID_diagnoses_LPR3$PERSON_ID %in% common_ids, ]
immunosuppression_diagnoses_LPR3_subset <- immunosuppression_diagnoses_LPR3[immunosuppression_diagnoses_LPR3$PERSON_ID %in% common_ids, ]
diabetes_diagnoses_LPR3_subset <- diabetes_diagnoses_LPR3[diabetes_diagnoses_LPR3$PERSON_ID %in% common_ids, ]
autoimmune_diagnoses_LPR3_subset <- autoimmune_diagnoses_LPR3[autoimmune_diagnoses_LPR3$PERSON_ID %in% common_ids, ]

IMID_diagnoses_LPR2_subset <- IMID_diagnoses_LPR2[IMID_diagnoses_LPR2$PERSON_ID %in% common_ids, ]
immunosuppression_diagnoses_LPR2_subset <- immunosuppression_diagnoses_LPR2[immunosuppression_diagnoses_LPR2$PERSON_ID %in% common_ids, ]
diabetes_diagnoses_LPR2_subset <- diabetes_diagnoses_LPR2[diabetes_diagnoses_LPR2$PERSON_ID %in% common_ids, ]
autoimmune_diagnoses_LPR2_subset <- autoimmune_diagnoses_LPR2[autoimmune_diagnoses_LPR2$PERSON_ID %in% common_ids, ]


# Function to add diagnosis columns with unique values -------------------
add_diagnosis_columns_combined <- function(df, LPR2_df, LPR3_df, diagnosis_type) {
  # Convert input dataframes to data.table for faster processing
  df <- as.data.table(df)
  if (!is.null(LPR2_df)) LPR2_df <- as.data.table(LPR2_df)
  if (!is.null(LPR3_df)) LPR3_df <- as.data.table(LPR3_df)
  
  # Initialize the diagnosis columns if they do not exist yet
  df[, (paste0(diagnosis_type, "_diagnosis")) := 0]
  df[, (paste0(diagnosis_type, "_codes")) := ""]
  df[, (paste0(diagnosis_type, "_dates")) := ""]
  df[, (paste0(diagnosis_type, "_source")) := ""]
  
  # Initialize the prepandemic diagnosis column
  df[, (paste0(diagnosis_type, "_diagnosis_prepandemic")) := 0]
  
  # ---- Handle LPR3 (with DATO_START) ----
  if (!is.null(LPR3_df)) {
    # Ensure dates are correctly formatted
    LPR3_df[, DATO_START := as.Date(DATO_START)]
    df[, DATESAMPLING := as.Date(DATESAMPLING)]
    
    # Merge on PERSON_ID and filter by DATO_START < DATESAMPLING
    merged_LPR3 <- merge(df, LPR3_df, by = "PERSON_ID", allow.cartesian = TRUE)
    merged_LPR3 <- merged_LPR3[DATO_START < DATESAMPLING]
    
    # Collect unique AKTIONSDIAGNOSE and update diagnosis, codes, dates, and source columns
    if (nrow(merged_LPR3) > 0) {
      result_LPR3 <- merged_LPR3[, .(
        unique_codes = paste(unique(AKTIONSDIAGNOSE), collapse = ", "),
        unique_dates = paste(unique(DATO_START), collapse = ", "),
        source = "LPR3"
      ), by = .(PERSON_ID)]
      
      df[result_LPR3, on = "PERSON_ID", (paste0(diagnosis_type, "_diagnosis")) := 1]
      
      # Combine values from LPR3 with existing values
      df[result_LPR3, on = "PERSON_ID", (paste0(diagnosis_type, "_codes")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_codes")]] == "", 
                  unique_codes, 
                  paste(.SD[[paste0(diagnosis_type, "_codes")]], unique_codes, sep = ", "))]
      
      df[result_LPR3, on = "PERSON_ID", (paste0(diagnosis_type, "_dates")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_dates")]] == "", 
                  unique_dates, 
                  paste(.SD[[paste0(diagnosis_type, "_dates")]], unique_dates, sep = ", "))]
      
      df[result_LPR3, on = "PERSON_ID", (paste0(diagnosis_type, "_source")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_source")]] == "", 
                  source, 
                  paste(.SD[[paste0(diagnosis_type, "_source")]], source, sep = ", "))]
      
      # Handle prepandemic: Filter for dates before 2020-01-01
      result_LPR3_prepandemic <- merged_LPR3[DATO_START < "2020-01-01", .(
        unique_codes = paste(unique(AKTIONSDIAGNOSE), collapse = ", "),
        unique_dates = paste(unique(DATO_START), collapse = ", "),
        source = "LPR3"
      ), by = .(PERSON_ID)]
      
      df[result_LPR3_prepandemic, on = "PERSON_ID", (paste0(diagnosis_type, "_diagnosis_prepandemic")) := 1]
    }
  }
  
  # ---- Handle LPR2 (with AAR for year only) ----
  if (!is.null(LPR2_df)) {
    # Merge on PERSON_ID for LPR2
    merged_LPR2 <- merge(df, LPR2_df, by = "PERSON_ID", allow.cartesian = TRUE)
    
    # Collect unique AKTIONSDIAGNOSE, AAR (year), and update columns
    if (nrow(merged_LPR2) > 0) {
      result_LPR2 <- merged_LPR2[, .(
        unique_codes = paste(unique(AKTIONSDIAGNOSE), collapse = ", "),
        unique_dates = paste(unique(AAR), collapse = ", "),  # Using AAR (year) here
        source = "LPR2"
      ), by = .(PERSON_ID)]
      
      df[result_LPR2, on = "PERSON_ID", (paste0(diagnosis_type, "_diagnosis")) := 1]
      
      # Combine values from LPR2 with existing values
      df[result_LPR2, on = "PERSON_ID", (paste0(diagnosis_type, "_codes")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_codes")]] == "", 
                  unique_codes, 
                  paste(.SD[[paste0(diagnosis_type, "_codes")]], unique_codes, sep = ", "))]
      
      df[result_LPR2, on = "PERSON_ID", (paste0(diagnosis_type, "_dates")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_dates")]] == "", 
                  unique_dates, 
                  paste(.SD[[paste0(diagnosis_type, "_dates")]], unique_dates, sep = ", "))]
      
      df[result_LPR2, on = "PERSON_ID", (paste0(diagnosis_type, "_source")) := 
           ifelse(.SD[[paste0(diagnosis_type, "_source")]] == "", 
                  source, 
                  paste(.SD[[paste0(diagnosis_type, "_source")]], source, sep = ", "))]
      
      # Handle prepandemic: Filter for records before 2020
      result_LPR2_prepandemic <- merged_LPR2[AAR < 2020, .(
        unique_codes = paste(unique(AKTIONSDIAGNOSE), collapse = ", "),
        unique_dates = paste(unique(AAR), collapse = ", "),  # Using AAR (year) here
        source = "LPR2"
      ), by = .(PERSON_ID)]
      
      df[result_LPR2_prepandemic, on = "PERSON_ID", (paste0(diagnosis_type, "_diagnosis_prepandemic")) := 1]
    }
  }
  
  # Convert back to data.frame for compatibility with other parts of the code
  return(as.data.frame(df))
}

#Adding data
subset_cases_and_controls_all_with_metadata_geo_diag <- add_diagnosis_columns_combined(
  subset_cases_and_controls_all_with_metadata_geo, 
  LPR2_df = IMID_diagnoses_LPR2_subset, 
  LPR3_df = IMID_diagnoses_LPR3_subset, 
  diagnosis_type = "IMID"
)

subset_cases_and_controls_all_with_metadata_geo_diag <- add_diagnosis_columns_combined(
  subset_cases_and_controls_all_with_metadata_geo_diag, 
  LPR2_df = diabetes_diagnoses_LPR2_subset, 
  LPR3_df = diabetes_diagnoses_LPR3_subset, 
  diagnosis_type = "diabetes"
)

subset_cases_and_controls_all_with_metadata_geo_diag <- add_diagnosis_columns_combined(
  subset_cases_and_controls_all_with_metadata_geo_diag, 
  LPR2_df = immunosuppression_diagnoses_LPR2_subset, 
  LPR3_df = immunosuppression_diagnoses_LPR3_subset, 
  diagnosis_type = "immunosuppression"
)

subset_cases_and_controls_all_with_metadata_geo_diag <- add_diagnosis_columns_combined(
  subset_cases_and_controls_all_with_metadata_geo_diag, 
  LPR2_df = autoimmune_diagnoses_LPR2_subset, 
  LPR3_df = autoimmune_diagnoses_LPR3_subset, 
  diagnosis_type = "autoimmune"
)

# Convert relevant variables to factors
subset_cases_and_controls_all_with_metadata_geo_diag$age_group <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$age_group, levels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"))
subset_cases_and_controls_all_with_metadata_geo_diag$KOEN <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$KOEN, levels = c(1, 2), labels = c("Male", "Female"))
subset_cases_and_controls_all_with_metadata_geo_diag$major_variant <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$major_variant)
subset_cases_and_controls_all_with_metadata_geo_diag$major_scorpio_call <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$major_scorpio_call)
subset_cases_and_controls_all_with_metadata_geo_diag$IMID_diagnosis <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$IMID_diagnosis, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$immunosuppression_diagnosis <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$immunosuppression_diagnosis, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$diabetes_diagnosis <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$diabetes_diagnosis, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$autoimmune_diagnosis <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$autoimmune_diagnosis, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$IMID_diagnosis_prepandemic <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$IMID_diagnosis_prepandemic, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$immunosuppression_diagnosis_prepandemic <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$immunosuppression_diagnosis_prepandemic, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$diabetes_diagnosis_prepandemic <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$diabetes_diagnosis_prepandemic, ordered = TRUE)
subset_cases_and_controls_all_with_metadata_geo_diag$autoimmune_diagnosis_prepandemic <- factor(subset_cases_and_controls_all_with_metadata_geo_diag$autoimmune_diagnosis_prepandemic, ordered = TRUE)








