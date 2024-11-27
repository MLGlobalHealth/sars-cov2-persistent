# Pre-processing to find potential persistent infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach); 
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb); 
library(DBI); library(parallel); library(ggsignif); library(Rcpp)

#Establish connection
drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej('', dbuser = '')
con2 <- DSTora::OraGenvej('', dbuser = '')

# To access data tables -------------
COVID_TEST <- dbReadTable(conn = conn,
                          name = "",
                          schema = '')

lifelines <- dbReadTable(conn = conn,
                         name = "",
                         schema = '')

lifelines_koen <- dbReadTable(conn = conn,
                              name = "",
                              schema = '')

sequence_metadata <- dbReadTable(conn = conn,
                                 name = "",
                                 schema = '')

# Loading metadata -------------
metadata_processed_multiple_rows_missing_CT <- read.csv(file="")

# Finding individuals where they have a potential prolonged infection, where at least ONE of their pairs of 
# sequences are >= 26 days and < 120 days apart. -------
metadata_processed_multiple_rows_missing_CT <- metadata_processed_multiple_rows_missing_CT %>%
  mutate(DATESAMPLING = ymd(DATESAMPLING))

# Create an empty data frame to store the result
valid_possible_persistent_rows <- data.frame()

# Loop through unique PERSON_ID values
for (pid in unique(metadata_processed_multiple_rows_missing_CT$PERSON_ID)) {
  # Subset the data for the current PERSON_ID
  subset_data <- metadata_processed_multiple_rows_missing_CT[metadata_processed_multiple_rows_missing_CT$PERSON_ID == pid, ]
  
  # Calculate pairwise time differences
  pairwise_diff <- outer(subset_data$DATESAMPLING, subset_data$DATESAMPLING, `-`)
  
  # Check if there are any missing values in the pairwise differences
  if (!any(is.na(pairwise_diff))) {
    # Check if any pair satisfies the condition
    if (any(pairwise_diff >= 26 & pairwise_diff < 120)) {
      # If satisfied, add all rows for this PERSON_ID to the result
      valid_possible_persistent_rows <- rbind(valid_possible_persistent_rows, subset_data)
    }
  }
}

# Reset row names
rownames(valid_possible_persistent_rows) <- NULL

# View the result
print(valid_possible_persistent_rows) # If you set threshold to 120 days, 10000ish rows


# Connecting to sequence metadata with sequences ---------
full_sequence_metadata <- read.csv(file="")
subset_persistent_sequence_metadata <- full_sequence_metadata %>%
  filter(File %in% valid_possible_persistent_rows$CONSENSUSFILENAME)

variant_mapping <- c(
  "Alpha (B.1.1.7-like)" = "Alpha",
  "Delta (AY.4-like)" = "Delta",
  "Delta (AY.4.2-like)" = "Delta",
  "Delta (B.1.617.2-like)" = "Delta",
  "Omicron (BA.1-like)" = "Omicron BA.1",
  "Omicron (BA.2-like)" = "Omicron BA.2",
  "Omicron (BA.5-like)" = "Omicron BA.5"
)

subset_persistent_sequence_metadata <- subset_persistent_sequence_metadata %>%
  mutate(major_scorpio_call = recode(scorpio_call, !!!variant_mapping))

subset_persistent_sequence_metadata <- subset_persistent_sequence_metadata %>%
  filter(major_scorpio_call %in% unique(variant_mapping))

subset_persistent_sequence_metadata <- subset_persistent_sequence_metadata %>%
  filter(!is.na(DATESAMPLING), !is.na(major_scorpio_call))

subset_persistent_sequence_metadata <- subset_persistent_sequence_metadata %>%
  filter(!is.na(DATESAMPLING))

subset_persistent_sequence_metadata$DATESAMPLING <- as.Date(subset_persistent_sequence_metadata$DATESAMPLING, format = "%Y-%m-%d")

valid_persons <- subset_persistent_sequence_metadata %>%
  group_by(PERSON_ID) %>%
  filter(n() > 1) %>%
  pull(PERSON_ID) %>%
  unique()

subset_persistent_sequence_metadata <- subset_persistent_sequence_metadata %>%
  filter(PERSON_ID %in% valid_persons)

# Creating infection episodes by date and if they have the same lineage -------
infection_episode_rows <- data.frame()

for (pid in unique(subset_persistent_sequence_metadata$PERSON_ID)) {
  # Subset the data for the current PERSON_ID
  subset_data <- subset_persistent_sequence_metadata %>% 
    filter(PERSON_ID == pid) %>% 
    arrange(DATESAMPLING)
  
  # Initialize infection episode counter
  episode_counter <- 1
  subset_data$infection_episode <- episode_counter
  
  # Iterate through the rows and assign infection episodes
  for (i in 2:nrow(subset_data)) {
    previous_date <- subset_data$DATESAMPLING[i - 1]
    current_date <- subset_data$DATESAMPLING[i]
    previous_lineage <- subset_data$major_scorpio_call[i - 1]
    current_lineage <- subset_data$major_scorpio_call[i]
    
    # Check the time difference and major_scorpio_call
    if (difftime(current_date, previous_date, units = "days") <= 120 && current_lineage == previous_lineage) {
      subset_data$infection_episode[i] <- episode_counter
    } else {
      episode_counter <- episode_counter + 1
      subset_data$infection_episode[i] <- episode_counter
    }
  }
  
  # Add the processed subset to the result dataframe
  infection_episode_rows <- rbind(infection_episode_rows, subset_data)
}

# Reset row names
rownames(infection_episode_rows) <- NULL
saveRDS(infection_episode_rows, file="")







