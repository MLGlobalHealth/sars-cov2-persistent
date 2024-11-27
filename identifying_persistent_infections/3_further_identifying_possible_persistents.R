# Cleaning up potential prolonged infections; preprocesseing step for identifying persistent infections

library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach); 
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb); 
library(DBI); library(parallel); library(ggsignif); library(Rcpp)

# Loading prolonged infections, but not definitively persistent ------
infection_episode_rows <- readRDS(file="")
nrow(unique(infection_episode_rows[c("PERSON_ID", "infection_episode")]))

# Removing infection episodes with only one sequence, as they cannot be classified as persistent infections --------
# Step 1: Group by PERSON_ID and infection_episode, then count the occurrences
counted_episodes <- infection_episode_rows %>%
  group_by(PERSON_ID, infection_episode) %>%
  summarise(count = n(), .groups = 'drop')

# Step 2: Filter out episodes that occur only once
episodes_to_keep <- counted_episodes %>%
  filter(count > 1) %>%
  select(PERSON_ID, infection_episode)

# Step 3: Merge the original data with the filtered episodes to keep
filtered_persistent_infections <- infection_episode_rows %>%
  inner_join(episodes_to_keep, by = c("PERSON_ID", "infection_episode"))

# Now ensuring again that the infection episodes are long enough --------------
episode_duration <- filtered_persistent_infections %>%
  group_by(PERSON_ID, infection_episode) %>%
  summarise(
    first_date = min(DATESAMPLING),
    last_date = max(DATESAMPLING),
    .groups = 'drop'
  ) %>%
  mutate(duration_days = as.numeric(difftime(last_date, first_date, units = "days")))

episode_duration <- episode_duration %>%
  filter(first_date < as.Date("2023-01-01"))

#  Filter episodes where the duration is at least 26 days
episodes_to_keep <- episode_duration %>%
  filter(duration_days >= 26) %>%
  select(PERSON_ID, infection_episode)

# Merge the original data with the filtered episodes to keep
new_filtered_persistent_infections <- filtered_persistent_infections %>%
  inner_join(episodes_to_keep, by = c("PERSON_ID", "infection_episode"))

#View
selected_columns <- new_filtered_persistent_infections %>%
  select(PERSON_ID, DATESAMPLING, CT, lineage, scorpio_call, major_variant, major_scorpio_call, infection_episode)

# Now removing people where we don't have CT values for all rows in the infection_episode
filtered_variants <- new_filtered_persistent_infections %>%
  filter(grepl("^(Alpha|Delta|Omicron)", major_scorpio_call))

episodes_with_null_ct <- filtered_variants %>%
  filter(is.na(CT) | CT == "NULL") %>%
  select(PERSON_ID, infection_episode) %>%
  distinct()

# Step 2: Exclude those infection episodes from the dataset
cleaned_infection_data_without_CT <- filtered_variants

# Status update: Now we have definite potential persistent infections by scorpio call
saveRDS(cleaned_infection_data_without_CT, file="")



