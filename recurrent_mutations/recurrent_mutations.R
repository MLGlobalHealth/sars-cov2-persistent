# Finding recurrent mutations

library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate);
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle);
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(doParallel);
library(foreach); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(furrr); library(future); library(future.apply); library(lubridate); library(seqinr);
library(adegenet); library(ggplot2); library(viridis)

#Load data
all_without_ct_changes_first_last_placements <- readRDS("")
all_without_ct_changes_df <- readRDS("")

# Comparing first and last sequences only -------------

# Ensure mutations are split and unnested properly
mutations_df <- all_without_ct_changes_first_last_placements %>%
  mutate(mutation_list = strsplit(nucleotide_changes, "; ")) %>%  
  unnest(cols = c(mutation_list)) %>%                             
  filter(!is.na(mutation_list))                                   

# Extract position and mutation type, verify as numeric
mutations_df <- mutations_df %>%
  mutate(
    position = as.numeric(sub("[A-Z](\\d+)[A-Z]", "\\1", mutation_list)),  
    mutation_type = mutation_list                                          
  )

# Count mutations by position
mutation_counts_by_position <- mutations_df %>%
  group_by(position) %>%                             
  summarise(mutation_count = n()) %>%                
  arrange(desc(mutation_count))                      

# Count unique mutations
mutation_counts_by_type <- mutations_df %>%
  distinct(PERSON_ID, infection_episode, mutation_type) %>%                        
  group_by(mutation_type) %>%                                    
  summarise(count = n(), .groups = 'drop') %>%                 
  arrange(desc(count))                                 

write.csv(mutation_counts_by_position, "", row.names = FALSE)
write.csv(mutation_counts_by_type, "", row.names = FALSE)
library(openxlsx)
#write.xlsx(mutation_counts_by_position, "", rowNames = FALSE)
#write.xlsx(mutation_counts_by_type, "", rowNames = FALSE)


# Recurrent mutations for all consecutive sequences -------------

mutations_df_all <- all_without_ct_changes_df %>%
  mutate(mutation_list = strsplit(nucleotide_changes, "; ")) %>%  
  unnest(cols = c(mutation_list)) %>%                             
  filter(!is.na(mutation_list))                                   

mutations_df_all <- mutations_df_all %>%
  mutate(
    # Convert to G12345T-like format
    mutation_list = sub("Position_(\\d+):([A-Z])->([A-Z])", "\\2\\1\\3", mutation_list)
  ) %>%
  mutate(
    # Extract position as a number
    position = as.numeric(sub("[A-Z](\\d+)[A-Z]", "\\1", mutation_list)),
    # Keep the full mutation type
    mutation_type = mutation_list
  )

mutations_df_all <- mutations_df_all %>%
  mutate(
    position = as.numeric(sub("[A-Z](\\d+)[A-Z]", "\\1", mutation_list)),  
    mutation_type = mutation_list                                          
  )

# Count mutations by position
mutation_counts_by_position_all <- mutations_df_all %>%
  group_by(position) %>%                             
  summarise(mutation_count = n()) %>%                
  arrange(desc(mutation_count))                      

# Count unique mutations
mutation_counts_by_type_all <- mutations_df_all %>%
  distinct(PERSON_ID, infection_episode, mutation_type) %>%       
  group_by(mutation_type) %>%                                    
  summarise(count = n(), .groups = 'drop') %>%                 
  arrange(desc(count))                                

# Save results as CSV files if needed
write.csv(mutation_counts_by_position_all, "", row.names = FALSE)
write.csv(mutation_counts_by_type_all, "", row.names = FALSE)
library(openxlsx)
write.xlsx(mutation_counts_by_position_all, "", rowNames = FALSE)
write.xlsx(mutation_counts_by_type_all, "", rowNames = FALSE)






