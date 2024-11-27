# Finding the final persistent infections using the rare SNP threshold

library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach); 
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb); 
library(DBI); library(parallel); library(ggsignif); library(Rcpp)
library(future);library(future.apply);library(seqinr)

#Loading data ------
conn <- DSTora::OraGenvej('', dbuser = '')
lifelines <- dbReadTable(conn = conn,
                         name = "",
                         schema = '')

lifelines_filtered <- lifelines %>%
  filter(!is.na(BIRTHDAY))

lifelines_koen <- dbReadTable(conn = conn,
                              name = "",
                              schema = '')
lifelines_koen_filtered <- lifelines_koen %>%
  filter(!is.na(KOEN))

all_possible_persistent_infections_without_CT <- readRDS(file="")

# Removing duplicate rows and ensuring that individuals have metadata --------
all_possible_persistent_infections_without_CT$has_birthday <- all_possible_persistent_infections_without_CT$PERSON_ID %in% lifelines_filtered$PERSON_ID
all_possible_persistent_infections_without_CT$has_koen <- all_possible_persistent_infections_without_CT$PERSON_ID %in% lifelines_koen_filtered$PERSON_ID
all_possible_persistent_infections_without_CT <- all_possible_persistent_infections_without_CT %>%
  filter(has_birthday | has_koen) %>%
  distinct(PERSON_ID, DATESAMPLING, .keep_all = TRUE)


#Finding persistent infections with rare SNPs and loading data -----
all_possible_persistent_infections_without_CT <- all_possible_persistent_infections_without_CT %>%
  group_by(strain) %>%
  slice(1) %>%
  ungroup()

#Getting relevant FASTA files -----------
all_possible_persistent_infections_strain_list_without_CT <- unique(all_possible_persistent_infections_without_CT$strain)
writeLines(all_possible_persistent_infections_strain_list_without_CT, "")

# run subset_all_potential_infections_persistent_without_CT_FASTA.py python code to extract sequences

# Load the relevant sequences
all_possible_persistent_fasta_seqs_without_CT <- read.fasta(file = "all_possible_persistent_infections_persistent_sequences_without_CT.fasta", as.string = TRUE)

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
  
  # Optionally remove the sequence column
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

fasta_df_all_possible_persistent_with_CT <- process_sequences(all_possible_persistent_fasta_seqs_with_CT)
fasta_df_all_possible_persistent_without_CT <- process_sequences(all_possible_persistent_fasta_seqs_without_CT)


# Load nucleotide proportions for different variants --------------
omicronBA1_proportions <- fread("omicronBA1_nucleotide_proportions.csv")
omicronBA2_proportions <- fread("omicronBA2_nucleotide_proportions.csv")
omicronBA5_proportions <- fread("omicronBA5_nucleotide_proportions.csv")
alpha_nucleotide_proportions <- fread("alpha_nucleotide_proportions.csv")
delta_nucleotide_proportions <- fread("delta_nucleotide_proportions.csv")

# Function to normalize selected columns and ensure only relevant ones are kept
normalize_and_select_columns <- function(dt) {
  # Keep only the relevant columns
  dt <- dt[, .(Position, a, t, c, g)]
  
  # Calculate row-wise sums of 'a', 't', 'c', and 'g'
  dt[, sum_cols := a + t + c + g]
  
  # Normalize 'a', 't', 'c', and 'g'
  dt[, `:=`(
    a = a / sum_cols,
    t = t / sum_cols,
    c = c / sum_cols,
    g = g / sum_cols
  )]
  
  # Remove the temporary sum column
  dt[, sum_cols := NULL]
  
  return(dt)
}

# Apply normalization to each dataset
omicronBA1_proportions <- normalize_and_select_columns(omicronBA1_proportions)
omicronBA2_proportions <- normalize_and_select_columns(omicronBA2_proportions)
omicronBA5_proportions <- normalize_and_select_columns(omicronBA5_proportions)
alpha_proportions <- normalize_and_select_columns(alpha_nucleotide_proportions)
delta_proportions <- normalize_and_select_columns(delta_nucleotide_proportions)


# Setting up dataframes to find rare SNPs ----------------------------
options(future.globals.maxSize = 1 * 1024 * 1024^2) # 1 GB

# Define the threshold
threshold <- 0.1090691

# Map scorpio_call to the relevant nucleotide proportions
scorpio_thresholds <- list(
  "Omicron BA.1" = omicronBA1_proportions,
  "Omicron BA.2" = omicronBA2_proportions,
  "Omicron BA.5" = omicronBA5_proportions,
  "Delta" = delta_nucleotide_proportions,
  "Alpha" = alpha_nucleotide_proportions
)

# Filter the dataframe to include only the specified variants
valid_variants <- names(scorpio_thresholds)
filtered_df_all_possible_persistent_with_CT <- all_possible_persistent_infections_with_CT[all_possible_persistent_infections_with_CT$major_scorpio_call %in% valid_variants, ]
filtered_df_all_possible_persistent_without_CT <- all_possible_persistent_infections_without_CT[all_possible_persistent_infections_without_CT$major_scorpio_call %in% valid_variants, ]


identify_rare_snps <- function(df, fasta_df, scorpio_thresholds, threshold) {
  # Convert df and fasta_df to data.table if not already
  plan(multisession, workers = 16)
  df <- as.data.table(df)
  fasta_df <- as.data.table(fasta_df)
  
  # Prepare a data.table to store results
  results <- data.table(
    PERSON_ID = integer(),
    major_scorpio_call = character(),
    infection_episode = integer(),
    rare_snp_count = integer(),
    rare_snp_positions = character(),
    rare_snp_proportions = character()
  )
  
  # Extract unique person IDs
  person_ids <- unique(df$PERSON_ID)
  
  # Define the parallelized operation for each person
  person_results <- future.apply::future_lapply(person_ids, function(person_id) {
    # Filter data for the current person
    print(person_id)
    person_sequences <- df[PERSON_ID == person_id]
    scorpio_calls <- unique(person_sequences$major_scorpio_call)
    infection_episodes <- unique(person_sequences$infection_episode)
    
    # Prepare a list to store the results for this person
    person_results_list <- list()
    
    for (episode in infection_episodes) {
      episode_sequences <- person_sequences[infection_episode == episode]
      major_scorpio_call <- unique(episode_sequences$major_scorpio_call)
      
      # Get the nucleotide proportions for the current major_scorpio_call
      nucleotide_proportions <- scorpio_thresholds[[major_scorpio_call]]
      
      # Initialize variables for rare SNP results
      rare_snp_count <- 0
      rare_snp_positions <- c()
      rare_snp_proportions <- c()
      
      # Extract strains for the current person
      strains <- episode_sequences$strain
      
      # Prepare a table for nucleotides at each position
      nucleotide_data <- fasta_df[strain %in% strains, .SD, .SDcols = names(fasta_df)[grepl("Position_", names(fasta_df))]]
      
      # Ensure nucleotide_data is a data.table with only relevant columns
      if (nrow(nucleotide_data) > 0) {
        # Iterate over each position
        for (position in names(nucleotide_data)) {  # Now position is directly from nucleotide_data
          position_number <- as.numeric(sub("Position_", "", position))
          
          # Extract nucleotides at this position for the person's strains
          nucleotides <- nucleotide_data[[position]]
          
          # Check if all nucleotides at this position are A, T, C, or G and are the same
          if (length(unique(nucleotides)) == 1 && nucleotides[1] %in% c("A", "T", "C", "G")) {
            nucleotide <- nucleotides[1]
            nucleotide_lower <- tolower(nucleotide)  # Convert to lowercase for comparison
            
            # Get the proportion of this nucleotide in the total population
            proportion_row <- nucleotide_proportions[Position == position_number]
            
            if (nrow(proportion_row) > 0) {
              # Check that the proportion is being looked up in the correct column
              proportion_column <- match(nucleotide_lower, names(proportion_row))
              proportion <- proportion_row[[proportion_column]]
              
              # Handle NA values in proportion
              if (!is.na(proportion)) {
                
                if (proportion < threshold) {
                  rare_snp_count <- rare_snp_count + 1
                  rare_snp_positions <- c(rare_snp_positions, position_number)
                  rare_snp_proportions <- c(rare_snp_proportions, proportion)
                }
              }
            }
          }
        }
      }
      
      # Add results to the person_results_list
      person_results_list[[length(person_results_list) + 1]] <- data.table(
        PERSON_ID = person_id,
        major_scorpio_call = major_scorpio_call,
        infection_episode = episode,
        rare_snp_count = rare_snp_count,
        rare_snp_positions = paste(rare_snp_positions, collapse = ", "),
        rare_snp_proportions = paste(rare_snp_proportions, collapse = ", ")
      )
    }
    
    return(rbindlist(person_results_list))
  })
  
  # Combine all results into a single data.table
  results <- rbindlist(person_results)
  
  return(results)
}

results_all_possible_persistent_without_CT <- identify_rare_snps(filtered_df_all_possible_persistent_without_CT, fasta_df_all_possible_persistent_without_CT, scorpio_thresholds, threshold)
saveRDS(results_all_possible_persistent_without_CT, file="results_all_possible_persistent_without_CT.RDS")

# Reset the future plan to default
plan(sequential)




















