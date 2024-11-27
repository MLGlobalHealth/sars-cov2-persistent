#Analyzing and visualizing persistent infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate);
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle);
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(doParallel);
library(foreach); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(furrr); library(future); library(future.apply); library(lubridate); library(seqinr);
library(adegenet); library(ggplot2); library(viridis)

# Testing, vaccination and general metadata
COVID_TEST_filtered <- readRDS(file = "")
COVID_VACC_filtered <- readRDS(file = "")
lifelines_koen_filtered <- readRDS(file = "")
lifelines_filtered <- readRDS(file = "")

# Persistent infection metadata
all_without_ct_metadata <- readRDS(file = "all_without_ct_metadata.rds")
all_without_ct_metadata_unique <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  slice_min(DATESAMPLING) %>%
  ungroup()


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



# Importing fasta data  ---------------------------------------------

alignment_all_persistent_infection_fasta_seqs_without_CT <- read.alignment(file = "all_possible_persistent_infections_persistent_sequences_without_CT.fasta", format = "fasta", forceToLower = TRUE)
all_possible_persistent_fasta_seqs_without_CT <- read.fasta(file = "all_possible_persistent_infections_persistent_sequences_without_CT.fasta", as.string = TRUE)
fasta_df_all_possible_persistent_without_CT <- process_sequences(all_possible_persistent_fasta_seqs_without_CT)


# Mutation locations -----------------------

# Function to count sequence changes
count_sequence_changes <- function(person_episode_pair, metadata_df, fasta_df) {
  
  # Extract PERSON_ID and infection_episode from the input pair
  person_id <- person_episode_pair$PERSON_ID
  episode <- person_episode_pair$infection_episode
  
  print(person_id)
  
  # Filter metadata for this person and episode
  person_data <- metadata_df %>% 
    filter(PERSON_ID == person_id, infection_episode == episode) %>%
    arrange(DATESAMPLING)
  
  # Extract strain names, dates, and major_scorpio_call
  strains <- person_data$strain
  dates <- person_data$DATESAMPLING
  major_scorpio_call <- person_data$major_scorpio_call
  
  # Retrieve the sequences corresponding to the strains
  sequences <- fasta_df %>% 
    filter(strain %in% strains)
  
  # Ensure sequences are ordered by DATESAMPLING
  sequences <- sequences[match(strains, sequences$strain), ]
  
  # Remove the 'strain' column to work with nucleotide positions
  nucleotide_positions <- sequences[, -1, with = FALSE]
  
  # Get position column names and their numeric indices
  position_columns <- names(nucleotide_positions)[grepl("Position_", names(nucleotide_positions))]
  position_numbers <- as.numeric(sub("Position_", "", position_columns))
  
  # Initialize vectors to store change counts, positions, nucleotide changes, and days between sequences
  change_counts <- integer()
  change_positions <- list()
  nucleotide_changes <- list()
  days_between <- numeric()
  
  # Compare consecutive sequences
  for (j in 1:(nrow(nucleotide_positions) - 1)) {
    seq1 <- nucleotide_positions[j, ..position_columns]
    seq2 <- nucleotide_positions[j + 1, ..position_columns]
    
    # Initialize counter for changes, a list for positions, and nucleotide change details
    count_changes <- 0
    positions <- integer()
    details <- character()
    
    # Loop through each position and compare
    for (k in seq_len(ncol(seq1))) {
      pos <- position_columns[k]
      # Check if both positions have valid nucleotide values
      if (seq1[[pos]] %in% c("A", "T", "C", "G") & seq2[[pos]] %in% c("A", "T", "C", "G")) {
        if (seq1[[pos]] != seq2[[pos]]) {
          count_changes <- count_changes + 1
          positions <- c(positions, position_numbers[k])
          details <- c(details, paste0("Position_", position_numbers[k], ":", seq1[[pos]], "->", seq2[[pos]]))
        }
      }
    }
    
    # Store the number of changes, positions, nucleotide change details, and major_scorpio_call
    change_counts <- c(change_counts, count_changes)
    change_positions[[j]] <- positions
    nucleotide_changes[[j]] <- paste(details, collapse = "; ")
    major_scorpio_call_this_pair <- major_scorpio_call[j]
    
    # Calculate days between the current and next sequence
    days_between <- c(days_between, as.numeric(difftime(dates[j + 1], dates[j], units = "days")))
  }
  
  # Get the first sampling date for the infection episode
  first_sampling_date <- person_data$DATESAMPLING[1]
  # Get the first date from the relevant sequences
  first_sequence_date <- dates[1]  # The first date in the sequence list
  
  # Combine results into a data frame
  results <- data.table(
    PERSON_ID = person_id,
    infection_episode = episode,
    infection_episode_start_date = first_sampling_date,
    first_sequence_date = dates[1:(nrow(nucleotide_positions) - 1)],  # Date of the first sequence in the comparison
    comparison = paste0("Sequence_", 1:(nrow(nucleotide_positions) - 1), "_to_Sequence_", 2:nrow(nucleotide_positions)),
    change_count = change_counts,
    change_positions = sapply(change_positions, function(pos) paste(pos, collapse = ", ")),
    nucleotide_changes = sapply(nucleotide_changes, function(details) paste(details, collapse = "; ")),
    days_between_sequences = days_between,
    major_scorpio_call = major_scorpio_call[1:(nrow(nucleotide_positions) - 1)]  # Add major_scorpio_call to the results
  )
  
  return(results)
}


# Function to parallelize
parallel_count_sequence_changes <- function(metadata_df, fasta_df) {
  # Get unique PERSON_ID and infection_episode combinations
  person_episodes <- metadata_df %>% 
    distinct(PERSON_ID, infection_episode) %>%
    as.data.table()  # Convert to data.table for compatibility with future_lapply
  
  # Set up parallel processing
  plan(multisession, workers = 16)
  
  # Apply the count_sequence_changes function in parallel
  results_list <- future.apply::future_lapply(seq_len(nrow(person_episodes)), function(i) {
    count_sequence_changes(person_episodes[i], metadata_df, fasta_df)
  })
  
  # Combine all results into a single data table
  final_results <- rbindlist(results_list)
  
  return(final_results)
}

options(future.globals.maxSize = 2 * 1024 * 1024^2) # 2 GB


# Run the parallelized function 
all_without_ct_changes_df <- parallel_count_sequence_changes(all_without_ct_metadata, fasta_df_all_possible_persistent_without_CT)
saveRDS(all_without_ct_changes_df, file= "all_without_ct_changes_df.RDS")
