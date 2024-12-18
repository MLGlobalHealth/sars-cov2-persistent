
# Synonymous vs. non-synonymous analysis for control / acute infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(ggpubr)

final_merged_new_id_QC_unique_metadata <- read.csv(file="")
conditional_subset <- readRDS(file="")
first_last_dnds_all_without_ct <- readRDS(file= "")

# Finding controls with multiple sequences -------

# Convert to data.table and filter for outcome == 0
subset_data <- as.data.table(conditional_subset)[outcome == 0]

# Keep the first occurrence of each unique combination of PERSON_ID and DATESAMPLING
subset_data <- subset_data[order(PERSON_ID, DATESAMPLING)][, .SD[1], by = .(PERSON_ID, DATESAMPLING)]

metadata <- as.data.table(final_merged_new_id_QC_unique_metadata)

# Ensure DATESAMPLING is in Date format
subset_data[, DATESAMPLING := as.Date(DATESAMPLING)]
metadata[, DATESAMPLING := as.Date(DATESAMPLING)]

# Create a unique identifier for PERSON_ID and event_id combinations
subset_data[, unique_id := paste(PERSON_ID, event_id, sep = "_")]

# Get unique strain values from the subset
unique_strains <- unique(subset_data$strain)
write.table(unique_strains, "unique_strains.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Initialize an empty list to store results
results_list <- list()
row_counter <- 0

# Loop through each unique combination of PERSON_ID and event_id
for (unique_row in subset_data[, unique(unique_id)]) {
  
  # Increment the counter
  row_counter <- row_counter + 1
  
  # Print progress every 5000 rows
  if (row_counter %% 5000 == 0) {
    cat("Processed", row_counter, "rows...\n")
  }
  
  # Get PERSON_ID and event_id from the unique identifier
  ids <- unlist(strsplit(unique_row, "_"))
  PERSON_ID <- ids[1]
  event_id <- ids[2]
  
  # Filter subset_data for the current unique_id
  current_sample <- subset_data[unique_id == unique_row]
  
  # Ensure scalar values for filtering
  current_date <- current_sample$DATESAMPLING[1]
  
  # Filter metadata for rows within the date range and matching PERSON_ID
  linked_rows <- metadata[
    PERSON_ID == current_sample$PERSON_ID[1] & 
      DATESAMPLING > (current_date + 1) & 
      DATESAMPLING < (current_date + 26)
  ]
  
  # If there are linked rows in metadata, process them
  if (nrow(linked_rows) > 0) {
    # Create a dummy variable for infection episode link
    linked_rows[, infection_episode := unique_row]  # Mark linked rows with infection episode ID
    
    # Add event_id to the rows from subset_data and set infection_episode variable
    current_sample[, infection_episode := unique_row]
    current_sample[, event_id := event_id]  # Add event_id from subset_data
    
    # Keep only relevant columns: strain, DATESAMPLING, PERSON_ID, event_id (for current_sample only)
    linked_rows <- linked_rows[, .(PERSON_ID, DATESAMPLING, strain, infection_episode)]
    current_sample <- current_sample[, .(PERSON_ID, DATESAMPLING, strain, event_id, infection_episode)]
    
    # Add the reference sample (from subset_data) to the linked rows
    linked_rows <- rbind(linked_rows, current_sample, fill = TRUE)
    
    # Append to results list
    results_list[[length(results_list) + 1]] <- linked_rows
  }
}

# Combine all results into a single data.table
if (length(results_list) > 0) {
  matched_data <- rbindlist(results_list, fill = TRUE)
  print(matched_data)
} else {
  print("No matched rows found.")
}

# Display the first few rows of the final result
head(matched_data)

# Reduce matched_data to two rows per infection_episode: earliest and latest DATESAMPLING
reduced_data <- matched_data[
  , .SD[which.min(DATESAMPLING)], by = infection_episode][
    , .SD[which.max(DATESAMPLING)], by = infection_episode]

# Ensure only one row for the earliest and one for the latest per infection_episode
reduced_data <- matched_data[
  , .SD[c(which.min(DATESAMPLING), which.max(DATESAMPLING))], by = infection_episode]

# Add a column for time since the first strain (in days)
reduced_data[, time_since_first_strain := as.numeric(DATESAMPLING - min(DATESAMPLING)), by = infection_episode]

# Remove rows with NA in infection_episode
reduced_data <- reduced_data[!is.na(infection_episode)]

# Get strain values from the subset
reduced_control_strains <- unique(reduced_data$strain)
write.table(
  reduced_control_strains, 
  "reduced_control_strains.txt", 
  quote = FALSE,  # No quotes around values
  row.names = FALSE,  # Exclude row numbers
  col.names = FALSE,  # Exclude column names
  sep = "\n"  # Ensure one strain per line
)

# Run get_reduced_control_sequences.py

# Run dnds analysis -----------

#Import back FASTA files
alignment_control_strains <- read.alignment(file = "reduced_control_sequences.fasta", format = "fasta", forceToLower = TRUE)

gene_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
)

create_alignment_for_episodes <- function(alignment, reduced_data) {
  # Extract sequence names from the alignment
  sequence_names <- alignment$nam
  
  # Initialize a list to store alignments for each episode
  episode_alignments <- list()
  
  # Loop through each infection_episode
  for (episode in unique(reduced_data$infection_episode)) {
    # Subset data for the current episode
    episode_data <- reduced_data[infection_episode == episode]
    
    # Ensure exactly two rows (one for earliest and one for latest)
    if (nrow(episode_data) != 2) {
      warning(paste("Skipping episode:", episode, "due to incomplete data"))
      next
    }
    
    # Identify strains for earliest and latest dates
    seq1 <- episode_data[which.min(DATESAMPLING), strain]
    seq2 <- episode_data[which.max(DATESAMPLING), strain]
    
    # Find indices of these sequences in the alignment
    indices_to_include <- which(sequence_names %in% c(seq1, seq2))
    
    # If both sequences are not found, skip this episode
    if (length(indices_to_include) < 2) {
      warning(paste("Sequences not found in alignment for episode:", episode))
      next
    }
    
    # Create a new alignment object with only the required sequences
    subset_alignment <- list(
      nam = sequence_names[indices_to_include],
      seq = as.character(alignment$seq[indices_to_include]),
      com = alignment$com[indices_to_include],
      nb = length(indices_to_include)
    )
    class(subset_alignment) <- "alignment"  # Ensure the class is set to 'alignment'
    
    # Store the alignment in the list using the episode ID
    episode_alignments[[as.character(episode)]] <- subset_alignment
  }
  
  return(episode_alignments)
}

split_alignment_by_gene <- function(alignment, gene_regions) {
  # Extract sequence names and sequences
  sequence_names <- alignment$nam
  sequences <- alignment$seq
  
  # Initialize a list to store alignments for each gene region
  gene_alignments <- list()
  
  # Loop through each gene region to split the alignment
  for (i in 1:nrow(gene_regions)) {
    gene_name <- gene_regions$region[i]
    start_pos <- gene_regions$start[i]
    end_pos <- gene_regions$end[i]
    
    # Create a new alignment object for this gene region
    gene_sequences <- lapply(sequences, function(seq) {
      if (nchar(seq) >= end_pos) {
        substring(seq, start_pos, end_pos)
      } else {
        # Handle shorter sequences by returning the original sequence
        warning(paste("Sequence length shorter than expected for gene", gene_name))
        seq
      }
    })
    
    # Create the gene alignment object
    gene_alignment <- list(
      nam = sequence_names,
      seq = gene_sequences,
      com = alignment$com,  # Assuming comments are the same
      nb = length(sequence_names)
    )
    class(gene_alignment) <- "alignment"  # Ensure the class is set to 'alignment'
    
    # Store the gene alignment in the list
    gene_alignments[[gene_name]] <- gene_alignment
  }
  
  return(gene_alignments)
}

run_kaks_on_gene_alignments <- function(gene_alignments) {
  results <- list()
  
  for (gene_name in names(gene_alignments)) {
    gene_alignment <- gene_alignments[[gene_name]]
    
    # Run Ka/Ks analysis and handle potential errors
    tryCatch({
      print(paste("Running Ka/Ks analysis for gene:", gene_name))
      kaks_result <- kaks(gene_alignment, verbose = TRUE, debug = FALSE, forceUpperCase = TRUE, rmgap = TRUE)
      results[[gene_name]] <- kaks_result
    }, error = function(e) {
      cat("Error during Ka/Ks analysis for gene", gene_name, ":", e$message, "\n")
    })
  }
  
  return(results)
}

extract_clean_dnds_counts_to_df <- function(results) {
  # Initialize an empty list to hold the rows for the dataframe
  rows <- list()
  
  # Loop through each pair of results
  for (pair_name in names(results)) {
    # Get the results for each gene in the pair
    gene_results <- results[[pair_name]]
    
    # Loop through each gene in the pair
    for (gene_name in names(gene_results)) {
      gene_result <- gene_results[[gene_name]]
      
      # Extract values from the lists or data frames for each field
      ka_value <- gene_result$ka[1]  # Getting the first value (since there's only one row per sequence)
      ks_value <- gene_result$ks[1]
      vka_value <- gene_result$vka[1]
      vks_value <- gene_result$vks[1]
      l0_value <- gene_result$l0[1]
      l2_value <- gene_result$l2[1]
      l4_value <- gene_result$l4[1]
      a0_value <- gene_result$a0[1]
      a2_value <- gene_result$a2[1]
      a4_value <- gene_result$a4[1]
      b0_value <- gene_result$b0[1]
      b2_value <- gene_result$b2[1]
      b4_value <- gene_result$b4[1]
      checksuml_value <- gene_result$checksuml[1]
      
      # Add a row for the current gene
      rows <- append(rows, list(data.frame(
        pair = pair_name,
        gene = gene_name,
        ka = ka_value,
        ks = ks_value,
        vka = vka_value,
        vks = vks_value,
        l0 = l0_value,
        l2 = l2_value,
        l4 = l4_value,
        a0 = a0_value,
        a2 = a2_value,
        a4 = a4_value,
        b0 = b0_value,
        b2 = b2_value,
        b4 = b4_value,
        checksuml = checksuml_value
      )))
    }
  }
  
  # Combine all the rows into a single dataframe
  df <- do.call(rbind, rows)
  return(df)
}



# Main analysis pipeline for all pairs (without CT)
class(alignment_control_strains) <- "alignment"

# Create alignments for each pair based on metadata
pair_alignments <- create_alignment_for_episodes(alignment_control_strains, reduced_data)

# Split each pair's alignment into gene regions
gene_alignments <- lapply(pair_alignments, function(example_alignment) {
  class(example_alignment) <- "alignment"  # Ensure the class is set to 'alignment'
  split_alignment_by_gene(example_alignment, gene_regions)
})

# Run Ka/Ks analysis for each gene region
kaks_results <- lapply(gene_alignments, run_kaks_on_gene_alignments)

# Extract and clean dN/dS counts
dnds_counts_controls <- extract_clean_dnds_counts_to_df(kaks_results)
saveRDS(dnds_counts_controls, file="dnds_counts_controls.RDS")









