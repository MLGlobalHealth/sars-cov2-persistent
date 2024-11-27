# Identifying mutation sites and (non)synonymous changes
# Comparing first and last sequences in each infection

# packages
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
all_without_ct_metadata <- readRDS(file = "")
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



# Importing fasta sequences ---------------------------------------------
alignment_all_persistent_infection_fasta_seqs_without_CT <- read.alignment(file = "", format = "fasta", forceToLower = TRUE)
all_possible_persistent_fasta_seqs_without_CT <- read.fasta(file = "", as.string = TRUE)
fasta_df_all_possible_persistent_without_CT <- process_sequences(all_possible_persistent_fasta_seqs_without_CT)


# Identifying placement of mutations -----------------------

count_sequence_changes <- function(person_episode_pair, metadata_df, fasta_df) {
  
  # Extract PERSON_ID and infection_episode from the input pair
  person_id <- person_episode_pair$PERSON_ID
  episode <- person_episode_pair$infection_episode
  
  print(person_id)
  
  # Filter metadata for this person and episode, removing duplicate dates
  person_data <- metadata_df %>% 
    filter(PERSON_ID == person_id, infection_episode == episode) %>%
    arrange(DATESAMPLING) %>%
    distinct(DATESAMPLING, .keep_all = TRUE)  # Keeps only one sample per date
  
  # Handle potential duplicate sequences on the last date
  last_date <- max(person_data$DATESAMPLING)  # Find the last date
  person_data <- person_data %>%
    group_by(DATESAMPLING) %>%
    mutate(is_last_date = (DATESAMPLING == last_date) & (n() > 1)) %>%
    filter(!(is_last_date & row_number() > 1)) %>%
    ungroup() %>%
    select(-is_last_date)  # Remove helper column
  
  # Exit if fewer than 2 sequences after handling duplicates
  if (nrow(person_data) < 2) {
    return(NULL)
  }
  
  # Extract strain names and dates
  strains <- person_data$strain
  dates <- person_data$DATESAMPLING
  
  # Retrieve sequences for the first and last strains in the episode
  sequences <- fasta_df %>% 
    filter(strain %in% c(strains[1], strains[nrow(person_data)])) %>%
    arrange(match(strain, c(strains[1], strains[nrow(person_data)])))  # Order: first and last sequence
  
  # Ensure we have exactly two sequences for comparison
  if (nrow(sequences) != 2) {
    return(NULL)
  }
  
  # Extract nucleotide position columns
  nucleotide_positions <- sequences[, -1, with = FALSE]
  position_columns <- names(nucleotide_positions)[grepl("Position_", names(nucleotide_positions))]
  position_numbers <- as.numeric(sub("Position_", "", position_columns))
  
  # Initialize lists to store change positions and detailed changes
  change_positions <- integer()
  nucleotide_changes <- character()
  
  # Compare first and last sequences
  seq1 <- nucleotide_positions[1, ..position_columns]
  seq2 <- nucleotide_positions[2, ..position_columns]
  
  for (k in seq_len(ncol(seq1))) {
    pos <- position_columns[k]
    if (seq1[[pos]] %in% c("A", "T", "C", "G") & seq2[[pos]] %in% c("A", "T", "C", "G")) {
      if (seq1[[pos]] != seq2[[pos]]) {
        # Store the position and format mutation as "G12345T"
        change_positions <- c(change_positions, position_numbers[k])
        nucleotide_changes <- c(nucleotide_changes, paste0(seq1[[pos]], position_numbers[k], seq2[[pos]]))
      }
    }
  }
  
  # Calculate duration of the infection episode
  total_days_between <- as.numeric(difftime(dates[nrow(person_data)], dates[1], units = "days"))
  
  # Return the results as a data table
  results <- data.table(
    PERSON_ID = person_id,
    infection_episode = episode,
    infection_episode_start_date = dates[1],
    infection_episode_end_date = dates[nrow(person_data)],
    total_days_between_sequences = total_days_between,
    total_change_count = length(nucleotide_changes),  # Total count of changes
    nucleotide_changes = paste(nucleotide_changes, collapse = "; "),  # Concatenated changes in "G12345T" format
    change_positions = paste(change_positions, collapse = ", ")  # Positions of changes
  )
  
  return(results)
}


# Apply the updated function with parallel processing
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
  final_results <- rbindlist(results_list, fill = TRUE)  # Handle any NULL values returned by rbindlist
  
  return(final_results)
}

# Execute the function and save results
all_without_ct_changes_first_last_placements <- parallel_count_sequence_changes(all_without_ct_metadata, fasta_df_all_possible_persistent_without_CT)
saveRDS(all_without_ct_changes_first_last_placements, file = "")

# Parse and count mutations by position ----------------------------------------

# Extract position and mutation details for each row
# Step 1: Ensure mutations are split and unnested properly
mutations_df <- all_without_ct_changes_first_last_placements %>%
  mutate(mutation_list = strsplit(nucleotide_changes, "; ")) %>%  # Split into individual mutations
  unnest(cols = c(mutation_list)) %>%                             # Unnest into individual rows
  filter(!is.na(mutation_list))                                   # Filter out any empty entries

# Step 2: Extract position and mutation type, verify as numeric
mutations_df <- mutations_df %>%
  mutate(
    position = as.numeric(sub("[A-Z](\\d+)[A-Z]", "\\1", mutation_list)),  # Extract position as a number
    mutation_type = mutation_list                                          # Keep full mutation type (e.g., "G12345T")
  )

# Confirm that position and mutation_type columns are correctly formatted
str(mutations_df)

# Step 3: Count mutations by position
mutation_counts_by_position <- mutations_df %>%
  group_by(position) %>%                             # Group by position
  summarise(mutation_count = n()) %>%                # Count occurrences per position
  arrange(desc(mutation_count))                      # Order by mutation frequency

# Step 4: Count unique mutations
mutation_counts_by_type <- mutations_df %>%
  distinct(PERSON_ID, infection_episode, mutation_type) %>%                        # Ensure only unique individuals contribute to counts
  group_by(mutation_type) %>%                                    # Group by unique mutation
  summarise(count = n(), .groups = 'drop') %>%                 # Count each unique mutation
  arrange(desc(count))                                 # Order by frequency

# Display result previews for both tables
head(mutation_counts_by_position)
head(mutation_counts_by_type)

# Save results as CSV files if needed
write.csv(mutation_counts_by_position, "", row.names = FALSE)
write.csv(mutation_counts_by_type, "", row.names = FALSE)




# Calculating rate  -------
num_sites <- 29891
# Add 'rate' and 'rate_per_site' columns to all_without_ct_changes_first_last_placements
all_without_ct_changes_first_last_placements$rate <- (all_without_ct_changes_first_last_placements$total_change_count / all_without_ct_changes_first_last_placements$total_days_between_sequences) * 365.25
all_without_ct_changes_first_last_placements$rate_per_site <- (all_without_ct_changes_first_last_placements$total_change_count / (all_without_ct_changes_first_last_placements$total_days_between_sequences * num_sites)) * 365.25
saveRDS(all_without_ct_changes_first_last_placements, file = "")
# Reading in metadata for these individuals
subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson <- readRDS(file="")

# Select the relevant columns from the metadata dataframe
columns_to_keep <- c("PERSON_ID", "infection_episode", "num_vacc_at_infection", 
                     "num_vacc_at_infection_cat", "age_at_infection", 
                     "age_group", "IMID_diagnosis", "IMID_codes", 
                     "IMID_dates", "IMID_source", 
                     "IMID_diagnosis_prepandemic", "diabetes_diagnosis", 
                     "diabetes_codes", "diabetes_dates", 
                     "diabetes_source", "diabetes_diagnosis_prepandemic", 
                     "immunosuppression_diagnosis", "immunosuppression_codes", 
                     "immunosuppression_dates", "immunosuppression_source", 
                     "immunosuppression_diagnosis_prepandemic", "autoimmune_diagnosis", 
                     "autoimmune_codes", "autoimmune_dates", 
                     "autoimmune_source", "autoimmune_diagnosis_prepandemic", 
                     "outcome", "charlson.index.5yrs", 
                     "charlson.index.10yrs", "KOEN", "major_scorpio_call")

# Subset the metadata dataframe
metadata_subset <- subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson %>%
  select(all_of(columns_to_keep))

#Join the two data frames
all_without_ct_changes_first_last_placements <- all_without_ct_changes_first_last_placements %>%
  left_join(metadata_subset, by = c("PERSON_ID", "infection_episode"))


#Save for final time -----------------
saveRDS(all_without_ct_changes_first_last_placements, file = "")






# Finding synonymous and non-synonymous changes --------
# Importing data and FASTA files if necessary again
all_without_ct_metadata <- readRDS(file = "")
alignment_all_persistent_infection_fasta_seqs_without_CT <- read.alignment(file = "", format = "fasta", forceToLower = TRUE)
all_possible_persistent_fasta_seqs_without_CT <- read.fasta(file = "", as.string = TRUE)

# Relevant functions and baseline data frames 

gene_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
)


# Function to create a comparison dataframe with first and last episodes for each PERSON_ID and infection_episode
create_comparison_pairs <- function(metadata_df) {
  # Ensure metadata is ordered correctly by PERSON_ID, infection_episode, and DATESAMPLING
  metadata_df <- metadata_df %>%
    arrange(PERSON_ID, infection_episode, DATESAMPLING) %>%
    group_by(PERSON_ID, infection_episode) %>%
    summarize(
      first_strain = first(strain),
      first_DATESAMPLING = first(DATESAMPLING),
      last_strain = last(strain),
      last_DATESAMPLING = last(DATESAMPLING),
      days_between_samples = as.numeric(difftime(last_DATESAMPLING, first_DATESAMPLING, units = "days")),
      .groups = 'drop'  # Drop grouping to return a normal data frame
    ) %>%
    mutate(
      pair_id = paste(PERSON_ID, infection_episode, sep = "_"),
      comparison = "First_to_Last"
    ) %>%
    select(PERSON_ID, infection_episode, first_strain, first_DATESAMPLING, last_strain, last_DATESAMPLING, days_between_samples, pair_id, comparison)
  
  return(metadata_df)
}

# Function to create an alignment for each unique pair from comparison_pairs_df
create_alignment_for_pairs <- function(alignment, comparison_pairs_df) {
  # Extract sequence names from the alignment
  sequence_names <- alignment$nam
  
  # Initialize a list to store alignments for each pair
  pair_alignments <- list()
  
  # Loop through each row in the comparison_pairs_df
  for (i in 1:nrow(comparison_pairs_df)) {
    seq1 <- comparison_pairs_df$first_strain[i]  # Use 'first_strain' here
    seq2 <- comparison_pairs_df$last_strain[i]   # Use 'last_strain' here
    
    # Find indices of these sequences in the alignment
    indices_to_include <- which(sequence_names %in% c(seq1, seq2))
    
    # Check if both sequences are present
    if (length(indices_to_include) < 2) {
      warning(paste("Insufficient sequences found in alignment for pair:", seq1, "and", seq2))
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
    
    # Use the unique pair identifier from the comparison_pairs_df
    pair_id <- comparison_pairs_df$pair_id[i]
    pair_alignments[[pair_id]] <- subset_alignment
  }
  
  return(pair_alignments)
}


# Function to split alignment by gene for each pair
split_alignment_by_gene_for_pairs <- function(pair_alignments, gene_regions, comparison_pairs_df) {
  # Initialize a list to store gene-region-specific alignments for each pair
  all_gene_splits_by_pair <- list()
  
  # Loop through each pair alignment
  for (pair_id in names(pair_alignments)) {
    pair_alignment <- pair_alignments[[pair_id]]
    
    # Extract sequence names and sequences for this pair
    sequence_names <- pair_alignment$nam
    sequences <- pair_alignment$seq
    
    # Initialize a list to store alignments for each gene region in this pair
    gene_alignments_for_pair <- list()
    
    # Loop through each gene region
    for (i in 1:nrow(gene_regions)) {
      gene_name <- gene_regions$region[i]
      start_pos <- gene_regions$start[i]
      end_pos <- gene_regions$end[i]
      
      # Create a new alignment object for the gene region
      gene_sequences <- lapply(sequences, function(seq) {
        if (nchar(seq) >= end_pos) {
          substring(seq, start_pos, end_pos)
        } else {
          # Handle cases where sequence length is shorter than end_pos
          warning(paste("Sequence length shorter than expected for gene", gene_name, "in pair", pair_id))
          seq  # Return the original sequence if it's too short
        }
      })
      
      # Create the gene alignment object
      gene_alignment <- list(
        nam = sequence_names,
        seq = gene_sequences,
        com = pair_alignment$com,  # Assuming comments are the same
        nb = length(sequence_names)
      )
      class(gene_alignment) <- "alignment"  # Ensure the class is set to 'alignment'
      
      # Store the gene alignment in the list for this pair
      gene_alignments_for_pair[[gene_name]] <- gene_alignment
    }
    
    # Add the gene-region-specific alignments for this pair to the overall list
    all_gene_splits_by_pair[[pair_id]] <- gene_alignments_for_pair
  }
  
  return(all_gene_splits_by_pair)
}


# Helper function to clean sequence and convert non-ATCG nucleotides to "N"
clean_sequence <- function(dna_seq) {
  dna_seq <- toupper(dna_seq)
  dna_seq[!dna_seq %in% c("A", "T", "C", "G")] <- "N"
  return(dna_seq)
}

# Function to count synonymous and non-synonymous changes and log details
count_syn_nonsyn_changes <- function(dna_seq1, dna_seq2) {
  # Clean sequences to replace non-ATCG values with "N"
  dna_seq1 <- clean_sequence(dna_seq1)
  dna_seq2 <- clean_sequence(dna_seq2)
  
  # Ensure lengths are divisible by 3 for codon-based translation
  dna_seq1 <- dna_seq1[1:(floor(length(dna_seq1) / 3) * 3)]
  dna_seq2 <- dna_seq2[1:(floor(length(dna_seq2) / 3) * 3)]
  
  # Translate DNA sequences to amino acid sequences
  protein_seq1 <- translate(s2c(paste(dna_seq1, collapse = "")), numcode = 1)  # Standard genetic code
  protein_seq2 <- translate(s2c(paste(dna_seq2, collapse = "")), numcode = 1)
  
  # Initialize counters and logs
  non_synonymous_changes <- 0
  synonymous_changes <- 0
  change_log <- list(synonymous = list(), non_synonymous = list())
  
  # Loop through codons
  for (i in seq(1, length(dna_seq1), by = 3)) {
    # Get codons for each sequence
    codon1 <- dna_seq1[i:(i+2)]
    codon2 <- dna_seq2[i:(i+2)]
    
    # Skip codons with any N or gaps
    if (any(codon1 == "N" | codon2 == "N" | codon1 == "-" | codon2 == "-")) {
      next
    }
    
    # Translate codons to amino acids
    aa1 <- translate(s2c(paste(codon1, collapse = "")), numcode = 1)
    aa2 <- translate(s2c(paste(codon2, collapse = "")), numcode = 1)
    
    # Check if there's a change in any nucleotide position
    nucleotide_changes <- codon1 != codon2
    positions_changed <- which(nucleotide_changes)
    
    # Prepare change details
    change_details <- list(
      codon_position = i,
      original_codon = paste(codon1, collapse = ""),
      modified_codon = paste(codon2, collapse = ""),
      positions_changed = positions_changed,
      nucleotides_changed = list()
    )
    
    # Log each nucleotide change within the codon
    for (pos in positions_changed) {
      change_details$nucleotides_changed[[length(change_details$nucleotides_changed) + 1]] <- list(
        position = i + pos - 1,  # Absolute position in the sequence
        original_base = codon1[pos],
        modified_base = codon2[pos]
      )
    }
    
    # Check for non-synonymous (amino acid change)
    if (aa1 != aa2) {
      non_synonymous_changes <- non_synonymous_changes + 1
      change_log$non_synonymous[[length(change_log$non_synonymous) + 1]] <- change_details
    } else {
      # Check for synonymous (nucleotide change without amino acid change)
      if (any(nucleotide_changes)) {
        synonymous_changes <- synonymous_changes + 1
        change_log$synonymous[[length(change_log$synonymous) + 1]] <- change_details
      }
    }
  }
  
  return(list(
    total_non_synonymous = non_synonymous_changes,
    total_synonymous = synonymous_changes,
    change_log = change_log
  ))
}


# Function to create a summary table of synonymous and non-synonymous counts for each gene region
create_gene_region_summary <- function(gene_splits_by_pair_all) {
  # Initialize an empty summary data frame
  gene_summary <- data.frame(
    pair_id = character(),
    region = character(),
    synonymous_count = integer(),
    nonsynonymous_count = integer(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each gene region
  for (i in 1:nrow(gene_regions)) {
    region_name <- gene_regions$region[i]
    print(region_name)
    
    # Loop through all pairs in the gene_splits_by_pair_all list
    for (pair_id in names(gene_splits_by_pair_all)) {
      pair_data <- gene_splits_by_pair_all[[pair_id]]
      
      # Check if the gene region exists for this pair
      if (region_name %in% names(pair_data)) {
        seq1 <- strsplit(pair_data[[region_name]]$seq[[1]], "")[[1]]
        seq2 <- strsplit(pair_data[[region_name]]$seq[[2]], "")[[1]]
        
        # Count the synonymous and non-synonymous changes for this pair
        counts <- count_syn_nonsyn_changes(seq1, seq2)
        
        # Add the counts to the summary table, along with pair_id and region
        gene_summary <- rbind(gene_summary, data.frame(
          pair_id = pair_id,
          region = region_name,
          synonymous_count = counts$total_synonymous,
          nonsynonymous_count = counts$total_non_synonymous,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(gene_summary)
}

append_counts_to_metadata <- function(pairs_all_metadata, gene_splits_by_pair_all) {
  # Create the gene region summary
  gene_summary_table <- create_gene_region_summary(gene_splits_by_pair_all)
  
  # Merge the summary table with pairs_all_metadata based on pair_id
  pairs_all_with_counts <- merge(pairs_all_metadata, gene_summary_table, by = "pair_id", all.x = TRUE)
  
  # Return the updated metadata with the appended counts
  return(pairs_all_with_counts)
}


# Call the functions in the appropriate order as needed in your analysis
comparison_pairs_df <- create_comparison_pairs(all_without_ct_metadata)
pair_alignments <- create_alignment_for_pairs(alignment_all_persistent_infection_fasta_seqs_without_CT, comparison_pairs_df)
gene_splits_by_pair <- split_alignment_by_gene_for_pairs(pair_alignments, gene_regions)
gene_summary_table_first_last <- create_gene_region_summary(gene_splits_by_pair)
updated_metadata_first_last <- append_counts_to_metadata(
  comparison_pairs_df, 
  gene_splits_by_pair
)
print(updated_metadata_first_last)

#Merge with metadata
updated_metadata_first_last <- updated_metadata_first_last %>%
  left_join(metadata_subset, by = c("PERSON_ID", "infection_episode"))

saveRDS(updated_metadata_first_last, file= "")
saveRDS(gene_summary_table_first_last, file= "")





