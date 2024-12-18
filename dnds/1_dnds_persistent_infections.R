
# Synonymous vs. non-synonymous analysis for persistent infections

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(ggpubr)


all_without_ct_changes_first_last_placements <- readRDS("")
dnds_summary_dataframe_all_without_ct <- readRDS("")
first_last <- dnds_summary_dataframe_all_without_ct %>%
  group_by(pair_id) %>%
  slice(1) %>%
  ungroup() 


# Importing fasta data from each folder ---------------------------------------------
alignment_all_persistent_infection_fasta_seqs_without_CT <- read.alignment(file = "", format = "fasta", forceToLower = TRUE)


# Sequential gene-specific synonymous vs nonsynonymous analysis ------------
gene_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
)



# Function to create an alignment for each unique pair
create_alignment_for_pairs <- function(alignment, metadata_df) {
  # Extract sequence names from the alignment
  sequence_names <- alignment$nam
  
  # Initialize a list to store alignments for each pair
  pair_alignments <- list()
  
  # Loop through each row in the metadata to get the strain pair (first_strain, last_strain)
  for (i in 1:nrow(metadata_df)) {
    seq1 <- metadata_df$first_strain[i]
    seq2 <- metadata_df$last_strain[i]
    
    # Skip if there's no second strain
    if (is.na(seq2)) next
    
    # Find indices of these sequences in the alignment
    indices_to_include <- which(sequence_names %in% c(seq1, seq2))
    
    # If both sequences are not found, skip this pair
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
    
    # Store the alignment in the list using a unique pair identifier
    pair_id <- paste(seq1, seq2, sep = "_")
    pair_alignments[[pair_id]] <- subset_alignment
  }
  
  return(pair_alignments)
}


# Function to split alignment into gene regions based on gene_regions dataframe
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
    
    # Debugging: print which gene is being processed
    print(paste("Processing gene region:", gene_name))
    
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

# Run Ka/Ks analysis on each gene alignment and capture results
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
class(alignment_all_persistent_infection_fasta_seqs_without_CT) <- "alignment"

# Create alignments for each pair based on metadata (first_last)
pair_alignments_all_without_ct <- create_alignment_for_pairs(alignment_all_persistent_infection_fasta_seqs_without_CT, first_last)

# Split each pair's alignment into gene regions
gene_alignments_all_without_ct <- lapply(pair_alignments_all_without_ct, function(example_alignment) {
  class(example_alignment) <- "alignment"  # Ensure the class is set to 'alignment'
  split_alignment_by_gene(example_alignment, gene_regions)
})

# Run Ka/Ks analysis for each gene region
kaks_results_all_without_ct <- lapply(gene_alignments_all_without_ct, run_kaks_on_gene_alignments)

# Extract and clean dN/dS counts
dnds_counts_all_without_ct <- extract_clean_dnds_counts_to_df(kaks_results_all_without_ct)
saveRDS(dnds_counts_all_without_ct, file="dnds_counts_all_without_ct.RDS")


