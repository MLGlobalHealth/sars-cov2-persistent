# Testing and identifying rare SNP thresholds for all variants 

library(ape);library(pracma);library(phytools);library(TreeTools);library(dplyr);library(tidyverse)
library(data.table);library(dbplyr);library(lubridate);library(rlang)
library(foreach);library(doParallel);library(DSTora);library(ROracle)
library(DSTcolectica);library(DSTdb);library(DBI);library(parallel);
library(foreach);library(ggsignif);library(Rcpp);library(ggplot2)
library(data.table);library(stringr);library(seqinr);library(parallel)
library(foreach);library(future);library(future.apply)
library(viridis);library(ggpubr)

#read in metadata
final_merged_new_id_QC_unique_metadata <- read.csv(file="")

# Generate random pairs function
generate_random_pairs <- function(df, n_pairs = 1000) {
  set.seed(42)  # For reproducibility
  
  # Initialize an empty dataframe to store the pairs
  pairs <- data.frame()
  
  while (nrow(pairs) < n_pairs) {
    # Sample two distinct rows
    sampled_df <- df %>%
      sample_n(2, replace = FALSE)
    
    # Check if PERSON_IDs are different and time difference is >= 26 days
    if (sampled_df$PERSON_ID[1] != sampled_df$PERSON_ID[2]) {
      date_diff <- abs(difftime(sampled_df$DATESAMPLING[1], sampled_df$DATESAMPLING[2], units = "days"))
      if (date_diff >= 26) {
        # Append to pairs dataframe
        new_pair <- data.frame(
          PERSON_ID_1 = sampled_df$PERSON_ID[1],
          DATESAMPLING_1 = sampled_df$DATESAMPLING[1],
          strain_1 = sampled_df$strain[1],
          PERSON_ID_2 = sampled_df$PERSON_ID[2],
          DATESAMPLING_2 = sampled_df$DATESAMPLING[2],
          strain_2 = sampled_df$strain[2]
        )
        pairs <- rbind(pairs, new_pair)
      }
    }
  }
  
  return(pairs)
}


# Create separate dataframes for each variant and generate random pairs ---------
omicron_BA1_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.1-like)") %>%
  generate_random_pairs()

omicron_BA2_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.2-like)") %>%
  generate_random_pairs()

omicron_BA5_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.5-like)") %>%
  generate_random_pairs()

delta_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(str_starts(scorpio_call, "Delta")) %>%
  generate_random_pairs()

alpha_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Alpha (B.1.1.7-like)") %>%
  generate_random_pairs()

omicron_BA1_sequences <- omicron_BA1_sequences %>% mutate(variant = "Omicron_BA1")
omicron_BA2_sequences <- omicron_BA2_sequences %>% mutate(variant = "Omicron_BA2")
omicron_BA5_sequences <- omicron_BA5_sequences %>% mutate(variant = "Omicron_BA5")
delta_sequences <- delta_sequences %>% mutate(variant = "Delta")
alpha_sequences <- alpha_sequences %>% mutate(variant = "Alpha")

# Combine the strain pairs into one dataframe
all_sequences <- bind_rows(omicron_BA1_sequences, omicron_BA2_sequences, omicron_BA5_sequences, delta_sequences, alpha_sequences)

# Export the strain values to a CSV for use in Python
write.table(all_sequences %>% select(strain_1, strain_2, variant), 
            "strain_pairs.txt", 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

#Use python code find_strain_sequences.py to get the relevant FASTA sequences




# Importing sequences back into R ----------

# Read the FASTA file using seqinr
omicron_BA1_seqs <- read.fasta(file = "omicron_BA1_sequences.fasta", as.string = TRUE)
omicron_BA2_seqs <- read.fasta(file = "omicron_BA2_sequences.fasta", as.string = TRUE)
omicron_BA5_seqs <- read.fasta(file = "omicron_BA5_sequences.fasta", as.string = TRUE)
delta_seqs <- read.fasta(file = "delta_sequences.fasta", as.string = TRUE)
alpha_seqs <- read.fasta(file = "alpha_sequences.fasta", as.string = TRUE)

# Function to read FASTA and process sequences in chunks
process_sequences <- function(file_path, chunk_size = 1000) {
  # Read sequences
  sequences <- read.fasta(file = file_path, as.string = TRUE)
  
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
    cat("Processed columns:", start_pos, "to", end_pos, "for", basename(file_path), "\n")
  }
  
  df[, sequence := NULL]  # Optionally remove the sequence column
  setcolorder(df, c("strain", setdiff(names(df), "strain")))  # Ensure 'strain' is the first column
  return(df)
}

# Function to add position columns
add_position_columns <- function(df, start_pos, end_pos) {
  for (i in start_pos:end_pos) {
    df[[paste0("Position_", i)]] <- substring(df$sequence, i, i)
  }
  return(df)
}

# List of FASTA files and corresponding output names
fasta_files <- list(
  "omicron_BA1_sequences.fasta" = "omicron_BA1_df",
  "omicron_BA2_sequences.fasta" = "omicron_BA2_df",
  "omicron_BA5_sequences.fasta" = "omicron_BA5_df",
  "delta_sequences.fasta" = "delta_df",
  "alpha_sequences.fasta" = "alpha_df"
)

# Process each FASTA file and store the output in a named list
results <- lapply(names(fasta_files), function(file) {
  df_name <- fasta_files[[file]]
  df <- process_sequences(file)
  assign(df_name, df, envir = .GlobalEnv)
  df
})

# Optionally, you can manually print the first few rows of each data.table for verification
print(omicron_BA1_df)
print(omicron_BA2_df)
print(omicron_BA5_df)
print(delta_df)
print(alpha_df)


# Define the FASTA data frames for each variant
fasta_dfs <- list(
  "Omicron_BA1" = omicron_BA1_df,
  "Omicron_BA2" = omicron_BA2_df,
  "Omicron_BA5" = omicron_BA5_df,
  "Delta" = delta_df,
  "Alpha" = alpha_df
)

saveRDS(fasta_dfs, file="combined_fasta_dataframes_all_variants.RDS")

# Defining arbitrary thresholds to be tested ----------
num_values <- 50
# Generate values between log10(0.001) and log10(0.1)
thresholds <- logspace(log10(0.00001), log10(0.20), num_values)
thresholds_df <- data.frame(threshold = thresholds)
write.csv(thresholds_df, file = "thresholds.csv", row.names = FALSE)

num_values_experiment <- 2
thresholds_experiment <- logspace(log10(0.00001), log10(0.05), num_values_experiment)
thresholds_df_experiment <- data.frame(threshold = thresholds_experiment)
write.csv(thresholds_df_experiment, file = "thresholds_experiment.csv", row.names = FALSE)


# Generating per site proportions for each of the variants ----------

#Alpha
alpha_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Alpha (B.1.1.7-like)")
write.csv(alpha_sequences, file="alpha_sequences.csv")
#Run alpha_sequence_dataframe.py python script
alpha_nucleotide_counts <- read.csv(file="alpha_nucleotide_counts_per_position.csv")
alpha_nucleotide_counts <- alpha_nucleotide_counts[, -which(names(alpha_nucleotide_counts) == "X.")]
alpha_nucleotide_counts$total <- rowSums(alpha_nucleotide_counts[, 2:7])
alpha_nucleotide_proportions <- alpha_nucleotide_counts
alpha_nucleotide_proportions[, 2:7] <- alpha_nucleotide_counts[, 2:7] / alpha_nucleotide_counts$total
alpha_nucleotide_proportions$Position <- alpha_nucleotide_proportions$Position + 1
write.csv(alpha_nucleotide_proportions, file="alpha_nucleotide_proportions.csv")


#Delta
all_delta_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(str_starts(scorpio_call, "Delta"))
write.csv(all_delta_sequences, file="delta_sequences.csv")
#Run delta_sequence_dataframe.py python script
delta_nucleotide_counts <- read.csv(file="delta_nucleotide_counts_per_position.csv")
delta_nucleotide_counts <- delta_nucleotide_counts[, -which(names(delta_nucleotide_counts) == "X.")]
#Calculate proportions
delta_nucleotide_counts$total <- rowSums(delta_nucleotide_counts[, 2:7])
delta_nucleotide_proportions <- delta_nucleotide_counts
delta_nucleotide_proportions[, 2:7] <- delta_nucleotide_counts[, 2:7] / delta_nucleotide_counts$total
delta_nucleotide_proportions$Position <- delta_nucleotide_proportions$Position + 1
write.csv(delta_nucleotide_proportions, file="delta_nucleotide_proportions.csv")

#Omicron BA.1
omicron_BA1_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.1-like)")
write.csv(omicron_BA1_sequences, file="omicron_BA1_sequences.csv")
# Run omicronBA1_sequence_dataframe.py python script
omicronBA1_nucleotide_counts <- read.csv(file="omicron_BA1_nucleotide_counts_per_position.csv")
omicronBA1_nucleotide_counts <- omicronBA1_nucleotide_counts[, -which(names(omicronBA1_nucleotide_counts) == "X.")]
#Calculate proportions
omicronBA1_nucleotide_counts$total <- rowSums(omicronBA1_nucleotide_counts[, 2:7])
omicronBA1_nucleotide_proportions <- omicronBA1_nucleotide_counts
omicronBA1_nucleotide_proportions[, 2:7] <- omicronBA1_nucleotide_counts[, 2:7] / omicronBA1_nucleotide_counts$total
omicronBA1_nucleotide_proportions$Position <- omicronBA1_nucleotide_proportions$Position + 1
write.csv(omicronBA1_nucleotide_proportions, file="omicronBA1_nucleotide_proportions.csv")

#Omicron BA.2
omicron_BA2_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.2-like)")
write.csv(omicron_BA2_sequences, file="omicron_BA2_sequences.csv")
#Run omicronBA2_sequence_dataframe.py python script
omicronBA2_nucleotide_counts <- read.csv(file="omicron_BA2_nucleotide_counts_per_position.csv")
omicronBA2_nucleotide_counts <- omicronBA2_nucleotide_counts[, -which(names(omicronBA2_nucleotide_counts) == "X.")]
#Calculate proportions
omicronBA2_nucleotide_counts$total <- rowSums(omicronBA2_nucleotide_counts[, 2:7])
omicronBA2_nucleotide_proportions <- omicronBA2_nucleotide_counts
omicronBA2_nucleotide_proportions[, 2:7] <- omicronBA2_nucleotide_counts[, 2:7] / omicronBA2_nucleotide_counts$total
omicronBA2_nucleotide_proportions$Position <- omicronBA2_nucleotide_proportions$Position + 1
write.csv(omicronBA2_nucleotide_proportions, file="omicronBA2_nucleotide_proportions.csv")
omicronBA2_nucleotide_proportions <- read.csv(file="omicronBA2_nucleotide_proportions.csv")
write.csv(omicronBA2_nucleotide_proportions, file="omicronBA2_nucleotide_proportions.csv")


#Omicron BA.5
omicron_BA5_sequences <- final_merged_new_id_QC_unique_metadata %>%
  filter(scorpio_call == "Omicron (BA.5-like)")
write.csv(omicron_BA5_sequences, file="omicron_BA5_sequences.csv")
#Run omicronBA5_sequence_dataframe.py python script
omicronBA5_nucleotide_counts <- read.csv(file="omicron_BA5_nucleotide_counts_per_position.csv")
omicronBA5_nucleotide_counts <- omicronBA5_nucleotide_counts[, -which(names(omicronBA5_nucleotide_counts) == "X.")]
#Calculate proportions
omicronBA5_nucleotide_counts$total <- rowSums(omicronBA5_nucleotide_counts[, 2:7])
omicronBA5_nucleotide_proportions <- omicronBA5_nucleotide_counts
omicronBA5_nucleotide_proportions[, 2:7] <- omicronBA5_nucleotide_counts[, 2:7] / omicronBA5_nucleotide_counts$total
omicronBA5_nucleotide_proportions$Position <- omicronBA5_nucleotide_proportions$Position + 1
write.csv(omicronBA5_nucleotide_proportions, file="omicronBA5_nucleotide_proportions.csv")






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
omicronBA1_proportions <- normalize_and_select_columns(omicronBA1_nucleotide_proportions)

omicronBA2_proportions <- normalize_and_select_columns(omicronBA2_nucleotide_proportions)

omicronBA5_proportions <- normalize_and_select_columns(omicronBA5_nucleotide_proportions)

alpha_proportions <- normalize_and_select_columns(alpha_nucleotide_nucleotide_proportions)

delta_proportions <- normalize_and_select_columns(delta_nucleotide_nucleotide_proportions)



#Finding false-positive rates for different thresholds -----------

# Define the list of thresholds
thresholds <- fread("thresholds.csv")$threshold
thresholds_experiment <- fread("thresholds_experiment.csv")$threshold

# Define the nucleotide proportions data
scorpio_thresholds <- list(
  "Omicron_BA1" = omicronBA1_proportions,
  "Omicron_BA2" = omicronBA2_proportions,
  "Omicron_BA5" = omicronBA5_proportions,
  "Delta" = delta_nucleotide_proportions,
  "Alpha" = alpha_nucleotide_proportions
)

# Convert each object to data.table if they are not already
all_sequences <- as.data.table(all_sequences)
fasta_dfs <- lapply(fasta_dfs, as.data.table)
scorpio_thresholds <- lapply(scorpio_thresholds, as.data.table)
thresholds <- as.data.table(thresholds)
thresholds_experiment <- as.data.table(thresholds_experiment)




# Running optimized version ---------

# Define the function to process a single variant with parallel processing
process_variant_parallel <- function(variant_name, all_sequences, fasta_dfs, scorpio_thresholds, thresholds) {
  cat("Processing variant:", variant_name, "\n")
  plan(multisession, workers = 16)
  
  # Filter pairs for the current variant
  variant_pairs <- all_sequences[variant == variant_name, ]
  
  # Get the FASTA data and nucleotide proportions for the current variant
  fasta_df <- fasta_dfs[[variant_name]]
  nucleotide_proportions <- scorpio_thresholds[[variant_name]]
  
  # Initialize a result data.table for this variant with fixed columns
  variant_result <- data.table(PERSON_ID_1 = character(), PERSON_ID_2 = character(), threshold = numeric(), has_rare_snp = logical())
  
  # Define a function to process each pair
  process_pair <- function(pair, index) {
    
    cat("Processing pair", index, "\n")
    
    strain_1 <- pair$strain_1
    strain_2 <- pair$strain_2

    
    # Filter the FASTA data for the current pair
    pair_fasta_df <- fasta_df[strain %in% c(strain_1, strain_2), ]
    
    # Initialize a list to store whether a rare SNP was found for each threshold
    rare_snp_found <- rep(FALSE, length(thresholds))
    
    # Iterate over each position in the sequences
    for (position in names(pair_fasta_df)[-1][grepl("Position_", names(pair_fasta_df)[-1])]) {
      position_number <- as.numeric(sub("Position_", "", position))
      
      nucleotides <- pair_fasta_df[[position]]
      
      # Check if both sequences have the same nucleotide (A, T, C, G)
      if (length(unique(nucleotides)) == 1 && nucleotides[1] %in% c("A", "T", "C", "G")) {
        nucleotide <- nucleotides[1]
        nucleotide_lower <- tolower(nucleotide)
        
        # Get the corresponding proportion from nucleotide_proportions
        proportion_row <- nucleotide_proportions[Position == position_number]
        
        if (!is.na(proportion_row[,2]) && nrow(proportion_row) == 1) {
          proportion_column <- match(nucleotide_lower, names(proportion_row))
          proportion <- proportion_row[[proportion_column]]
          
          # Loop through each threshold
          for (j in seq_along(thresholds)) {
            threshold <- thresholds[j]
            
            # Check if the proportion is less than the threshold
            if (!is.na(proportion) && proportion < threshold) {
              rare_snp_found[j] <- TRUE
            }
          }
          
          # If any rare SNP is found, break out of the position loop for efficiency
          if (any(rare_snp_found)) {
            break
          }
        }
      }
    }
    
    # Add results for this pair to the variant_result data.table
    results <- data.table(
      PERSON_ID_1 = strain_1,
      PERSON_ID_2 = strain_2,
      threshold = thresholds,
      has_rare_snp = rare_snp_found
    )
    
    return(results)
  }
  
  # Apply processing function in parallel
  result_list <- future.apply::future_lapply(seq_len(nrow(variant_pairs)), function(i) {
    process_pair(variant_pairs[i, ], i)
  })
  
  # Combine results into a single data.table
  variant_result <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
  
  # Return the result for this variant
  return(variant_result)
}

# Assuming `all_sequences`, `fasta_dfs`, `scorpio_thresholds`, and `thresholds_experiment` are already defined
result_Omicron_BA1 <- process_variant_parallel("Omicron_BA1", all_sequences, fasta_dfs, scorpio_thresholds, thresholds)

result_Omicron_BA2 <- process_variant_parallel("Omicron_BA2", all_sequences, fasta_dfs, scorpio_thresholds, thresholds)

result_Omicron_BA5 <- process_variant_parallel("Omicron_BA5", all_sequences, fasta_dfs, scorpio_thresholds, thresholds)

result_Delta <- process_variant_parallel("Delta", all_sequences, fasta_dfs, scorpio_thresholds, thresholds)

result_Alpha <- process_variant_parallel("Alpha", all_sequences, fasta_dfs, scorpio_thresholds, thresholds)


#Summarizing false-positive rates ------

#Omicron BA1
summary_Omicron_BA1 <- result_Omicron_BA1 %>%
  group_by(threshold) %>%
  summarise(
    count_TRUE = sum(has_rare_snp == TRUE),
    total = n(),
    percentage_TRUE = (count_TRUE / total) * 100
  )

#Omicron BA2
summary_Omicron_BA2 <- result_Omicron_BA2 %>%
  group_by(threshold) %>%
  summarise(
    count_TRUE = sum(has_rare_snp == TRUE),
    total = n(),
    percentage_TRUE = (count_TRUE / total) * 100
  )

#Omicron BA5
summary_Omicron_BA5 <- result_Omicron_BA5 %>%
  group_by(threshold) %>%
  summarise(
    count_TRUE = sum(has_rare_snp == TRUE),
    total = n(),
    percentage_TRUE = (count_TRUE / total) * 100
  )

#Delta
summary_Delta <- result_Delta %>%
  group_by(threshold) %>%
  summarise(
    count_TRUE = sum(has_rare_snp == TRUE),
    total = n(),
    percentage_TRUE = (count_TRUE / total) * 100
  )

#Alpha
summary_Alpha <- result_Alpha %>%
  group_by(threshold) %>%
  summarise(
    count_TRUE = sum(has_rare_snp == TRUE),
    total = n(),
    percentage_TRUE = (count_TRUE / total) * 100
  )


# Based on these: 1.0907e-01 is the threshold we can set
# Specific variant false-positive rates: Omicron BA.1: 3.3, Omicron BA.2: 2.1, Omicron BA.5: 2.2, Delta: 3.0, Alpha: 3.1

# Plotting the false-positive rates vs. the thresholds
vline_x <- log10(1.090691e-01)

# Function to plot the data
plot_percentage_vs_log_threshold <- function(data, title) {
  ggplot(data, aes(x = log10(threshold), y = percentage_TRUE)) +
    geom_line(color = viridis(1)) +
    geom_point(color = viridis(1)) +
    geom_vline(xintercept = vline_x, linetype = "dashed", color = "black") +
    labs(
      title = title,
      x = "Log10(Threshold)",
      y = "False-Positive Rate, %"
    ) +
    theme_minimal()
}

# Plot for Omicron BA1
plot_Omicron_BA1 <- plot_percentage_vs_log_threshold(summary_Omicron_BA1, "Omicron BA.1")

# Plot for Omicron BA2
plot_Omicron_BA2 <- plot_percentage_vs_log_threshold(summary_Omicron_BA2, "Omicron BA.2")

# Plot for Omicron BA5
plot_Omicron_BA5 <- plot_percentage_vs_log_threshold(summary_Omicron_BA5, "Omicron BA.5")

# Plot for Delta
plot_Delta <- plot_percentage_vs_log_threshold(summary_Delta, "Delta")

# Plot for Alpha
plot_Alpha <- plot_percentage_vs_log_threshold(summary_Alpha, "Alpha")

# Display the plots
print(plot_Omicron_BA1)
print(plot_Omicron_BA2)
print(plot_Omicron_BA5)
print(plot_Delta)
print(plot_Alpha)

# Combine the plots using ggarrange
combined_threshold_plot <- ggarrange(
  plot_Omicron_BA1, plot_Omicron_BA2, plot_Omicron_BA5, plot_Delta, plot_Alpha,
  ncol = 2, nrow = 3
)

# Export the combined plot to a PDF file in portrait mode
ggsave(
  filename = "",     # Output file name
  plot = combined_threshold_plot,               # Plot object to save
  device = "pdf",                               # File format
  width = 8.5,                            # Set the width of the PDF
  height = 11,                          # Set the height of the PDF
  units = "in",                                 # Units for width and height
  dpi = 300                                     # Resolution for high quality
)






