#Figure 4, Plotting Fitness of Mutations

library('ape');library('adegenet');library('TreeTools')
library('tidyverse');library(phangorn);library(patchwork);library(ggplot2);library(readr)
library(phytools);library(openxlsx);library(readxl);library(lubridate);library(dplyr)
library(data.table);library(magrittr);library(tidyr);library(parallel);library(reshape2);library(pegas)
library(stats);library(zoo);library(ggplot2);library(dplyr);library(brms);library(ggpubr)

#Read in data
file_path <- ""
node_data <- read.delim(file_path, sep = "\t")
node_data$nt_mutations <- as.character(node_data$nt_mutations)
node_data$aa_mutations <- as.character(node_data$aa_mutations)

file_path_new <- ""
node_data_new <- read.delim(file_path_new, sep = "\t")
node_data_new$nt_mutations <- as.character(node_data_new$nt_mutations)
node_data_new$aa_mutations <- as.character(node_data_new$aa_mutations)

combined_data <- rbind(node_data, node_data_new)

#Getting mutations and amino acid changes
mutation_dict <- data.frame(
  nt_mutation = character(),
  aa_mutation = character(),
  stringsAsFactors = FALSE
)

# Iterate over each row in the original data frame
for (i in 1:nrow(combined_data)) {
  # Split the nt_mutations and aa_mutations fields by ";"
  nt_mut_list <- unlist(strsplit(combined_data$nt_mutations[i], ";"))
  aa_mut_list <- unlist(strsplit(combined_data$aa_mutations[i], ";"))
  
  # Ensure both lists have the same length to avoid misalignment
  if (length(nt_mut_list) == length(aa_mut_list)) {
    # Create a temporary data frame for the current row's mappings
    temp_df <- data.frame(
      nt_mutation = nt_mut_list,
      aa_mutation = aa_mut_list,
      stringsAsFactors = FALSE
    )
    # Append to the main mutation_dict data frame
    mutation_dict <- rbind(mutation_dict, temp_df)
  } else {
    warning(paste("Row", i, "has mismatched counts of nt_mutations and aa_mutations"))
  }
}
head(mutation_dict)


# Reading in recurrent mutations ---------
recurrent_file_path <- ""
recurrent_mutation_data <- read.delim(recurrent_file_path, sep = ",")

filtered_mutation_dict <- mutation_dict %>%
  group_by(nt_mutation) %>%
  filter(grepl("ORF1ab", aa_mutation) | !any(grepl("ORF1ab", aa_mutation))) %>%
  slice(1) %>%  # keep only the first occurrence if there are multiple
  ungroup()

# Next, join this with recurrent_mutation_data based on mutation_type and nt_mutation
recurrent_mutation_data <- recurrent_mutation_data %>%
  left_join(filtered_mutation_dict, by = c("mutation_type" = "nt_mutation"))

# Display the resulting data frame with the added amino acid mutations
head(recurrent_mutation_data)

recurrent_mutation_data_clean <- recurrent_mutation_data[!is.na(recurrent_mutation_data$aa_mutation) & recurrent_mutation_data$aa_mutation != "", ]
# Function to extract aa_site and aa
extract_aa_info <- function(aa_mutation) {
  # Extract the last character as the amino acid
  aa <- substr(aa_mutation, nchar(aa_mutation), nchar(aa_mutation))
  
  # Remove everything before the colon, the colon, and the first letter after the colon
  site_str <- sub(".*:[A-Za-z](\\d+)[A-Za-z]$", "\\1", aa_mutation)
  
  # Convert site to integer
  site <- as.integer(site_str)
  
  # Return the site and the amino acid
  return(list(aa_site = site, aa = aa))
}

# Apply the function to extract aa_site and aa for the clean data
aa_info <- apply(recurrent_mutation_data_clean, 1, function(x) extract_aa_info(x["aa_mutation"]))

# Add the extracted values back to the dataframe
recurrent_mutation_data_clean$aa_site <- sapply(aa_info, function(x) x$aa_site)
recurrent_mutation_data_clean$aa <- sapply(aa_info, function(x) x$aa)

# Check the result
head(recurrent_mutation_data_clean)


# Set working directory and read the fitness data
fitness_file_path <- "aa_fitness.csv"
fitness_mutation_data <- read.delim(fitness_file_path, sep = ",")

# Filter fitness data to remove genes that start with 'nsp'
fitness_mutation_data <- fitness_mutation_data[!grepl("^nsp", fitness_mutation_data$gene), ]

# Aggregate fitness data by gene, aa_site, and aa
fitness_mutation_data_agg <- fitness_mutation_data %>% 
  group_by(gene, aa_site, aa) %>% 
  summarise(
    fitness = mean(fitness, na.rm = TRUE),  # Average fitness values
    expected_count = mean(expected_count, na.rm = TRUE),  # Average expected counts
    .groups = 'drop'  # Drop the grouping after summarizing
  )

# Function to extract gene, amino acid (aa), and site information from aa_mutation
extract_aa_info <- function(aa_mutation) {
  # Extract gene (everything before the colon)
  gene <- sub(":.*", "", aa_mutation)
  
  # Extract the amino acid (last character)
  aa <- substr(aa_mutation, nchar(aa_mutation), nchar(aa_mutation))
  
  # Extract the site number (remove gene and amino acid, keep the numeric value)
  site_str <- sub(".*:[A-Za-z](\\d+)[A-Za-z]$", "\\1", aa_mutation)
  site <- as.integer(site_str)
  
  # Return a list containing gene, aa_site, and aa
  return(list(gene = gene, aa_site = site, aa = aa))
}

# Apply the function to extract aa_site, aa, and gene for the mutations
recurrent_mutation_data_clean$aa_info <- lapply(recurrent_mutation_data_clean$aa_mutation, extract_aa_info)

# Add the extracted information to the dataframe
recurrent_mutation_data_clean$gene <- sapply(recurrent_mutation_data_clean$aa_info, function(x) x$gene)
recurrent_mutation_data_clean$aa_site <- sapply(recurrent_mutation_data_clean$aa_info, function(x) x$aa_site)
recurrent_mutation_data_clean$aa <- sapply(recurrent_mutation_data_clean$aa_info, function(x) x$aa)

# Split the data into recurrent and non-recurrent mutations based on the count
recurrent_mutation_data_recurrent_clean <- recurrent_mutation_data_clean %>%
  filter(count >= 2)  # Recurrent mutations (count >= 2)

recurrent_mutation_data_non_recurrent_clean <- recurrent_mutation_data_clean %>%
  filter(count < 2)  # Non-recurrent mutations (count < 2)

# Remove rows with NA or empty aa_mutation for both recurrent and non-recurrent mutations
recurrent_mutation_data_recurrent_clean <- recurrent_mutation_data_recurrent_clean %>%
  filter(!is.na(aa_mutation) & aa_mutation != "")

recurrent_mutation_data_non_recurrent_clean <- recurrent_mutation_data_non_recurrent_clean %>%
  filter(!is.na(aa_mutation) & aa_mutation != "")

# Remove duplicates based on aa_site, aa, and gene
recurrent_mutation_data_recurrent_unique <- recurrent_mutation_data_recurrent_clean %>%
  distinct(gene, aa_site, aa, .keep_all = TRUE)

recurrent_mutation_data_non_recurrent_unique <- recurrent_mutation_data_non_recurrent_clean %>%
  distinct(gene, aa_site, aa, .keep_all = TRUE)

# Merge fitness values based on gene, aa_site, and aa for both recurrent and non-recurrent mutations
recurrent_mutation_data_recurrent_updated <- recurrent_mutation_data_recurrent_unique %>%
  left_join(fitness_mutation_data_agg, by = c("gene", "aa_site", "aa"))

recurrent_mutation_data_non_recurrent_updated <- recurrent_mutation_data_non_recurrent_unique %>%
  left_join(fitness_mutation_data_agg, by = c("gene", "aa_site", "aa"))

# View the final updated data for recurrent mutations
head(recurrent_mutation_data_recurrent_updated)

# View the final updated data for non-recurrent mutations
head(recurrent_mutation_data_non_recurrent_updated)


#Getting the other corresponding sites and amino acids (i.e., the inverse of the recurrents)
unique_recurrent <- unique(recurrent_mutation_data_non_recurrent_updated[, c("gene", "aa_site", "aa")])
filtered_complementary_aa_data <- fitness_mutation_data[
  fitness_mutation_data$gene %in% unique_recurrent$gene &
    fitness_mutation_data$aa_site %in% unique_recurrent$aa_site &
    !(paste(fitness_mutation_data$gene, fitness_mutation_data$aa_site, fitness_mutation_data$aa) %in% 
        paste(unique_recurrent$gene, unique_recurrent$aa_site, unique_recurrent$aa)),
]
head(filtered_complementary_aa_data)





# Plotting fitness --------

#Plot 1, showing fitness for single recurrents vs multiple recurrents
recurrent_mutation_data_recurrent_updated$mutation_type <- "Recurrent mutations"
recurrent_mutation_data_non_recurrent_updated$mutation_type <- "Single mutations"

# ggplot
fitness_density_histogram <- ggplot() +
  # Histogram for the ">1" mutation type
  geom_histogram(
    data = recurrent_mutation_data_recurrent_updated,
    aes(x = fitness, y = ..count.. / sum(..count..), fill = "Recurrent mutations"),
    position = "identity",
    alpha = 0.6,
    bins = 30,
    color = "black"
  ) +
  # Histogram for the "1" mutation type
  geom_histogram(
    data = recurrent_mutation_data_non_recurrent_updated,
    aes(x = fitness, y = ..count.. / sum(..count..), fill = "Single mutations"),
    position = "identity",
    alpha = 0.6,
    bins = 30,
    color = "black"
  ) +
  # Customizing the color palette
  scale_fill_manual(values = c("Recurrent mutations" = "#009E73", "Single mutations" = "#E69F00")) +
  # Vertical line at fitness = 0
  geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
  labs(
    title = "",
    x = "Fitness",
    y = "Proportion",
    fill = "# Appearances"
  ) +
  theme_minimal(base_size = 14) +
  # Customizing the legend
  theme(
    legend.position = "top",                 # Position legend at the top
    legend.title = element_blank(),  # Smaller title
    legend.text = element_text(size = 9),    # Smaller text
    legend.key.size = unit(0.5, "cm")        # Adjust the size of the legend keys
  )

# Display the proportion histogram
print(fitness_density_histogram)


# Plot 2 - comparing recurrent vs all other potnetial amino acids from the same location
filtered_complementary_aa_data$mutation_type <- "Matched sites, unrealized mutation"
recurrent_mutation_data_recurrent_updated$mutation_type <- "Recurrent mutations"

# ggplot
fitness_recurrent_vs_others_histogram <- ggplot() +
  # Histogram for the "Recurrent" mutation type
  geom_histogram(
    data = recurrent_mutation_data_recurrent_updated,
    aes(x = fitness, y = ..count.. / sum(..count..), fill = "Recurrent mutations"),
    position = "identity",
    alpha = 0.6,
    bins = 30,
    color = "black"
  ) +
  # Histogram for the "Matched sites, unrealized mutations" mutation type
  geom_histogram(
    data = filtered_complementary_aa_data,
    aes(x = fitness, y = ..count.. / sum(..count..), fill = "Matched sites, unrealized mutation"),
    position = "identity",
    alpha = 0.6,
    bins = 30,
    color = "black"
  ) +
  # Customizing the color palette
  scale_fill_manual(values = c("Recurrent mutations" = "#009E73", "Matched sites, unrealized mutation" = "#888888")) +
  # Vertical line at fitness = 0
  geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
  labs(
    title = "",
    x = "Fitness",
    y = "Proportion",
    fill = "Mutation Type"
  ) +
  theme_minimal(base_size = 14) +
  # Customizing the legend
  theme(
    legend.position = "top",                 # Position legend at the top
    legend.title = element_blank(),  # Smaller title
    legend.text = element_text(size = 9),    # Smaller text
    legend.key.size = unit(0.5, "cm")        # Adjust the size of the legend keys
  )

print(fitness_recurrent_vs_others_histogram)

# Combine into one plot ------
combined_fitness_plot <- ggarrange(
  fitness_density_histogram, fitness_recurrent_vs_others_histogram,
  labels = c("a)", "b)"),
  ncol = 2, nrow = 1,
  common.legend = FALSE,
  legend = "bottom"
)

ggsave("", 
       plot = combined_fitness_plot, 
       device = "pdf", 
       width = 8,      
       height = 4,     
       dpi = 300)




# Statistical tests ------

#Single mutation vs recurrents
fitness_data_combined <- rbind(recurrent_mutation_data_recurrent_updated, recurrent_mutation_data_non_recurrent_updated)
fitness_data_combined$mutation_type <- as.factor(fitness_data_combined$mutation_type)

brm_model <- brm(
  formula = fitness ~ mutation_type,
  data = fitness_data_combined,
  family = gaussian(),  
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 1), class = "b")           
  ),
  iter = 2000,  
  chains = 4,   
  cores = 4,    
  seed = 123    
)
summary(brm_model)
# Extract the posterior samples for the mutation_type effect
posterior_summary <- posterior_summary(brm_model, variable = "b_mutation_typeRecurrent")
print(posterior_summary)
print(cr_interval)# Calculate the 95% CrI for the effect of "Recurrent" mutation_type
cr_interval <- quantile(posterior_samples$b_mutation_typeRecurrent, probs = c(0.025, 0.975))
print(cr_interval)
hypothesis_test <- hypothesis(brm_model, "mutation_typeRecurrent > 0")
print(hypothesis_test)


#Recurrent vs others
recurrent_subset <- recurrent_mutation_data_recurrent_updated[, c("fitness", "mutation_type")]
unrealized_subset <- filtered_complementary_aa_data[, c("fitness", "mutation_type")]
recurrent_subset$mutation_type <- "Recurrent mutations"
unrealized_subset$mutation_type <- "Matched sites, unrealized mutation"
comparison_data <- rbind(recurrent_subset, unrealized_subset)
comparison_data$mutation_type <- as.factor(comparison_data$mutation_type)
comparison_data$mutation_type <- relevel(comparison_data$mutation_type, ref = "Recurrent mutations")


# Fit the Bayesian model with the updated reference level
brm_model_comparison_recurrent_others <- brm(
  formula = fitness ~ mutation_type,  # Linear model with mutation_type as predictor
  data = comparison_data,
  family = gaussian(),  # Gaussian family for continuous fitness
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 1), class = "b")           
  ),
  iter = 2000,  
  chains = 4,   
  cores = 4,    
  seed = 123    
)

# Summarize the model
summary(brm_model_comparison_recurrent_others)
posterior_summary <- posterior_summary(brm_model_comparison_recurrent_others, variable = "b_mutation_typeMatchedsitesunrealizedmutation")
print(posterior_summary)
posterior_samples <- as_draws_df(brm_model_comparison_recurrent_others)
cr_interval <- quantile(posterior_samples$b_mutation_typeMatchedsitesunrealizedmutation, probs = c(0.025, 0.975))
print(cr_interval)


