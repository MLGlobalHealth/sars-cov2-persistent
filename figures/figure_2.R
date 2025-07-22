# Figure 2: Visualizing the locations of mutations and changes

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(viridis); library(patchwork); library(openxlsx); library(dplyr); library(gridExtra); library(ggpubr)


# Reading in data -----------
all_without_ct_changes_df <- readRDS("")
updated_metadata_with_dnds_counts_and_totals_all_without_ct <- readRDS("")

# ALL: Summarizing placements of mutations ---------------
full_coding_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
)

# Function to map mutation positions to coding regions and summarize by major_scorpio_call
map_mutations_to_regions <- function(mutations_df, coding_regions_df) {
  
  # Convert position data from string to numeric vectors
  mutations_df$change_positions <- lapply(strsplit(mutations_df$change_positions, ", "), as.numeric)
  
  # Initialize an empty list to store the results
  region_summary <- list()
  
  # Loop over each unique major_scorpio_call
  for (call in unique(mutations_df$major_scorpio_call)) {
    # Filter the mutations for the current major_scorpio_call
    mutations_for_call <- mutations_df[mutations_df$major_scorpio_call == call, ]
    
    # Loop over each coding region
    for (i in seq_len(nrow(coding_regions_df))) {
      region <- coding_regions_df$region[i]
      start <- coding_regions_df$start[i]
      end <- coding_regions_df$end[i]
      
      # Initialize a count for this region and call
      region_count <- 0
      
      # Loop through each mutation's positions
      for (positions in mutations_for_call$change_positions) {
        # Count how many positions fall within the current region
        region_count <- region_count + sum(positions >= start & positions <= end)
      }
      
      # Store the result for this region and call
      scorpio_summary <- data.frame(
        region = region,
        major_scorpio_call = call,
        mutation_count = region_count
      )
      
      # Append the result to the list
      region_summary[[length(region_summary) + 1]] <- scorpio_summary
    }
  }
  
  # Combine all the region summaries into a single data frame
  final_summary <- bind_rows(region_summary)
  
  # Ensure that all regions and major_scorpio_call combinations are represented
  all_combinations <- expand.grid(
    region = coding_regions_df$region,
    major_scorpio_call = unique(mutations_df$major_scorpio_call),
    stringsAsFactors = FALSE
  )
  
  # Merge the summary with all combinations to ensure no missing combinations
  final_summary <- full_join(all_combinations, final_summary, by = c("region", "major_scorpio_call"))
  
  # Replace NA mutation counts with 0
  final_summary$mutation_count[is.na(final_summary$mutation_count)] <- 0
  
  # Set factor levels to ensure correct order
  final_summary$region <- factor(final_summary$region, levels = coding_regions_df$region)
  
  return(final_summary)
}

all_without_ct_summary <- map_mutations_to_regions(all_without_ct_changes_df, full_coding_regions)

# Updated plotting function to stack bars and respect the region order
plot_histogram_by_region <- function(region_summary) {
  ggplot(region_summary, aes(x = region, y = mutation_count, fill = major_scorpio_call)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis_d(option = "D") +  # Applying Viridis colors
    theme_classic() +
    labs(
      title = "",
      x = "Coding Region",
      y = "# Changes",
      fill = "Major Variant"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.8, 1),  # Inset legend at the top right
          legend.justification = c("right", "top"))
}

mutation_all_without_plot <- plot_histogram_by_region(all_without_ct_summary)
mutation_all_without_plot


# Doing version per site -----
# Normalize mutation counts by coding region length
normalize_mutation_counts <- function(region_summary, coding_regions_df) {
  # Add lengths of coding regions
  coding_regions_df <- coding_regions_df %>%
    mutate(region_length = end - start + 1)
  
  # Join region summary with coding regions length
  normalized_summary <- region_summary %>%
    left_join(coding_regions_df, by = "region") %>%
    mutate(mutations_per_site = mutation_count / region_length)  # Calculate mutation rate
  
  return(normalized_summary)
}
normalized_summary <- normalize_mutation_counts(all_without_ct_summary, full_coding_regions)

normalized_summary$region <- factor(
  normalized_summary$region,
  levels = full_coding_regions$region
)

plot_grouped_bars_by_region <- function(normalized_summary) {
  ggplot(normalized_summary, aes(x = region, y = mutations_per_site, fill = major_scorpio_call)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_viridis_d(option = "D") +  # Applying Viridis colors
    theme_classic() +
    labs(
      title = "",
      x = "Coding Region",
      y = "# Changes per Site",
      fill = "Major Variant"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

# Generate the updated plot
grouped_mutation_plot <- plot_grouped_bars_by_region(normalized_summary)
grouped_mutation_plot


# Collapsed for each gene ----
# Collapse all variants into one by summarizing mutations per site for each gene
collapsed_summary <- normalized_summary %>%
  group_by(region) %>%
  summarize(mutations_per_site = sum(mutations_per_site)) %>%
  ungroup()

# Set the order of regions to match the original gene dataframe
collapsed_summary$region <- factor(
  collapsed_summary$region,
  levels = full_coding_regions$region
)

# Plot collapsed mutation rates
plot_collapsed_bars_by_region <- function(collapsed_summary) {
  ggplot(collapsed_summary, aes(x = region, y = mutations_per_site, fill = region)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_viridis_d(option = "D") +  # Applying Viridis colors
    theme_classic() +
    labs(
      title = "",
      x = "Coding Region",
      y = "# Changes per Site"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Generate the collapsed plot
collapsed_mutation_plot <- plot_collapsed_bars_by_region(collapsed_summary)
collapsed_mutation_plot








# Synonymous vs Nonsynonymous plot for all ----------------

prepare_summary_data <- function(df) {
  # Check if necessary columns exist
  if (!"region" %in% colnames(df)) {
    stop("Column 'region' not found in the dataframe.")
  }
  
  if (!"synonymous_count" %in% colnames(df)) {
    stop("Column 'synonymous_count' not found in the dataframe.")
  }
  
  if (!"nonsynonymous_count" %in% colnames(df)) {
    stop("Column 'nonsynonymous_count' not found in the dataframe.")
  }
  
  # Summarize the data by region, aggregating synonymous and nonsynonymous counts
  summary_data <- df %>%
    group_by(region) %>%
    summarise(
      Synonymous = sum(synonymous_count, na.rm = TRUE),
      Non_Synonymous = sum(nonsynonymous_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_data)
}

# Function to prepare plot data
prepare_plot_data <- function(summary_data) {
  # Reshape the data into a long format for plotting
  plot_data <- summary_data %>%
    pivot_longer(
      cols = c("Synonymous", "Non_Synonymous"),
      names_to = "mutation_type",
      values_to = "count"
    )
  
  return(plot_data)
}

# Full coding region definitions
full_coding_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
)

# Function to prepare histogram data with ordered regions
prepare_histogram_data <- function(plot_data) {
  # Ensure the regions are ordered according to the specified levels
  plot_data <- plot_data %>%
    mutate(region = factor(region, levels = full_coding_regions$region)) %>%
    # Remove rows where 'region' is NA
    filter(!is.na(region))
  
  return(plot_data)
}

# Example preparation of data
prepared_data_all_without_ct <- prepare_summary_data(updated_metadata_with_dnds_counts_and_totals_all_without_ct)
plot_data_all_without_ct <- prepare_plot_data(prepared_data_all_without_ct)
cleaned_data <- prepare_histogram_data(plot_data_all_without_ct)

# Plot the data
plot_all_without_ct <- ggplot(cleaned_data, aes(x = region, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "D") +
  theme_minimal(base_size = 10) +
  labs(
    title = "",
    x = "Coding Region",
    y = "# Changes",
    fill = "Mutation Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.9, 1.0),
        legend.justification = c("right", "top"))

# Print the plot
print(plot_all_without_ct)

total_counts <- cleaned_data %>%
  group_by(mutation_type) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
print(total_counts)





# FIRST AND LAST ANALYSES -----------
all_without_ct_changes_first_last_placements <- readRDS("")
dnds_summary_dataframe_all_without_ct <- readRDS("")

# FIRST AND LAST: Mutation locations ----------
coding_regions_df <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),
  start = c(266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
)

# Function to map mutation positions to coding regions and summarize by major_scorpio_call
map_mutations_to_regions <- function(mutations_df, coding_regions_df) {
  
  # Convert position data from string to numeric vectors
  mutations_df$change_positions <- lapply(strsplit(mutations_df$change_positions, ", "), as.numeric)
  
  # Initialize an empty list to store the results
  region_summary <- list()
  
  # Loop over each unique major_scorpio_call
  for (call in unique(mutations_df$major_scorpio_call)) {
    # Filter the mutations for the current major_scorpio_call
    mutations_for_call <- mutations_df[mutations_df$major_scorpio_call == call, ]
    
    # Loop over each coding region
    for (i in seq_len(nrow(coding_regions_df))) {
      region <- coding_regions_df$region[i]
      start <- coding_regions_df$start[i]
      end <- coding_regions_df$end[i]
      
      # Initialize a count for this region and call
      region_count <- 0
      
      # Loop through each mutation's positions
      for (positions in mutations_for_call$change_positions) {
        # Count how many positions fall within the current region
        region_count <- region_count + sum(positions >= start & positions <= end)
      }
      
      # Store the result for this region and call
      scorpio_summary <- data.frame(
        region = region,
        major_scorpio_call = call,
        mutation_count = region_count
      )
      
      # Append the result to the list
      region_summary[[length(region_summary) + 1]] <- scorpio_summary
    }
  }
  
  # Combine all the region summaries into a single data frame
  final_summary <- bind_rows(region_summary)
  
  # Ensure that all regions and major_scorpio_call combinations are represented
  all_combinations <- expand.grid(
    region = coding_regions_df$region,
    major_scorpio_call = unique(mutations_df$major_scorpio_call),
    stringsAsFactors = FALSE
  )
  
  # Merge the summary with all combinations to ensure no missing combinations
  final_summary <- full_join(all_combinations, final_summary, by = c("region", "major_scorpio_call"))
  
  # Replace NA mutation counts with 0
  final_summary$mutation_count[is.na(final_summary$mutation_count)] <- 0
  
  # Set factor levels to ensure correct order
  final_summary$region <- factor(final_summary$region, levels = coding_regions_df$region)
  
  return(final_summary)
}


plot_histogram_by_region <- function(region_summary) {
  ggplot(region_summary, aes(x = region, y = mutation_count, fill = major_scorpio_call)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis_d(option = "D") +  # Applying Viridis colors
    theme_classic() +
    labs(
      title = "",
      x = "Coding Region",
      y = "# Changes",
      fill = "Major Variant"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = c(0.9, 1.0),  # Inset legend at the top right
          legend.justification = c("right", "top"))
}

all_without_ct_first_last_summary <- map_mutations_to_regions(all_without_ct_changes_first_last_placements, coding_regions_df)
mutation_first_last_all_without_plot <- plot_histogram_by_region(all_without_ct_first_last_summary)
print(mutation_first_last_all_without_plot)




# Commonly occurring mutations ---------

mutation_aa_counts <- read.xlsx("") #From recurrent_mutations

mutation_aa_counts <- mutation_aa_counts %>%
  mutate(Region = sub(":.*", "", Amino.acid.change),  # Extract region before the colon
         Amino.acid.change = sub(".*:", "", Amino.acid.change))  # Extract change after the colon

coding_regions_df <- data.frame(
  region = c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
)

mutation_aa_counts$Region <- factor(mutation_aa_counts$Region, levels = coding_regions_df$region)

mutation_aa_counts$Amino.acid.change <- factor(mutation_aa_counts$Amino.acid.change, 
                                               levels = mutation_aa_counts$Amino.acid.change[order(-mutation_aa_counts$Count)])

aa_change_plot <- ggplot(mutation_aa_counts, aes(x = Amino.acid.change, y = Count, fill = Region)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +  # Add count labels on top of bars
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8, direction = 1) +  # Green to purple range
  labs(y = "Counts", x = "AA Change") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = c(0.5, 0.9),  # Place legend inside the plotting area
    legend.justification = c("center", "top"),  # Center the legend horizontally
    legend.direction = "horizontal",  # Make the legend horizontal
    legend.box = "horizontal",  # Use horizontal box layout for the legend
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  )





# FIRST AND LAST: Synonymous and non-synonymous locations ----------
prepare_summary_data <- function(df) {
  # Check if necessary columns exist
  if (!"region" %in% colnames(df)) {
    stop("Column 'region' not found in the dataframe.")
  }
  
  if (!"synonymous_count" %in% colnames(df)) {
    stop("Column 'synonymous_count' not found in the dataframe.")
  }
  
  if (!"nonsynonymous_count" %in% colnames(df)) {
    stop("Column 'nonsynonymous_count' not found in the dataframe.")
  }
  
  # Summarize the data by region, aggregating synonymous and nonsynonymous counts
  summary_data <- df %>%
    group_by(region) %>%
    summarise(
      Synonymous = sum(synonymous_count, na.rm = TRUE),
      Non_Synonymous = sum(nonsynonymous_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_data)
}

# Function to prepare plot data
prepare_plot_data <- function(summary_data) {
  # Reshape the data into a long format for plotting
  plot_data <- summary_data %>%
    pivot_longer(
      cols = c("Synonymous", "Non_Synonymous"),
      names_to = "mutation_type",
      values_to = "count"
    )
  
  return(plot_data)
}

# Full coding region definitions
full_coding_regions <- data.frame(
  region = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
)

# Function to prepare histogram data with ordered regions
prepare_histogram_data <- function(plot_data) {
  # Ensure the regions are ordered according to the specified levels
  plot_data <- plot_data %>%
    mutate(region = factor(region, levels = full_coding_regions$region)) %>%
    # Remove rows where 'region' is NA
    filter(!is.na(region))
  
  return(plot_data)
}

# Prepare and plot data for dnds_summary_dataframe_all_without_ct
prepared_data <- prepare_summary_data(dnds_summary_dataframe_all_without_ct)
plot_data_first_last <- prepare_plot_data(prepared_data)
first_last_plot_synonymous_nonsynonymous <- prepare_histogram_data(plot_data_first_last)
print(first_last_plot_synonymous_nonsynonymous)

# Plot the data
first_last_plot_synonymous_nonsynonymous <- ggplot(first_last_plot_synonymous_nonsynonymous, aes(x = region, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "D") +
  theme_classic() +
  labs(
    title = "",
    x = "Coding Region",
    y = "# Changes",
    fill = "Mutation Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.9, 1.0),
        legend.justification = c("right", "top"))














# FIRST AND LAST: Average rate vs infection duration ----------
all_without_ct_changes_first_last_placements <- readRDS("")

# Step 2: Filter out rows with total_change_count = 0
filtered_data <- all_without_ct_changes_first_last_placements %>%
  filter(total_change_count > 0)

# Step 3: Prepare data with log transformations
filtered_data <- filtered_data %>%
  mutate(
    log_days_between_sequences = log(total_days_between_sequences),  # Using natural log
    log_change_count = log(total_change_count)  # Using natural log
  )

correlation_result <- cor(filtered_data$log_days_between_sequences, 
                          filtered_data$log_change_count, 
                          use = "complete.obs")  # Use complete observations (ignore NAs)

# Print the correlation result
print(correlation_result)

# Step 4: Fit a linear-log model without zeros
linear_log_model_no_zeros <- brm(
  formula = log_change_count ~ log_days_between_sequences,
  data = filtered_data,
  family = gaussian(),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 2000,
  seed = 123
)

# Step 5: Fit a quadratic-log model without zeros
quadratic_log_model_no_zeros <- brm(
  formula = log_change_count ~ log_days_between_sequences + I(log_days_between_sequences^2),
  data = filtered_data,
  family = gaussian(),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 2000,
  seed = 123
)

# Step 6: Compute PSIS-LOO with moment matching if needed
loo_linear_log_no_zeros <- loo(linear_log_model_no_zeros, moment_match = TRUE)
loo_quadratic_log_no_zeros <- loo(quadratic_log_model_no_zeros, moment_match = TRUE)

# Step 7: Compare models
comparison_no_zeros <- loo_compare(loo_linear_log_no_zeros, loo_quadratic_log_no_zeros)
print(comparison_no_zeros)

# Step 8: Create a new data frame for predictions based on the log-transformed variable
new_data_log_no_zeros <- data.frame(
  log_days_between_sequences = seq(min(filtered_data$log_days_between_sequences, na.rm = TRUE),
                                   max(filtered_data$log_days_between_sequences, na.rm = TRUE),
                                   length.out = 100)
)

# Step 9: Get predictions from the linear log model without zeros
new_data_log_no_zeros$predicted_log_linear_changes <- fitted(linear_log_model_no_zeros, newdata = new_data_log_no_zeros, re_formula = NA, scale = "response")[, "Estimate"]

# Step 10: Get predictions from the quadratic log model without zeros
new_data_log_no_zeros$predicted_log_quadratic_changes <- fitted(quadratic_log_model_no_zeros, newdata = new_data_log_no_zeros, re_formula = NA, scale = "response")[, "Estimate"]

pred_intervals <- fitted(quadratic_log_model_no_zeros, 
                         newdata = new_data_log_no_zeros, 
                         re_formula = NA, 
                         scale = "response")

# Step 11: Plot the log-transformed data with the model fits
log_changes_time_plot_excluding_zeros <- ggplot(filtered_data, aes(x = log_days_between_sequences, y = log_change_count)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = new_data_log_no_zeros, 
              aes(x = log_days_between_sequences, ymin = lwr, ymax = upr),
              fill = 'darkslateblue', alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line(data = new_data_log_no_zeros, aes(x = log_days_between_sequences, y = predicted_log_quadratic_changes),
            color = 'darkslateblue', linewidth = 1,
            inherit.aes = FALSE) +
  labs(x = expression(Ln~(Total~Days~Between~Sequences)),
       y = expression(Ln~(Total~Change~Count))) +
  theme_minimal()

correlation_result <- cor(filtered_data$log_days_between_sequences, 
                          filtered_data$log_change_count, 
                          use = "complete.obs")  # Use complete observations (ignore NAs)

# Print the correlation result
print(correlation_result)









# Tip length analysis -----------
tip_lengths_df <- readRDS(file= "") #From tip_length_analysis folder and file

filtered_data <- tip_lengths_df %>%
  filter(Type %in% c("Case", "Control"))

# Step 2: Create bins for 'Length * 29891' (0-1, 1-2, ..., 49-50)
filtered_data$bin <- cut(filtered_data$Length * 29891, breaks = seq(0, 50, by = 1), right = FALSE, labels = paste0(seq(0, 49)))

# Step 3: Count observations in each bin for 'Case' and 'Control'
bin_counts <- filtered_data %>%
  group_by(bin, Type) %>%
  summarise(count = n(), .groups = 'drop')

# Step 4: Normalize the counts by dividing by the total counts for each 'Type' group
bin_counts <- bin_counts %>%
  group_by(Type) %>%
  mutate(proportion = count / sum(count))

# View the binned data and proportions
print(bin_counts)

# Step 5: Plot the results
histogram <- ggplot(bin_counts, aes(x = bin, y = proportion, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  labs(
    title = "",
    x = "Tip Length (Nucleotide Substitutions)",
    y = "Proportion",
    fill = "Group"
  ) +
  scale_x_discrete(name = "Tip Length (Nucleotide Substitutions)") +
  scale_fill_viridis_d(name = "Group", option = "viridis") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.85),  # Position the legend inside the plot area (top-right corner)
    legend.background = element_rect(fill = alpha("white", 0.5), color = "black"),  # Optional: background for the legend
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )





# Combined  plot --------------------

# Combine the plots
mutation_and_rate_combined_plot <- ggarrange(
  # First row (a and b)
  ggarrange(
    mutation_first_last_all_without_plot,                
    first_last_plot_synonymous_nonsynonymous,            
    ncol = 2,                                            
    labels = c("a", "b")
  ),
  # Second row (c and d)
  ggarrange(
    collapsed_mutation_plot,
    aa_change_plot,
    ncol = 2,                                           
    labels = c("c", "d")
  ),
  # Third row (e and f)
  ggarrange(
    log_changes_time_plot_excluding_zeros,
    histogram,
    ncol = 2,                                            
    labels = c("e", "f")
  ),
  ncol = 1,  # Align the three rows vertically
  heights = c(1, 0.8, 1)  # Adjust the heights of the rows
)

# Print the combined plot
print(mutation_and_rate_combined_plot)

ggsave("mutation_and_rate_combined_plot.pdf", 
       plot = mutation_and_rate_combined_plot, 
       device = "pdf", 
       width = 9,        # Width in inches
       height = 12,      # Height in inches
       dpi = 300)        # Resolution

supplementary_mutation_and_rate_combined_plot <- ggarrange(
  mutation_all_without_plot,  # First plot
  plot_all_without_ct,        # Second plot
  labels = c("a)", "b)"),
  ncol = 2,                   # Number of columns
  nrow = 1                    # Number of rows
)

ggsave("supplementary_mutation_and_rate_combined_plot.pdf", 
       plot = supplementary_mutation_and_rate_combined_plot, 
       device = "pdf", 
       width = 8,      
       height = 4,     
       dpi = 300)




