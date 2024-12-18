# Figure 3: Visualization of dnds results

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr);
library(tidyverse); library(data.table); library(dbplyr); library(lubridate); 
library(rlang); library(foreach); library(doParallel); library(DSTora); library(ROracle); 
library(DSTcolectica); library(DSTdb); library(DBI); library(parallel); library(ggsignif); 
library(Rcpp); library(purrr); library(tidyr); library(furrr); library(future); library(future.apply); 
library(seqinr); library(adegenet); library(ggplot2); library(viridis); library(lme4); library(broom.mixed); 
library(brms); library(ggpubr); library(scales)

dnds_counts_all_without_ct <- readRDS(file="")
dnds_counts_controls <- readRDS(file = "")


# Persistent infection analysis ----------

# Step 1: Filter the data to only include rows where ka and ks are non-negative and finite
filtered_data <- dnds_counts_all_without_ct %>%
  filter(is.finite(ka), is.finite(ks), ka >= 0, ka < 1, ks >= 0, ks < 1)

# Step 2: Collapse (sum) the values for each pair_id
collapsed_data <- filtered_data %>%
  group_by(pair) %>%
  summarise(
    total_ka = sum(ka, na.rm = TRUE),
    total_ks = sum(ks, na.rm = TRUE)) %>%
  filter(total_ka > 0 | total_ks > 0)

# Step 3: Recalculate ka/ks ratio for each collapsed pair
collapsed_data$collapsed_ka_ks <- collapsed_data$total_ka - collapsed_data$total_ks

collapsed_data_clean <- collapsed_data %>%
  filter(!is.na(collapsed_ka_ks))

# View the clean collapsed data
head(collapsed_data_clean)
nrow(collapsed_data_clean)
hist(collapsed_data_clean$collapsed_ka_ks)
mean(collapsed_data_clean$collapsed_ka_ks)
median(collapsed_data_clean$collapsed_ka_ks)
quantile(collapsed_data_clean$collapsed_ka_ks, probs = c(0.25, 0.75), na.rm = TRUE)

# Calculate the median of collapsed_ka_ks
median_ka_ks <- median(collapsed_data_clean$collapsed_ka_ks, na.rm = TRUE)

# Plot with histogram and median line
individual_dnds <- ggplot(collapsed_data_clean, aes(x = collapsed_ka_ks)) +
  geom_histogram(bins = 30, fill = viridis(1, option = "D"), color = "black", alpha = 0.7) +  # Dark viridis color
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  geom_vline(xintercept = median_ka_ks, color = "#1F968BFF", linetype = "solid", linewidth = 1) +  # Median line with viridis color
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = "Counts"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
nrow(collapsed_data_clean)

# Statistical testing ------
model <- brm(
  formula = collapsed_ka_ks ~ 1,   # Intercept-only model
  data = collapsed_data_clean,      # Use your data
  family = gaussian(),              # Normal distribution (Gaussian family)
  prior = c(
    prior(normal(0, 0.01), class = "Intercept")  # Prior for the intercept (mean)
  ),
  chains = 4,                      # Number of MCMC chains
  iter = 2000,                     # Number of iterations
  warmup = 1000,                   # Number of warm-up iterations
  control = list(adapt_delta = 0.95)  # Control for sampler stability
)

# Check the posterior samples
summary(model)

# Extract the posterior distribution of the intercept (mean)
posterior <- as.data.frame(model)

posterior_summary <- summary(model)$fixed
posterior_summary["Intercept", "Estimate"]
posterior_summary["Intercept", "l-95% CI"]
posterior_summary["Intercept", "u-95% CI"]


# Same but for each gene separately
collapsed_gene_data <- filtered_data %>%
  filter(ka > 0 | ks > 0) %>%
  mutate(collapsed_ka_ks = ka - ks)

collapsed_gene_data$gene <- factor(collapsed_gene_data$gene, 
                                   levels = gene_regions$region)

# Compute the median for each gene
medians <- collapsed_gene_data %>%
  group_by(gene) %>%
  summarise(median_ka_ks = median(collapsed_ka_ks, na.rm = TRUE))


# Plot with two vertical lines: one at 0 and one at the median
gene_dnds <- ggplot(collapsed_gene_data, aes(x = collapsed_ka_ks)) +
  geom_histogram(bins = 40, fill = viridis(1, option = "D"), color = "black", alpha = 0.7) +  # Dark viridis color
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  geom_vline(data = medians, aes(xintercept = median_ka_ks), 
             linetype = "solid", color="#1F968BFF", linewidth = 1) +  # Line at median with viridis color
  facet_wrap(~ gene, scales = "free", nrow = 4, ncol = 3) +  
  theme_bw() +
  labs(
    title = "",
    x = "Difference between Ka and Ks Values (Ka-Ks)",
    y = "Counts"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Statistical testing for each gene

results_list <- list()

# Loop through each gene
for (gene_name in levels(collapsed_gene_data$gene)) {
  # Subset data for this gene
  gene_data <- collapsed_gene_data %>%
    filter(gene == gene_name)
  
  # Fit Bayesian model for this gene
  model <- brm(
    formula = collapsed_ka_ks ~ 1,   # Intercept-only model
    data = gene_data,
    family = gaussian(),             # Normal distribution (Gaussian family)
    prior = c(
      prior(normal(0, 0.01), class = "Intercept")  # Prior for the intercept (mean)
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.95)
  )
  
  # Extract posterior summary
  posterior_summary <- summary(model)$fixed
  estimate <- posterior_summary["Intercept", "Estimate"]
  lower_95 <- posterior_summary["Intercept", "l-95% CI"]
  upper_95 <- posterior_summary["Intercept", "u-95% CI"]
  
  # Store results
  results_list[[gene_name]] <- data.frame(
    gene = gene_name,
    estimate = estimate,
    lower_95 = lower_95,
    upper_95 = upper_95
  )
}

# Combine results into a single table
results_table_persistent <- do.call(rbind, results_list)
saveRDS(results_table_persistent, file="results_table_persistent.RDS")








# Controls / acute infection analysis ----------

#Collapsing by person

# Step 1: Filter the data to only include rows where ka and ks are non-negative and finite
filtered_data_controls <- dnds_counts_controls %>%
  filter(is.finite(ka), is.finite(ks), ka >= 0, ka < 1, ks >= 0, ks < 1) %>%
  mutate(person_id = sub("_.*", "", pair)) %>% # Extract the first part of 'pair'
  group_by(person_id, gene) %>%
  slice(1) %>% # Keep the first occurrence
  ungroup() %>%
  select(-person_id) # Remove the temporary 'person_id' column if no longer needed


# Step 2: Collapse (sum) the values for each pair_id
collapsed_data_controls <- filtered_data_controls %>%
  group_by(pair) %>%
  summarise(
    total_ka = sum(ka, na.rm = TRUE),
    total_ks = sum(ks, na.rm = TRUE)) %>%
  filter(total_ka > 0 | total_ks > 0)

# Step 3: Recalculate ka/ks ratio for each collapsed pair
collapsed_data_controls$collapsed_ka_ks <- collapsed_data_controls$total_ka - collapsed_data_controls$total_ks

collapsed_data_clean_controls <- collapsed_data_controls %>%
  filter(!is.na(collapsed_ka_ks))

# View the clean collapsed data
head(collapsed_data_clean_controls)
nrow(collapsed_data_clean_controls)
hist(collapsed_data_clean_controls$collapsed_ka_ks)
mean(collapsed_data_clean_controls$collapsed_ka_ks)
median(collapsed_data_clean_controls$collapsed_ka_ks)
quantile(collapsed_data_clean_controls$collapsed_ka_ks, probs = c(0.25, 0.75), na.rm = TRUE)

# Calculate the median of collapsed_ka_ks
median_ka_ks <- median(collapsed_data_clean_controls$collapsed_ka_ks, na.rm = TRUE)

# Plot with histogram and median line
individual_dnds_controls <- ggplot(collapsed_data_clean_controls, aes(x = collapsed_ka_ks)) +
  geom_histogram(bins = 30, fill = viridis(1, option = "D"), color = "black", alpha = 0.7) +  # Dark viridis color
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  geom_vline(xintercept = median_ka_ks, color = "#1F968BFF", linetype = "solid", linewidth = 1) +  # Median line with viridis color
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
nrow(collapsed_data_clean_controls) #167

# Statistical testing ------
model_controls <- brm(
  formula = collapsed_ka_ks ~ 1,   # Intercept-only model
  data = collapsed_data_clean_controls,      # Use your data
  family = gaussian(),              # Normal distribution (Gaussian family)
  prior = c(
    prior(normal(0, 0.01), class = "Intercept")  # Prior for the intercept (mean)
  ),
  chains = 4,                      # Number of MCMC chains
  iter = 2000,                     # Number of iterations
  warmup = 1000,                   # Number of warm-up iterations
  control = list(adapt_delta = 0.95)  # Control for sampler stability
)

# Check the posterior samples
summary(model_controls)

# Extract the posterior distribution of the intercept (mean)
posterior <- as.data.frame(model_controls)

posterior_summary <- summary(model_controls)$fixed
posterior_summary["Intercept", "Estimate"]
posterior_summary["Intercept", "l-95% CI"]
posterior_summary["Intercept", "u-95% CI"]




# Same but for each gene separately

collapsed_gene_data_controls <- filtered_data_controls %>%
  filter(ka > 0 | ks > 0) %>%
  mutate(collapsed_ka_ks = ka - ks)

collapsed_gene_data_controls$gene <- factor(collapsed_gene_data_controls$gene, 
                                   levels = gene_regions$region)

# Compute the median for each gene
medians_controls <- collapsed_gene_data_controls %>%
  group_by(gene) %>%
  summarise(median_ka_ks = median(collapsed_ka_ks, na.rm = TRUE))

# Plot with two vertical lines: one at 0 and one at the median
gene_dnds_controls <- ggplot(collapsed_gene_data_controls, aes(x = collapsed_ka_ks)) +
  geom_histogram(bins = 40, fill = viridis(1, option = "D"), color = "black", alpha = 0.7) +  # Dark viridis color
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  geom_vline(data = medians_controls, aes(xintercept = median_ka_ks), 
             linetype = "solid", color="#1F968BFF", linewidth = 1) +  # Line at median with viridis color
  facet_wrap(~ gene, scales = "free", nrow = 4, ncol = 3) +  # Facet by gene
  theme_bw() +
  labs(
    title = "",
    x = "Difference between Ka and Ks Values (Ka-Ks)",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Statistical testing for each gene ---------

results_list_controls <- list()

# Loop through each gene
for (gene_name in levels(collapsed_gene_data_controls$gene)) {
  # Subset data for this gene
  gene_data <- collapsed_gene_data_controls %>%
    filter(gene == gene_name)
  
  # Fit Bayesian model for this gene
  model <- brm(
    formula = collapsed_ka_ks ~ 1,   # Intercept-only model
    data = gene_data,
    family = gaussian(),             # Normal distribution (Gaussian family)
    prior = c(
      prior(normal(0, 0.01), class = "Intercept")  # Prior for the intercept (mean)
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.95)
  )
  
  # Extract posterior summary
  posterior_summary <- summary(model)$fixed
  estimate <- posterior_summary["Intercept", "Estimate"]
  lower_95 <- posterior_summary["Intercept", "l-95% CI"]
  upper_95 <- posterior_summary["Intercept", "u-95% CI"]
  
  # Store results
  results_list_controls[[gene_name]] <- data.frame(
    gene = gene_name,
    estimate = estimate,
    lower_95 = lower_95,
    upper_95 = upper_95
  )
}

# Combine results into a single table
results_table_controls <- do.call(rbind, results_list_controls)



# Combined plots per gene ------------

collapsed_gene_data_filtered <- collapsed_gene_data %>% 
  filter(gene %in% c("ORF1a", "ORF1b", "S", "N"))

collapsed_gene_data_controls_filtered <- collapsed_gene_data_controls %>% 
  filter(gene %in% c("ORF1a", "ORF1b", "S", "N"))

collapsed_gene_data_filtered <- collapsed_gene_data_filtered %>% 
  mutate(Group = "Persistent")

collapsed_gene_data_controls_filtered <- collapsed_gene_data_controls_filtered %>% 
  mutate(Group = "Controls")

combined_data <- bind_rows(collapsed_gene_data_filtered, collapsed_gene_data_controls_filtered)

medians_filtered <- medians %>% 
  filter(gene %in% c("ORF1a", "ORF1b", "S", "N")) %>% 
  mutate(Group = "Persistent")

medians_controls_filtered <- medians_controls %>% 
  filter(gene %in% c("ORF1a", "ORF1b", "S", "N")) %>% 
  mutate(Group = "Controls")

combined_medians <- bind_rows(medians_filtered, medians_controls_filtered)

descriptive_stats <- combined_data %>%
  group_by(gene, Group) %>%
  summarise(
    mean_value = mean(collapsed_ka_ks, na.rm = TRUE),
    median_value = median(collapsed_ka_ks, na.rm = TRUE),
    lower_quartile = quantile(collapsed_ka_ks, 0.25, na.rm = TRUE),
    upper_quartile = quantile(collapsed_ka_ks, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Print descriptive statistics
print(descriptive_stats)

ks_results <- combined_data %>%
  group_by(gene) %>%
  summarise(
    ks_statistic = ks.test(
      x = collapsed_ka_ks[Group == "Persistent"],
      y = collapsed_ka_ks[Group == "Controls"]
    )$statistic,
    p_value = ks.test(
      x = collapsed_ka_ks[Group == "Persistent"],
      y = collapsed_ka_ks[Group == "Controls"]
    )$p.value
  ) %>%
  mutate(p_label = paste0("p = ", format(p_value, digits = 3, scientific = TRUE)))

# Merge KS results with combined_data for plotting
combined_data <- combined_data %>%
  left_join(ks_results %>% select(gene, p_label), by = "gene")

# Create the plot with p-values displayed
combined_plot <- ggplot(combined_data, aes(x = collapsed_ka_ks, fill = Group)) +
  geom_histogram(bins = 40, color = "black", alpha = 0.7, position = "identity") +  # Overlapping histograms
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  geom_vline(data = combined_medians, aes(xintercept = median_ka_ks, color = Group),
             linetype = "solid", linewidth = 1) +  # Median lines
  facet_wrap(~ gene, scales = "free", nrow = 2, ncol = 2) +  # Automatically adjusts layout
  scale_fill_viridis_d(option = "D") +  # Use viridis colors for fill
  scale_color_viridis_d(option = "D") +  # Use viridis colors for lines
  geom_text(
    data = ks_results, 
    aes(x = -Inf, y = Inf, label = p_label), 
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.5, size = 2.7 # Smaller text size
  ) +  # Add p-value text to top-left corner of each facet
  theme_bw() +
  labs(
    title = "",
    x = "Difference between Ka and Ks Values (Ka-Ks)",
    y = "Counts",
    fill = "Group",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(combined_plot)




# Combining plots ------------

collapsed_data_clean$Group <- "Persistent"
collapsed_data_clean_controls$Group <- "Controls"

# Combine Persistent and Controls data into a new dataframe
merged_data <- bind_rows(collapsed_data_clean, collapsed_data_clean_controls)

# Calculate the median for each group
median_persistent <- median(collapsed_data_clean$collapsed_ka_ks, na.rm = TRUE)
median_controls <- median(collapsed_data_clean_controls$collapsed_ka_ks, na.rm = TRUE)

# Perform KS test comparing the two groups
ks_test <- ks.test(
  collapsed_data_clean$collapsed_ka_ks,
  collapsed_data_clean_controls$collapsed_ka_ks
)

# Extract and format p-value for scientific notation
p_value <- format(ks_test$p.value, scientific = TRUE, digits = 2)
p_label <- paste0("p = ", p_value)

merged_plot <- ggplot(merged_data, aes(x = collapsed_ka_ks, fill = Group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7, position = "identity") +  # Overlapping histograms
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  # Line at 0
  scale_fill_viridis_d(option = "D") +  # Use viridis color palette for fill
  annotate("text", x = -Inf, y = Inf, label = p_label, hjust = -0.1, vjust = 1.5, size = 3.5) +  # Add p-value on top-left corner
  theme_minimal() +
  labs(
    title = "",
    x = "Difference between Ka and Ks Values (Ka-Ks)",
    y = "Counts",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )

print(merged_plot)

combined_dnds_plots <- ggarrange(
  individual_dnds,
  individual_dnds_controls,
  gene_dnds,
  gene_dnds_controls,
  labels = c("a", "b", "c", "d"),
  ncol = 2,                   
  nrow = 2,   
  heights = c(0.4, 1.5)
)
combined_dnds_plots

ggsave("combined_dnds_plots.pdf", 
       plot = combined_dnds_plots, 
       device = "pdf", 
       width = 9,      
       height = 12,     
       dpi = 300)

streamlined_dnds_plots <- ggarrange(
  merged_plot,
  combined_plot,
  labels = c("a", "b"),
  ncol = 1,                   
  nrow = 2,   
  heights = c(0.6, 1.5), 
  common.legend = TRUE,   
  legend = "right"      
)

ggsave("streamlined_dnds_plots.pdf", 
       plot = streamlined_dnds_plots, 
       device = "pdf", 
       width = 6.75,      
       height = 9,     
       dpi = 300)




