#Figure 1, Baseline characteristics

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr);
library(tidyr); library(patchwork); library(viridis); library(gridExtra); library(ggpubr); library(brms)

all_without_ct_changes_first_last_placements <- readRDS("")
all_without_ct_changes_first_last_placements$num_vacc_at_infection <- as.numeric(as.character(all_without_ct_changes_first_last_placements$num_vacc_at_infection))
all_without_ct_metadata <- readRDS(file = "")
likely_without_ct_metadata <- readRDS(file = "")
confirmed_without_ct_metadata <- readRDS(file = "")

#Slicing data to get first row
all_without_ct_metadata_unique <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  slice_min(DATESAMPLING, with_ties = FALSE) %>%  # Automatically resolves ties by choosing one row
  ungroup() %>%
  mutate(DATESAMPLING = as.Date(DATESAMPLING),    # Convert DATESAMPLING to Date type if needed
         age_group = cut(age_at_infection,        # Group age_at_infection into brackets
                         breaks = c(0, 15, 30, 45, 60, 75, Inf),  # Specify age group boundaries
                         labels = c("0-15", "15-30", "30-45", "45-60", "60-75", "75+"),  # Labels for groups
                         right = FALSE))

#Information for Supplementary Tables -------
table(all_without_ct_metadata_unique$major_scorpio_call)
table(all_without_ct_metadata_unique$age_group)
table(all_without_ct_metadata_unique$KOEN)
table(all_without_ct_changes_first_last_placements$num_vacc_at_infection)
#Naive status
likely_without_ct_metadata_unique <- likely_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  slice_min(DATESAMPLING, with_ties = FALSE)
nrow(likely_without_ct_metadata_unique)
confirmed_without_ct_metadata_unique <- confirmed_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  slice_min(DATESAMPLING, with_ties = FALSE)
nrow(confirmed_without_ct_metadata_unique)


# Plot A) Time of infection by variant --------
all_without_ct_metadata_unique$DATESAMPLING <- as.Date(all_without_ct_metadata_unique$DATESAMPLING)

obs_count <- all_without_ct_metadata_unique %>%
  group_by(major_scorpio_call) %>%
  summarise(count = n(), .groups = "drop")

# Plot with box plot, colored boxes, and observation counts
infection_duration_plot <- ggplot(all_without_ct_metadata_unique, aes(x = DATESAMPLING, y = major_scorpio_call, fill = major_scorpio_call)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Box plot without outliers
  scale_fill_viridis(discrete = TRUE) +  # Color the boxes using viridis palette
  geom_text(data = obs_count, 
            aes(x = as.Date("2020-11-01"), y = major_scorpio_call, label = paste("n=", count)), # Adjust 'x' to control label position
            color = "black", size = 2.5, nudge_x = 20) +  # Position the count slightly above the boxes
  labs(title = "",
       x = "Infection Start Date",
       y = "Variant") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")


# Plot B) Duration of infection episode--------
infection_duration_df <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  summarize(
    first_date = min(DATESAMPLING),
    last_date = max(DATESAMPLING),
    duration_days = as.numeric(last_date - first_date),
    major_scorpio_call = first(major_scorpio_call),  # Get the first major_scorpio_call for each group
    .groups = 'drop'
  ) %>%
  select(PERSON_ID, infection_episode, duration_days, major_scorpio_call)

duration_plot <- ggplot(infection_duration_df, aes(x = major_scorpio_call, y = duration_days, fill = major_scorpio_call)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Box plot without outliers shown
  scale_fill_viridis(discrete = TRUE) +  # Color the boxes using viridis palette
  labs(
    title = "",
    x = "",
    y = "Infection Duration (Days)"
  ) +
  coord_flip() +  # Flip the axes
  ylim(10, 100) +  # Set the max value for the x-axis to 150 days
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = "none")  # Hide the legend






# Plot C and D) Age and sex distribution--------
age_bins <- c(0, 15, 30, 45, 60, 75, Inf)
age_labels <- c("0-15", "15-30", "30-45", "45-60", "60-75", "75+")
age_sex_distribution <- all_without_ct_metadata_unique %>%
  mutate(
    age_group = cut(age_at_infection, breaks = age_bins, labels = age_labels, right = FALSE),
    Sex = ifelse(KOEN == 1, "Male", "Female")
  ) %>%
  group_by(age_group, Sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(!is.na(Sex)) %>%  # Remove rows with NA in sex
  mutate(Sex = factor(Sex, levels = c("Male", "Female")))  # Convert sex to a factor

# Age Distribution Plot (No Sex Information)
age_plot <- ggplot(age_sex_distribution, aes(x = age_group, y = count)) +
  geom_bar(stat = "identity", fill = "#A3DAD2FF") +  # No fill color for sex
  labs(x = "Age Group", y = "Counts") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Hide legend
  )

# Sex Distribution Plot (No Age Information)
sex_plot <- ggplot(age_sex_distribution, aes(x = Sex, y = count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +  # Side-by-side bars for each sex
  scale_fill_manual(values = c("Male" = "#A3DAD2FF", "Female" = "#2A788DFF")) +  # Color for male and female
  labs( x = "Sex", y = "Percent") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend
  )





# Plot E) Number of sequences per person--------
sequence_counts <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  summarise(num_sequences = n()) %>%
  ungroup()
sequence_counts_grouped <- sequence_counts %>%
  mutate(
    num_sequences_grouped = ifelse(num_sequences >= 8, "8+", as.character(num_sequences))
  )

sequence_count_plot <- ggplot(sequence_counts_grouped, aes(x = num_sequences_grouped)) +
  geom_bar(
    fill = "#004B87FF", 
    color = "black", 
    alpha = 0.7
  ) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    size=3.5,
    vjust = -0.5  # Position text slightly above the bars
  ) +
  labs(
    title = "",
    x = "# Sequences Per Person",
    y = "Counts"
  ) +
  scale_x_discrete(
    limits = c(as.character(2:7), "8+")  # Ensure 8 or more appears at the end
  ) +
  ylim(0, 200) +  # Set the y-axis limit
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text.x = element_text(angle = 30, hjust = 1)  # Rotate x-axis labels for readability
  )





# Plot F) Histograms of First - Last values ---------
two_or_more_all_without_ct_metadata <- all_without_ct_metadata %>%
  group_by(PERSON_ID, infection_episode) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  mutate(CT = as.numeric(CT)) %>% # Ensure CT is numeric
  arrange(PERSON_ID, infection_episode, DATESAMPLING) %>%
  group_by(PERSON_ID, infection_episode) %>%
  mutate(days_since_first = as.numeric(difftime(DATESAMPLING, min(DATESAMPLING), units = "days"))) %>%
  ungroup() %>%
  filter(!is.na(CT)) %>%
  group_by(PERSON_ID, infection_episode) %>%
  filter(n() >= 2) %>%
  ungroup()

# Function to calculate the CT difference between the first and last days
calculate_ct_difference <- function(df) {
  df %>%
    group_by(PERSON_ID, infection_episode) %>%
    summarise(
      first_ct = CT[which.min(days_since_first)],
      last_ct = CT[which.max(days_since_first)],
      ct_difference = first_ct - last_ct
    ) %>%
    ungroup()
}

# Function to create histogram of CT differences
plot_ct_difference_histogram <- function(df) {
  ggplot(df, aes(x = ct_difference)) +
    geom_histogram(
      binwidth = 1,
      fill = viridis::viridis(1, option = "D")[1],  # Use purple from Viridis palette
      color = "black",
      alpha = 0.6  # Set transparency level
    ) +
    labs(
      x = "Ct Value Difference (First - Last)",
      y = "Counts"
    ) +
    theme_minimal() +
    geom_vline(
      xintercept = 0,
      color = "black",
      linewidth = 1.5  # Set the thickness of the vertical line
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks())  # Ensure y-axis shows only integers
}

# Calculate differences and plot histograms
all_without_ct_diff <- calculate_ct_difference(two_or_more_all_without_ct_metadata)
# Generate histograms
hist_all_without_ct_plot <- plot_ct_difference_histogram(all_without_ct_diff)



# Combining plots --------
combined_plot <- ggarrange(
  infection_duration_plot, duration_plot,
  age_plot, sex_plot, sequence_count_plot, hist_all_without_ct_plot,
  ncol = 2, nrow = 3,
  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
  common.legend = FALSE,
  align = "h",
  widths = c(1, 1),  # Give column with infection_duration_plot slightly more space
  heights = c(1, 1, 1)  # Adjust row heights as needed
)

ggsave(
  filename = "combined_baseline_plot.pdf",         
  plot = combined_plot,                     
  device = "pdf",                          
  width = 7,                               
  height = 7,                              
  units = "in",                            
  dpi = 300                                
)


# Ct differences, Statistical testing ---------
all_without_ct_diff
# Compute summary statistics for the ct_difference column
summary_stats <- all_without_ct_diff %>%
  summarise(
    median_ct_diff = median(ct_difference, na.rm = TRUE),
    q1_ct_diff = quantile(ct_difference, 0.25, na.rm = TRUE),  # 25th percentile
    q3_ct_diff = quantile(ct_difference, 0.75, na.rm = TRUE),  # 75th percentile
    mean_ct_diff = mean(ct_difference, na.rm = TRUE),
    sd_ct_diff = sd(ct_difference, na.rm = TRUE),
    min_ct_diff = min(ct_difference, na.rm = TRUE),
    max_ct_diff = max(ct_difference, na.rm = TRUE)
  )
# Define the prior: Normally distributed with mean 0, wide SD to reflect uncertainty
priors <- c(
  prior(normal(0, 10), class = "Intercept")  # Prior on the intercept (mean ct_difference)
)

# Fit the Bayesian model
bayesian_model <- brm(
  ct_difference ~ 1,                     # Model for ct_difference with only intercept
  data = all_without_ct_diff,            # Data frame with ct_difference column
  prior = priors,                         # Specified priors
  iter = 4000,                            # Increase iterations for better convergence
  warmup = 1000,                          # Warmup samples
  chains = 4,                             # Number of chains for MCMC
  seed = 1234                             # Set seed for reproducibility
)

# Print a summary of the model
summary(bayesian_model)

# Extract posterior samples using as_draws_df
posterior_samples <- as_draws_df(bayesian_model)

# Calculate the posterior probability that ct_difference mean is less than 0
prob_mean_less_than_zero <- mean(posterior_samples$b_Intercept < 0)
cat("Posterior probability that the mean ct_difference is less than 0:", prob_mean_less_than_zero, "\n")

# Extract 95% credible interval for the mean (intercept) of ct_difference
posterior_summary <- posterior_summary(bayesian_model, pars = "b_Intercept")
mean_ct_diff <- posterior_summary[1, "Estimate"]
ci_lower <- posterior_summary[1, "Q2.5"]
ci_upper <- posterior_summary[1, "Q97.5"]
posterior_summary
mean_ct_diff
# Print 95% credible interval
cat("95% Credible Interval for mean ct_difference:", ci_lower, "to", ci_upper, "\n")
