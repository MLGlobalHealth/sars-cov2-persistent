# Scripts for comparing tip lengths for persistent infections  ---------
# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(broom); library(mediation); library(brms); library(RcppEigen); library(ggplot2); library(brms);
library(bayestestR); library(mediation); library(rstanarm); library(survival); library(arm)


# Read in relevant data
conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- readRDS(file="")
conditional_subset <- conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson
first_last_dnds_all_without_ct <- readRDS(file= "")

# Filter event_id with at least 11 unique rows
filtered_event_ids <- conditional_subset %>%
  group_by(event_id) %>%
  filter(n() >= 11) %>%
  ungroup() %>%
  distinct(event_id) #293 unique

# Initialize a list to store sequence names per event_id
sequence_lists <- list()

# Process each event_id
for (eid in filtered_event_ids$event_id) {
  # Subset data for the current event_id
  event_data <- conditional_subset %>% filter(event_id == eid)
  
  # Separate by outcome
  outcome_1 <- event_data %>% filter(outcome == 1)
  outcome_0 <- event_data %>% filter(outcome == 0)
  
  # Get sequence names for outcome == 1
  seq_names_1 <- outcome_1 %>%
    inner_join(first_last_dnds_all_without_ct, by = c("PERSON_ID", "infection_episode")) %>%
    pull(last_strain)
  
  # Get sequence names for outcome == 0
  seq_names_0 <- outcome_0 %>% pull(strain)
  
  # Combine sequence names
  all_seq_names <- unique(c(seq_names_1, seq_names_0))
  
  # Save sequence names in the list
  sequence_lists[[eid]] <- all_seq_names
  
  # Write to a text file
  output_dir <- ""
  writeLines(all_seq_names, file.path(output_dir, paste0(eid, "_sequences.txt")))
}

# FASTA processing will be handled by Python script, get_eventid_sequences.py

# Then run run_IQTREE_all.py to run IQTREE for each set of sequences

# Importing trees from IQTREE back
# List all .treefile files in the directory
tree_files <- list.files(pattern = "\\.treefile$")

# Initialize an empty dataframe to store results
tip_lengths_df <- data.frame(Tree = character(),
                             Tip = character(),
                             Length = numeric(),
                             Type = character(),
                             EventID = character(),
                             stringsAsFactors = FALSE)


# Process each treefile
for (tree_file in tree_files) {
  # Read the tree
  tree <- read.tree(tree_file)
  
  # Extract event_id from tree filename
  event_id <- sub("\\.fasta.treefile$", "", tree_file)
  
  # Extract tip labels and branch lengths
  tips <- tree$tip.label
  lengths <- tree$edge.length[tree$edge[, 2] <= length(tips)]  # Only tip edges
  
  # Classify each tip
  classifications <- sapply(tips, function(tip_name) {
    if (tip_name == "hCoV-19_Wuhan_WIV04_2019") {
      return("Outgroup")
    } else if (tip_name %in% first_last_dnds_all_without_ct$last_strain) {
      return("Case")
    } else if (tip_name %in% conditional_subset$strain) {
      outcome <- conditional_subset %>%
        filter(strain == tip_name) %>%
        pull(outcome)
      if (length(outcome) > 0 && outcome[1] == 0) {
        return("Control")
      }
    }
    return("unknown")
  })
  
  # Create a dataframe for this tree
  tree_df <- data.frame(
    Tree = tree_file,
    Tip = tips,
    Length = lengths,
    Type = classifications,
    EventID = event_id,
    stringsAsFactors = FALSE
  )
  
  # Append to the main dataframe
  tip_lengths_df <- bind_rows(tip_lengths_df, tree_df)
}


# -----------------------------------
# Plot the density of tip lengths
# -----------------------------------
# Step 1: Filter out only 'Case' and 'Control' types
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
    title = "Histogram of Tip Lengths",
    x = "# Nucleotide Substitutions",
    y = "Proportion",
    fill = "Type"
  ) +
  scale_x_discrete(name = "# Nucleotide Substitutions") +
  scale_fill_viridis_d(name = "Type", option = "viridis") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.85),  # Position the legend inside the plot area (top-right corner)
    legend.background = element_rect(fill = alpha("white", 0.5), color = "black"),  # Optional: background for the legend
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Statistical testing ----------
filtered_data$Length_adjusted <- filtered_data$Length * 29891 + 1e-6

model <- brm(
  formula = Length_adjusted ~ Type,  # 'Type' as predictor
  data = filtered_data,
  family = lognormal(),  # Log-normal distribution (log-transformed response variable)
  prior = c(
    prior(normal(0, 1), class = "b"),  # Prior for the coefficients
    prior(normal(0, 1), class = "Intercept"),  # Prior for the intercept
    prior(exponential(1), class = "sigma")  # Prior for the standard deviation
  ),
  chains = 4,  # Number of MCMC chains
  iter = 2000,  # Number of iterations per chain
  warmup = 1000,  # Number of warmup iterations
  control = list(adapt_delta = 0.95)  # Control the convergence criteria
)

# View the model summary
summary(model)