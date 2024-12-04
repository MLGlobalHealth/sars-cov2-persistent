
# Preparing controls and adding Charlson comorbidity index values ----------------------

# read packages
library(ape); library(phytools); library(TreeTools); library(dplyr); library(tidyverse);
library(data.table); library(dbplyr); library(lubridate); library(rlang); library(foreach);
library(doParallel); library(DSTora); library(ROracle); library(DSTcolectica); library(DSTdb);
library(DBI); library(parallel); library(ggsignif); library(Rcpp); library(purrr); library(tidyr);
library(MatchIt);library(dplyr);

# Case-control with only sociodemographic information ------------------------

# Loading persistent infections and controls
subset_cases_and_controls_all_with_metadata_geo_diag <- readRDS(file="subset_cases_and_controls_all_with_metadata_geo_diag.RDS")

# Step 1: Keep only the first row for each PERSON_ID and DATESAMPLING where type == "case"
cases_and_controls_geo_diag_filtered <- subset_cases_and_controls_all_with_metadata_geo_diag %>%
  group_by(PERSON_ID, DATESAMPLING) %>%
  filter(type == "case") %>%
  slice(1) %>%  # Keep only the first occurrence for each group
  ungroup()

# Step 2: Identify the event_id values from the filtered rows
event_ids_to_keep <- cases_and_controls_geo_diag_filtered$event_id

# Step 3: Keep all rows (cases or controls) with the identified event_id values
subset_cases_and_controls_all_with_metadata_geo_diag <- subset_cases_and_controls_all_with_metadata_geo_diag %>%
  filter(event_id %in% event_ids_to_keep)

#Filtering out controls without sufficient metadata (birthday and sex)
subset_cases_and_controls_all_with_metadata_geo_diag <- subset_cases_and_controls_all_with_metadata_geo_diag %>%
  filter(!is.na(BIRTHDAY) & !is.na(KOEN))

subset_cases_and_controls_all_with_metadata_geo_diag_limited <- subset_cases_and_controls_all_with_metadata_geo_diag

subset_cases_and_controls_all_with_metadata_geo_diag_limited$outcome <- ifelse(
  subset_cases_and_controls_all_with_metadata_geo_diag_limited$type == "case", 1, 0
)
subset_cases_and_controls_all_with_metadata_geo_diag_limited$age_group <- relevel(
  subset_cases_and_controls_all_with_metadata_geo_diag_limited$age_group, ref = "15-30"
)


# Adding Charlson comorbidity index values ----------------
library(comorbidity)
library(heaven)
drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej('STATPROD', dbuser = 'XP2')

# Define charlson.codes with codes starting with "D"
charlson.codes <- list(
  myocardial.infarction = c('410', 'DI21', 'DI22'),
  heart.failure = c('42709', '42710', '42711', '42719', '42899', '78249', 'DI099', 'DI110', 'DI130', 'DI132', 'DI255', 'DI425', 'DI426', 'DI427', 'DI429', 'DI428A', 'DP290', 'DI43', 'DI50', 'DE105', 'DE115', 'DE125', 'DE135', 'DE145'),
  peripheral.vascular.disease = c('440', '441', '442', '443', '444', '445', 'DI70', 'DI71', 'DI72', 'DI731', 'DI738', 'DI739', 'DI77', 'DI790', 'DI792', 'DK551', 'DK558', 'DK559', 'DZ958', 'DZ959'),
  cerebrovascular.disease = c(paste0('43', 0:8), paste0('DI6', 0:9), 'DG45', 'DG46', 'DH340'),
  dementia = c('290', paste0('DF0', 0:3), 'DG30', 'DF051', 'DG311'),
  chronic.pulmonary.disease = c(paste0('51', 5:8), paste0('49', 0:3), paste0('DJ4', 0:7), paste0('DJ6', 0:7), 'DJ684', 'DI278', 'DI279', 'DJ84', 'DJ701', 'DJ703', 'DJ920', 'DJ953', 'DJ961', 'DJ982', 'DJ983'), 
  rheumatic.disease = c('712', '716', '734', '446', '13599', 'DM05', 'DM06', 'DM08', 'DM09', 'DM30', 'DM31', 'DM32', 'DM33', 'DM34', 'DM35', 'DM36', 'D86'),
  peptic.ulcer.disease = c('53091', '53098', paste0('53', 1:4), 'DK25', 'DK26', 'DK27', 'DK28', 'DK221'),
  mild.liver.disease = c('571', '57301', '57304', 'DB18', 'DK700', 'DK701', 'DK702', 'DK709', 'DK703', 'DK713', 'DK714', 'DK715', 'DK717', 'DK73', 'DK74', 'DK760', 'DK762', 'DK763', 'DK764', 'DK769', 'DZ944'),
  severe.liver.disease = c('07000', '07002', '07004', '07006', '07008', '57300', '45601', '45602', '45603', '45604', '45605', '45606', '45607', '45608', '45609', 'DB150', 'DB160', 'DB162', 'DB190', 'DI850', 'DI859', 'DI864', 'DI982', 'DK704', 'DK711', 'DK721', 'DK729', 'DK765', 'DK766', 'DK767'), 
  diabetes.without.complications = c('24900', '24906', '24907', '24909', '25000', '25006', '25007', '25009', 'DE100', 'DE101', 'DE108', 'DE109', 'DE110', 'DE111', 'DE119', 'DE120', 'DE121', 'DE129', 'DE130', 'DE131', 'DE139', 'DE140', 'DE141', 'DE149'),
  diabetes.with.complications = c(paste0('2490', 1:5), '24908', paste0('2500', 1:5), '25008', paste0('DE10', 2:7), paste0('DE11', 2:8), paste0('DE12', 2:8), paste0('DE13', 2:8), paste0('DE14', 2:8)),
  hemiplegia.paraplegia = c('344', paste0('DG83', 0:4), 'DG81', 'DG82', 'DG041', 'DG114', 'DG801', 'DG802', 'DG839'),
  renal.disease = c('403', '404', paste0('58', 0:4), '59009', '59319', paste0('7531', 0:9), '792', paste0('DN03', 2:7), paste0('DN05', 2:7), 'DZ490', 'DZ491', 'DZ492', 'DN18', 'DN19', 'DI120', 'DI131', 'DI132', 'DN250', 'DZ940', 'DZ992', 'DN26'),
  any.malignancy = c(paste0('1', 40:72), paste0(174:194), '27559', paste0('DC', 0:3), paste0('DC4', 0:9), 'DC5', 'DC6', paste0('DC7', 0:6), 'DC86', 'DC97'),
  metastatic.solid.tumor = c(paste0('19', 5:9), paste0('DC', 77:80)),
  AIDS.HIV = c('07983', 'DB20', 'DB21', 'DB22', 'DB23', 'DB24'),
  leukemia = c(paste0('20', 4:7), paste0('DC9', 1:5)),
  lymphoma = c(paste0('20', 0:3), '27559', paste0('DC8', 1:5), 'DC88', 'DC90', 'DC96')
)

add_trailing_digits <- function(codes) {
  extended_codes <- unlist(lapply(codes, function(code) {
    c(code, paste0(code, 0:9))  # Keep the original code and append 0-9
  }))
  return(extended_codes)
}

charlson.codes.extended <- lapply(charlson.codes, add_trailing_digits)

# Function to filter out the codes starting with "D"
filter_codes_starting_with_D <- function(codes) {
  filtered_codes <- lapply(codes, function(code_list) {
    code_list[startsWith(code_list, "D")]
  })
  return(filtered_codes)
}

# Convert charlson.codes into a data.table for efficient filtering
charlson_codes_filtered <- as.data.table(unlist(filter_codes_starting_with_D(charlson.codes.extended)))

# Function to filter diagnoses based on filtered charlson codes
filter_diagnoses <- function(diagnoses) {
  relevant_codes <- charlson_codes_filtered$V1
  diagnoses[AKTIONSDIAGNOSE %in% relevant_codes]  # Filter using data.table syntax
}

# Function to split large PERSON_ID lists into chunks of 1000 (due to SQL limit)
split_ids_into_chunks <- function(ids, chunk_size = 1000) {
  split(ids, ceiling(seq_along(ids) / chunk_size))
}

# Function to query SQL database for each chunk of PERSON_IDs
query_sql_database_all_icd10 <- function(chunk, conn, sampling_dates) {
  # Convert chunk into a string for SQL query
  person_ids_list <- paste0("'", chunk, "'", collapse = ", ")
  
  # Query LPR3KONTAKTER
  query_LPR3 <- paste0("SELECT PERSON_ID, AKTIONSDIAGNOSE, DATO_START FROM D222008.LPR3KONTAKTER WHERE PERSON_ID IN (", person_ids_list, ")")
  diagnoses_LPR3 <- as.data.table(dbGetQuery(conn, query_LPR3))
  
  # Query LPRDIAG_PID
  query_LPR2 <- paste0("SELECT PERSON_ID, C_DIAG, AAR FROM D222008.LPRDIAG_PID WHERE PERSON_ID IN (", person_ids_list, ")")
  diagnoses_LPR2 <- as.data.table(dbGetQuery(conn, query_LPR2))
  
  # Standardize column names and format dates
  diagnoses_LPR2[, AKTIONSDIAGNOSE := C_DIAG]
  diagnoses_LPR2[, diagnosis_date := as.Date(paste0(AAR, "-01-01"))]
  diagnoses_LPR3[, diagnosis_date := as.Date(DATO_START)]
  
  # Combine both datasets
  diagnoses_combined <- rbindlist(list(
    diagnoses_LPR3[, .(PERSON_ID, AKTIONSDIAGNOSE, diagnosis_date)],
    diagnoses_LPR2[, .(PERSON_ID, AKTIONSDIAGNOSE, diagnosis_date)]
  ))
  
  # Filter the diagnoses using the filtered charlson codes
  diagnoses_filtered <- filter_diagnoses(diagnoses_combined)
  
  # Join the filtered diagnoses with the sampling_dates
  setkey(diagnoses_filtered, PERSON_ID)
  setkey(sampling_dates, PERSON_ID)
  
  # Perform the join
  diagnoses_final <- diagnoses_filtered[sampling_dates, nomatch = 0]
  
  return(diagnoses_final)
}

# Main processing function
process_person_ids_all_icd10 <- function(person_id_data, conn) {
  # Convert person_id_data to a data.table
  person_id_data_dt <- as.data.table(person_id_data)
  sampling_dates <- person_id_data_dt[, .(PERSON_ID, DATESAMPLING)]
  
  # Split PERSON_IDs into chunks of 1000
  person_id_chunks <- split_ids_into_chunks(person_id_data_dt$PERSON_ID)
  
  # Query the SQL database for each chunk and combine the results
  all_diagnoses <- rbindlist(lapply(person_id_chunks, query_sql_database_all_icd10, conn = conn, sampling_dates = sampling_dates))
  
  # Rename the column to 'code' as expected by the comorbidity function
  setnames(all_diagnoses, "AKTIONSDIAGNOSE", "code")
  
  return(all_diagnoses)
}
# Example of running the function with the provided data
clean_icd10_df_all <- process_person_ids_all_icd10(subset_cases_and_controls_all_with_metadata_geo_diag_limited, conn)
# Add pandemic_start_date to the clean_icd10_df_all dataset
clean_icd10_df_all$pandemic_start_date <- as.Date("2020-01-01")
subset_cases_and_controls_all_with_metadata_geo_diag_limited$pandemic_start_date <- as.Date("2020-01-01")

# Generate Charlson results for 5 years and 10 years
charlson_results_all_5yrs <- charlsonIndex(
  data = clean_icd10_df_all,
  ptid = "PERSON_ID",
  vars = "code",
  data.date = "diagnosis_date",
  charlson.date = "DATESAMPLING",
  look.back = 5,
  ccodes = charlson.codes.extended
)

charlson_results_all_10yrs <- charlsonIndex(
  data = clean_icd10_df_all,
  ptid = "PERSON_ID",
  vars = "code",
  data.date = "diagnosis_date",
  charlson.date = "DATESAMPLING",
  look.back = 10,
  ccodes = charlson.codes.extended
)

# Generate prepandemic results
charlson_results_all_5yrs_prepandemic <- charlsonIndex(
  data = clean_icd10_df_all,
  ptid = "PERSON_ID",
  vars = "code",
  data.date = "diagnosis_date",
  charlson.date = "pandemic_start_date",
  look.back = 5,
  ccodes = charlson.codes.extended
)

charlson_results_all_10yrs_prepandemic <- charlsonIndex(
  data = clean_icd10_df_all,
  ptid = "PERSON_ID",
  vars = "code",
  data.date = "diagnosis_date",
  charlson.date = "pandemic_start_date",
  look.back = 10,
  ccodes = charlson.codes.extended
)

# Check the structure of the prepandemic results to ensure 'pandemic_start_date' is included
print(str(charlson_results_all_5yrs_prepandemic[[1]]))
print(str(charlson_results_all_10yrs_prepandemic[[1]]))

# Merging the Charlson results with metadata for the first dataset (all cases)
subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson <- subset_cases_and_controls_all_with_metadata_geo_diag_limited %>%
  left_join(charlson_results_all_5yrs[[1]] %>% rename(charlson.index.5yrs = charlson.index), by = c("PERSON_ID", "DATESAMPLING")) %>%
  mutate(charlson.index.5yrs = ifelse(is.na(charlson.index.5yrs), 0, charlson.index.5yrs)) %>%
  
  left_join(charlson_results_all_10yrs[[1]] %>% rename(charlson.index.10yrs = charlson.index), by = c("PERSON_ID", "DATESAMPLING")) %>%
  mutate(charlson.index.10yrs = ifelse(is.na(charlson.index.10yrs), 0, charlson.index.10yrs)) %>%
  
  left_join(charlson_results_all_5yrs_prepandemic[[1]] %>% rename(charlson.index.5yrs_prepandemic = charlson.index), by = c("PERSON_ID", "pandemic_start_date")) %>%
  mutate(charlson.index.5yrs_prepandemic = ifelse(is.na(charlson.index.5yrs_prepandemic), 0, charlson.index.5yrs_prepandemic)) %>%
  
  left_join(charlson_results_all_10yrs_prepandemic[[1]] %>% rename(charlson.index.10yrs_prepandemic = charlson.index), by = c("PERSON_ID", "pandemic_start_date")) %>%
  mutate(charlson.index.10yrs_prepandemic = ifelse(is.na(charlson.index.10yrs_prepandemic), 0, charlson.index.10yrs_prepandemic))

conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson <- subset_cases_and_controls_all_with_metadata_geo_diag_limited_charlson

# Saving dataframes -----------
saveRDS(conditional_subset_cases_and_controls_all_with_metadata_geo_diag_charlson, file="")





