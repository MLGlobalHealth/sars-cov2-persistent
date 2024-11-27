# First pre-processing step to identify persistent infections

# read packages
library(ape);library(pracma);library(phytools);library(TreeTools);library(dplyr);library(tidyverse)
library(data.table);library(dbplyr);library(lubridate);library(rlang)
library(foreach);library(doParallel);library(DSTora);library(ROracle)
library(DSTcolectica);library(DSTdb);library(DBI);library(parallel);
library(foreach);library(ggsignif);library(Rcpp);library(ggplot2)
library(data.table);library(stringr);library(seqinr);library(parallel)
library(foreach);library(future);library(future.apply)

#Establish connection
drv <- dbDriver('Oracle')
conn <- DSTora::OraGenvej('', dbuser = '')
con2 <- DSTora::OraGenvej('', dbuser = '')

COVID_TEST <- dbReadTable(conn = conn,
                          name = "",
                          schema = '')

lifelines <- dbReadTable(conn = conn,
                         name = "",
                         schema = '')

lifelines_koen <- dbReadTable(conn = conn,
                              name = "",
                              schema = '')

sequence_metadata <- dbReadTable(conn = conn,
                                 name = "",
                                 schema = '')

sequence_metadata_v2 <- dbReadTable(conn = conn,
                                 name = "",
                                 schema = '')


# Cleaning up metadata ---------
# Getting all those with a person ID, and keeping only those that passed QC:
metadata_passed_QC_valid_ID <- subset(sequence_metadata, !is.na(PERSON_ID) & GENOMESTATUS == "Genome") # 738944 rows

metadata_multiple_rows <- metadata_passed_QC_valid_ID %>%
  group_by(PERSON_ID) %>%
  filter(n() > 1) %>%
  ungroup()

metadata_processed_missing_CT <- metadata_multiple_rows %>%
  group_by(PERSON_ID, DATESAMPLING) %>%
  mutate(same_CT_same_sequencestatus = ifelse(n() > 1 & all(CT == first(CT)) & all(SEQUENCESTATUS == "Main sequence"), TRUE, FALSE)) %>%
  mutate(keep_row = case_when(
    n() == 1 ~ TRUE,
    same_CT_same_sequencestatus & SEQUENCESTATUS == "Main sequence" ~ TRUE,
    TRUE ~ row_number() == 1
  )) %>%
  filter(keep_row) %>%
  ungroup() %>%
  select(-same_CT_same_sequencestatus, -keep_row)

updated_metadata_processed_multiple_rows <- metadata_processed_missing_CT %>%
  group_by(PERSON_ID) %>%
  filter(n() > 1) %>%
  ungroup() #81326 rows, 36601 unique individuals

write.csv(updated_metadata_processed_multiple_rows, file="")








