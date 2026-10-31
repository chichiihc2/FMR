# Phenotype table of the gene expression analysis (Section 5.1): baseline FDG, PTAU, ABETA,
# ADAS13 and MMSE from ADNIMERGE, with age and sex from the ADNI demographics.
#
# Input (data/raw/): ADNI_adnimerge_20170629_QT-freeze.csv, ADNINIGHTINGALE_20210219_04Jan2024.csv,
#                    Demo.csv
# Output: data/raw/data_cleaned_demo.csv (and data/processed/data_cleaned.csv)
#   Rscript realdata/prep/R/01_clean_phenotypes.R
library(tidyverse)
library(readr)
library(here)
here::i_am("realdata/prep/R/01_clean_phenotypes.R")
raw_dir       <- here::here("data", "raw")
processed_dir <- here::here("data", "processed")
if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)

# baseline outcomes
df1 <- read_csv(file.path(raw_dir, "ADNI_adnimerge_20170629_QT-freeze.csv"))
df1 <- df1 %>%
  filter(VISCODE == "bl") %>%
  select(RID, PTID, FDG, PTAU, ABETA, ADAS13, MMSE)
df.meta <- read_csv(file.path(raw_dir, "ADNINIGHTINGALE_20210219_04Jan2024.csv"))
df.meta <- df.meta %>%
  filter(VISCODE == "bl") %>%
  select(RID, XL.HDL.P, XL.HDL.FC, XL.HDL.PL, XL.HDL.L, GP)
df1 <- df1 %>% left_join(df.meta, by = "RID")
# ABETA values outside the assay range
df1$ABETA[df1$ABETA == ">1700"] <- "1700"
df1$ABETA[df1$ABETA == "<200"]  <- "100"
df1$ABETA <- as.numeric(df1$ABETA)
write_csv(df1, file.path(processed_dir, "data_cleaned.csv"))

# demographics: one row per subject
df <- read_csv(file.path(processed_dir, "data_cleaned.csv"))
df.demo <- read_csv(file.path(raw_dir, "Demo.csv"))
df.demo <- df.demo %>%
  rename(PTID = subject_id) %>%
  filter(!is.na(research_group))
df.demo <- df.demo %>%
  group_by(PTID) %>%
  summarise(
    Age            = mean(subject_age, rm.na = T),
    Gender         = mean(PTGENDER),
    research_group = paste0(unique(research_group), collapse = ",")
  ) %>%
  ungroup()
df <- df %>% left_join(df.demo, by = "PTID")
df <- df %>%
  mutate(grouped_category = case_when(
    research_group %in% c("CN", "SMC")          ~ "CN",
    research_group %in% c("EMCI", "MCI", "LMCI") ~ "MCI",
    research_group == "AD"                       ~ "AD"
  ))
write_csv(df, file.path(raw_dir, "data_cleaned_demo.csv"))
message("01_clean_phenotypes.R done -> data/raw/data_cleaned_demo.csv")
