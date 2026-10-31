# Gene expression (Section 5.1): ADNI blood expression profile (Affymetrix U219), probe sets
# averaged within gene symbol; one value per subject and gene.
#
# Input: data/raw/ADNI_Gene_Expression_Profile.csv (as downloaded: 8 header rows of sample
#        information, then one row per probe set)
# Output: data/processed/ADNI_Gene_Expression.csv (PTID, gene.name, value) and
#         data/processed/ADNI_Gene_Expression_wide.csv
#   Rscript realdata/prep/R/02_clean_gene_expression.R
library(tidyverse)
library(readr)
library(here)
here::i_am("realdata/prep/R/02_clean_gene_expression.R")
raw_dir       <- here::here("data", "raw")
processed_dir <- here::here("data", "processed")
if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)
profile_csv <- file.path(raw_dir, "ADNI_Gene_Expression_Profile.csv")

# probe-level matrix without the sample-information rows and annotation columns
d1 <- read_csv(profile_csv, skip = 8)
d1 <- d1 %>% select(-ProbeSet, -LocusLink) %>% drop_na()
d1 <- d1[, -746]

# mean over the probe sets of each gene symbol
dd1 <- d1 %>%
  mutate(Symbol = as.factor(Symbol)) %>%
  group_by(Symbol) %>%
  summarise_all(mean)

# sample identifiers: row 1 = visit code, row 2 = PTID
ID <- read_csv(profile_csv)
ID <- tibble(PTID = unlist(ID[2, ]), visit = unlist(ID[1, ])) %>% drop_na()
ID <- ID[-1, ]

# one row per sample, one column per gene
d2 <- dd1[, -1]
d2 <- t(d2)
d2 <- as_tibble(d2)
colnames(d2) <- dd1[, 1]$Symbol
d3 <- tibble(ID, as_tibble(d2))
write_csv(d3, file.path(processed_dir, "ADNI_Gene_Expression_wide.csv"))

# long format: PTID, gene.name, value
d4 <- d3 %>%
  select(-visit) %>%
  gather(key = "gene.name", value = "value", -PTID) %>%
  drop_na()
write_csv(d4, file.path(processed_dir, "ADNI_Gene_Expression.csv"))
message("02_clean_gene_expression.R done -> data/processed/ADNI_Gene_Expression.csv")
