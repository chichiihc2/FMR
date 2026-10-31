# The 50 genotype PCs of 03_pca_plink.sh as a CSV (PTID, PC1, ..., PC50).
# Input: data/processed/filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.eigenvec
# Output: data/processed/filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.csv
#   Rscript realdata/prep/R/04_convert_pca_eigenvec.R
library(readr)
library(here)

here::i_am("realdata/prep/R/04_convert_pca_eigenvec.R")

raw_dir <- here::here("data", "raw")
out_dir <- here::here("data", "processed", "filtered_ADNI_808_plink2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file <- "pca50_ADNI_808_plink2"

# the "#FID IID ..." header line is skipped as a comment
df.pca <- read.table(file.path(out_dir, paste0(file, ".eigenvec")), header = FALSE)

# Drop FID; keep IID (renamed PTID) + the PC columns.
df.pca <- df.pca[, -1]
colnames(df.pca) <- c("PTID", paste0("PC", seq_len(ncol(df.pca) - 1)))

write_csv(df.pca, file.path(out_dir, paste0(file, ".csv")))

message("04_convert_pca_eigenvec.R done -> ",
        "data/processed/filtered_ADNI_808_plink2/", file, ".csv")
