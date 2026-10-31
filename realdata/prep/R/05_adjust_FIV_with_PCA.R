# The functional instrument of every gene, adjusted for population structure (Section 5).
# data/processed/Genes_pca6/<gene>.csv: one column per variant assigned to the gene (header =
# position), one row per subject in .fam order; each column is the allele count of the variant
# (PLINK A1) minus its least-squares fit on PC1-PC6 (with intercept). Genes with >= 5 variants.
#
# Input: data/raw/filtered_ADNI_808_plink2.{bed,bim,fam},
#        data/processed/filtered_ADNI_808_plink2/{pca50_ADNI_808_plink2.csv, bim.csv, gene_snps.csv}
# Output: data/processed/Genes_pca6/ (about 56 GB)
#   Rscript realdata/prep/R/05_adjust_FIV_with_PCA.R          (FMR_CORES workers; default 8)
# Genes already written are skipped, so an interrupted run resumes.
suppressPackageStartupMessages({ library(readr); library(genio); library(dplyr); library(parallel); library(here) })
here::i_am("realdata/prep/R/05_adjust_FIV_with_PCA.R")
file <- "filtered_ADNI_808_plink2"
path <- here::here("data", "raw", file)
path2 <- here::here("data", "processed", file)
npca <- 6L
out_dir <- here::here("data", "processed", paste0("Genes_pca", npca))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
NC  <- as.integer(Sys.getenv("FMR_CORES", "8"))

bim <- read_bim(paste0(path, ".bim"))
fam <- read_fam(paste0(path, ".fam"))
df.pca <- read_csv(file.path(path2, "pca50_ADNI_808_plink2.csv"), show_col_types = FALSE)
stopifnot(identical(df.pca$PTID, fam$id))                    # PCs in .fam order
df.gene.snps <- read_csv(file.path(path2, "gene_snps.csv"), show_col_types = FALSE) %>% filter(n >= 5)
df.gene.snp  <- read_csv(file.path(path2, "bim.csv"), show_col_types = FALSE)
df.gene.snp$row <- match(df.gene.snp$id, bim$id)            # .bed row of each variant, by variant ID

# --- the .bed, one variant (row) at a time: variant-major, 2 bits per subject -------------
n  <- nrow(fam); nb <- ceiling(n / 4)
bed <- paste0(path, ".bed")
local({ con <- file(bed, "rb"); on.exit(close(con))
  stopifnot(identical(readBin(con, "raw", 3), as.raw(c(0x6c, 0x1b, 0x01)))) })
DECODE <- c(2, NA, 1, 0)                                      # 00 = A1/A1, 01 = missing, 10 = het, 11 = A2/A2
read_rows <- function(rows) {
  con <- file(bed, "rb"); on.exit(close(con))
  vapply(rows, function(r) {
    seek(con, 3 + (r - 1) * nb)
    b <- as.integer(readBin(con, "raw", nb))
    DECODE[bitwAnd(bitwShiftR(rep(b, each = 4), c(0L, 2L, 4L, 6L)), 3L)[seq_len(n)] + 1L]
  }, numeric(n))
}

run <- function(i) {
  gene <- df.gene.snps$gene.name[i]
  df.gene.info <- df.gene.snp %>% filter(gene.name == gene)
  Z <- read_rows(df.gene.info$row)
  Zt <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  for (j in 1:ncol(Z)) {
    fit <- lm(Z[, j] ~ as.matrix(df.pca[, 2:(npca + 1)]))
    Zt[, j] <- Z[, j] - predict(fit, df.pca[, 2:(npca + 1)])
  }
  colnames(Zt) <- df.gene.info$pos
  tmp <- file.path(out_dir, paste0(".", gene, ".csv.part"))     # renamed when complete
  write_csv(as.data.frame(Zt), tmp)
  file.rename(tmp, file.path(out_dir, paste0(gene, ".csv")))
  gene
}

todo <- which(!file.exists(file.path(out_dir, paste0(df.gene.snps$gene.name, ".csv"))))
message(sprintf("05_adjust_FIV_with_PCA: %d of %d genes to write, %d workers", length(todo), nrow(df.gene.snps), NC))
t0 <- Sys.time()
done <- unlist(mclapply(todo, function(i) tryCatch(run(i), error = function(e) NA_character_),
                        mc.cores = NC, mc.preschedule = FALSE))
message(sprintf("done: %d genes in %.1f min, %d failed", sum(!is.na(done)),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(is.na(done))))
