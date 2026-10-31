# Assignment of variants to genes (Section 5): a variant belongs to a gene if its position is
# within the gene's GRCh37 coordinates (the first gene listed in gene_locations37.csv if
# several). At positions carrying several variants, only those with an rs ID are kept.
#
# Input: data/raw/filtered_ADNI_808_plink2.bim, data/external/gene_locations37.csv
# Output: data/processed/filtered_ADNI_808_plink2/bim.csv (chr, id, posg, pos, alt, ref,
#         gene.name, geneID) and gene_snps.csv (gene.name, number of variants)
#   Rscript realdata/prep/R/03_assign_snps_to_genes.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(genio)
  library(data.table)
  library(here)
})

here::i_am("realdata/prep/R/03_assign_snps_to_genes.R")

raw_dir      <- here::here("data", "raw")
external_dir <- here::here("data", "external")
out_dir      <- here::here("data", "processed", "filtered_ADNI_808_plink2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

t0 <- Sys.time()

# --- inputs -----------------------------------------------------------------
bim <- read_bim(file.path(raw_dir, "filtered_ADNI_808_plink2.bim"))
DT_snp <- as.data.table(bim)                 # cols: chr, id, posg, pos, ref, alt
DT_snp[, snp_row := .I]
DT_snp[, `:=`(start = pos, end = pos)]       # SNP "interval" is a point

gene_locations <- fread(file.path(external_dir, "gene_locations37.csv"))
# cols: GeneID, chr, start, end, p, gene.name
gene_locations[, gl_row := .I]               # preserve original row order
gene_locations[, chr := as.character(chr)]
DT_snp[, chr := as.character(chr)]

cat(sprintf("  read inputs: %d SNPs, %d genes (%.1fs)\n",
            nrow(DT_snp), nrow(gene_locations),
            as.numeric(Sys.time() - t0, units = "secs")))

# --- per-SNP gene annotation via vectorized interval overlap ----------------
# foverlaps requires both tables to be keyed on (chr, start, end).
setkey(gene_locations, chr, start, end)
ov <- foverlaps(DT_snp[, .(snp_row, chr, start, end)],
                gene_locations[, .(chr, start, end, GeneID, gene.name, gl_row)],
                type   = "any",
                nomatch = NULL)
# Per variant, the first overlapping gene in gene_locations37.csv order.
setorder(ov, snp_row, gl_row)
ov_first <- ov[, .SD[1L], by = snp_row, .SDcols = c("GeneID", "gene.name")]

cat(sprintf("  overlap join: %d SNP-gene pairs, %d SNPs assigned (%.1fs)\n",
            nrow(ov), nrow(ov_first),
            as.numeric(Sys.time() - t0, units = "secs")))

# Attach gene.name / geneID back to the SNP table. Use typed NAs so that
# data.table's reference assignment does not coerce the character / integer
# values to logical when filling them in.
DT_snp[, gene.name := NA_character_]
DT_snp[, geneID    := NA_integer_]
DT_snp[ov_first, on = "snp_row", `:=`(gene.name = i.gene.name,
                                      geneID    = i.GeneID)]

# --- drop SNPs with no gene assignment --------------------------------------
df2 <- DT_snp[!is.na(gene.name)]

# --- de-duplicate multi-SNP positions --------------------------------------
# Where several variants share a (chr, pos), keep only those with an rsID
# (id starting with "rs"); WGS variants without an rsID at such positions are
# dropped. Positions with a single variant are untouched.
multi_pos <- df2[, .N, by = .(chr, pos)][N > 1]
bad_ids   <- df2[multi_pos, on = c("chr", "pos")][!grepl("^rs", id), id]
df2 <- df2[!id %in% bad_ids]

cat(sprintf("  after dedup: %d SNPs retained (%.1fs)\n",
            nrow(df2),
            as.numeric(Sys.time() - t0, units = "secs")))

# --- write per-SNP table ----------------------------------------------------
df2_out <- df2[, .(chr, id, posg, pos, alt, ref, gene.name, geneID)]
fwrite(df2_out, file.path(out_dir, "bim.csv"))

# --- write gene -> SNP-count table (unfiltered) -----------------------------
df.gene.snps <- df2_out[, .(n = .N), by = gene.name][order(gene.name)]
fwrite(df.gene.snps, file.path(out_dir, "gene_snps.csv"))

cat(sprintf("03_assign_snps_to_genes.R done in %.1fs -> ",
            as.numeric(Sys.time() - t0, units = "secs")),
    "data/processed/filtered_ADNI_808_plink2/{bim.csv, gene_snps.csv}\n",
    sep = "")
