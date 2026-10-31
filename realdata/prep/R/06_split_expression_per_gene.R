# Expression of each gene in its own file, data/processed/Gene_expression_per_gene/<gene>.csv
# (PTID, value), for the genes with >= 5 variants in gene_snps.csv. Gene symbols are matched
# exactly, else ignoring case (e.g. C10orf107 in the gene annotation, C10ORF107 in the
# expression data).
# Input: data/processed/ADNI_Gene_Expression.csv, data/processed/filtered_ADNI_808_plink2/gene_snps.csv
#   Rscript realdata/prep/R/06_split_expression_per_gene.R
suppressPackageStartupMessages({ library(data.table); library(here) })
here::i_am("realdata/prep/R/06_split_expression_per_gene.R")
expr_dir <- here::here("data", "processed", "Gene_expression_per_gene")
dir.create(expr_dir, showWarnings = FALSE, recursive = TRUE)
gs <- fread(here::here("data", "processed", "filtered_ADNI_808_plink2", "gene_snps.csv"))[n >= 5]
dge <- fread(here::here("data", "processed", "ADNI_Gene_Expression.csv"))
sym <- unique(dge$gene.name)
stopifnot(!anyDuplicated(tolower(sym)))
src <- ifelse(gs$gene.name %in% sym, gs$gene.name, sym[match(tolower(gs$gene.name), tolower(sym))])
setkey(dge, gene.name)
k <- 0L
for (j in which(!is.na(src))) {
  sub <- dge[J(src[j]), .(PTID, value), nomatch = NULL]
  if (nrow(sub) > 0L) { fwrite(sub, file.path(expr_dir, paste0(gs$gene.name[j], ".csv"))); k <- k + 1L }
}
message(sprintf("06_split_expression_per_gene: %d gene files (%d matched ignoring case) in %s",
                k, sum(!is.na(src) & src != gs$gene.name), expr_dir))
