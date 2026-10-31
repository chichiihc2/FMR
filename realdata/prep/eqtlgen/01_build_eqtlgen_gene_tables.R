# eQTLGen cis-eQTL z-scores of the variants assigned to each gene in ADNI, harmonized to the
# allele counted in the ADNI genotypes (PLINK A1); strand-ambiguous SNPs are dropped. Used to
# rank the cis-eQTL instruments of the conventional 2SLS (realdata/gene_expression/cis_2sls.R).
# Both data sets are on GRCh37.
# Input: data/external/eqtlgen/ (00_download_eqtlgen.sh), data/processed/filtered_ADNI_808_plink2/bim.csv
# Output: data/processed/eqtlgen/eqtlgen_adni_cis.csv.gz
#   Rscript realdata/prep/eqtlgen/01_build_eqtlgen_gene_tables.R

library(here)
here::i_am("realdata/prep/eqtlgen/01_build_eqtlgen_gene_tables.R")

library(data.table)

ext_dir <- here::here("data", "external", "eqtlgen")
out_dir <- here::here("data", "processed", "eqtlgen")
res_dir <- here::here("data", "processed", "eqtlgen")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

cis_file <- file.path(
  ext_dir, "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
af_file <- file.path(
  ext_dir, "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")

for (f in c(cis_file, af_file)) {
  if (!file.exists(f))
    stop("Missing ", basename(f), ". Run 00_download_eqtlgen.sh first.")
}

msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                     paste0(...)))

# --- ADNI SNP-gene map --------------------------------------------------------
msg("Reading ADNI bim.csv ...")
bim <- fread(here::here("data", "processed", "filtered_ADNI_808_plink2", "bim.csv"))
# chr, id, posg, pos, alt, ref, gene.name, geneID; `alt` is PLINK A1, the allele counted
setnames(bim, "gene.name", "gene")
bim <- bim[!is.na(chr) & !is.na(pos) & gene != ""]
msg(sprintf("  %s SNP-gene rows, %s genes, %s SNPs without rsID",
            format(nrow(bim), big.mark = ","),
            format(uniqueN(bim$gene), big.mark = ","),
            format(sum(bim$id == "."), big.mark = ",")))

# --- Stream-filter the multi-GB cis file --------------------------------------
# The full cis file is ~10^8 rows. Never read it whole. Prefilter with awk on a
# hash of ADNI's (gene, chr, pos) keys, so at most nrow(bim) rows survive.
hdr <- names(fread(cmd = sprintf("gzip -dc %s | head -1", shQuote(cis_file))))
idx <- function(nm) {
  j <- match(nm, hdr)
  if (is.na(j)) stop("Column '", nm, "' not in cis file header: ",
                     paste(hdr, collapse = ", "))
  j
}
j_sym <- idx("GeneSymbol"); j_chr <- idx("SNPChr"); j_pos <- idx("SNPPos")

key_file <- tempfile(fileext = ".txt")
fwrite(unique(bim[, .(gene, chr, pos)]), key_file, sep = "\t", col.names = FALSE)
on.exit(unlink(key_file), add = TRUE)

msg("Stream-filtering cis file against ADNI keys (this is the slow step) ...")
# awk reads the key file first (NR==FNR), then the decompressed stream as "-".
# Piping rather than using <(...) keeps this working under /bin/sh, which is what
# fread(cmd=) invokes -- process substitution is a bash-ism and is not available there.
awk_prog <- sprintf(
  'NR==FNR{k[$1 FS $2 FS $3]; next} FNR==1{print; next} ($%d FS $%d FS $%d) in k',
  j_sym, j_chr, j_pos)
cis <- fread(cmd = sprintf("gzip -dc %s | awk -F'\\t' %s %s -",
                           shQuote(cis_file), shQuote(awk_prog), shQuote(key_file)),
             sep = "\t")
msg(sprintf("  %s eQTLGen rows matched ADNI (gene, chr, pos)",
            format(nrow(cis), big.mark = ",")))
if (nrow(cis) == 0)
  stop("Zero rows matched. Check genome build and gene-symbol conventions.")

setnames(cis, "GeneSymbol", "gene")
cis <- cis[, .(gene, chr = SNPChr, pos = SNPPos, rsid = SNP,
               assessed = toupper(AssessedAllele), other = toupper(OtherAllele),
               z = Zscore, n_ext = NrSamples, ensg = Gene)]

# --- Allele frequency (needed for Z -> beta) ----------------------------------
msg("Reading allele-frequency file ...")
af <- fread(cmd = sprintf("gzip -dc %s", shQuote(af_file)))
af <- af[, .(rsid = SNP, af_alleleB = AlleleB_all,
             alleleA = toupper(AlleleA), alleleB = toupper(AlleleB))]
af <- unique(af, by = "rsid")

cis <- merge(cis, af, by = "rsid", all.x = FALSE)
msg(sprintf("  %s rows have allele frequency", format(nrow(cis), big.mark = ",")))

# AlleleB is the assessed allele in the AF file, but align explicitly rather
# than assume the two files agree row-for-row.
cis[, af_assessed := fifelse(alleleB == assessed, af_alleleB,
                     fifelse(alleleA == assessed, 1 - af_alleleB, NA_real_))]
n_af_bad <- sum(is.na(cis$af_assessed))
cis <- cis[!is.na(af_assessed)]
msg(sprintf("  dropped %s rows where AF alleles disagreed with the cis file",
            format(n_af_bad, big.mark = ",")))

# --- Join to ADNI and harmonize to ADNI's effect allele -----------------------
adni <- bim[, .(gene, chr, pos, rsid_adni = id, alt = toupper(alt), ref = toupper(ref))]
dt <- merge(cis, adni, by = c("gene", "chr", "pos"), all.x = FALSE)

# Strand-ambiguous SNPs cannot be harmonized without strand information. Drop
# rather than guess -- a wrong guess silently inverts the instrument.
amb <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") |
                      (a == "C" & b == "G") | (a == "G" & b == "C")
n_amb <- sum(amb(dt$alt, dt$ref))
dt <- dt[!amb(alt, ref)]
msg(sprintf("  dropped %s strand-ambiguous (A/T, G/C) SNPs",
            format(n_amb, big.mark = ",")))

# Flip Z (and AF) when eQTLGen's assessed allele is ADNI's ref rather than alt,
# so every z below refers to the SAME allele ADNI's dosages count.
dt[, flip := fcase(assessed == alt & other == ref, FALSE,
                   assessed == ref & other == alt, TRUE,
                   default = NA)]
n_allele_mismatch <- sum(is.na(dt$flip))
dt <- dt[!is.na(flip)]
dt[flip == TRUE, `:=`(z = -z, af_assessed = 1 - af_assessed)]
msg(sprintf("  flipped %s SNPs to ADNI's alt allele; dropped %s with allele mismatch",
            format(sum(dt$flip), big.mark = ","),
            format(n_allele_mismatch, big.mark = ",")))

# --- Z -> beta (Zhu et al. 2016, SMR), expression on a standardized scale ------
# beta = z / sqrt(2p(1-p)(N + z^2)),  se = 1 / sqrt(2p(1-p)(N + z^2))
dt[, denom := sqrt(2 * af_assessed * (1 - af_assessed) * (n_ext + z^2))]
dt <- dt[is.finite(denom) & denom > 0]
dt[, `:=`(beta = z / denom, se = 1 / denom)]
dt[, pval := 2 * pnorm(-abs(z))]
dt[, denom := NULL]

setorder(dt, gene, chr, pos)
out <- dt[, .(gene, ensg, chr, pos, rsid, rsid_adni, alt, ref,
              af = af_assessed, z, beta, se, pval, n_ext, flipped = flip)]

out_file <- file.path(out_dir, "eqtlgen_adni_cis.csv.gz")
fwrite(out, out_file)
msg("Wrote ", out_file)

# --- Coverage diagnostics -----------------------------------------------------
cov_by_gene <- merge(
  bim[, .(n_snps_adni = .N), by = gene],
  out[, .(n_snps_matched = .N, max_abs_z = max(abs(z))), by = gene],
  by = "gene", all.x = TRUE)
cov_by_gene[is.na(n_snps_matched), `:=`(n_snps_matched = 0L, max_abs_z = NA_real_)]
cov_by_gene[, frac_matched := n_snps_matched / n_snps_adni]
fwrite(cov_by_gene, file.path(res_dir, "eqtlgen_join_coverage.csv"))

cat("\n=== Join summary ===\n")
cat(sprintf("ADNI SNP-gene rows          : %s\n", format(nrow(bim), big.mark = ",")))
cat(sprintf("Matched to eQTLGen          : %s (%.1f%%)\n",
            format(nrow(out), big.mark = ","), 100 * nrow(out) / nrow(bim)))
cat(sprintf("ADNI genes                  : %s\n", format(uniqueN(bim$gene), big.mark = ",")))
cat(sprintf("Genes with >=1 matched SNP  : %s\n",
            format(sum(cov_by_gene$n_snps_matched > 0), big.mark = ",")))
cat(sprintf("Genes with ZERO matched SNPs: %s  (gene-symbol alias or not expressed in blood)\n",
            format(sum(cov_by_gene$n_snps_matched == 0), big.mark = ",")))

# check: the APOE/TOMM40/APOC1 locus has strong cis-eQTLs in blood
cat("\n=== Positive control (chr19 APOE locus) ===\n")
print(cov_by_gene[gene %chin% c("APOE", "TOMM40", "APOC1")])
cat("\nExpect large max_abs_z (tens to hundreds). If these are small or NA,\n")
cat("STOP -- the join, the genome build, or the gene symbols are wrong.\n")
