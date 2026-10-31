#!/usr/bin/env Rscript
# =============================================================================
# run_realdata.R -- the ADNI analysis of Section 5: Table 2, the numbers of Sections 5.1 and
# 5.2, and Figure S.2.
#
#   1. Put the ADNI files in data/raw/ as described in data/README.md.
#   2. Rscript run_realdata.R
#
#   prep             genotype PCs, variant -> gene map, instruments of all genes, expression,
#                    eQTLGen (downloaded)                              -> data/processed/
#   gene_expression  Section 5.1                         -> realdata/gene_expression/results/
#   apoe             Section 5.2                                    -> realdata/apoe/results/
#   paper            Table 2, Figure S.2, numbers of the text                     -> output/
# The run ends by comparing the outputs with the paper's (reference/).
#
# A step whose outputs exist is skipped, so an interrupted run resumes; --force reruns all.
# Options: --cores=N (default: all), --only=prep|gene_expression|apoe|paper, --force.
# Needs plink2 and PLINK 1.9 (as `plink`) on the PATH, curl, and about 70 GB of disk.
# Run from the repository root.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
opt  <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  v <- sub(paste0("^--", key, "=?"), "", hit[1]); if (nzchar(v)) v else TRUE
}
CORES <- opt("cores", as.character(parallel::detectCores()))
ONLY  <- opt("only", "all")
FORCE <- isTRUE(opt("force"))
stopifnot(file.exists("paper/table2.R"))                       # run from the repository root
LOG <- "output/logs"; dir.create(LOG, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
part <- function(p) ONLY %in% c("all", p)

# ---- inputs ------------------------------------------------------------------------------
R0 <- "data/raw"
fam <- file.path(R0, "filtered_ADNI_808_plink2.fam")
vcf <- file.path(R0, paste0("ADNI.808_indiv.minGQ_21.pass_combined.vcf", c("", ".gz")))
if (!file.exists(fam) && any(file.exists(vcf))) {             # genotype QC from the ADNI VCF
  for (s in c("01_filter_vcf.sh", "02_qc_plink.sh")) {
    log_msg("run  genotype QC ", s)
    if (system2("bash", file.path("realdata/prep/wgs_qc", s), env = sprintf("VCF_IN=%s", vcf[file.exists(vcf)][1]),
                stdout = file.path(LOG, paste0(s, ".log")), stderr = file.path(LOG, paste0(s, ".log"))) != 0)
      stop(s, " failed; see ", LOG)
  }
}
if (!file.exists(file.path(R0, "PTID.csv")) && file.exists(fam))          # subject IDs in .fam order
  data.table::fwrite(data.frame(PTID = utils::read.table(fam)$V2), file.path(R0, "PTID.csv"))
pheno_in <- file.path(R0, c("ADNI_adnimerge_20170629_QT-freeze.csv", "ADNINIGHTINGALE_20210219_04Jan2024.csv", "Demo.csv"))
if (!file.exists(file.path(R0, "data_cleaned_demo.csv")) && all(file.exists(pheno_in))) {
  log_msg("run  phenotype table (realdata/prep/R/01_clean_phenotypes.R)")
  if (system2("Rscript", "realdata/prep/R/01_clean_phenotypes.R", stdout = file.path(LOG, "phenotypes.log"),
              stderr = file.path(LOG, "phenotypes.log")) != 0) stop("phenotype table failed; see ", LOG)
}
if (!dir.exists(file.path(R0, "ADNIMERGE2", "data"))) {                    # unpack ADNIMERGE2
  tgz <- list.files(R0, "^ADNIMERGE2.*\\.tar\\.gz$", full.names = TRUE)
  zip <- list.files(R0, "^adnimerge.*\\.zip$", full.names = TRUE, ignore.case = TRUE)
  if (!length(tgz) && length(zip)) { utils::unzip(zip[1], exdir = file.path(R0, "unz"))
    tgz <- list.files(file.path(R0, "unz"), "^ADNIMERGE2.*\\.tar\\.gz$", full.names = TRUE) }
  if (length(tgz)) { log_msg("unpacking ", tgz[1]); utils::untar(tgz[1], exdir = R0) }
}
need <- c(file.path(R0, paste0("filtered_ADNI_808_plink2.", c("bed", "bim", "fam"))),
          file.path(R0, c("PTID.csv", "data_cleaned_demo.csv", "ADNI_Gene_Expression_Profile.csv")),
          file.path(R0, "ADNIMERGE2", "data", paste0(c("ADAS", "ADSP_PHC_PET_Amyloid_Simple", "FAQ", "PTDEMOG", "VITALS"), ".rda")))
miss <- need[!file.exists(need)]
if (length(miss) && ONLY != "paper") stop("missing ADNI input(s):\n  ", paste(miss, collapse = "\n  "), "\nSee data/README.md.")
if (part("prep") && !nzchar(Sys.which("plink2"))) stop("plink2 must be on the PATH (www.cog-genomics.org/plink/2.0)")
if (part("apoe") && !nzchar(Sys.which("plink"))) stop("PLINK 1.9 must be on the PATH as `plink` (www.cog-genomics.org/plink/1.9)")

# data release: compare with the files the paper used (a mismatch is a warning)
ck <- read.table("data/checksums.md5", col.names = c("md5", "file"), stringsAsFactors = FALSE)
ck <- ck[file.exists(ck$file), ]
if (nrow(ck)) {
  bad <- ck$file[unname(tools::md5sum(ck$file)) != ck$md5]
  if (length(bad)) log_msg("WARNING: these inputs differ from the paper's (another ADNI release?): ",
                           paste(basename(bad), collapse = ", "), " -- the results may differ.")
  else log_msg("all ", nrow(ck), " input files found match the paper's")
}

# ---- steps: name, part, command, outputs ----------------------------------------------------
P <- "data/processed"; FA <- file.path(P, "filtered_ADNI_808_plink2")
GE <- "realdata/gene_expression/results"; AP <- "realdata/apoe/results"
S <- function(name, part, cmd, out) list(name = name, part = part, cmd = cmd, out = out)
STEPS <- list(
  S("genotype PCs",            "prep", c("bash", "realdata/prep/wgs_qc/03_pca_plink.sh"), file.path(FA, "pca50_ADNI_808_plink2.eigenvec")),
  S("PCs to csv",              "prep", c("Rscript", "realdata/prep/R/04_convert_pca_eigenvec.R"), file.path(FA, "pca50_ADNI_808_plink2.csv")),
  S("variant -> gene map",     "prep", c("Rscript", "realdata/prep/R/03_assign_snps_to_genes.R"), file.path(FA, c("bim.csv", "gene_snps.csv"))),
  S("instruments (Genes_pca6)", "prep", c("Rscript", "realdata/prep/R/05_adjust_FIV_with_PCA.R"), file.path(P, "Genes_pca6", "TTN.csv")),
  S("expression",              "prep", c("Rscript", "realdata/prep/R/02_clean_gene_expression.R"), file.path(P, "ADNI_Gene_Expression.csv")),
  S("expression per gene",     "prep", c("Rscript", "realdata/prep/R/06_split_expression_per_gene.R"), file.path(P, "Gene_expression_per_gene", "TTN.csv")),
  S("eQTLGen download",        "prep", c("bash", "realdata/prep/eqtlgen/00_download_eqtlgen.sh"),
    file.path("data/external/eqtlgen", "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")),
  S("eQTLGen table",           "prep", c("Rscript", "realdata/prep/eqtlgen/01_build_eqtlgen_gene_tables.R"), file.path(P, "eqtlgen", "eqtlgen_adni_cis.csv.gz")),
  S("FMR, all genes",          "gene_expression", c("Rscript", "realdata/gene_expression/fmr_genomewide.R"), file.path(GE, "fmr", "master.rds")),
  S("2SLS, cis-eQTLs",         "gene_expression", c("Rscript", "realdata/gene_expression/cis_2sls.R"), file.path(GE, "cis_2sls", "master.rds")),
  S("APOE genotypes",          "apoe", c("bash", "realdata/apoe/genotypes.sh"), file.path(P, "apoe", c("win.raw", "gpc.eigenvec", "rel.genome"))),
  S("APOE phenotypes",         "apoe", c("Rscript", "realdata/apoe/phenotypes.R"), file.path(P, "apoe", "pheno.csv")),
  S("APOE analysis",           "apoe", c("Rscript", "realdata/apoe/fmr_apoe.R"), file.path(AP, c("amy_cl_faq.rds", "bmi_adas13.rds"))))

for (s in STEPS) {
  if (!part(s$part)) next
  if (!FORCE && all(file.exists(s$out))) { log_msg("skip ", s$name, " (done)"); next }
  lf <- file.path(LOG, paste0(gsub("[^A-Za-z0-9]+", "_", s$name), ".log"))
  log_msg("run  ", s$name); t0 <- Sys.time()
  st <- system2(s$cmd[1], s$cmd[-1], env = sprintf("FMR_CORES=%s", CORES), stdout = lf, stderr = lf)
  if (st != 0 || !all(file.exists(s$out))) stop(s$name, " failed; see ", lf)
  log_msg("done ", s$name, sprintf(" in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- Table 2, Figure S.2, numbers of the text ------------------------------------------------
if (part("paper")) {
  for (p in c("table2", "section5_2", "figureS2")) {
    lf <- file.path(LOG, paste0(p, ".log"))
    if (system2("Rscript", c("--vanilla", file.path("paper", paste0(p, ".R"))), stdout = lf, stderr = lf) != 0)
      stop("paper/", p, ".R failed; see ", lf)
    log_msg("made ", p)
  }
  ok <- vapply(c("table2.tex", "section5_1.txt", "section5_2.txt"), function(f) {
    same <- identical(readLines(file.path("reference", f)), readLines(file.path("output", f)))
    cat(sprintf("  %-15s %s\n", f, if (same) "same as the paper" else "DIFFERS from the paper"))
    same
  }, TRUE)
  log_msg(if (all(ok)) "DONE: all results as in the paper" else "DONE: SOME RESULTS DIFFER")
  quit(status = if (all(ok)) 0 else 1)
}
