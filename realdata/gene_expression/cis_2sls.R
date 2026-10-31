# Section 5.1, conventional 2SLS of Table 2: for each gene, its variants that are cis-eQTLs in
# eQTLGen (whole blood) ranked by |z|, LD-clumped at r^2 < 0.5; the top k (k as for FMR) are
# the instruments. Y, X and the variants are adjusted for age and sex.
#
# Usage: Rscript realdata/gene_expression/cis_2sls.R      (after fmr_genomewide.R)
# Output: realdata/gene_expression/results/cis_2sls/master.rds
suppressPackageStartupMessages({ library(here); library(data.table); library(parallel) })
here::i_am("realdata/gene_expression/cis_2sls.R")
source(here::here("lib", "adni_gene_data.R"))

OUTCOMES <- c("ADAS13", "MMSE", "PTAU", "FDG", "ABETA")
MIN_SNPS <- 10L; R2 <- 0.5
NCORES <- as.integer(Sys.getenv("FMR_CORES", max(1L, min(20L, detectCores() - 4L))))
out_dir <- here::here("realdata", "gene_expression", "results", "cis_2sls")
dir.create(file.path(out_dir, "genes"), recursive = TRUE, showWarnings = FALSE)

demo <- fread(here::here("data", "raw", "data_cleaned_demo.csv"))[!duplicated(PTID)]
ptid <- fread(here::here("data", "raw", "PTID.csv"))
expr_dir <- here::here("data", "processed", "Gene_expression_per_gene")
demo_out <- as.data.frame(demo)[, c("PTID", OUTCOMES)]
ext <- fread(here::here("data", "processed", "eqtlgen", "eqtlgen_adni_cis.csv.gz"),
             select = c("gene", "pos", "rsid", "z", "pval", "n_ext"))
ext <- ext[order(gene, pos, -abs(z))][!duplicated(ext[order(gene, pos, -abs(z)), .(gene, pos)])]
fmr <- readRDS(here::here("realdata", "gene_expression", "results", "fmr", "master.rds"))[status == "ok", .(gene, outcome, k)]
genes <- intersect(unique(fmr$gene), unique(ext$gene))

# Greedy clumping: the first k candidates (in order of |z|) with r^2 < r2max to those kept.
clump <- function(Z, cand, k, r2max) {
  sel <- integer(0)
  for (j in cand) {
    if (stats::sd(Z[, j]) == 0) next
    if (!length(sel) || max(drop(cor(Z[, sel, drop = FALSE], Z[, j]))^2) < r2max) sel <- c(sel, j)
    if (length(sel) == k) break
  }
  sel
}
tsls <- function(Y, X, Z) {
  Zq <- qr.Q(qr(Z)); Xh <- drop(Zq %*% crossprod(Zq, X))
  b  <- sum(Xh * Y) / sum(Xh * X)
  n  <- length(Y); s2 <- sum((Y - X * b)^2) / (n - 1)
  se <- sqrt(s2 / sum(Xh * X))
  r2 <- sum(Xh^2) / sum(X^2); q <- ncol(Z)
  c(beta = b, se = se, p = 2 * pt(-abs(b / se), n - 1), k_cis = q, F = (r2 / q) / ((1 - r2) / (n - q - 3L)))
}

run_gene <- function(g) {
  gd <- tryCatch(load_gene_frame(g, demo, ptid, expr_dir, min_snps = MIN_SNPS, min_obs = 50L), error = function(e) NULL)
  if (is.null(gd)) return(data.table(gene = g, status = "no_data"))
  eg <- ext[gene == g]
  mi <- match(gd$POS, eg$pos)
  cand <- which(!is.na(mi))
  if (!length(cand)) return(data.table(gene = g, status = "no_cis_eqtl"))
  cand <- cand[order(-abs(eg$z[mi[cand]]))]
  cells <- fmr[gene == g]; rows <- list()
  for (r in seq_len(nrow(cells))) {
    y <- demo_out[[cells$outcome[r]]][match(gd$df0$PTID, demo_out$PTID)]
    keep <- which(is.finite(y) & y > 0)
    d  <- gd$df0[keep, , drop = FALSE]
    cv <- stats::model.matrix(~ d$Age + factor(d$Gender))
    Y  <- drop(stats::resid(stats::lm(scale(log(y[keep]), TRUE, TRUE) ~ cv - 1)))
    X  <- drop(stats::resid(stats::lm(gd$Xexp[keep] ~ cv - 1)))
    Z  <- stats::resid(stats::lm(gd$Z[keep, , drop = FALSE] ~ cv - 1))
    s  <- clump(Z, cand, cells$k[r], R2)
    if (!length(s)) next
    rows[[r]] <- as.data.table(c(list(gene = g, outcome = cells$outcome[r], n_obs = length(keep)),
                                 as.list(tsls(Y, X, Z[, s, drop = FALSE]))))
  }
  fwrite(rbindlist(rows), file.path(out_dir, "genes", paste0(g, ".csv")))
  data.table(gene = g, status = "ok")
}

message(sprintf("%d genes on %d cores", length(genes), NCORES)); t0 <- Sys.time()
st <- rbindlist(mclapply(genes, function(g)
  tryCatch(run_gene(g), error = function(e) data.table(gene = g, status = paste("error:", conditionMessage(e)))),
  mc.cores = NCORES, mc.preschedule = FALSE), fill = TRUE)
message(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
fs <- file.path(out_dir, "genes", paste0(st[status == "ok"]$gene, ".csv"))
saveRDS(rbindlist(lapply(fs, fread), fill = TRUE), file.path(out_dir, "master.rds"))
print(st[, .N, by = status])
