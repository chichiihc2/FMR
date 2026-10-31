# Section 5.1: genome-wide MR of blood gene expression on five AD outcomes (Table 2).
#
# For each gene (at least 10 variants, with expression data) and outcome (ADAS13, MMSE, PTAU,
# FDG, ABETA; log-transformed and standardized):
#   Y, X  outcome and expression, adjusted for age and sex;
#   FPCA of the gene's genotype function; k by first-stage BIC over 3..10 (pairs with fewer
#         than three FPCA scores are not analysed);
#   first stage: partial R^2 and F of X on the k scores;
#   FMR: lambda by GCV, beta-hat, sandwich SE, bootstrap-t p-value (B = 199);
#   T_alpha test of H0: alpha = 0;
#   F2SLS: 2SLS with the k scores as instruments.
#
# Usage: Rscript realdata/gene_expression/fmr_genomewide.R [genes=A,B,...]
# Output: realdata/gene_expression/results/fmr/{master.rds, gene_status.csv}
suppressPackageStartupMessages({
  library(here); library(data.table); library(parallel); library(fda); library(refund); library(MASS)
})
here::i_am("realdata/gene_expression/fmr_genomewide.R")
source(here::here("lib", "estimators.R"))
source(here::here("lib", "adni_gene_data.R"))
source(here::here("lib", "gene_fpca.R"))
source(here::here("lib", "fmr_talpha.R"))

OUTCOMES <- c("ADAS13", "MMSE", "PTAU", "FDG", "ABETA")
MIN_SNPS <- 10L; NBOOT <- 199L; TOL <- 1e-15; LGRID <- 2^(-10:10)
NCORES <- as.integer(Sys.getenv("FMR_CORES", max(1L, min(20L, detectCores() - 4L))))
args <- commandArgs(trailingOnly = TRUE)
GSEL <- sub("^genes=", "", grep("^genes=", args, value = TRUE))

out_dir <- here::here("realdata", "gene_expression", "results", "fmr")
dir.create(file.path(out_dir, "genes"), recursive = TRUE, showWarnings = FALSE)
pen  <- make_pen(30L)
demo <- fread(here::here("data", "raw", "data_cleaned_demo.csv"))[!duplicated(PTID)]
ptid <- fread(here::here("data", "raw", "PTID.csv"))
expr_dir <- here::here("data", "processed", "Gene_expression_per_gene")
demo_out <- as.data.frame(demo)[, c("PTID", OUTCOMES)]

gs <- fread(here::here("data", "processed", "filtered_ADNI_808_plink2", "gene_snps.csv"))[n >= MIN_SNPS]
genes <- gs$gene.name[file.exists(file.path(expr_dir, paste0(gs$gene.name, ".csv"))) &
                      file.exists(genes_pca6_path(gs$gene.name))]
seed_of <- function(g, oi) 7000000L + match(g, genes) * 10L + oi      # bootstrap seed of a pair
run_genes <- if (length(GSEL)) intersect(genes, strsplit(GSEL, ",")[[1]]) else genes

# FMR at lambda_GCV with the bootstrap-t p-value, and T_alpha. Y, X, Xi centred.
fmr_pair <- function(Y, X, Xi, M2, seed) {
  n <- length(Y); k <- ncol(Xi)
  A <- Xi %*% M2
  cp0 <- xprod(Y, cbind(X, A), cbind(1, Xi), TOL)
  lam <- gcv_lambda(cp0, pen$PP_rough, LGRID, k)
  if (!is.finite(lam)) return(list(status = "no_solution"))
  f0 <- se_sandwich(cp0, lam, pen$PP_rough)
  if (is.null(f0) || !is.finite(f0$se) || f0$se <= 0) return(list(status = "no_se"))
  set.seed(seed)
  ts <- vapply(seq_len(NBOOT), function(b) {
    id  <- sample.int(n, n, replace = TRUE)
    cpb <- tryCatch(xprod(Y[id], cbind(X[id], A[id, , drop = FALSE]), cbind(1, Xi[id, , drop = FALSE]), TOL),
                    error = function(e) NULL)
    if (is.null(cpb)) return(NA_real_)
    lb <- tryCatch(gcv_lambda(cpb, pen$PP_rough, LGRID, k), error = function(e) NA_real_)
    if (!is.finite(lb)) lb <- lam
    sb <- tryCatch(se_sandwich(cpb, lb, pen$PP_rough), error = function(e) NULL)
    if (!is.null(sb) && is.finite(sb$se) && sb$se > 0) (sb$beta - f0$beta) / sb$se else NA_real_
  }, numeric(1))
  ts <- ts[is.finite(ts)]
  if (length(ts) < 50L) return(list(status = "too_few_draws"))
  c_hat <- unname(stats::quantile(abs(ts), 0.95))
  ta <- fmr_talpha(Y, X, Xi, A, lam, pen$PP_rough, pen$BO, tau = 0.1, tol = TOL)
  list(status = "ok", lambda = lam, beta = f0$beta, se = f0$se, c_hat = c_hat,
       p = 2 * pnorm(-qnorm(0.975) * abs(f0$beta) / (c_hat * f0$se)),      # bootstrap-t p-value
       T_alpha = ta[["T"]], p_talpha = ta[["p"]])
}

f2sls <- function(Y, X, Xi) tryCatch({
  sm <- summary(ivreg::ivreg(Y ~ -1 + X | Xi, data = tibble::tibble(Y = as.numeric(Y), X = as.numeric(X),
                                                                   Xi = unname(as.matrix(Xi)))))$coefficients
  c(beta = sm["X", 1], se = sm["X", 2], p = sm["X", 4])
}, error = function(e) c(beta = NA_real_, se = NA_real_, p = NA_real_))

run_gene <- function(g) {
  gd <- tryCatch(load_gene_frame(g, demo, ptid, expr_dir, min_snps = MIN_SNPS, min_obs = 50L),
                 error = function(e) NULL)
  if (is.null(gd)) return(data.table(gene = g, status = "no_data"))
  fa <- gene_fpca(gd$Z, gd$t01)
  if (is.null(fa)) return(data.table(gene = g, status = "fpca_failed"))
  wq <- gd$t01 - dplyr::lag(gd$t01); wq[1] <- 0
  Psiw <- fda::eval.basis(gd$t01, pen$BO) * wq
  rows <- list()
  for (oi in seq_along(OUTCOMES)) {
    y <- demo_out[[OUTCOMES[oi]]][match(gd$df0$PTID, demo_out$PTID)]
    keep <- which(is.finite(y) & y > 0)
    if (length(keep) < 50L) next
    d  <- gd$df0[keep, , drop = FALSE]; n <- length(keep)
    cv <- stats::model.matrix(~ d$Age + factor(d$Gender))
    Y  <- stats::resid(stats::lm(scale(log(y[keep]), TRUE, TRUE) ~ cv - 1))
    X  <- stats::resid(stats::lm(gd$Xexp[keep] ~ cv - 1))
    S  <- matrix(stats::resid(stats::lm(fa$scores[keep, , drop = FALSE] ~ cv - 1)), nrow = n)
    k  <- select_k_bic(X, S)
    r2 <- summary(stats::lm(X ~ S[, seq_len(k), drop = FALSE]))$r.squared
    Xi <- fa$scores[keep, seq_len(k), drop = FALSE]
    res <- if (k < 3L) list(status = "fewer_than_3_scores") else
      tryCatch(fmr_pair(drop(scale(Y, TRUE, FALSE)), drop(scale(X, TRUE, FALSE)), scale(Xi, TRUE, FALSE),
                        t(fa$efunctions[, seq_len(k), drop = FALSE]) %*% Psiw, seed_of(g, oi)),
               error = function(e) list(status = paste("error:", conditionMessage(e))))
    f2 <- f2sls(Y, X, Xi)
    rows[[oi]] <- as.data.table(c(list(gene = g, outcome = OUTCOMES[oi], n_obs = n, n_snp = ncol(gd$Z),
      K_pool = fa$K_pool, k = k, partial_R2 = r2, F = (r2 / k) / ((1 - r2) / (n - k - 3L)),
      f2sls_beta = f2[["beta"]], f2sls_se = f2[["se"]], f2sls_p = f2[["p"]]), res))
  }
  fwrite(rbindlist(rows, fill = TRUE), file.path(out_dir, "genes", paste0(g, ".csv")))
  data.table(gene = g, status = "ok")
}

message(sprintf("%d genes x %d outcomes on %d cores", length(run_genes), length(OUTCOMES), NCORES))
t0 <- Sys.time()
st <- rbindlist(mclapply(run_genes, function(g)
  tryCatch(run_gene(g), error = function(e) data.table(gene = g, status = paste("error:", conditionMessage(e)))),
  mc.cores = NCORES, mc.preschedule = FALSE), fill = TRUE)
message(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (!length(GSEL)) {
  fwrite(st, file.path(out_dir, "gene_status.csv"))
  fs <- file.path(out_dir, "genes", paste0(st[status == "ok"]$gene, ".csv"))
  saveRDS(rbindlist(lapply(fs, fread), fill = TRUE), file.path(out_dir, "master.rds"))
}
print(st[, .N, by = status])
