# Scenario S1 (Section 4.1): the instrument is the ADNI genotype of the DLGAP2 gene (its
# variants ordered by position, adjusted for six genotype PCs), restricted to the subjects
# with gene expression data. The FPCA of Z(t) is computed once; each replication draws
# new errors, X and Y on this fixed design.
#
# Needs the ADNI-derived files of data/README.md (Genes_pca6/DLGAP2.csv,
# Gene_expression_per_gene/DLGAP2.csv, data_cleaned_demo.csv, PTID.csv).

s1_build_design <- function(gene = "DLGAP2", pve = 0.80, knots = 35, n_scores = 11L) {
  source(here::here("lib", "adni_gene_data.R"), local = TRUE)
  demo <- data.table::fread(here::here("data", "raw", "data_cleaned_demo.csv"), data.table = FALSE)
  demo <- demo[!duplicated(demo$PTID), , drop = FALSE]
  ptid <- data.table::fread(here::here("data", "raw", "PTID.csv"), data.table = FALSE)
  fr <- load_gene_frame(gene, demo, ptid, here::here("data", "processed", "Gene_expression_per_gene"))
  if (is.null(fr)) stop("no analysis data for ", gene)
  fp <- refund::fpca.face(fr$Z, argvals = fr$t01, knots = knots, pve = pve)
  if (ncol(fp$scores) < n_scores) stop("fpca.face retained fewer than ", n_scores, " components")
  list(gene = gene, Z = fr$Z, POS = fr$t01, n = fr$n_obs, p = fr$n_snps, fpca = fp,
       scores = fp$scores, efunctions = fp$efunctions, n_scores = n_scores)
}

# Direct effect alpha(t) of Cases 1-6, with size s_alpha = Dsignal.
alpha_S1 <- function(case, Dsignal) {
  switch(as.character(case),
    `1` = function(t) 2 * Dsignal * t,
    `2` = function(t) Dsignal * (t + sin(pi * t)),
    `3` = function(t) 2 * Dsignal * (1 / (1 + exp(-8 * (t - 0.5)))),
    `4` = function(t) 5 * Dsignal * ifelse(t < 0.5, t, 1 - t),
    `5` = function(t) 1.5 * ifelse((t - 0.2) * (t - 0.8) < 0, Dsignal, 0),
    `6` = function(t) 2.9420 * Dsignal * ifelse(t > 0.25 & t < 0.75, sin(pi * (t - 0.25) / 0.5), 0),
    stop("case must be 1-6"))
}

# Returns a data generator on the fixed design `des`.
make_generate_S1 <- function(des) {
  POS <- des$POS
  wq  <- c(0, diff(POS))
  Xi  <- des$scores
  efs <- des$efunctions
  function(effect = 1, signal = 0.02, Dsignal = 0, s2 = 1, case = 1, coef_PC1 = 2, rho = 0.8) {
    sf <- alpha_S1(case, Dsignal)
    e  <- mvtnorm::rmvnorm(des$n, mean = c(0, 0), sigma = s2 * matrix(c(1, rho, rho, 1), 2, 2))
    gamma <- signal * c(coef_PC1, 2, -1, 1, -1, 1, -1, 1, -1, 1, -1)
    X  <- as.numeric(Xi[, seq_len(des$n_scores), drop = FALSE] %*% gamma + e[, 1])
    ft <- as.numeric(Xi %*% (t(efs) %*% (sf(POS) * wq)))
    list(Y = 2 + effect * X + ft + e[, 2], X = X, Xi = Xi, Z = des$Z, POS = POS,
         fpca = des$fpca, alpha_t = sf)
  }
}
