# Table S.1: rejection rates of the T_alpha test and the Sargan test under H0: alpha(t) = 0
# (Case 1, s_alpha = 0), 1000 replications per cell:
#   S2: Z(t) simulated on 101 points, s = 0.015, n = 500, 1500, 2500, 3500;
#       data of replication i from set.seed(i + 1000); lambda grid 2^(-10:10);
#   S1: DLGAP2 genotype of all 808 subjects (needs the ADNI files), s = 0.01, ..., 0.04;
#       data of replication i from set.seed(i); lambda grid 2^(-16:16); the
#       genotype is adjusted for PC1-PC10 (missing genotypes left missing), the FPCA uses
#       pve 0.95, and the signs of the 11 eigenfunctions are fixed (S1_SIGNS below).
# Both use the first 11 FPCA scores as instruments (S2: pve 0.99) and q = 30 B-splines.
#
# T_alpha: penalized 2SLS with lambda chosen by GCV; the sandwich variance V of the spline
# coefficients delta at 0.1 * lambda; the statistic ||alpha-hat||^2 = delta' M delta referred
# to its weighted chi-square null distribution (Wood's approximation, dr::dr.pvalue).
# Sargan: over-identification test of the 2SLS fit with the 11 scores as instruments.
#
# Usage: Rscript simulations/talpha_sargan.R S2      (or S1)
# Output: simulations/results/talpha_sargan/<S2|S1>/master.rds, one row per replication.
here::i_am("simulations/talpha_sargan.R")
source(here::here("simulations", "helpers", "setup.R"))
suppressPackageStartupMessages({ library(mgcv); library(dr) })
SC <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(SC %in% c("S1", "S2"))
K <- 11L; NBASIS <- 30L

# p-value of T_alpha. Y, X: centred; Xi: centred scores (n x K); efunctions: nT x K on POS.
talpha_pvalue <- function(Y, X, Xi, POS, efunctions, L, theta = 0.1) {
  POS <- (POS - min(POS)) / (max(POS) - min(POS))
  bo <- create.bspline.basis(c(0, 1), NBASIS)
  w  <- c(0, diff(POS))
  A  <- Xi %*% t(efunctions) %*% (eval.basis(POS, bo) * w)
  Z1 <- cbind(1, Xi)
  Ps <- Z1 %*% ginv(t(Z1) %*% Z1, tol = 1e-15) %*% t(Z1)
  B  <- cbind(X, A)
  pen <- rbind(0, cbind(0, bsplinepen(bo)))
  max_df <- min(ncol(Xi), NBASIS + 1)
  dat <- list(Y.tilde = Ps %*% Y, B.tilde = Ps %*% B)
  fit_at <- function(lam) {
    f <- gam(data = dat, Y.tilde ~ -1 + B.tilde, paraPen = list(B.tilde = list(pen, sp = c(lam))))
    list(coef = unname(f$coefficients), df = f$df.null - f$df.residual)
  }
  lambdas <- 2^(-L:L)
  gcv <- vapply(lambdas, function(lam) {
    f <- fit_at(lam); e <- Y - B %*% f$coef
    max_df * abs((t(e) %*% Ps %*% e)[1]) / (max_df - f$df)^2
  }, numeric(1))
  lam <- lambdas[which.min(gcv)]
  f <- fit_at(lam); e <- Y - B %*% f$coef
  G <- ginv(t(B) %*% Ps %*% B + theta * lam * pen, tol = 1e-15)
  KB <- t(B) %*% Xi; U <- ginv(t(Xi) %*% Xi)
  V <- G %*% KB %*% U %*% (t(Xi) %*% diag(c(e^2)) %*% Xi) %*% U %*% t(KB) %*% G
  M <- bsplinepen(bo, Lfdobj = 0)
  delta <- f$coef[-1]
  D  <- M %*% V[-1, -1]
  ll <- eigen(D %*% t(D), symmetric = TRUE)$values
  ll[ll < 0] <- 0
  p <- dr.pvalue(ll^0.5, c(t(delta) %*% M %*% delta), chi2approx = "wood")$pval.adj
  if (is.na(p)) 1 else p
}

sargan_pvalue <- function(Y, X, Xi) {
  fit <- ivreg(Y ~ -1 + X | Xi, data = tibble(Y = Y, X = X, Xi = Xi))
  summary(fit)$diagnostics["Sargan", "p-value"]
}

# Under H0 the outcome is Y = 2 + X + eps, with X loading on the first 11 scores.
null_data <- function(Xi, signal) {
  e <- rmvnorm(nrow(Xi), mean = c(0, 0), sigma = SET$s2 * matrix(c(1, SET$rho, SET$rho, 1), 2, 2))
  X <- Xi[, 1:K] %*% (signal * c(SET$coef_PC1, 2, -1, 1, -1, 1, -1, 1, -1, 1, -1)) + e[, 1]
  list(Y = 2 + X * SET$beta + e[, 2], X = X)
}

one_rep <- function(Y, X, scores, POS, efunctions, L) {
  Y <- scale(Y, TRUE, FALSE); X <- scale(X, TRUE, FALSE); Xi <- scale(scores[, 1:K], TRUE, FALSE)
  c(p_Talpha = talpha_pvalue(Y, X, Xi, POS, efunctions[, 1:K], L), p_Sargan = sargan_pvalue(Y, X, Xi))
}

if (SC == "S2") {
  grid <- tibble(n = c(500L, 1500L, 2500L, 3500L), case = 1L, Dsignal = 0, signal = 0.015)
  rep_fun <- function(i, cell) {
    set.seed(i + 1000)
    D <- generate_S2(n = cell$n, nT = 101, effect = SET$beta, signal = cell$signal, Dsignal = 0,
                     s2 = SET$s2, case = 1, coef_PC1 = SET$coef_PC1, rho = SET$rho)
    one_rep(D$Y, D$X, D$Xi, D$POS, D$fpca$efunctions, L = 10)
  }
} else {
  source(here::here("lib", "adni_gene_data.R"))
  f_Z  <- genes_pca6_path("DLGAP2")
  f_pc <- here::here("data", "processed", "filtered_ADNI_808_plink2", "pca50_ADNI_808_plink2.csv")
  if (!all(file.exists(f_Z, f_pc))) stop("needs ", f_Z, " and ", f_pc, " (see data/README.md)")
  # Genes_pca6 is adjusted for PC1-PC6, so residualizing it on PC1-PC10 gives the genotype
  # adjusted for PC1-PC10. Missing genotypes stay missing (fpca.face handles them).
  Z  <- as.matrix(data.table::fread(f_Z, header = TRUE, data.table = FALSE))
  PC <- as.matrix(data.table::fread(f_pc, data.table = FALSE)[, paste0("PC", 1:10)])
  for (j in seq_len(ncol(Z))) {
    ok <- !is.na(Z[, j])
    Z[ok, j] <- lm.fit(cbind(1, PC[ok, ]), Z[ok, j])$residuals
  }
  pos <- as.numeric(colnames(Z)); POS1 <- (pos - min(pos)) / (max(pos) - min(pos))
  fp <- refund::fpca.face(Z, argvals = POS1, knots = 35, pve = 0.95)
  # The sign of an eigenfunction is arbitrary, but X loads on the scores with alternating
  # signs, so the signs are part of the design. Each of the first 11 eigenfunctions is
  # oriented to a positive integral, then multiplied by S1_SIGNS.
  S1_SIGNS <- c(1, 1, 1, -1, 1, 1, -1, 1, -1, -1, 1)
  sg <- sign(colSums(fp$efunctions[, 1:K] * c(0, diff(POS1)))) * S1_SIGNS
  fp$scores[, 1:K]     <- sweep(fp$scores[, 1:K], 2, sg, "*")
  fp$efunctions[, 1:K] <- sweep(fp$efunctions[, 1:K], 2, sg, "*")
  grid <- tibble(n = nrow(Z), case = 1L, Dsignal = 0, signal = c(0.01, 0.02, 0.03, 0.04))
  rep_fun <- function(i, cell) {
    set.seed(i)
    d <- null_data(fp$scores, cell$signal)
    one_rep(d$Y, d$X, fp$scores, POS1, fp$efunctions, L = 16)
  }
}

out <- results_dir("talpha_sargan", SC)
master <- run_cells(grid, function(i, cell) {
  row <- tibble(i = i, n = cell$n, case = cell$case, Dsignal = cell$Dsignal, signal = cell$signal)
  tryCatch(bind_cols(row, as_tibble_row(rep_fun(i, cell)), error = NA_character_),
           error = function(e) mutate(row, p_Talpha = NA_real_, p_Sargan = NA_real_, error = conditionMessage(e)))
}, out, n_reps(1000))
saveRDS(master, file.path(out, "master.rds"))
msg("wrote ", file.path(out, "master.rds"))
