# Section 5.2 and Figure S.2: the APOE locus (APOE +- 50 kb, chr19, GRCh37) as a functional
# instrument for amyloid PET -> FAQ and BMI -> ADAS13.
#
# Exposure, outcome and genotypes are standardized and adjusted for age, sex, education and
# PC1-PC6. FPCA of the locus (fpca.face, 35 knots), k by first-stage BIC (k >= 3). For each pair:
#   FMR: lambda by GCV, beta-hat, sandwich SE, bootstrap-t interval (B = 1999);
#   F2SLS on the k scores; 2SLS and sisVIVE with the APOE e4 and e2 SNPs (rs429358, rs7412);
#   paired bootstrap tests of FMR - F2SLS and FMR - 2SLS (B = 1999);
#   T_alpha test of H0: alpha = 0; Sargan test of the e4/e2 2SLS;
#   alpha-hat(t) with a 95% bootstrap simultaneous band (B = 999).
#
# Usage: Rscript realdata/apoe/fmr_apoe.R       (after genotypes.sh and phenotypes.R)
# Output: realdata/apoe/results/{amy_cl_faq, bmi_adas13}.rds
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(fda); library(refund); library(MASS)
  library(ivreg); library(sisVIVE); library(parallel)
})
source("lib/estimators.R"); source("lib/fmr_talpha.R")
NCORES <- as.integer(Sys.getenv("FMR_CORES", max(1L, parallel::detectCores() - 1L)))
OUT <- "realdata/apoe/results"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PAIRS <- list(c("amy_cl", "faq"), c("bmi", "adas13"))
APOE_BP <- 45411941                                     # APOE (rs429358)

# ---- data ------------------------------------------------------------------------------
P <- "data/processed/apoe"
W <- read_csv(file.path(P, "pheno.csv"), show_col_types = FALSE)
raw <- read.table(file.path(P, "win.raw"), header = TRUE, check.names = FALSE)
G0 <- as.matrix(raw[, -(1:6)]); rownames(G0) <- raw$IID
bim <- read.table(file.path(P, "apoe_region.bim"), col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))
bim <- bim[bim$bp >= 45.0e6 & bim$bp <= 45.8e6, ]
o <- order(bim$bp); G0 <- G0[, o]; bim <- bim[o, ]
k <- !duplicated(bim$bp); G0 <- G0[, k]; bim <- bim[k, ]
for (j in seq_len(ncol(G0))) { m <- is.na(G0[, j]); if (any(m)) G0[m, j] <- mean(G0[, j], na.rm = TRUE) }
jw <- which(abs(bim$bp - APOE_BP) <= 50e3)

# ---- functions -------------------------------------------------------------------------
# k minimizing the first-stage BIC over k >= 3, and the first-stage F of each k.
select_k <- function(X, scores, k_max = 20L, k_min = 3L) {
  n <- length(X); ks <- seq_len(min(k_max, ncol(scores)))
  rss <- vapply(ks, function(k) sum(lm.fit(cbind(1, scores[, 1:k, drop = FALSE]), X)$residuals^2), 0)
  bic <- n * log(rss / n) + (ks + 1) * log(n)
  F <- vapply(ks, function(k) ((sum((X - mean(X))^2) - rss[k]) / k) / (rss[k] / (n - k - 1)), 0)
  list(k = ks[ks >= k_min][which.min(bic[ks >= k_min])], F = F)
}

design_real <- function(Y, X, P01, fpca, k, BO) {
  Xi <- scale(fpca$scores[, seq_len(k), drop = FALSE], TRUE, FALSE)
  M2 <- t(fpca$efunctions[, seq_len(k), drop = FALSE]) %*% (eval.basis(P01, BO) * c(0, diff(P01)))
  list(Y = drop(scale(Y, TRUE, FALSE)), X = drop(scale(X, TRUE, FALSE)), Xi = Xi, A = Xi %*% M2)
}

# FMR at lambda_GCV with the bootstrap-t critical value c-hat (pairs bootstrap, lambda re-chosen).
fmr_fit <- function(d, pen, B, seed, LG = 2^(-10:10), tol = 1e-15) {
  n <- length(d$Y); k <- ncol(d$Xi)
  cp0 <- xprod(d$Y, cbind(d$X, d$A), cbind(1, d$Xi), tol)
  lam <- gcv_lambda(cp0, pen$PP_rough, LG, k)
  s0  <- se_sandwich(cp0, lam, pen$PP_rough)
  set.seed(seed)
  tt <- vapply(seq_len(B), function(b) {
    id  <- sample.int(n, n, replace = TRUE)
    cpb <- tryCatch(xprod(d$Y[id], cbind(d$X[id], d$A[id, , drop = FALSE]), cbind(1, d$Xi[id, , drop = FALSE]), tol),
                    error = function(e) NULL)
    if (is.null(cpb)) return(NA_real_)
    lb <- gcv_lambda(cpb, pen$PP_rough, LG, k)
    if (!is.finite(lb)) lb <- lam
    sb <- tryCatch(se_sandwich(cpb, lb, pen$PP_rough), error = function(e) NULL)
    if (is.null(sb) || !is.finite(sb$beta) || !is.finite(sb$se) || sb$se <= 0) NA_real_ else (sb$beta - s0$beta) / sb$se
  }, numeric(1))
  c_hat <- unname(quantile(abs(tt[is.finite(tt)]), 0.95))
  se_bt <- c_hat * s0$se / qnorm(0.975)                 # half-width of the interval / 1.96
  list(beta = s0$beta, lambda = lam, se = s0$se, c_hat = c_hat, se_bt = se_bt,
       p = 2 * pnorm(-abs(s0$beta) / se_bt), ci = s0$beta + c(-1, 1) * c_hat * s0$se)
}

# 2SLS (ivreg, with intercept): beta-hat, SE, p, first-stage F, Sargan p.
iv_fit <- function(Y, X, Zm) {
  s <- summary(ivreg(Y ~ X | Zm), diagnostics = TRUE)
  b <- unname(s$coefficients["X", 1]); se <- unname(s$coefficients["X", 2])
  list(beta = b, se = se, p = 2 * pnorm(-abs(b / se)), F = unname(summary(lm(X ~ Zm))$fstatistic[1]),
       sargan_p = if ("Sargan" %in% rownames(s$diagnostics)) unname(s$diagnostics["Sargan", "p-value"]) else NA)
}
b2sls <- function(Y, X, Z) { Z <- cbind(1, Z); P <- Z %*% solve(crossprod(Z), t(Z))
  drop(crossprod(X, P %*% Y) / crossprod(X, P %*% X)) }

# alpha-hat(t) of the efficient GMM fit and its bootstrap simultaneous band: each draw
# redoes the FPCA, lambda and the fit; band = alpha-hat(t) +- q sd*(t), q the 95% quantile of
# max_t |alpha*(t) - alpha-hat(t)| / sd*(t).
alpha_band <- function(Y, X, Z, POS, n_scores, B = 999L, seed = 2L, nbasis = 30L, knots = 35, pve = 0.99,
                       grid = 0:100 / 100, LG = 2^(-16:16), tol = 1e-15) {
  n <- length(Y); POS01 <- (POS - min(POS)) / (max(POS) - min(POS))
  BO <- create.bspline.basis(c(0, 1), nbasis)
  PP <- rbind(0, cbind(0, bsplinepen(BO)))
  Psiw <- eval.basis(POS01, BO) * c(0, diff(POS01))
  Psi_g <- eval.basis(grid, BO)
  fit <- function(id) {
    fp <- fpca.face(Z[id, , drop = FALSE], argvals = POS01, knots = knots, pve = pve)
    if (ncol(fp$scores) < n_scores) return(NULL)
    Xi <- scale(fp$scores[, seq_len(n_scores), drop = FALSE], TRUE, FALSE)
    Yc <- drop(scale(Y[id], TRUE, FALSE)); Xc <- drop(scale(X[id], TRUE, FALSE))
    Bm <- cbind(Xc, Xi %*% (t(fp$efunctions[, seq_len(n_scores), drop = FALSE]) %*% Psiw))
    cp <- xprod(Yc, Bm, cbind(1, Xi), tol)
    lam <- gcv_lambda(cp, PP, LG, n_scores)
    if (!is.finite(lam)) return(NULL)
    th1 <- drop(solve(cp$BPsB + lam * PP, cp$BPsY))
    Wi <- ginv(crossprod(Xi * drop(Yc - Bm %*% th1)), tol = tol)
    G1 <- crossprod(Bm, Xi)
    th2 <- tryCatch(drop(solve(G1 %*% Wi %*% crossprod(Xi, Bm) + lam * PP, G1 %*% Wi %*% crossprod(Xi, Yc))),
                    error = function(e) rep(NA_real_, ncol(Bm)))
    list(lambda = lam, curve = drop(Psi_g %*% th2[-1]))
  }
  f0 <- fit(seq_len(n))
  set.seed(seed)
  idx <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
  Cb <- do.call(rbind, lapply(idx, function(id) {
    fb <- tryCatch(fit(id), error = function(e) NULL)
    if (is.null(fb)) rep(NA_real_, length(grid)) else fb$curve
  }))
  Cb <- Cb[stats::complete.cases(Cb), , drop = FALSE]
  sdt <- apply(Cb, 2, stats::sd)
  sdt[!is.finite(sdt) | sdt <= 0] <- NA_real_
  flo <- stats::median(sdt, na.rm = TRUE) * 1e-6
  sdt[is.na(sdt)] <- flo; sdt <- pmax(sdt, flo)
  Tsup <- sort(apply(abs(sweep(sweep(Cb, 2, f0$curve, "-"), 2, sdt, "/")), 1, max))
  q <- Tsup[min(length(Tsup), max(1L, ceiling(0.95 * (length(Tsup) + 1L))))]
  list(bp = min(POS) + grid * diff(range(POS)), alpha = f0$curve, lo = f0$curve - q * sdt, hi = f0$curve + q * sdt,
       q = q, n_boot = nrow(Cb))
}

# ---- analysis of one exposure-outcome pair ----------------------------------------------
fit_pair <- function(x, y, B = 1999L) {
  D <- W %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  C <- model.matrix(~ age_exact + sex + educ + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = D)
  rs <- function(M) M - C %*% qr.coef(qr(C), M)
  X <- drop(rs(matrix(as.numeric(scale(D[[x]])))))
  Y <- drop(rs(matrix(as.numeric(scale(D[[y]])))))
  Gr <- rs(G0[D$PTID, jw]); ok <- apply(Gr, 2, sd) > 1e-8; Gr <- Gr[, ok]
  bp <- bim$bp[jw][ok]; snp <- bim$snp[jw][ok]; P01 <- (bp - min(bp)) / diff(range(bp))
  fp <- fpca.face(Gr, argvals = P01, knots = 35, npc = 20)
  ks <- select_k(X, fp$scores); K <- ks$k; n <- length(Y)
  pen <- make_pen(30L); d <- design_real(Y, X, P01, fp, K, pen$BO)
  E <- scale(cbind(e4 = Gr[, snp == "rs429358"], e2 = Gr[, snp == "rs7412"]), TRUE, FALSE)

  fm <- fmr_fit(d, pen, B = B, seed = 20261007)
  f2 <- iv_fit(d$Y, d$X, d$Xi)
  ee <- iv_fit(d$Y, d$X, E)
  set.seed(9e6)
  sv <- cv.sisVIVE(Y = as.numeric(d$Y), D = as.numeric(d$X), Z = unname(E), K = 10, intercept = FALSE)
  ta <- fmr_talpha(d$Y, d$X, d$Xi, d$A, fm$lambda, pen$PP_rough, pen$BO, tau = 0.1)

  # paired bootstrap of the differences FMR - F2SLS and FMR - 2SLS(e4, e2)
  dd <- t(simplify2array(mclapply(seq_len(B), function(b) {
    set.seed(7e6 + b); id <- sample.int(n, n, TRUE)
    cpb <- xprod(d$Y[id], cbind(d$X[id], d$A[id, , drop = FALSE]), cbind(1, d$Xi[id, , drop = FALSE]))
    lb <- gcv_lambda(cpb, pen$PP_rough, 2^(-10:10), K); if (!is.finite(lb)) lb <- fm$lambda
    bf <- tryCatch(se_sandwich(cpb, lb, pen$PP_rough)$beta, error = function(e) NA_real_)
    c(bf - b2sls(d$Y[id], d$X[id], d$Xi[id, , drop = FALSE]), bf - b2sls(d$Y[id], d$X[id], E[id, , drop = FALSE]))
  }, mc.cores = NCORES)))
  d0 <- c(fm$beta - f2$beta, fm$beta - ee$beta)
  p_diff <- vapply(1:2, function(j) { z <- dd[is.finite(dd[, j]), j]
    (1 + sum(abs(z - mean(z)) >= abs(d0[j]))) / (length(z) + 1) }, 0)

  F <- ks$F[K]
  list(x = x, y = y, n = n, n_variants = ncol(Gr), k = K, F = F, partial_R2 = F * K / (F * K + n - K - 1),
       fmr = fm, f2sls = f2, tsls_e4e2 = ee, sisvive = list(beta = sv$beta, n_invalid = sum(sv$alpha != 0)),
       talpha = list(T = ta[["T"]], p = ta[["p"]]), diff = list(est = d0, p = p_diff),
       band = alpha_band(d$Y, d$X, Gr, bp, n_scores = K))
}

for (pr in PAIRS) {
  r <- fit_pair(pr[1], pr[2])
  saveRDS(r, file.path(OUT, paste0(pr[1], "_", pr[2], ".rds")))
  cat(sprintf("%s -> %s: n = %d, k = %d, F = %.1f, FMR beta = %.3f (p = %.3g), F2SLS %.3f, 2SLS %.3f, T_alpha p = %.3f\n",
              pr[1], pr[2], r$n, r$k, r$F, r$fmr$beta, r$fmr$p, r$f2sls$beta, r$tsls_e4e2$beta, r$talpha$p))
}
