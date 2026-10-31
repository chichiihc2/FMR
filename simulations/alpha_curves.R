# Figure S.1: the estimated direct effect alpha-hat(t) over 1000 replications, S2 at n = 1500
# (s = 0.015) and S1 at s = 0.01; Cases 1-6 at a low and a high s_alpha, and s_alpha = 0.
# For each cell: the mean of alpha-hat(t), its empirical pointwise 2.5% and 97.5% quantiles,
# and the mean of the estimated pointwise 95% bands alpha-hat(t) +- 1.96 sd(t).
#
# Usage: Rscript simulations/alpha_curves.R S2      (or S1)
# Output: simulations/results/alpha_curves/curves_S2.rds
here::i_am("simulations/alpha_curves.R")
source(here::here("simulations", "helpers", "setup.R"))
SC <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(SC %in% c("S1", "S2"))

pen  <- make_pen(SET$nbasis)
TG   <- seq(0, 1, length.out = 201)               # evaluation grid of alpha-hat(t)
PSI  <- eval.basis(TG, pen$BO)
REPS <- n_reps(1000)

# alpha-hat(t) and its pointwise sd, sd(t)^2 = psi(t)' V_delta psi(t), V the sandwich of theta.
one_curve <- function(i, case, Dsignal) {
  if (SC == "S2") {
    set.seed(i)
    D <- generate_S2(n = 1500L, nT = SET$nT, effect = SET$beta, signal = 0.015, Dsignal = Dsignal,
                     s2 = SET$s2, case = case, coef_PC1 = SET$coef_PC1, rho = SET$rho)
    d <- design(D, pen$BO, SET$n_scores); LG <- SET$lambda_S2
  } else {
    d <- make_S1(i, S1$des$n, case, Dsignal, 0.01, pen); LG <- SET$lambda_S1
  }
  cp <- xprod(d$Y, cbind(d$X, d$A), cbind(1, d$Xi), SET$tol)
  lam <- gcv_lambda(cp, pen$PP_rough, LG, SET$max_df)
  if (!is.finite(lam)) return(NULL)
  Gi <- solve(cp$BPsB + lam * pen$PP_rough)
  th <- drop(Gi %*% cp$BPsY)
  e  <- drop(cp$Y - cp$B %*% th)
  s2 <- sum(e^2) / (length(e) - sum(diag(Gi %*% cp$BPsB)))
  Vd <- (s2 * Gi %*% cp$BPsB %*% Gi)[-1, -1]
  list(i = i, lambda = lam, beta = th[1], alpha = drop(PSI %*% th[-1]),
       sd = sqrt(pmax(rowSums((PSI %*% Vd) * PSI), 0)))
}

run_cell <- function(case, Dsignal) {
  fits <- Filter(Negate(is.null), mclapply(seq_len(REPS), function(i)
    tryCatch(one_curve(i, case, Dsignal), error = function(e) NULL), mc.cores = max(1L, n_cores() - 1L)))
  A  <- do.call(rbind, lapply(fits, `[[`, "alpha"))
  SD <- do.call(rbind, lapply(fits, `[[`, "sd"))
  tt <- seq(0, 1, length.out = SET$nT)
  alpha_t <- if (SC == "S2") alpha_S2(case, Dsignal) else alpha_S1(case, Dsignal)
  list(case = case, Dsignal = Dsignal, reps = nrow(A),
       betas = tibble(i = vapply(fits, `[[`, 0, "i"), beta = vapply(fits, `[[`, 0, "beta"),
                      lambda = vapply(fits, `[[`, 0, "lambda")),
       t = TG, mean = colMeans(A),
       q_lo = apply(A, 2, quantile, 0.025, names = FALSE),
       q_hi = apply(A, 2, quantile, 0.975, names = FALSE),
       w_lo = colMeans(A - ZQ * SD), w_hi = colMeans(A + ZQ * SD),
       truth_t = tt, truth = vapply(tt, alpha_t, 0))
}

if (SC == "S1") s1_init(pen)
Ds <- if (SC == "S2") c(0.3, 0.6) else c(10, 20)
cells <- bind_rows(tibble(case = 1L, Dsignal = 0), expand_grid(case = 1:6, Dsignal = Ds))
out <- results_dir("alpha_curves"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
res <- lapply(seq_len(nrow(cells)), function(k) {
  msg(SC, " case ", cells$case[k], ", s_alpha = ", cells$Dsignal[k])
  run_cell(cells$case[k], cells$Dsignal[k])
})
saveRDS(res, file.path(out, paste0("curves_", SC, ".rds")))
msg("wrote ", file.path(out, paste0("curves_", SC, ".rds")))
