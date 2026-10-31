# Table S.4: sisVIVE with the genetic variants as candidate instruments (after LD pruning:
# sliding window of 50 variants, r^2 < 0.1) and with the FPCA scores, and FMR, on the same
# replications. S2 with s = 0.015, n = 1500, 2500; Cases 1-4 at s_alpha = 0.6, and s_alpha = 0.
# 500 replications per cell.
#
# Usage: Rscript simulations/sisvive_snp.R
# Output: simulations/results/sisvive_snp/master.rds
here::i_am("simulations/sisvive_snp.R")
source(here::here("simulations", "helpers", "setup.R"))

# Greedy LD pruning: keep variant j unless r^2 >= r2 with a kept variant among the previous
# `window` positions; then drop linearly dependent columns.
ld_prune <- function(Z, window = 50L, r2 = 0.1) {
  Z <- scale(Z, TRUE, FALSE); s <- sqrt(colSums(Z^2))
  keep <- integer(0)
  for (j in which(s > 0)) {
    near <- keep[keep >= j - window]
    if (length(near)) {
      r <- drop(crossprod(Z[, near, drop = FALSE], Z[, j])) / (s[near] * s[j])
      if (any(r^2 >= r2)) next
    }
    keep <- c(keep, j)
  }
  q <- qr(Z[, keep, drop = FALSE])
  sort(keep[q$pivot[seq_len(q$rank)]])
}

pen <- make_pen(SET$nbasis)
rep_snp <- function(i, cell) {
  row <- tibble(i = i, n = cell$n, case = cell$case, Dsignal = cell$Dsignal, signal = cell$signal)
  tryCatch({
    set.seed(i)
    D <- generate_S2(n = cell$n, nT = SET$nT, effect = SET$beta, signal = cell$signal,
                     Dsignal = cell$Dsignal, s2 = SET$s2, case = cell$case,
                     coef_PC1 = SET$coef_PC1, rho = SET$rho)
    d <- design(D, pen$BO, SET$n_scores)
    kp <- ld_prune(D$Z)
    cp <- xprod(d$Y, cbind(d$X, d$A), cbind(1, d$Xi), SET$tol)
    fmr <- se_sandwich(cp, gcv_lambda(cp, pen$PP_rough, SET$lambda_S2, SET$max_df), pen$PP_rough)
    mutate(row, L = length(kp), fmr_beta = fmr$beta,
           sisvive_pc_beta = sisvive_of(d$Y, d$X, d$Xi, i),
           sisvive_snp_beta = sisvive_of(d$Y, d$X, scale(D$Z[, kp], TRUE, FALSE), i),
           error = NA_character_)
  }, error = function(e) mutate(row, error = conditionMessage(e)))
}

grid <- expand_grid(n = c(1500L, 2500L),
                    bind_rows(tibble(case = 1L, Dsignal = 0), tibble(case = 1:4, Dsignal = 0.6))) %>%
  mutate(signal = 0.015)
out <- results_dir("sisvive_snp")
master <- run_cells(grid, rep_snp, out, n_reps(500))
saveRDS(master, file.path(out, "master.rds"))
msg("wrote ", file.path(out, "master.rds"))
