# Table S.2: sensitivity of FMR to the number q of B-spline functions of alpha(t). S2 with
# s = 0.015, n = 500, 1500 and q = 15, 20, 30, 60; Cases 1-5 at s_alpha = 0.6, and s_alpha = 0.
# 1000 replications per cell.
#
# Usage: Rscript simulations/q_sensitivity.R
# Output: simulations/results/q_sensitivity/master.rds
here::i_am("simulations/q_sensitivity.R")
source(here::here("simulations", "helpers", "setup.R"))

QS <- c(15L, 20L, 30L, 60L)
pens <- setNames(lapply(QS, make_pen), QS)
grid <- expand_grid(nbasis = QS, n = c(500L, 1500L),
                    bind_rows(tibble(case = 1L, Dsignal = 0), tibble(case = 1:5, Dsignal = 0.6))) %>%
  mutate(signal = 0.015)

out <- results_dir("q_sensitivity")
master <- run_cells(grid, function(i, cell)
  sim_rep(i, cell$n, cell$case, cell$Dsignal, cell$signal, make_S2, pens[[as.character(cell$nbasis)]],
          SET$lambda_S2, ridge = FALSE, others = FALSE), out, n_reps(1000)) %>%
  dplyr::select(i, n, case, Dsignal, signal, nbasis, lambda, beta, se, c_hat, error)
saveRDS(master, file.path(out, "master.rds"))
msg("wrote ", file.path(out, "master.rds"))
