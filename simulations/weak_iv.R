# Table S.3: coverage of the FMR intervals under weak instruments. S2 with s = 0.002, 0.004,
# 0.007 and n = 500, 1500, 2500; Cases 1-4 at s_alpha = 0.3, 0.6, and s_alpha = 0.
# (The s = 0.015 rows of the table come from main.R S2.) 1000 replications per cell.
#
# Usage: Rscript simulations/weak_iv.R
# Output: simulations/results/weak_iv/master.rds
here::i_am("simulations/weak_iv.R")
source(here::here("simulations", "helpers", "setup.R"))

pen <- make_pen(SET$nbasis)
grid <- expand_grid(signal = c(0.002, 0.004, 0.007), n = c(500L, 1500L, 2500L),
                    bind_rows(tibble(case = 1L, Dsignal = 0), expand_grid(case = 1:4, Dsignal = c(0.3, 0.6))))

out <- results_dir("weak_iv")
master <- run_cells(grid, function(i, cell)
  sim_rep(i, cell$n, cell$case, cell$Dsignal, cell$signal, make_S2, pen, SET$lambda_S2,
          ridge = FALSE, others = FALSE), out, n_reps(1000))
saveRDS(master, file.path(out, "master.rds"))
msg("wrote ", file.path(out, "master.rds"))
