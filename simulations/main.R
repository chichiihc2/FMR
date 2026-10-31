# Figures 2 and 3 (and the s = 0.015 rows of Table S.3): FMR, Ridge, F2SLS, sisVIVE and LS in
#   S2: simulated Z(t), s = 0.015, n = 500, 1500, 2500, 3500;
#   S1: DLGAP2 genotype (needs the ADNI files), n = 712, s = 0.01, 0.02, 0.03, 0.04;
# Cases 1-6 at a low and a high s_alpha, and s_alpha = 0 (the same model for all cases).
# 1000 replications per cell.
#
# Usage: Rscript simulations/main.R S2      (or S1)
# Output: simulations/results/main_S2/master.rds, one row per replication.
here::i_am("simulations/main.R")
source(here::here("simulations", "helpers", "setup.R"))
SC <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(SC %in% c("S1", "S2"))

pen <- make_pen(SET$nbasis)
cases <- function(Ds) bind_rows(tibble(case = 1L, Dsignal = 0), expand_grid(case = 1:6, Dsignal = Ds))
if (SC == "S2") {
  grid <- expand_grid(n = c(500L, 1500L, 2500L, 3500L), cases(c(0.3, 0.6))) %>% mutate(signal = 0.015)
  maker <- make_S2; LG <- SET$lambda_S2
} else {
  s1_init(pen)
  grid <- expand_grid(signal = c(0.01, 0.02, 0.03, 0.04), cases(c(10, 20))) %>% mutate(n = S1$des$n)
  maker <- make_S1; LG <- SET$lambda_S1
}

out <- results_dir(paste0("main_", SC))
master <- run_cells(grid, function(i, cell)
  sim_rep(i, cell$n, cell$case, cell$Dsignal, cell$signal, maker, pen, LG), out, n_reps(1000))
saveRDS(master, file.path(out, "master.rds"))
msg("wrote ", file.path(out, "master.rds"))
