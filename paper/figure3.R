# Figure 3: coverage rates of the 95% confidence intervals for beta in (a) S1 and (b) S2.
# FMR: bootstrap-t interval; FMR (tau = 0.1): Wald interval with the SE at 0.1 lambda_GCV.
# Input:  simulations/results/main_S1, main_S2 (simulations/main.R)
# Output: output/figure3a.pdf (S1), output/figure3b.pdf (S2)
here::i_am("paper/figure3.R")
source(here::here("paper", "utils.R"))
M <- c("FMR", "FMR (tau=0.1)", "Ridge", "F2SLS", "LS")
method_figure(fig23_data("S1"), "signal", "CP", "Signal strength s", "CR", M, out_file("figure3a.pdf"), hline = 0.95)
method_figure(fig23_data("S2"), "n", "CP", "n", "CR", M, out_file("figure3b.pdf"), hline = 0.95)
