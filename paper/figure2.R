# Figure 2: MSE of the estimators of the causal effect beta in (a) S1 and (b) S2.
# Input:  simulations/results/main_S1, main_S2 (simulations/main.R)
# Output: output/figure2a.pdf (S1), output/figure2b.pdf (S2)
here::i_am("paper/figure2.R")
source(here::here("paper", "utils.R"))
M <- c("FMR", "Ridge", "F2SLS", "sisVIVE", "LS")
method_figure(fig23_data("S1"), "signal", "MSE", "Signal strength s", "MSE", M, out_file("figure2a.pdf"), logy = TRUE)
method_figure(fig23_data("S2"), "n", "MSE", "n", "MSE", M, out_file("figure2b.pdf"), logy = TRUE)
