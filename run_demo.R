#!/usr/bin/env Rscript
# =============================================================================
# run_demo.R -- one-click demo of the FMR package (a few seconds).
#
#   Rscript run_demo.R
#
# Simulates one data set from the paper's Scenario S2 (a direct effect alpha(t) of the
# functional instrument on the outcome, i.e. pleiotropy, and unmeasured confounding of X and
# Y), fits FMR, and compares it with F2SLS and OLS, which ignore alpha(t). Writes the plot of
# alpha-hat(t) to output/demo/alpha_hat.pdf. Installs the package first if needed.
# =============================================================================
if (!requireNamespace("FMR", quietly = TRUE)) source("install.R")
library(FMR)

set.seed(2026)
d <- simulate_fmr_data(n = 1500,        # subjects
                       nT = 501,        # positions of the functional instrument Z(t)
                       beta = 1,        # true causal effect
                       s_alpha = 0.6,   # size of the direct (pleiotropic) effect
                       case = 2)        # alpha(t) = s_alpha (t + sin(pi t))
str(d[c("y", "x", "Z", "pos")])

fit <- fmr(y = d$y, x = d$x, Z = d$Z, pos = d$pos,
           k = NULL,          # number of FPCA scores: first-stage BIC over 3..10
           B = 199,           # bootstrap-t draws
           tau = 0.1,         # Wald SE evaluated at tau * lambda
           seed = 1)
summary(fit)
cat("\nTrue beta = 1. F2SLS and OLS ignore the direct effect alpha(t) and are biased;",
    "FMR estimates it.\n")

dir.create("output/demo", recursive = TRUE, showWarnings = FALSE)
pdf("output/demo/alpha_hat.pdf", width = 6, height = 4.2)
plot(fit, xlab = "position t")
tt <- fit$alpha$pos; lines(tt, d$alpha(tt), lty = 2, lwd = 2)
legend("bottomright", "true alpha(t)", lty = 2, lwd = 2, bty = "n")
invisible(dev.off())
cat("plot of alpha-hat(t): output/demo/alpha_hat.pdf\n")

# With your own data:
#   fit <- fmr(y, x, Z, pos = bp_positions, covariates = data.frame(age, sex, PC1, PC2, ...))
