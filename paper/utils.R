# Shared code of the figure and table scripts.
#   FMR_SIM_RESULTS  folder of the simulation results (default simulations/results)
#   FMR_OUTPUT       folder for the figures and tables (default output)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

sim_results <- function(...) file.path(Sys.getenv("FMR_SIM_RESULTS", here::here("simulations", "results")), ...)
read_master <- function(name) readRDS(sim_results(name, "master.rds"))
out_file <- function(f) {
  d <- Sys.getenv("FMR_OUTPUT", here::here("output"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, f)
}
ZQ <- qnorm(0.975)
BETA <- 1                                   # true causal effect in the simulations

# Bias, SD, SE, MSE and coverage of the intervals beta-hat +- mult * se, per cell.
# SE is the half-width / 1.96 (for bootstrap-t intervals, mult = c-hat).
cell_stats <- function(d, keys) {
  d %>% group_by(across(all_of(keys))) %>%
    summarise(reps = dplyr::n(), Bias = mean(beta - BETA), SD = sd(beta),
              SE = mean(mult * se) / ZQ, MSE = mean((beta - BETA)^2),
              CP = mean(abs(beta - BETA) < mult * se), .groups = "drop")
}

# Per-cell results of the methods of Figures 2 and 3 in scenario "S1" or "S2". At s_alpha = 0
# all cases are the same model; that cell (computed once, as Case 1) is shown for Cases 1-6.
fig23_data <- function(scenario) {
  m <- read_master(paste0("main_", scenario)) %>% filter(is.na(error))
  keys <- c("n", "signal", "case", "Dsignal")
  meth <- list(
    "FMR"           = m %>% transmute(n, signal, case, Dsignal, beta, se, mult = c_hat),
    "FMR (tau=0.1)" = m %>% transmute(n, signal, case, Dsignal, beta, se = se_tau01, mult = ZQ),
    "Ridge"         = m %>% transmute(n, signal, case, Dsignal, beta = ridge_beta, se = ridge_se, mult = ridge_c_hat),
    "F2SLS"         = m %>% transmute(n, signal, case, Dsignal, beta = f2sls_beta, se = f2sls_se, mult = ZQ),
    "sisVIVE"       = m %>% transmute(n, signal, case, Dsignal, beta = sisvive_beta, se = NA_real_, mult = ZQ),
    "LS"            = m %>% transmute(n, signal, case, Dsignal, beta = ls_beta, se = ls_se, mult = ZQ))
  s <- bind_rows(lapply(names(meth), function(k) cell_stats(meth[[k]], keys) %>% mutate(method = k)))
  bind_rows(s %>% filter(Dsignal == 0) %>% select(-case) %>% crossing(case = 1:6), s %>% filter(Dsignal > 0)) %>%
    arrange(method)
}

# Figures 2, 3: one panel per Case (rows) and size of the direct effect (columns).
COL <- c("FMR" = "#F8766D", "FMR (tau=0.1)" = "#A50F15", "Ridge" = "#E69F00",
         "F2SLS" = "#00BA38", "sisVIVE" = "#C77CFF", "LS" = "#619CFF")
SHP <- c("FMR" = 16, "FMR (tau=0.1)" = 1, "Ridge" = 18, "F2SLS" = 17, "sisVIVE" = 3, "LS" = 15)
LTY <- c("FMR" = "solid", "FMR (tau=0.1)" = "dashed", "Ridge" = "dotdash", "F2SLS" = "dotted",
         "sisVIVE" = "longdash", "LS" = "dashed")
LAB <- c("FMR" = "FMR", "FMR (tau=0.1)" = "FMR (τ = 0.1)", "Ridge" = "Ridge", "F2SLS" = "F2SLS",
         "sisVIVE" = "sisVIVE", "LS" = "LS")
alpha_level <- function(Ds) factor(ifelse(Ds == 0, "Zero", ifelse(Ds %in% c(0.3, 10), "Low", "High")),
                                   levels = c("Zero", "Low", "High"))
panel_labels <- labeller(case = function(v) paste("Case:", v), sa = function(v) paste("Alpha signal:", v))

method_figure <- function(d, x, y, xlab, ylab, methods, file, logy = FALSE, hline = NULL) {
  d <- d %>% filter(method %in% methods) %>%
    mutate(method = factor(method, levels = methods), sa = alpha_level(Dsignal))
  g <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], colour = method, shape = method, linetype = method)) +
    { if (!is.null(hline)) geom_hline(yintercept = hline, linetype = "dashed", linewidth = 0.3) } +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    facet_grid(case ~ sa, labeller = panel_labels) +
    scale_colour_manual(values = COL, labels = LAB, name = "Method") +
    scale_shape_manual(values = SHP, labels = LAB, name = "Method") +
    scale_linetype_manual(values = LTY, labels = LAB, name = "Method") +
    labs(x = xlab, y = ylab) + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  if (logy) g <- g + scale_y_log10(breaks = c(1e-3, 1e-2, 1e-1, 1), labels = c("0.001", "0.010", "0.100", "1.000"))
  if (!is.null(hline)) g <- g + scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1))
  ggsave(file, g, width = 8.5, height = 7.6, device = cairo_pdf)
  message("wrote ", file)
}

# LaTeX tables: numbers without the leading zero, as in the paper.
f3 <- function(x) ifelse(is.na(x), "--", sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x)))
f4 <- function(x) ifelse(is.na(x), "--", sub("^(-?)0\\.", "\\1.", sprintf("%.4f", x)))
fp <- function(x) ifelse(is.na(x), "--", ifelse(x >= 0.9995, "1", f3(x)))
show_first <- function(label, key = label) ifelse(duplicated(key), "", as.character(label))   # label once per group
write_table <- function(file, caption, label, align, head, rows) {
  writeLines(c("\\begin{table}[H]\\centering\\small",
               sprintf("\\caption{%s}\\label{%s}", caption, label),
               sprintf("\\begin{tabular}{%s}", align), "\\toprule", head, "\\midrule", rows,
               "\\bottomrule", "\\end{tabular}", "\\end{table}"), file)
  message("wrote ", file)
}
