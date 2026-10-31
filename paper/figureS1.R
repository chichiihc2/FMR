# Figure S.1: estimated direct effects alpha-hat(t) over 1000 replications, (a) S1 at s = 0.01
# and (b) S2 at n = 1500. Dotted: true alpha(t); solid: mean of alpha-hat(t) and its empirical
# pointwise 95% band; dashed: mean of the estimated pointwise 95% bands.
# Input:  simulations/results/alpha_curves (simulations/alpha_curves.R)
# Output: output/figureS1a.pdf (S1), output/figureS1b.pdf (S2)
here::i_am("paper/figureS1.R")
source(here::here("paper", "utils.R"))

curves <- function(scenario) {
  res <- readRDS(sim_results("alpha_curves", paste0("curves_", scenario, ".rds")))
  null <- Filter(function(z) z$Dsignal == 0, res)[[1]]       # s_alpha = 0: one model for all cases
  Ds <- sort(unique(vapply(res, `[[`, 0, "Dsignal")))
  bind_rows(lapply(1:6, function(cs) bind_rows(lapply(Ds, function(D) {
    r <- if (D == 0) null else Filter(function(z) z$case == cs && z$Dsignal == D, res)[[1]]
    tibble(t = r$t, mean = r$mean, q_lo = r$q_lo, q_hi = r$q_hi, w_lo = r$w_lo, w_hi = r$w_hi,
           truth = approx(r$truth_t, r$truth, xout = r$t)$y, case = cs, Dsignal = D)
  }))))
}

alpha_figure <- function(d, file) {
  d <- d %>% mutate(sa = alpha_level(Dsignal))
  g <- ggplot(d, aes(x = t)) +
    geom_line(aes(y = truth), linetype = "dotted", linewidth = 0.7) +
    geom_line(aes(y = mean, colour = "FMR"), linewidth = 0.5) +
    geom_line(aes(y = q_lo, colour = "FMR"), linewidth = 0.3) + geom_line(aes(y = q_hi, colour = "FMR"), linewidth = 0.3) +
    geom_line(aes(y = w_lo, colour = "FMR"), linewidth = 0.3, linetype = "dashed") +
    geom_line(aes(y = w_hi, colour = "FMR"), linewidth = 0.3, linetype = "dashed") +
    facet_grid(case ~ sa, scales = "free_y", labeller = panel_labels) +
    scale_colour_manual(values = c(FMR = "#F8766D"), name = "Method") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "Position t", y = "Alpha") + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  ggsave(file, g, width = 8.5, height = 7.6, device = cairo_pdf)
  message("wrote ", file)
}
alpha_figure(curves("S1"), out_file("figureS1a.pdf"))
alpha_figure(curves("S2"), out_file("figureS1b.pdf"))
