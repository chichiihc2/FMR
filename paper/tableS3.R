# Table S.3: coverage of the 95% intervals for beta, S2, averaged over nine settings (Cases 1-4
# at s_alpha = 0.3, 0.6, and s_alpha = 0). Wald, tau: beta-hat(lambda_GCV) +- 1.96 SE(tau lambda_GCV);
# boot: bootstrap-t interval. F: mean first-stage F statistic.
# Input:  simulations/results/weak_iv (simulations/weak_iv.R) and main_S2 (s = 0.015)
# Output: output/tableS3.tex
here::i_am("paper/tableS3.R")
source(here::here("paper", "utils.R"))

m <- bind_rows(read_master("weak_iv"), read_master("main_S2") %>% filter(case <= 4)) %>% filter(is.na(error))
iv <- list("0.05" = m %>% mutate(se = se_tau005, mult = ZQ), "0.1" = m %>% mutate(se = se_tau01, mult = ZQ),
           "0.2" = m %>% mutate(se = se_tau02, mult = ZQ), "1" = m %>% mutate(mult = ZQ),
           "boot" = m %>% mutate(mult = c_hat))
cp <- bind_rows(lapply(names(iv), function(k)
  cell_stats(iv[[k]], c("signal", "n", "case", "Dsignal")) %>% mutate(interval = k))) %>%
  group_by(signal, n, interval) %>% summarise(CP = mean(CP), cells = dplyr::n(), .groups = "drop") %>%
  pivot_wider(names_from = interval, values_from = CP)
Fm <- m %>% group_by(signal, n, case, Dsignal) %>% summarise(F = mean(Fstat), .groups = "drop") %>%
  group_by(signal, n) %>% summarise(F = mean(F), .groups = "drop")
t <- cp %>% left_join(Fm, by = c("signal", "n")) %>% arrange(signal, n)
stopifnot(all(t$cells == 9))
rows <- sprintf("%s & %d & %.1f & %s & %s & %s & %s & %s \\\\", show_first(sub("^0", "", sprintf("%.3f", t$signal))),
                t$n, t$F, fp(t$`0.05`), fp(t$`0.1`), fp(t$`0.2`), fp(t$`1`), fp(t$boot))
write_table(out_file("tableS3.tex"), "Sensitivity of confidence interval coverage to the smoothing parameter in S2.",
  "tab:S3", "ccrrrrrr",
  c("& & & \\multicolumn{4}{c}{Wald, $\\tau$} & \\\\ \\cmidrule(lr){4-7}",
    "$s$ & $n$ & $F$ & 0.05 & 0.1 & 0.2 & 1 & boot \\\\"), rows)
