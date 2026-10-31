# Table S.2: sensitivity of FMR to the number q of basis functions, S2, s = 0.015.
# SE (half-width / 1.96) and CP of the bootstrap-t interval. At s_alpha = 0, Cases 1-5 are one model.
# Input:  simulations/results/q_sensitivity (simulations/q_sensitivity.R)
# Output: output/tableS2.tex
here::i_am("paper/tableS2.R")
source(here::here("paper", "utils.R"))

s <- read_master("q_sensitivity") %>% filter(is.na(error)) %>%
  transmute(Dsignal, case, n, nbasis, beta, se, mult = c_hat) %>%
  cell_stats(c("Dsignal", "case", "n", "nbasis")) %>% arrange(Dsignal, case, n, nbasis)
sa <- sprintf("%.1f", s$Dsignal); cs <- ifelse(s$Dsignal == 0, "1--5", s$case)
rows <- sprintf("%s & %s & %s & %d & %s & %s & %s & %s & %s \\\\", show_first(sa), show_first(cs, paste(sa, cs)),
                show_first(s$n, paste(sa, cs, s$n)), s$nbasis, f3(s$Bias), f3(s$SD), f3(s$SE), f4(s$MSE), fp(s$CP))
write_table(out_file("tableS2.tex"), "Sensitivity of FMR to the number of basis functions in S2.", "tab:S2",
  "ccccrrrrr", "$s_\\alpha$ & Case & $n$ & $q$ & Bias & SD & SE & MSE & CP \\\\", rows)
