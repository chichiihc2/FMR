# Table S.1: rejection rates at level 0.05 of the T_alpha test and the Sargan test under
# H0: alpha(t) = 0 (Case 1, s_alpha = 0), S1 by s and S2 by n, 1000 replications.
# Input:  simulations/results/talpha_sargan/{S1,S2}/master.rds (simulations/talpha_sargan.R)
# Output: output/tableS1.tex
here::i_am("paper/tableS1.R")
source(here::here("paper", "utils.R"))

r <- bind_rows(lapply(c("S1", "S2"), function(s)
  mutate(readRDS(sim_results("talpha_sargan", s, "master.rds")), scenario = s))) %>%
  group_by(scenario, signal, n) %>%
  summarise(Ta = mean(p_Talpha < 0.05), Sg = mean(p_Sargan < 0.05), .groups = "drop")
s1 <- r %>% filter(scenario == "S1") %>% arrange(signal)
s2 <- r %>% filter(scenario == "S2") %>% arrange(n)
rows <- sprintf("%.3f & %.3f & %.3f & %d & %.3f & %.3f \\\\", s1$signal, s1$Ta, s1$Sg, s2$n, s2$Ta, s2$Sg)
write_table(out_file("tableS1.tex"),
  "Rejection rates under the null hypothesis $\\alpha^*(t)=0$ across two scenarios.", "tab:S1", "cccccc",
  c("\\multicolumn{3}{c}{Scenario 1} & \\multicolumn{3}{c}{Scenario 2} \\\\ \\cmidrule(lr){1-3} \\cmidrule(lr){4-6}",
    "$s$ & $T_\\alpha$ & $T_{\\rm Sargan}$ & $n$ & $T_\\alpha$ & $T_{\\rm Sargan}$ \\\\"), rows)
