# Table S.4: FMR and sisVIVE with the FPCA scores (sisVIVE-PC) and with the LD-pruned genetic
# variants (sisVIVE-SNP) as instruments, S2, s = 0.015, the same 500 replications.
# Input:  simulations/results/sisvive_snp (simulations/sisvive_snp.R)
# Output: output/tableS4.tex
here::i_am("paper/tableS4.R")
source(here::here("paper", "utils.R"))

m <- read_master("sisvive_snp") %>% filter(is.na(error))
M <- c(FMR = "fmr_beta", "sisVIVE-PC" = "sisvive_pc_beta", "sisVIVE-SNP" = "sisvive_snp_beta")
s <- bind_rows(lapply(names(M), function(k)
  m %>% transmute(n, Dsignal, case, beta = .data[[M[[k]]]], se = NA_real_, mult = NA_real_) %>%
    cell_stats(c("n", "Dsignal", "case")) %>% mutate(method = k)))
w <- s %>% select(n, Dsignal, case, method, Bias, SD, MSE) %>%
  pivot_wider(names_from = method, values_from = c(Bias, SD, MSE)) %>% arrange(n, Dsignal, case)
cols <- c(paste0("Bias_", names(M)), paste0("SD_", names(M)), paste0("MSE_", names(M)))
rows <- vapply(seq_len(nrow(w)), function(r) paste(c(show_first(w$n)[r], sprintf("%g", w$Dsignal[r]),
  if (w$Dsignal[r] == 0) "1--4" else w$case[r],
  vapply(cols, function(cc) if (grepl("^MSE", cc)) f4(w[[cc]][r]) else f3(w[[cc]][r]), "")), collapse = " & "), "")
write_table(out_file("tableS4.tex"),
  sprintf("Comparison of FMR and sisVIVE using original genetic variants and FPCA scores in S2 (on average %.0f variants retained).", mean(m$L)),
  "tab:S4", "cccrrrrrrrrr",
  c("& & & \\multicolumn{3}{c}{Bias} & \\multicolumn{3}{c}{SD} & \\multicolumn{3}{c}{MSE} \\\\",
    "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9} \\cmidrule(lr){10-12}",
    paste0("$n$ & $s_\\alpha$ & Case & ", paste(rep(names(M), 3), collapse = " & "), " \\\\")),
  paste0(rows, " \\\\"))
