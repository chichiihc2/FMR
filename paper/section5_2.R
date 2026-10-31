# Section 5.2: the numbers of the APOE locus analysis quoted in the text.
# Input:  realdata/apoe/results/{amy_cl_faq, bmi_adas13}.rds (realdata/apoe/fmr_apoe.R)
# Output: output/section5_2.txt
here::i_am("paper/section5_2.R")
source(here::here("paper", "utils.R"))

pf <- function(p) if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)
one <- function(pair, title) {
  r <- readRDS(here::here("realdata", "apoe", "results", paste0(pair, ".rds")))
  c(sprintf("%s (n = %d, %d variants)", title, r$n, r$n_variants),
    sprintf("  FPCA scores k = %d, first-stage F = %.1f, partial R^2 = %.3f", r$k, r$F, r$partial_R2),
    sprintf("  FMR        beta = %5.2f, bootstrap-t P = %s, 95%% CI [%.2f, %.2f]", r$fmr$beta, pf(r$fmr$p), r$fmr$ci[1], r$fmr$ci[2]),
    sprintf("  F2SLS      beta = %5.2f, P = %s", r$f2sls$beta, pf(r$f2sls$p)),
    sprintf("  2SLS       beta = %5.2f, P = %s (rs429358, rs7412; first-stage F = %.1f)", r$tsls_e4e2$beta, pf(r$tsls_e4e2$p), r$tsls_e4e2$F),
    sprintf("  sisVIVE    beta = %5.2f, invalid instruments selected: %d", r$sisvive$beta, r$sisvive$n_invalid),
    sprintf("  paired bootstrap: FMR - F2SLS P = %s, FMR - 2SLS P = %s", pf(r$diff$p[1]), pf(r$diff$p[2])),
    sprintf("  T_alpha P = %s; Sargan P = %s", pf(r$talpha$p), pf(r$tsls_e4e2$sargan_p)), "")
}
txt <- c(one("amy_cl_faq", "Amyloid PET -> FAQ"), one("bmi_adas13", "BMI -> ADAS13"))
writeLines(txt, out_file("section5_2.txt")); cat(txt, sep = "\n")
