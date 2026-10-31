# Table 2: the ten gene expression-outcome pairs with the smallest FMR bootstrap p-values,
# and the numbers of Section 5.1 quoted in the text (output/section5_1.txt).
# Input:  realdata/gene_expression/results/{fmr,cis_2sls} (realdata/gene_expression/*.R)
# Output: output/table2.tex, output/section5_1.txt
here::i_am("paper/table2.R")
source(here::here("paper", "utils.R"))
suppressPackageStartupMessages(library(data.table))
GE <- function(...) file.path(Sys.getenv("FMR_GE_RESULTS", here::here("realdata", "gene_expression", "results")), ...)

fmr <- readRDS(GE("fmr", "master.rds"))[status == "ok"]
cis <- readRDS(GE("cis_2sls", "master.rds"))[, .(gene, outcome, beta_2sls = beta, p_2sls = p)]
A <- merge(fmr, cis, by = c("gene", "outcome"), all.x = TRUE)
top <- A[order(p)][1:10]

pv <- function(p) vapply(p, function(x) {
  if (!is.finite(x)) return("--")
  if (x >= 1e-3) return(sprintf("%.3f", x))
  e <- floor(log10(x)); sprintf("$%.1f\\times10^{%d}$", x / 10^e, e)
}, "")
b2 <- function(x) ifelse(is.finite(x), sub("^-", "$-$", sprintf("%.2f", x)), "--")
rows <- top[, sprintf("\\textit{%s} & %s & %d & %.1f & %.3f & %s & %s & %s & %s & %s & %s & %s \\\\",
                      gene, outcome, k, F, partial_R2, b2(beta), pv(p), pv(p_talpha), b2(f2sls_beta), pv(f2sls_p),
                      b2(beta_2sls), pv(p_2sls))]
write_table(out_file("table2.tex"),
  "Top estimated causal effects of gene expression on AD-related outcomes ranked by the FMR bootstrap-based $P$-values.",
  "tab:2", "llrrrrrrrrrr",
  c(" & & \\multicolumn{3}{c}{FIV} & \\multicolumn{3}{c}{FMR} & \\multicolumn{2}{c}{F2SLS} & \\multicolumn{2}{c}{2SLS} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
    "Gene & Outcome & $k$ & $F$ & Partial $R^2$ & $\\hat\\beta$ & $P$ & $P_\\alpha$ & $\\hat\\beta$ & $P$ & $\\hat\\beta$ & $P$ \\\\"),
  rows)

# numbers quoted in Section 5.1
st <- fread(GE("fmr", "gene_status.csv"))
N <- nrow(A); bonf <- 0.05 / N
txt <- c(
  sprintf("genes with expression data and at least 10 variants: %d", nrow(st)),
  sprintf("genes with at least three FPCA scores (analysed):     %d", uniqueN(A$gene)),
  sprintf("gene expression-outcome tests:                        %d", N),
  sprintf("Bonferroni threshold 0.05 / %d = %.2e", N, bonf),
  sprintf("FMR bootstrap p below the threshold:  %d (smallest p = %.1e)", sum(A$p < bonf), min(A$p)),
  sprintf("T_alpha p below the threshold:        %s",
          paste(A[p_talpha < bonf, sprintf("%s-%s (p = %.1e)", gene, outcome, p_talpha)], collapse = "; ")),
  sprintf("ten leading pairs: F %.1f to %.1f, partial R^2 %.3f to %.3f",
          min(top$F), max(top$F), min(top$partial_R2), max(top$partial_R2)))
writeLines(txt, out_file("section5_1.txt")); cat(txt, sep = "\n")
