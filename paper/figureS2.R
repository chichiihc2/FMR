# Figure S.2: alpha-hat(t) of the APOE locus with 95% bootstrap simultaneous bands;
# (a) amyloid PET -> FAQ, (b) BMI -> ADAS13.
# Input:  realdata/apoe/results/{amy_cl_faq, bmi_adas13}.rds (realdata/apoe/fmr_apoe.R)
# Output: output/figureS2a.pdf, output/figureS2b.pdf
here::i_am("paper/figureS2.R")
source(here::here("paper", "utils.R"))

GENES <- data.frame(gene = c("NECTIN2", "TOMM40", "APOE", "APOC1"),               # GRCh37
                    s = c(45349432, 45394477, 45409039, 45417504) / 1e6,
                    e = c(45392485, 45406946, 45412650, 45422606) / 1e6)
plot_alpha <- function(pair, file) {
  b <- readRDS(here::here("realdata", "apoe", "results", paste0(pair, ".rds")))$band
  ad <- data.frame(bp = b$bp / 1e6, a = b$alpha, lo = b$lo, hi = b$hi)
  genes <- GENES %>% mutate(s = pmax(s, min(ad$bp)), e = pmin(e, max(ad$bp)))
  g <- ggplot(ad, aes(bp, a)) +
    geom_rect(data = genes, aes(xmin = s, xmax = e, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
              fill = "steelblue", alpha = .12) +
    geom_text(data = genes, aes(x = (s + e) / 2, y = max(ad$hi), label = gene), inherit.aes = FALSE,
              size = 3, vjust = 1) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") + geom_line(linewidth = .8) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "chr19 position (Mb, hg19)", y = expression(hat(alpha)(t))) + theme_bw(10)
  ggsave(file, g, width = 7, height = 3.4)
  message("wrote ", file)
}
plot_alpha("amy_cl_faq", out_file("figureS2a.pdf"))
plot_alpha("bmi_adas13", out_file("figureS2b.pdf"))
