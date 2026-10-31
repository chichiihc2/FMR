# Reading the ADNI data of one gene: the genotype matrix of its variants (adjusted for six
# genotype PCs; data/processed/Genes_pca6/<gene>.csv, header = variant positions, rows in
# PTID.csv order), its expression, and the subjects' age and sex.
suppressPackageStartupMessages(library(data.table))

genes_pca6_path <- function(gene, dir = here::here("data", "processed", "Genes_pca6")) {
  file.path(dir, paste0(gene, ".csv"))
}

# n x p genotype matrix with the positions as column names; missing values are set to the
# column mean. NULL if the file is missing or has fewer than two distinct positions.
read_genes_pca6 <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  dt <- tryCatch(data.table::fread(file_path, header = TRUE, showProgress = FALSE, data.table = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || !NCOL(dt) || !NROW(dt)) return(NULL)
  Z <- as.matrix(dt)
  storage.mode(Z) <- "double"
  na_idx <- which(is.na(Z), arr.ind = TRUE)
  if (NROW(na_idx) > 0) {
    cm <- colMeans(Z, na.rm = TRUE)
    cm[!is.finite(cm)] <- 0
    Z[na_idx] <- cm[na_idx[, "col"]]
  }
  pos <- suppressWarnings(as.numeric(colnames(Z)))
  if (anyNA(pos) || length(pos) < 2L || diff(range(pos)) == 0) return(NULL)
  Z
}

# Complete-case data of one gene: Z (n x p), positions POS and t01 (scaled to [0, 1]),
# expression Xexp, and df0 (PTID, Age, Gender). NULL if fewer than min_snps variants or
# min_obs subjects.
load_gene_frame <- function(gene, demo, ptid, expr_dir, min_snps = 5L, min_obs = 50L) {
  G <- read_genes_pca6(genes_pca6_path(gene))
  if (is.null(G)) return(NULL)
  POS <- as.numeric(colnames(G))
  expr_file <- file.path(expr_dir, paste0(gene, ".csv"))
  if (!file.exists(expr_file)) return(NULL)
  df0  <- as.data.frame(demo)[, c("PTID", "Age", "Gender")]
  expr <- as.data.frame(data.table::fread(expr_file))
  df0$value <- expr$value[match(df0$PTID, expr$PTID)]
  zi <- match(df0$PTID, ptid$PTID)
  df0$Z <- G[zi, , drop = FALSE]
  keep <- stats::complete.cases(df0[, setdiff(names(df0), "Z"), drop = FALSE]) &
          !is.na(zi) & stats::complete.cases(df0$Z)
  df0 <- df0[keep, , drop = FALSE]
  n_obs <- nrow(df0); n_snps <- ncol(G)
  if (n_snps < min_snps || n_obs < min_obs) return(NULL)
  Z <- as.matrix(df0$Z)
  if (!is.matrix(Z) || ncol(Z) != n_snps) Z <- matrix(unlist(df0$Z), nrow = n_obs)
  list(Z = Z, POS = POS, Xexp = as.numeric(df0$value), df0 = df0,
       t01 = (POS - min(POS)) / diff(range(POS)), n_obs = n_obs, n_snps = n_snps)
}
