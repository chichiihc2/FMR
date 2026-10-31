# FPCA of a gene's genotype function and the choice of the number k of FPCA scores
# (Section 5): fpca.face with min(p - 4, 35) knots, components up to 99.9% of the variance
# (at most 10); k by the BIC of the first-stage regression over k = 3, ..., 10.
suppressPackageStartupMessages(library(refund))

gene_fpca <- function(Z, t01, pve = 0.999, cap = 10L) {
  p <- ncol(Z); n <- nrow(Z)
  fp <- tryCatch(refund::fpca.face(Z, argvals = t01, knots = min(p - 4L, 35L), pve = pve),
                 error = function(e) NULL)
  if (is.null(fp)) return(NULL)
  K <- fp$npc
  sc  <- matrix(fp$scores,     nrow = n, ncol = K)
  efn <- matrix(fp$efunctions, nrow = p, ncol = K)
  Kk <- min(K, cap)
  list(scores = sc[, seq_len(Kk), drop = FALSE], efunctions = efn[, seq_len(Kk), drop = FALSE],
       evalues = fp$evalues, K_pool = K, K = Kk)
}

# k minimizing n log(RSS_k / n) + log(n) (k + 1) of the regression of x on the first k
# columns of S, over k = floor_k, ..., min(cap_k, ncol(S)).
select_k_bic <- function(x, S, floor_k = 3L, cap_k = 10L) {
  n <- length(x); Kmax <- min(ncol(S), cap_k, n - 5L)
  if (Kmax < floor_k) return(max(1L, Kmax))
  bic <- rep(NA_real_, Kmax)
  for (m in floor_k:Kmax) {
    rss <- sum(stats::lm.fit(cbind(1, S[, seq_len(m), drop = FALSE]), x)$residuals^2)
    bic[m] <- n * log(rss / n) + log(n) * (m + 1)
  }
  which.min(replace(bic, is.na(bic), Inf))
}
