# Numerical core of FMR: the same computations as the paper's simulation and data-analysis code
# (lib/estimators.R, lib/fmr_talpha.R, lib/gene_fpca.R in the repository);
# tests/testthat/test-reproduce.R checks the package against stored replications of the
# simulation study.

# Cubic B-spline sieve on [0, 1] for alpha(t) and its roughness penalty (beta unpenalized).
.make_pen <- function(nbasis) {
  BO <- fda::create.bspline.basis(c(0, 1), nbasis)
  P0 <- fda::bsplinepen(BO)
  list(BO = BO, nbasis = nbasis, P0 = P0, PP_rough = rbind(0, cbind(0, P0)))
}

# Cross-products of the second stage, Q = [1, Xi], B = [X, A]; every lambda then costs one solve.
.xprod <- function(Y, B, Q, tol = 1e-15) {
  Si <- MASS::ginv(crossprod(Q), tol = tol)
  M  <- crossprod(Q, B); v <- crossprod(Q, Y)
  list(Si = Si, M = M, v = v, BPsB = t(M) %*% Si %*% M, BPsY = t(M) %*% Si %*% v,
       Y = Y, B = B, Q = Q)
}

# lambda minimizing GCV(lambda) = MD ||P_Q(Y - B theta)||^2 / (MD - df)^2 over the grid LG.
.gcv_lambda <- function(cp, PP, LG, MD) {
  sc <- vapply(LG, function(lam) {
    G  <- cp$BPsB + lam * PP
    th <- tryCatch(drop(solve(G, cp$BPsY)), error = function(e) NULL)
    if (is.null(th)) return(Inf)
    df <- sum(diag(solve(G, cp$BPsB)))
    Qe <- cp$v - cp$M %*% th
    mse <- abs(drop(crossprod(Qe, cp$Si %*% Qe)))
    if (!is.finite(mse) || abs(MD - df) < 1e-8) return(Inf)
    MD * mse / (MD - df)^2
  }, numeric(1))
  if (all(!is.finite(sc))) NA_real_ else LG[which.min(sc)]
}

# beta-hat(lambda) and its V-hat standard error.
.se_sandwich <- function(cp, lambda, PP) {
  G  <- cp$BPsB + lambda * PP
  th <- drop(solve(G, cp$BPsY))
  e  <- drop(cp$Y - cp$B %*% th)
  df <- sum(diag(solve(G, cp$BPsB)))
  nn <- length(cp$Y)
  if (!is.finite(df) || nn - df <= 0) return(NULL)
  s2 <- sum(e^2) / (nn - df)
  u  <- drop(solve(G, c(1, rep(0, ncol(cp$B) - 1L))))
  v  <- drop(crossprod(u, cp$BPsB %*% u))
  list(beta = unname(th[1]), df = df, sigma2 = s2,
       se = if (is.finite(v) && v > 0 && is.finite(s2)) sqrt(s2 * v) else NA_real_)
}

# T_alpha for H0: alpha = 0 (step-2 GMM at lambda, sandwich at tau * lambda, Wood approximation).
.fmr_talpha <- function(Y, X, Xi, A, lambda, PP, basisobj, tau = 0.1, tol = 1e-15) {
  Bmat <- cbind(X, A); Q <- cbind(1, Xi)
  Si <- MASS::ginv(crossprod(Q), tol = tol)
  M  <- crossprod(Q, Bmat)
  th1 <- drop(solve(t(M) %*% Si %*% M + lambda * PP, t(M) %*% Si %*% crossprod(Q, Y)))
  e1 <- drop(Y - Bmat %*% th1)
  Wi <- MASS::ginv(crossprod(Xi * e1), tol = tol)
  G1 <- crossprod(Bmat, Xi); XtB <- t(G1)
  Gm <- G1 %*% Wi %*% XtB
  rhs <- G1 %*% Wi %*% crossprod(Xi, Y)
  th2 <- drop(tryCatch(solve(Gm + lambda * PP, rhs), error = function(e) MASS::ginv(Gm + lambda * PP, tol = tol) %*% rhs))
  e2 <- drop(Y - Bmat %*% th2)
  mid <- G1 %*% Wi %*% crossprod(Xi * e2) %*% Wi %*% XtB
  Ainv <- MASS::ginv(Gm + tau * lambda * PP, tol = tol)
  Vd <- (Ainv %*% mid %*% Ainv)[-1, -1, drop = FALSE]
  M0 <- fda::bsplinepen(basisobj, Lfdobj = 0)
  dl <- th2[-1]; Tst <- drop(crossprod(dl, M0 %*% dl))
  D2 <- M0 %*% Vd
  ll <- eigen(D2 %*% t(D2), symmetric = TRUE, only.values = TRUE)$values
  ll[ll < 0] <- 0
  p <- tryCatch(dr::dr.pvalue(sqrt(ll), Tst, chi2approx = "wood")$pval.adj, error = function(e) NA_real_)
  if (is.na(p)) p <- 1
  c(T = Tst, p = p)
}

# k by first-stage BIC along the eigenvalue-ordered prefix of the scores, floor k_min, cap k_max.
.select_k_bic <- function(x, S, floor_k = 3L, cap_k = 10L) {
  n <- length(x); Kmax <- min(ncol(S), cap_k, n - 5L)
  if (Kmax < floor_k) return(max(1L, Kmax))
  bic <- rep(NA_real_, Kmax)
  for (m in floor_k:Kmax) {
    rss <- sum(stats::lm.fit(cbind(1, S[, seq_len(m), drop = FALSE]), x)$residuals^2)
    bic[m] <- n * log(rss / n) + log(n) * (m + 1)
  }
  which.min(replace(bic, is.na(bic), Inf))
}
