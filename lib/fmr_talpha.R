# Test of H0: alpha(t) = 0 (Section 3, statistic T_alpha = ||alpha-hat||^2).
# delta-hat is the efficient (two-step) GMM estimate at the given lambda; under H0, T_alpha is
# approximately sum_j w_j chi^2_1 with w the eigenvalues of (int psi psi') V_delta, V_delta the
# robust sandwich with the penalty at tau * lambda; the p-value uses Wood's approximation.
#   Y, X: centred outcome and exposure; Xi: n x k centred FPCA scores; A = Xi M2;
#   PP: (1 + q) x (1 + q) penalty (0 for beta); basisobj: B-spline basis of alpha(t).
fmr_talpha <- function(Y, X, Xi, A, lambda, PP, basisobj, tau = 0.1, tol = 1e-15) {
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
