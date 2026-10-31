# Estimators shared by the simulation study and the real-data analysis.
#
# FMR (Section 3): with Xi the n x k centred FPCA scores of Z(t), psi a cubic B-spline basis of
# dimension q for alpha(t) and A = Xi M2 (M2 = int phi(t) psi(t)' dt), the estimator is
#   theta = (beta, delta) = argmin (Y - B theta)' P_Q (Y - B theta) + lambda delta' P delta,
# B = [X, A], Q = [1, Xi], P the roughness penalty int psi'' psi''' (P = I for Ridge).
# lambda is chosen by GCV; SE(beta) is the sandwich of the penalized 2SLS fit.

ZQ <- qnorm(0.975)
LG_RIDGE <- 2^(-30:10)          # lambda grid of Ridge

# B-spline basis of alpha(t) on [0, 1] with roughness (P0) and ridge (I0) penalties;
# PP_* pad the penalty with a zero row/column for beta.
make_pen <- function(nbasis) {
  BO <- fda::create.bspline.basis(c(0, 1), nbasis)
  P0 <- fda::bsplinepen(BO)
  I0 <- diag(1, nbasis)
  list(BO = BO, nbasis = nbasis, P0 = P0, I0 = I0,
       PP_rough = rbind(0, cbind(0, P0)), PP_ridge = rbind(0, cbind(0, I0)))
}

# Cross-products of the penalized 2SLS problem.
xprod <- function(Y, B, Q, tol = 1e-15) {
  Si <- MASS::ginv(crossprod(Q), tol = tol)
  M  <- crossprod(Q, B); v <- crossprod(Q, Y)
  list(Si = Si, M = M, v = v, BPsB = t(M) %*% Si %*% M, BPsY = t(M) %*% Si %*% v,
       Y = Y, B = B, Q = Q)
}

# lambda minimizing GCV(lambda) = MD * ||P_Q (Y - B theta)||^2 / (MD - df)^2 over the grid LG.
gcv_lambda <- function(cp, PP, LG, MD) {
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

# Fit at lambda and the sandwich SE of beta-hat.
se_sandwich <- function(cp, lambda, PP) {
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

# The same SE computed through a QR decomposition, numerically stable for the small
# lambda values that GCV selects for Ridge.
se_sandwich_qr <- function(cp, lambda, PP) {
  s <- se_sandwich(cp, lambda, PP)
  if (is.null(s)) return(NULL)
  es <- eigen(cp$Si, symmetric = TRUE)
  R  <- (es$vectors %*% (sqrt(pmax(es$values, 0)) * t(es$vectors))) %*% cp$M
  pr <- which(diag(PP) > 0)
  A  <- rbind(R, sqrt(lambda) * PP[pr, , drop = FALSE])
  qa <- qr(A, LAPACK = TRUE)
  Ra <- qr.R(qa); Qa <- qr.Q(qa)
  e1 <- c(1, rep(0, ncol(A) - 1L))[qa$pivot]
  w  <- forwardsolve(t(Ra), e1)
  v  <- sum((Qa[seq_len(nrow(R)), , drop = FALSE] %*% w)^2)
  s$se <- if (is.finite(v) && v > 0 && is.finite(s$sigma2)) sqrt(s$sigma2 * v) else NA_real_
  s
}

# Design of one simulated data set D (from generate_S2 / generate_S1): centred Y, X and
# scores, and A = Xi M2 with M2 = int phi(t) psi(t)' dt by a Riemann sum over the positions.
design <- function(D, BO, n_scores) {
  Y  <- drop(scale(D$Y, TRUE, FALSE)); X <- drop(scale(D$X, TRUE, FALSE))
  Xi <- scale(D$Xi[, seq_len(n_scores)], TRUE, FALSE)
  efs <- D$fpca$efunctions[, seq_len(n_scores)]
  POS01 <- (D$POS - min(D$POS)) / (max(D$POS) - min(D$POS))
  wq <- POS01 - dplyr::lag(POS01); wq[1] <- 0
  Psiw <- fda::eval.basis(POS01, BO) * wq
  M2 <- t(efs) %*% Psiw
  list(Y = Y, X = X, Xi = Xi, A = Xi %*% M2)
}

# First-stage regression of X on the scores: c(R^2, F).
first_stage_F <- function(d) {
  sm <- summary(lm(d$X ~ d$Xi))
  c(unname(sm$r.squared), unname(sm$fstatistic[1]))
}

# F2SLS: 2SLS with the FPCA scores as instruments (ivreg); c(beta-hat, SE).
f2sls_of <- function(d) {
  fit <- ivreg::ivreg(Y ~ -1 + X | Xi, data = tibble::tibble(Y = as.numeric(d$Y), X = as.numeric(d$X),
                                                            Xi = unname(as.matrix(d$Xi))))
  c(unname(fit$coefficients[["X"]]), unname(summary(fit)$coefficients["X", "Std. Error"]))
}

# sisVIVE with the FPCA scores (or other instruments Z) as candidate instruments; 10-fold CV.
sisvive_of <- function(Y, X, Z, i) {
  set.seed(9e6 + i)                       # cv.sisVIVE draws its folds with sample()
  unname(sisVIVE::cv.sisVIVE(Y = as.numeric(Y), D = as.numeric(X), Z = unname(as.matrix(Z)),
                             K = 10, intercept = FALSE)$beta)
}

# LS: regression of Y on X; c(beta-hat, SE).
ls_of <- function(d) {
  s <- summary(stats::lm(d$Y ~ -1 + d$X))$coefficients
  c(unname(s[1, 1]), unname(s[1, 2]))
}
