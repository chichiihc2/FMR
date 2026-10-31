# Scenario S2 (Section 4.1): the instrument Z(t) is simulated on nT equally spaced points of
# [0, 1] from a sine-cosine expansion with eigenvalues 1/l (l <= 10) and 1/l^2 (l > 10).
# The exposure loads on the first 11 FPCA scores with strength `signal` (s), and
#   Y = 2 + beta X + int Z(t) alpha(t) dt + eps,   corr(eps_X, eps) = rho.

generate_functional_data <- function(n, nT) {
  POS <- seq(0, 1, length.out = nT)
  sig2 <- numeric(100)
  for (l in 1:10) sig2[l] <- 1 / l
  for (l in 11:100) sig2[l] <- 1 / (l^2)
  Xi.Sin <- matrix(rnorm(n * 100), n, 100)
  Xi.Cos <- matrix(rnorm(n * 100), n, 100)
  for (l in 1:100) {
    Xi.Sin[, l] <- Xi.Sin[, l] * sqrt(sig2[l])
    Xi.Cos[, l] <- Xi.Cos[, l] * sqrt(sig2[l])
  }
  Phi.Sin <- sapply(1:100, function(l) sqrt(2) * sin(l * pi * POS))
  Phi.Cos <- sapply(1:100, function(l) sqrt(2) * cos(l * pi * POS))
  list(POS = POS, Z = Xi.Sin %*% t(Phi.Sin) + Xi.Cos %*% t(Phi.Cos))
}

# Direct effect alpha(t) of Cases 1-6, with size s_alpha = Dsignal.
alpha_S2 <- function(case, Dsignal) {
  switch(as.character(case),
    `1` = function(t) 2 * Dsignal * t,
    `2` = function(t) Dsignal * (t + sin(pi * t)),
    `3` = function(t) 2 * Dsignal * (1 / (1 + exp(-8 * (t - 0.5)))),
    `4` = function(t) 5 * Dsignal * ifelse(t < 0.5, t, 1 - t),
    `5` = function(t) 1.5 * ifelse((t - 0.2) * (t - 0.8) < 0, Dsignal, 0),
    `6` = function(t) 2.1967 * Dsignal * ifelse(t > 0.25 & t < 0.75, sin(pi * (t - 0.25) / 0.5), 0),
    stop("case must be 1-6"))
}

generate_S2 <- function(n, nT = 501, effect = 1, signal = 0.015, Dsignal = 0, s2 = 1,
                        case = 1, coef_PC1 = 2, rho = 0.8) {
  alpha_t <- alpha_S2(case, Dsignal)
  data.Z <- generate_functional_data(n = n, nT = nT)
  POS <- data.Z$POS
  fpca0 <- refund::fpca.face(data.Z$Z, argvals = POS, knots = 35, pve = 0.99)
  Xi <- fpca0$scores
  error <- mvtnorm::rmvnorm(n = n, mean = c(0, 0), sigma = s2 * matrix(c(1, rho, rho, 1), 2, 2))
  gamma <- signal * c(coef_PC1, 2, -1, 1, -1, 1, -1, 1, -1, 1, -1)
  X <- Xi[, 1:11] %*% gamma + error[, 1]
  weight <- POS - dplyr::lag(POS)
  weight[1] <- 0
  alpha_coef <- t(fpca0$efunctions) %*% (sapply(POS, alpha_t) * weight)
  Y <- 2 + X * effect + Xi %*% alpha_coef + error[, 2]
  list(Y = Y, X = X, Xi = Xi, Z = data.Z$Z, POS = POS, fpca = fpca0, alpha_t = alpha_t)
}
