#' Simulate data from the paper's Scenario S2
#'
#' The functional instrument is \eqn{Z(t)=\sum_l \xi_l\phi_l(t)} on a grid of \code{nT}
#' points (sine and cosine basis, eigenvalues \eqn{1/l} for \eqn{l\le10} and \eqn{1/l^2}
#' after); the exposure loads on its first 11 FPCA scores, and the outcome is
#' \deqn{Y = 2 + \beta X + \int Z(t)\alpha(t)\,dt + \varepsilon,}
#' with \eqn{\mathrm{corr}(\varepsilon_X,\varepsilon)=\rho} (unmeasured confounding).
#' The code is the paper's generator, so \code{set.seed(i)} followed by this call
#' reproduces replication \code{i} of the simulation study.
#'
#' @param n Sample size.
#' @param nT Number of positions of \eqn{Z(t)}.
#' @param beta Causal effect.
#' @param signal Instrument strength \eqn{s} (paper: 0.015).
#' @param s_alpha Size of the direct effect \eqn{s_\alpha} (0 = no pleiotropy).
#' @param case Shape of \eqn{\alpha(t)}, numbered as in the paper: 1 linear, 2 linear plus
#'   sine, 3 logistic, 4 triangle, 5 step, 6 localized half-wave on \eqn{[0.25, 0.75]}.
#' @param rho Correlation of the exposure and outcome errors.
#' @return A list with \code{y}, \code{x}, \code{Z} (n x nT), \code{pos} (positions in [0, 1]),
#'   \code{alpha} (the true \eqn{\alpha(t)} as a function) and \code{beta}.
#' @examples
#' set.seed(1)
#' d <- simulate_fmr_data(n = 500)
#' str(d[c("y", "x", "pos")])
#' @export
simulate_fmr_data <- function(n = 1000, nT = 501, beta = 1, signal = 0.015, s_alpha = 0.3,
                              case = 1, rho = 0.8) {
  if (!case %in% 1:6) stop("case must be 1-6")
  D <- .generate_data_S2(n = n, nT = nT, effect = beta, signal = signal, Dsignal = s_alpha,
                         s2 = 1, case = case, coef_PC1 = 2, rho = rho)
  list(y = drop(D$Y), x = drop(D$X), Z = D$Z, pos = D$POS, alpha = D$alpha_t, beta = beta)
}

# The generator of simulations/helpers/generate_S2.R in the repository.
.generate_functional_data <- function(n, nT) {
  POS <- seq(0, 1, length.out = nT)
  sig2 <- numeric(100)
  for (l in 1:10) sig2[l] <- 1 / l
  for (l in 11:100) sig2[l] <- 1 / (l^2)
  Xi.Sin <- matrix(stats::rnorm(n * 100), n, 100)
  Xi.Cos <- matrix(stats::rnorm(n * 100), n, 100)
  for (l in 1:100) {
    Xi.Sin[, l] <- Xi.Sin[, l] * sqrt(sig2[l])
    Xi.Cos[, l] <- Xi.Cos[, l] * sqrt(sig2[l])
  }
  Phi.Sin <- sapply(1:100, function(l) sqrt(2) * sin(l * pi * POS))
  Phi.Cos <- sapply(1:100, function(l) sqrt(2) * cos(l * pi * POS))
  Z <- Xi.Sin %*% t(Phi.Sin) + Xi.Cos %*% t(Phi.Cos)
  list(POS = POS, Xi = cbind(Xi.Sin, Xi.Cos), Z = Z, sig2 = sig2)
}

.generate_data_S2 <- function(n, nT, effect, signal, Dsignal, s2, case, coef_PC1, rho) {
  Dsignals <- switch(as.character(case),
    `1` = function(t) 2 * Dsignal * t,
    `2` = function(t) Dsignal * (t + sin(pi * t)),
    `3` = function(t) 2 * Dsignal * (1 / (1 + exp(-8 * (t - 0.5)))),
    `4` = function(t) 5 * Dsignal * ifelse(t < 0.5, t, 1 - t),
    `5` = function(t) 1.5 * ifelse((t - 0.2) * (t - 0.8) < 0, Dsignal, 0),
    `6` = function(t) 2.1967 * Dsignal * ifelse(t > 0.25 & t < 0.75, sin(pi * (t - 0.25) / 0.5), 0))
  data.Z <- .generate_functional_data(n = n, nT = nT)
  POS <- data.Z$POS
  fpca0 <- refund::fpca.face(data.Z$Z, argvals = POS, knots = 35, pve = 0.99)
  Xi <- fpca0$scores
  auto_corr_matrix <- matrix(c(1.0, rho, rho, 1.0), nrow = 2, byrow = TRUE)
  error <- mvtnorm::rmvnorm(n = n, mean = c(0, 0), sigma = s2 * auto_corr_matrix)
  gamma <- signal * c(coef_PC1, 2, -1, 1, -1, 1, -1, 1, -1, 1, -1)
  R11 <- diag(11)
  X <- Xi[, 1:11] %*% (R11 %*% gamma) + error[, 1]
  weight <- POS - c(NA, POS[-length(POS)])
  weight[1] <- 0
  alpha_coef <- t(fpca0$efunctions) %*% (sapply(POS, Dsignals) * weight)
  functional.term <- Xi %*% alpha_coef
  Y <- 2 + X * effect + functional.term + error[, 2]
  list(Y = Y, X = X, Xi = Xi, Z = data.Z$Z, POS = POS, fpca = fpca0, alpha_t = Dsignals)
}
