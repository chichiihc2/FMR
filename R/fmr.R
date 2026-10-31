#' Functional Mendelian randomization
#'
#' Estimates the causal effect \eqn{\beta} of an exposure \code{x} on an outcome \code{y}
#' with a functional instrument \eqn{Z(t)}, e.g. the genotype of the variants of a gene
#' indexed by their genomic position \eqn{t}, allowing a direct (pleiotropic) effect
#' \eqn{\alpha(t)} of the instrument on the outcome:
#' \deqn{Y = \beta X + \int Z(t)\alpha(t)\,dt + \varepsilon.}
#'
#' The instrument is reduced to its first \code{k} functional principal component scores
#' (\code{refund::fpca.face}); \eqn{\alpha(t)} is a cubic B-spline with a roughness penalty
#' whose weight \eqn{\lambda} is chosen by GCV. Inference on \eqn{\beta} uses the pairs
#' bootstrap-\eqn{t} interval studentized by the sandwich standard error
#' \eqn{\hat V^{1/2}}, with \eqn{\lambda} re-chosen in every draw; the Wald interval
#' \eqn{\hat\beta(\lambda)\pm z\,\mathrm{SE}(\tau\lambda)} is also reported.
#' \eqn{H_0:\alpha\equiv 0} is tested by \eqn{T_\alpha}.
#'
#' @param y Outcome, numeric vector of length n.
#' @param x Exposure, numeric vector of length n.
#' @param Z Functional instrument: n x p matrix, one column per position (e.g. allele counts
#'   of the p variants of a gene). Missing entries are replaced by the column mean.
#' @param pos Positions of the p columns of \code{Z} (e.g. base-pair positions); default
#'   \code{1:p}. Rescaled to [0, 1].
#' @param covariates Optional data frame or matrix of covariates (n rows). \code{y}, \code{x}
#'   (and \code{Z}, if \code{adjust_Z}) are replaced by their residuals on an intercept and the
#'   covariates.
#' @param adjust_Z Residualize \code{Z} on the covariates as well (default \code{TRUE}).
#' @param k Number of FPCA scores used as instruments. \code{NULL} (default) chooses it by the
#'   first-stage BIC over \code{k_min:k_max}.
#' @param k_min,k_max Range of \code{k} for the BIC choice (defaults 3 and 10).
#' @param pve Proportion of variance explained that \code{fpca.face} retains (default 0.999).
#' @param knots Number of knots of \code{fpca.face} (default 35, at most p - 4).
#' @param nbasis Number of cubic B-spline functions for \eqn{\alpha(t)} (default 30).
#' @param lambda_grid Grid for the GCV choice of \eqn{\lambda} (default \code{2^(-10:10)}).
#' @param B Bootstrap draws for the bootstrap-\eqn{t} interval (default 199; 0 skips it).
#' @param tau Fraction of \eqn{\lambda} at which the Wald standard error is evaluated
#'   (default 0.1).
#' @param level Confidence level (default 0.95).
#' @param seed Seed of the bootstrap draws; the caller's random-number stream is restored.
#' @param talpha Compute the test of \eqn{H_0:\alpha\equiv0} (default \code{TRUE}).
#' @param compare Also fit F2SLS (2SLS on the same scores; needs \pkg{ivreg}) and OLS.
#' @param n_grid Number of points at which \eqn{\hat\alpha(t)} is returned.
#' @return An object of class \code{"fmr"}; see \code{\link{summary.fmr}} and
#'   \code{\link{plot.fmr}}. Its main elements: \code{coef} (table of estimates, standard
#'   errors, intervals and p-values by method), \code{beta}, \code{lambda}, \code{k},
#'   \code{first_stage} (F and partial R^2), \code{talpha}, \code{alpha} (data frame of
#'   \eqn{\hat\alpha(t)} with a pointwise band).
#' @examples
#' set.seed(1)
#' d <- simulate_fmr_data(n = 1000, s_alpha = 0.3, case = 2)
#' fit <- fmr(d$y, d$x, d$Z, pos = d$pos, B = 49)
#' fit
#' plot(fit)
#' @export
fmr <- function(y, x, Z, pos = NULL, covariates = NULL, adjust_Z = TRUE,
                k = NULL, k_min = 3L, k_max = 10L, pve = 0.999, knots = 35L,
                nbasis = 30L, lambda_grid = 2^(-10:10), B = 199L, tau = 0.1,
                level = 0.95, seed = 1L, talpha = TRUE, compare = TRUE, n_grid = 101L) {
  cl <- match.call()
  y <- as.numeric(y); x <- as.numeric(x); Z <- as.matrix(Z)
  n0 <- length(y); p <- ncol(Z)
  if (length(x) != n0 || nrow(Z) != n0) stop("y, x and Z must have the same number of rows")
  if (is.null(pos)) pos <- seq_len(p)
  if (length(pos) != p) stop("pos must have one entry per column of Z")
  if (anyDuplicated(pos)) stop("pos has duplicated positions")
  if (is.unsorted(pos)) { o <- order(pos); Z <- Z[, o, drop = FALSE]; pos <- pos[o] }
  if (p < 5L) stop("the instrument needs at least 5 positions")

  # ---- complete cases, missing genotypes, covariates --------------------------------
  C <- NULL
  if (!is.null(covariates)) {
    C <- stats::model.matrix(~ ., data = as.data.frame(covariates), na.action = stats::na.pass)
    if (nrow(C) != n0) stop("covariates must have one row per observation")
  }
  ok <- is.finite(y) & is.finite(x) & (if (is.null(C)) TRUE else stats::complete.cases(C))
  y <- y[ok]; x <- x[ok]; Z <- Z[ok, , drop = FALSE]; if (!is.null(C)) C <- C[ok, , drop = FALSE]
  n <- length(y)
  miss <- is.na(Z)
  if (any(miss)) for (j in which(colSums(miss) > 0)) Z[miss[, j], j] <- mean(Z[, j], na.rm = TRUE)
  if (!is.null(C)) {
    rs <- function(M) M - C %*% qr.coef(qr(C), M)
    y <- drop(rs(matrix(y))); x <- drop(rs(matrix(x)))
    if (adjust_Z) Z <- rs(Z)
  }
  t01 <- (pos - min(pos)) / (max(pos) - min(pos))

  # ---- instrument: FPCA scores and k -------------------------------------------------
  fp <- refund::fpca.face(Z, argvals = t01, knots = min(p - 4L, knots), pve = pve)
  K <- fp$npc
  scores <- matrix(fp$scores, nrow = n, ncol = K)
  efun <- matrix(fp$efunctions, nrow = p, ncol = K)
  S <- if (is.null(C)) scores else matrix(stats::resid(stats::lm(scores ~ C - 1)), nrow = n)
  if (is.null(k)) {
    if (K < k_min) stop(sprintf("fpca.face retained %d component(s); FMR needs at least k_min = %d", K, k_min))
    k <- .select_k_bic(x, S[, seq_len(min(K, k_max)), drop = FALSE], floor_k = k_min, cap_k = k_max)
  }
  if (k > K) stop(sprintf("k = %d exceeds the %d components retained at pve = %g", k, K, pve))
  d0 <- if (is.null(C)) 1L else ncol(C)
  r2 <- summary(stats::lm(x ~ S[, seq_len(k), drop = FALSE]))$r.squared
  Fst <- (r2 / k) / ((1 - r2) / (n - k - d0))

  # ---- second stage -------------------------------------------------------------------
  pen <- .make_pen(nbasis)
  Y  <- drop(scale(y, TRUE, FALSE)); X <- drop(scale(x, TRUE, FALSE))
  Xi <- scale(scores[, seq_len(k), drop = FALSE], TRUE, FALSE)
  wq <- c(0, diff(t01))
  Psiw <- fda::eval.basis(t01, pen$BO) * wq
  M2 <- t(efun[, seq_len(k), drop = FALSE]) %*% Psiw
  A  <- Xi %*% M2
  cp0 <- .xprod(Y, cbind(X, A), cbind(1, Xi))
  lam <- .gcv_lambda(cp0, pen$PP_rough, lambda_grid, k)
  if (!is.finite(lam)) stop("GCV found no finite lambda on the grid")
  f0 <- .se_sandwich(cp0, lam, pen$PP_rough)
  if (is.null(f0) || !is.finite(f0$se)) stop("the sandwich variance is undefined at lambda-hat")
  ft <- .se_sandwich(cp0, tau * lam, pen$PP_rough)
  se_tau <- if (is.null(ft)) NA_real_ else ft$se
  z <- stats::qnorm(1 - (1 - level) / 2)
  b <- f0$beta

  # ---- bootstrap-t (pairs; lambda re-chosen by GCV in each draw) ------------------------
  chat <- NA_real_; nb_ok <- 0L
  if (B > 0) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
              assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    ts <- vapply(seq_len(B), function(bb) {
      id  <- sample.int(n, n, replace = TRUE)
      cpb <- tryCatch(.xprod(Y[id], cbind(X[id], A[id, , drop = FALSE]), cbind(1, Xi[id, , drop = FALSE])),
                      error = function(e) NULL)
      if (is.null(cpb)) return(NA_real_)
      lb <- tryCatch(.gcv_lambda(cpb, pen$PP_rough, lambda_grid, k), error = function(e) NA_real_)
      if (!is.finite(lb)) lb <- lam
      sb <- tryCatch(.se_sandwich(cpb, lb, pen$PP_rough), error = function(e) NULL)
      if (!is.null(sb) && is.finite(sb$se) && sb$se > 0) (sb$beta - b) / sb$se else NA_real_
    }, numeric(1))
    ts <- ts[is.finite(ts)]; nb_ok <- length(ts)
    if (nb_ok >= min(50L, B)) chat <- unname(stats::quantile(abs(ts), level))
    else warning("fewer than 50 usable bootstrap draws; the bootstrap-t interval is not reported")
  }

  # ---- alpha-hat(t) and its pointwise band (V-hat at lambda-hat) -------------------------
  G  <- cp0$BPsB + lam * pen$PP_rough; Gi <- solve(G)
  th <- drop(Gi %*% cp0$BPsY)
  V  <- f0$sigma2 * Gi %*% cp0$BPsB %*% Gi
  tg <- seq(0, 1, length.out = n_grid)
  PSI <- fda::eval.basis(tg, pen$BO)
  ah <- drop(PSI %*% th[-1])
  sd_a <- sqrt(pmax(rowSums((PSI %*% V[-1, -1]) * PSI), 0))
  alpha <- data.frame(t = tg, pos = min(pos) + tg * (max(pos) - min(pos)),
                      alpha = ah, lower = ah - z * sd_a, upper = ah + z * sd_a)

  ta <- if (talpha) .fmr_talpha(Y, X, Xi, A, lam, pen$PP_rough, pen$BO, tau) else c(T = NA_real_, p = NA_real_)

  # ---- table ----------------------------------------------------------------------------
  row <- function(method, est, se, lo, hi, pval) data.frame(method = method, estimate = est, se = se,
                                                         lower = lo, upper = hi, p = pval)
  tab <- rbind(
    row("FMR, bootstrap-t", b, chat * f0$se / z, b - chat * f0$se, b + chat * f0$se,
        2 * stats::pnorm(-z * abs(b) / (chat * f0$se))),
    row(sprintf("FMR, Wald (tau = %g)", tau), b, se_tau, b - z * se_tau, b + z * se_tau,
        2 * stats::pnorm(-abs(b) / se_tau)))
  if (compare) {
    if (requireNamespace("ivreg", quietly = TRUE)) {
      s2 <- summary(ivreg::ivreg(Y ~ -1 + X | Xi, data = data.frame(Y = Y, X = X, Xi = I(unname(Xi)))))$coefficients
      tab <- rbind(tab, row("F2SLS", s2["X", 1], s2["X", 2], s2["X", 1] - z * s2["X", 2],
                            s2["X", 1] + z * s2["X", 2], 2 * stats::pnorm(-abs(s2["X", 1] / s2["X", 2]))))
    }
    so <- summary(stats::lm(Y ~ -1 + X))$coefficients
    tab <- rbind(tab, row("OLS", so["X", 1], so["X", 2], so["X", 1] - z * so["X", 2],
                          so["X", 1] + z * so["X", 2], 2 * stats::pnorm(-abs(so["X", 1] / so["X", 2]))))
  }
  rownames(tab) <- NULL

  structure(list(call = cl, coef = tab, beta = b, se_vhat = f0$se, c_hat = chat, n_boot = nb_ok,
                 se_tau = se_tau, lambda = lam, df = f0$df, theta = th, k = k, n = n, p = p,
                 K_retained = K, first_stage = c(F = Fst, partial_R2 = r2),
                 talpha = c(T_alpha = unname(ta["T"]), p = unname(ta["p"])), alpha = alpha,
                 level = level, tau = tau, B = B),
            class = "fmr")
}
