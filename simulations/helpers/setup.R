# Common code of the simulation scripts: settings, data makers for S1 and S2, one
# replication of all methods, and a runner over the cells of a design.
#
# Seeds: replication i uses set.seed(i) for the data, set.seed(5e6 + i) for the bootstrap
# and set.seed(9e6 + i) for sisVIVE's folds, so a replication gives the same numbers in
# every script and on any number of cores.
#
# Environment variables: FMR_SIM_RESULTS (results folder, default simulations/results),
# FMR_NREPS (replications per cell, default as in the paper), FMR_CORES (default: all).
suppressPackageStartupMessages({
  library(parallel); library(dplyr); library(tidyr); library(tibble)
  library(fda); library(refund); library(mvtnorm); library(MASS)
  library(ivreg); library(sisVIVE); library(data.table)
})
source(here::here("lib", "estimators.R"))
source(here::here("simulations", "helpers", "generate_S2.R"))
source(here::here("simulations", "helpers", "generate_S1.R"))

results_dir <- function(...) file.path(Sys.getenv("FMR_SIM_RESULTS", here::here("simulations", "results")), ...)
n_reps  <- function(default) as.integer(Sys.getenv("FMR_NREPS", default))
n_cores <- function() as.integer(Sys.getenv("FMR_CORES", parallel::detectCores()))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "", file = stderr())

# Settings of Section 4.1.
SET <- list(beta = 1, nT = 501L, n_scores = 11L, nbasis = 30L, max_df = 11, s2 = 1, rho = 0.8,
            coef_PC1 = 2, tol = 1e-15, n_boot = 199L,
            lambda_S2 = 2^(-10:10), lambda_S1 = 2^(-16:16), taus = c(0.05, 0.1, 0.2))

# ---- data makers: replication i of a cell -> design list (Y, X, Xi, A) --------------------
make_S2 <- function(i, n, case, Dsignal, signal, pen) {
  set.seed(i)
  D <- generate_S2(n = n, nT = SET$nT, effect = SET$beta, signal = signal, Dsignal = Dsignal,
                   s2 = SET$s2, case = case, coef_PC1 = SET$coef_PC1, rho = SET$rho)
  design(D, pen$BO, SET$n_scores)
}

S1 <- NULL                      # the DLGAP2 design, built once by s1_init()
s1_init <- function(pen) {
  des <- s1_build_design(gene = "DLGAP2", pve = 0.80, knots = 35, n_scores = SET$n_scores)
  gen <- make_generate_S1(des)
  S1 <<- list(des = des, gen = gen, dfix = design(gen(), pen$BO, SET$n_scores))
  invisible(S1)
}
make_S1 <- function(i, n, case, Dsignal, signal, pen) {
  set.seed(i)
  D <- S1$gen(effect = SET$beta, signal = signal, Dsignal = Dsignal, s2 = SET$s2, case = case,
              coef_PC1 = SET$coef_PC1, rho = SET$rho)
  d <- S1$dfix
  d$Y <- drop(scale(D$Y, TRUE, FALSE)); d$X <- drop(scale(D$X, TRUE, FALSE))
  d
}

# ---- one replication --------------------------------------------------------------------
# FMR: lambda by GCV, beta-hat and its sandwich SE, the SE at tau * lambda (Wald intervals
# with tau = 0.05, 0.1, 0.2) and, with boot = TRUE, the bootstrap-t critical value c-hat
# (95% quantile of |t*| over B = 199 pairs-bootstrap draws, lambda re-chosen by GCV in each
# draw). The same for Ridge (ridge = TRUE); F2SLS, sisVIVE and LS (others = TRUE).
sim_rep <- function(i, n, case, Dsignal, signal, make_data, pen, LG,
                    boot = TRUE, ridge = TRUE, others = TRUE) {
  row <- tibble(i = i, n = n, case = case, Dsignal = Dsignal, signal = signal, nbasis = pen$nbasis)
  tryCatch({
    d   <- make_data(i, n, case, Dsignal, signal, pen)
    cp0 <- xprod(d$Y, cbind(d$X, d$A), cbind(1, d$Xi), SET$tol)
    lf  <- gcv_lambda(cp0, pen$PP_rough, LG, SET$max_df)
    if (!is.finite(lf)) stop("no finite FMR lambda")
    f0  <- se_sandwich(cp0, lf, pen$PP_rough)
    if (is.null(f0) || !is.finite(f0$se)) stop("FMR sandwich undefined")
    se_tau <- vapply(SET$taus, function(t) { s <- se_sandwich(cp0, t * lf, pen$PP_rough)
                                             if (is.null(s)) NA_real_ else s$se }, numeric(1))
    r0 <- NULL
    if (ridge) {
      lr <- gcv_lambda(cp0, pen$PP_ridge, LG_RIDGE, SET$max_df)
      r0 <- if (is.finite(lr)) se_sandwich_qr(cp0, lr, pen$PP_ridge) else NULL
    }
    cvf <- cvr <- NA_real_
    if (boot) {
      set.seed(5e6 + i)
      tt <- vapply(seq_len(SET$n_boot), function(b) {
        id  <- sample.int(n, n, replace = TRUE)
        cpb <- tryCatch(xprod(d$Y[id], cbind(d$X[id], d$A[id, , drop = FALSE]),
                              cbind(1, d$Xi[id, , drop = FALSE]), SET$tol), error = function(e) NULL)
        if (is.null(cpb)) return(c(NA_real_, NA_real_))
        lb <- tryCatch(gcv_lambda(cpb, pen$PP_rough, LG, SET$max_df), error = function(e) NA_real_)
        if (!is.finite(lb)) lb <- lf
        sb <- tryCatch(se_sandwich(cpb, lb, pen$PP_rough), error = function(e) NULL)
        tf <- if (!is.null(sb) && is.finite(sb$se) && sb$se > 0) (sb$beta - f0$beta) / sb$se else NA_real_
        tr <- NA_real_
        if (!is.null(r0) && is.finite(r0$se)) {
          lbr <- tryCatch(gcv_lambda(cpb, pen$PP_ridge, LG_RIDGE, SET$max_df), error = function(e) NA_real_)
          if (!is.finite(lbr)) lbr <- lr
          sr <- tryCatch(se_sandwich_qr(cpb, lbr, pen$PP_ridge), error = function(e) NULL)
          if (!is.null(sr) && is.finite(sr$se) && sr$se > 0) tr <- (sr$beta - r0$beta) / sr$se
        }
        c(tf, tr)
      }, numeric(2))
      q95 <- function(x) { x <- x[is.finite(x)]; if (length(x) < 50L) NA_real_ else unname(quantile(abs(x), 0.95)) }
      cvf <- q95(tt[1, ]); cvr <- q95(tt[2, ])
    }
    row <- mutate(row, Fstat = first_stage_F(d)[2], lambda = lf, beta = f0$beta, se = f0$se,
                  se_tau005 = se_tau[1], se_tau01 = se_tau[2], se_tau02 = se_tau[3], c_hat = cvf)
    if (ridge) row <- mutate(row, ridge_lambda = if (is.null(r0)) NA_real_ else lr,
                             ridge_beta = if (is.null(r0)) NA_real_ else r0$beta,
                             ridge_se = if (is.null(r0)) NA_real_ else r0$se, ridge_c_hat = cvr)
    if (others) {
      fq <- tryCatch(f2sls_of(d), error = function(e) c(NA_real_, NA_real_))
      sv <- tryCatch(sisvive_of(d$Y, d$X, d$Xi, i), error = function(e) NA_real_)
      ls <- ls_of(d)
      row <- mutate(row, f2sls_beta = fq[1], f2sls_se = fq[2], sisvive_beta = sv,
                    ls_beta = ls[1], ls_se = ls[2])
    }
    mutate(row, error = NA_character_)
  }, error = function(e) mutate(row, error = conditionMessage(e)))
}

# ---- runner ------------------------------------------------------------------------------
# Runs fun(i, cell) for i = 1..reps in every cell (row of grid) and returns all rows. Each
# cell is saved as it finishes, so an interrupted run resumes. With few replications per
# cell (a check run) all cells are run at once.
run_cells <- function(grid, fun, out_dir, reps) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  nc <- max(1L, n_cores() - 1L)
  msg(basename(out_dir), ": ", nrow(grid), " cells x ", reps, " replications on ", nc, " cores")
  if (reps < nc) {
    tasks <- expand_grid(k = seq_len(nrow(grid)), i = seq_len(reps))
    return(bind_rows(mclapply(seq_len(nrow(tasks)), function(r)
      fun(tasks$i[r], grid[tasks$k[r], ]), mc.cores = nc)))
  }
  bind_rows(lapply(seq_len(nrow(grid)), function(k) {
    cf <- file.path(out_dir, sprintf("cell_%03d.rds", k))
    if (file.exists(cf)) {
      prev <- readRDS(cf)
      if (nrow(prev) == reps) return(prev)
    }
    t0 <- Sys.time()
    out <- bind_rows(mclapply(seq_len(reps), fun, cell = grid[k, ], mc.cores = nc))
    saveRDS(out, cf)
    msg(sprintf("  cell %d/%d (%s): %.1f min, %d failed", k, nrow(grid),
                paste(names(grid), unlist(grid[k, ]), sep = " = ", collapse = ", "),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(!is.na(out$error))))
    out
  }))
}
