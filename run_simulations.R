#!/usr/bin/env Rscript
# =============================================================================
# run_simulations.R -- Figures 2, 3, S.1 and Tables S.1-S.4 (Section 4 and the supplement).
#
#   Rscript run_simulations.R            figures and tables from the stored replications in
#                                        simulations/results/ (about a minute)
#   Rscript run_simulations.R --check    re-simulates replications 1 and 2 of every cell and
#                                        compares them with the stored replications (a few minutes
#                                        on 30 cores)
#   Rscript run_simulations.R --rerun    re-runs every simulation (1000 replications; about 5 hours
#                                        on 30 cores) into simulations/results_rerun/, then makes
#                                        the figures and tables from it
#
# Options: --cores=N (default: all cores).
# Scenario S1 uses the ADNI genotype of DLGAP2 (see data/README.md). Without these files,
# --check and --rerun skip S1 and use its stored replications.
# Figures and tables are written to output/; each table is compared with the paper's
# (reference/). Run from the repository root.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
opt  <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  v <- sub(paste0("^--", key, "=?"), "", hit[1]); if (nzchar(v)) v else TRUE
}
MODE  <- if (isTRUE(opt("rerun"))) "rerun" else if (isTRUE(opt("check"))) "check" else "stored"
CORES <- opt("cores", as.character(parallel::detectCores()))
stopifnot(file.exists("paper/figure2.R"))        # run from the repository root
STORED <- "simulations/results"
HAVE_S1 <- all(file.exists(c("data/raw/data_cleaned_demo.csv", "data/raw/PTID.csv",
                             "data/processed/Genes_pca6/DLGAP2.csv",
                             "data/processed/Gene_expression_per_gene/DLGAP2.csv",
                             "data/processed/filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.csv")))
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# simulation scripts: results folder, script, argument, needs the ADNI files
SIMS <- list(
  list("main_S2",          "main.R",          "S2", FALSE),
  list("weak_iv",          "weak_iv.R",       NULL, FALSE),
  list("q_sensitivity",    "q_sensitivity.R", NULL, FALSE),
  list("sisvive_snp",      "sisvive_snp.R",   NULL, FALSE),
  list("alpha_curves/S2",  "alpha_curves.R",  "S2", FALSE),
  list("talpha_sargan/S2", "talpha_sargan.R", "S2", FALSE),
  list("main_S1",          "main.R",          "S1", TRUE),
  list("alpha_curves/S1",  "alpha_curves.R",  "S1", TRUE),
  list("talpha_sargan/S1", "talpha_sargan.R", "S1", TRUE))
PAPER <- c("figure2", "figure3", "figureS1", "tableS1", "tableS2", "tableS3", "tableS4")

run_sims <- function(root, nreps = NULL) {
  dir.create(file.path(root, "logs"), recursive = TRUE, showWarnings = FALSE)
  for (s in SIMS) {
    if (s[[4]] && !HAVE_S1) { log_msg("skip ", s[[1]], " (needs the ADNI files; stored results are used)"); next }
    log_msg("run  simulations/", s[[2]], " ", s[[3]])
    lf <- file.path(root, "logs", paste0(gsub("/", "_", s[[1]]), ".log"))
    env <- c(sprintf("FMR_SIM_RESULTS=%s", normalizePath(root)), sprintf("FMR_CORES=%s", CORES),
             if (!is.null(nreps)) sprintf("FMR_NREPS=%d", nreps))
    if (system2("Rscript", c("--vanilla", file.path("simulations", s[[2]]), s[[3]]), env = env,
                stdout = lf, stderr = lf) != 0) stop(s[[1]], " failed; see ", lf)
  }
}

make_outputs <- function(root, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (p in PAPER) {
    lf <- file.path(out, paste0(p, ".log"))
    st <- system2("Rscript", c("--vanilla", file.path("paper", paste0(p, ".R"))),
                  env = c(sprintf("FMR_SIM_RESULTS=%s", normalizePath(root)), sprintf("FMR_OUTPUT=%s", out)),
                  stdout = lf, stderr = lf)
    if (st != 0) stop("paper/", p, ".R failed; see ", lf)
    file.remove(lf)
    log_msg("made ", p)
  }
  ok <- vapply(paste0(c("tableS1", "tableS2", "tableS3", "tableS4"), ".tex"), function(f) {
    same <- identical(readLines(file.path("reference", f)), readLines(file.path(out, f)))
    cat(sprintf("  %-12s %s\n", f, if (same) "same as the paper" else "DIFFERS from the paper"))
    same
  }, TRUE)
  all(ok)
}

# Replications 1..R of every cell, new run vs stored.
compare <- function(root, R) {
  cat("\nReplications 1 to", R, "of every cell, re-simulated vs stored:\n")
  d <- vapply(SIMS, function(s) {
    if (s[[4]] && !HAVE_S1) return(0)
    if (grepl("^alpha_curves", s[[1]])) {           # curves: compare the per-replication beta-hats
      f <- paste0("curves_", s[[3]], ".rds")
      a <- readRDS(file.path(root, "alpha_curves", f)); b <- readRDS(file.path(STORED, "alpha_curves", f))
      dd <- max(mapply(function(x, y) { y <- y$betas[match(x$betas$i, y$betas$i), ]
        max(abs(c(x$betas$beta - y$beta, x$betas$lambda - y$lambda))) }, a, b))
      cat(sprintf("  %-16s beta-hat and lambda, max |difference| %.1e\n", s[[1]], dd)); return(dd)
    }
    a <- readRDS(file.path(root, s[[1]], "master.rds"))
    b <- readRDS(file.path(STORED, s[[1]], "master.rds"))
    b <- b[b$i <= R, ]
    key <- intersect(c("i", "n", "case", "Dsignal", "signal", "nbasis"), names(a))
    b <- b[match(do.call(paste, a[key]), do.call(paste, b[key])), ]
    if (anyNA(b$i)) stop(s[[1]], ": cells differ")
    num <- names(a)[vapply(a, is.numeric, TRUE)]
    dd <- max(vapply(num, function(v) { x <- a[[v]]; y <- b[[v]]
      if (!identical(is.na(x), is.na(y))) return(Inf)
      if (all(is.na(x))) 0 else max(abs(x - y) / pmax(abs(y), 1), na.rm = TRUE) }, 0))
    cat(sprintf("  %-16s %5d replications, %2d columns, max relative difference %.1e\n",
                s[[1]], nrow(a), length(num), dd))
    dd
  }, 0)
  cat(sprintf("\nlargest difference: %.1e -> %s\n", max(d),
              if (max(d) == 0) "identical" else if (max(d) < 1e-8) "identical up to rounding" else "DIFFERENT"))
  max(d) < 1e-8
}

log_msg("mode: ", MODE, " | cores: ", CORES, " | ADNI files for S1: ", if (HAVE_S1) "found" else "not found")
if (MODE == "stored") {
  ok <- make_outputs(STORED, "output")
} else if (MODE == "check") {
  root <- "simulations/results_check"                  # an interrupted check resumes here
  run_sims(root, 2L)
  ok <- compare(root, 2L)
} else {
  root <- "simulations/results_rerun"
  run_sims(root)
  if (!HAVE_S1) {                                      # S1 needs ADNI: use the stored replications
    for (f in c("main_S1/master.rds", "alpha_curves/curves_S1.rds", "talpha_sargan/S1/master.rds")) {
      dir.create(dirname(file.path(root, f)), recursive = TRUE, showWarnings = FALSE)
      file.copy(file.path(STORED, f), file.path(root, f), overwrite = TRUE)
    }
  }
  ok <- make_outputs(root, "output/rerun")
}
log_msg(if (ok) "DONE: all checks passed" else "DONE: SOME CHECKS FAILED")
quit(status = if (ok) 0 else 1)
