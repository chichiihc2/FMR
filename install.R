# Installs the R packages used by the FMR package and the reproduction scripts, then the FMR
# package from this repository:   Rscript install.R
# The paper's results were produced with R 4.4.2 and the package versions in renv.lock
# (renv::restore() installs exactly those versions).
pkgs <- c("dr", "fda", "MASS", "mvtnorm", "refund", "ivreg",                       # FMR package
          "data.table", "dplyr", "genio", "ggplot2", "here", "readr", "sisVIVE",     # scripts
          "tibble", "tidyr", "tidyverse")
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
install.packages(".", repos = NULL, type = "source")
cat("FMR", as.character(packageVersion("FMR")), "installed\n")
