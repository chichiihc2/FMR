# FMR: Functional Mendelian Randomization

R package and reproduction code for

> Chi-Shian Dai and Cai Li. *Functional Mendelian Randomization under Pleiotropy with
> Application to Alzheimer's Disease.*

FMR estimates the causal effect $\beta$ of an exposure $X$ on an outcome $Y$ when the
instrument is a function $Z(t)$, such as the genotypes of the variants of a gene ordered by
position, that may also act on the outcome directly (pleiotropy):

$$Y = \beta X + \int Z(t)\alpha(t)\,dt + \varepsilon .$$

The instrument enters through its leading functional principal component scores, the direct
effect $\alpha(t)$ is a penalized B-spline with the smoothing parameter chosen by GCV, and
$\beta$ comes with a bootstrap-$t$ confidence interval. The package also tests
$H_0:\alpha\equiv0$.

## Quick start

| Task | Command | Time |
|---|---|---|
| Install the package and the packages the scripts use | `Rscript install.R` (or `remotes::install_github("chichiihc2/FMR")` for the package only) | minutes |
| **Demo** | `Rscript run_demo.R` | seconds |
| **Simulation study** (Figures 2, 3, S.1; Tables S.1-S.4) | `Rscript run_simulations.R` | 1 minute (from stored results) |
| **Real-data analysis** (Table 2, Section 5, Figure S.2) | put the ADNI files in `data/raw/` ([data/README.md](data/README.md)), then `Rscript run_realdata.R` | about 2 hours |

Run the scripts from the repository root. Figures and tables are written to `output/`, and
each table is compared with the paper's (`reference/`): the scripts print `same as the paper`
or `DIFFERS from the paper`. R >= 4.3; the paper's results were produced with R 4.4.2 and
the package versions in `renv.lock` (`renv::restore()` installs them).

## Which script gives which table or figure

| Paper | Script | Results it reads (made by) | Output |
|---|---|---|---|
| Figure 2 (MSE) | `paper/figure2.R` | `simulations/results/main_S1`, `main_S2` (`simulations/main.R`) | `output/figure2a.pdf` (S1), `figure2b.pdf` (S2) |
| Figure 3 (coverage) | `paper/figure3.R` | as Figure 2 | `output/figure3a.pdf`, `figure3b.pdf` |
| Table 2 and the numbers of Section 5.1 | `paper/table2.R` | `realdata/gene_expression/results` (`realdata/gene_expression/fmr_genomewide.R`, `cis_2sls.R`) | `output/table2.tex`, `section5_1.txt` |
| Section 5.2 (APOE locus) | `paper/section5_2.R` | `realdata/apoe/results` (`realdata/apoe/fmr_apoe.R`) | `output/section5_2.txt` |
| Figure S.1 | `paper/figureS1.R` | `simulations/results/alpha_curves` (`simulations/alpha_curves.R`) | `output/figureS1a.pdf` (S1), `figureS1b.pdf` (S2) |
| Figure S.2 | `paper/figureS2.R` | `realdata/apoe/results` (`realdata/apoe/fmr_apoe.R`) | `output/figureS2a.pdf` (amyloid PET -> FAQ), `figureS2b.pdf` (BMI -> ADAS13) |
| Table S.1 | `paper/tableS1.R` | `simulations/results/talpha_sargan` (`simulations/talpha_sargan.R`) | `output/tableS1.tex` |
| Table S.2 | `paper/tableS2.R` | `simulations/results/q_sensitivity` (`simulations/q_sensitivity.R`) | `output/tableS2.tex` |
| Table S.3 | `paper/tableS3.R` | `simulations/results/weak_iv` (`simulations/weak_iv.R`), `main_S2` | `output/tableS3.tex` |
| Table S.4 | `paper/tableS4.R` | `simulations/results/sisvive_snp` (`simulations/sisvive_snp.R`) | `output/tableS4.tex` |

Table 1 (assumptions) and Figure 1 (geometric illustration) contain no computed results.
Each script in `paper/` runs on its own, e.g. `Rscript paper/figure2.R`.

## The FMR package

```r
library(FMR)
set.seed(1)
d <- simulate_fmr_data(n = 1500, s_alpha = 0.6, case = 2)   # Scenario S2 of the paper

fit <- fmr(y = d$y, x = d$x, Z = d$Z, pos = d$pos)
fit           # beta-hat of FMR (bootstrap-t and Wald intervals), F2SLS and OLS; first stage; T_alpha
plot(fit)     # alpha-hat(t) with a pointwise band
```

With real data, `Z` is the n x p matrix of allele counts of the variants of a gene, `pos`
their positions, and `covariates` a data frame (age, sex, genotype PCs, ...) for which `y`,
`x` and `Z` are adjusted:

```r
fit <- fmr(y, x, Z, pos = bp, covariates = data.frame(age, sex, PC1, PC2, PC3, PC4, PC5, PC6))
```

Main options: `k` (number of FPCA scores; default: first-stage BIC over 3 to 10), `nbasis`
(B-spline functions for $\alpha$, 30), `lambda_grid` (GCV grid, `2^(-10:10)`), `B`
(bootstrap draws, 199), `tau` (0.1), `seed`. See `?fmr`. The test
`tests/testthat/test-reproduce.R` checks that `fmr()` reproduces stored replications of the
simulation study.

## Simulation study (Section 4 and the supplement)

`Rscript run_simulations.R` makes Figures 2, 3, S.1 and Tables S.1-S.4 from the stored
replications in `simulations/results/` (one row per replication; seeds: `set.seed(i)` for
the data of replication `i`, `5e6 + i` for its bootstrap, `9e6 + i` for sisVIVE's folds).

* `Rscript run_simulations.R --check` re-simulates replications 1 and 2 of every cell of every
  simulation and compares them with the stored ones (a few minutes on 30 cores).
* `Rscript run_simulations.R --rerun` runs all simulations from scratch (1000 replications per
  cell; about 5 hours on 30 cores) into `simulations/results_rerun/` and makes the figures and
  tables from them.

| Script | Design | Results |
|---|---|---|
| `simulations/main.R S2` | S2: simulated $Z(t)$, $s = 0.015$, $n = 500, \dots, 3500$, Cases 1-6; FMR, Ridge, F2SLS, sisVIVE, LS | `main_S2` |
| `simulations/main.R S1` | S1: *DLGAP2* genotype, $n = 712$, $s = 0.01, \dots, 0.04$, Cases 1-6 | `main_S1` |
| `simulations/alpha_curves.R S2`, `S1` | $\hat\alpha(t)$ over 1000 replications | `alpha_curves` |
| `simulations/q_sensitivity.R` | $q = 15, 20, 30, 60$ | `q_sensitivity` |
| `simulations/weak_iv.R` | $s = 0.002, 0.004, 0.007$ | `weak_iv` |
| `simulations/sisvive_snp.R` | sisVIVE on LD-pruned variants and on FPCA scores | `sisvive_snp` |
| `simulations/talpha_sargan.R S2`, `S1` | $T_\alpha$ and Sargan tests under $\alpha \equiv 0$: S2 with $Z(t)$ on 101 points (data seed `i + 1000`); S1 with all 808 subjects, the *DLGAP2* genotype adjusted for 10 PCs, FPCA at pve 0.95 with fixed eigenfunction signs | `talpha_sargan` |

Shared code: `simulations/helpers/` (data generators, one replication of all methods) and
`lib/estimators.R`. Scenario S1 uses the ADNI genotype of *DLGAP2*; without the ADNI files,
`--check` and `--rerun` skip S1 and use its stored replications.

## Real-data analysis (Section 5)

ADNI data cannot be shared. With ADNI access, put the files listed in
[data/README.md](data/README.md) in `data/raw/` and run `Rscript run_realdata.R`. It needs
`plink2` and PLINK 1.9 (as `plink`) on the PATH, `curl` (it downloads eQTLGen, 4 GB) and about
70 GB of disk. Steps already done are skipped, so an interrupted run resumes.

| Step | Scripts | Output | Time (20 cores) |
|---|---|---|---|
| data preparation | `realdata/prep/` (genotype QC and PCs, variant to gene map, instrument of each gene, expression, eQTLGen) | `data/processed/` | about 1.5 h |
| Section 5.1, gene expression | `realdata/gene_expression/fmr_genomewide.R`, `cis_2sls.R` | `realdata/gene_expression/results/` | about 30 min |
| Section 5.2, APOE locus | `realdata/apoe/genotypes.sh`, `phenotypes.R`, `fmr_apoe.R` | `realdata/apoe/results/` | 5 min |
| Table 2, Figure S.2, numbers of the text | `paper/table2.R`, `section5_2.R`, `figureS2.R` | `output/` | seconds |

## Reproducibility

Checked on R 4.4.2 (x86_64, reference BLAS) on macOS:

* every table in `reference/` (Tables 2, S.1-S.4 and the numbers of Sections 5.1 and 5.2) has
  the same numbers as the paper, and Figures 2, 3, S.1 and S.2 match the paper's;
* `run_simulations.R --check` reproduces the stored replications exactly;
* `run_realdata.R` run from the ADNI files reproduces Table 2, Figure S.2 and the numbers of
  Section 5.

The smoothing parameter (GCV), the number of scores (BIC) and bootstrap quantiles are chosen
on grids, so near ties can resolve differently on another BLAS or platform and change the last
digits of some results.

## Repository layout

```
R/, man/, tests/, DESCRIPTION, NAMESPACE   the FMR package
run_demo.R, run_simulations.R, run_realdata.R, install.R
paper/                one script per table and figure of the paper
simulations/          simulation scripts; helpers/; results/ (stored replications)
realdata/prep/        ADNI data preparation
realdata/gene_expression/, realdata/apoe/    the analyses of Sections 5.1 and 5.2
lib/                  functions shared by the simulations and the real-data analysis
reference/            the paper's tables, against which the outputs are compared
data/                 input data (ADNI files are not included; see data/README.md)
output/               figures and tables made by the scripts
```

## License

MIT (see `LICENSE`).
