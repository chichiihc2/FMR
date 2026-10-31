# Section 5.2: phenotypes of the APOE locus analysis from the ADNIMERGE2 package. One row per
# genotyped subject at baseline (the first ADAS-Cog 13 visit); each variable is the
# measurement nearest to baseline within 365 days. Of each pair of related subjects
# (pi-hat > 0.2) one is dropped.
#   amy_cl  amyloid PET (centiloids)     faq     log(1 + FAQ total)
#   bmi     body mass index              adas13  ADAS-Cog 13
# Usage: Rscript realdata/apoe/phenotypes.R       (after genotypes.sh)
# Output: data/processed/apoe/pheno.csv
suppressPackageStartupMessages({ library(dplyr); library(readr) })
RAW <- "data/raw/ADNIMERGE2/data"; P <- "data/processed/apoe"
ld <- function(f) { e <- new.env(); load(file.path(RAW, paste0(f, ".rda")), envir = e)
  as.data.frame(get(ls(e)[1], e)) }
num <- function(v) suppressWarnings(as.numeric(as.character(v)))

fam <- read.table("data/raw/filtered_ADNI_808_plink2.fam")$V2
adas <- ld("ADAS") %>% transmute(PTID = as.character(PTID), d = as.Date(VISDATE), v = num(TOTAL13)) %>%
  filter(is.finite(v), !is.na(d))
anchor <- adas %>% group_by(PTID) %>% summarise(anchor = min(d), .groups = "drop") %>% filter(PTID %in% fam)
near <- function(tab, name, tol = 365) {
  tab %>% filter(is.finite(v), !is.na(d)) %>% inner_join(anchor, "PTID") %>%
    mutate(gap = abs(as.numeric(d - anchor))) %>% filter(gap <= tol) %>%
    group_by(PTID) %>% slice_min(gap, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(PTID, !!name := v)
}
pet <- ld("ADSP_PHC_PET_Amyloid_Simple") %>% mutate(PTID = as.character(PTID), d = as.Date(PHC_SCANDATE))
faq <- ld("FAQ") %>% mutate(PTID = as.character(PTID), d = as.Date(VISDATE))
vit <- ld("VITALS") %>% mutate(PTID = as.character(PTID), d = as.Date(VISDATE),
  kg = ifelse(VSWTUNIT == "pounds", num(VSWEIGHT) * 0.4536, num(VSWEIGHT)),
  m  = ifelse(VSHTUNIT == "inches", num(VSHEIGHT) * 0.0254, num(VSHEIGHT) / 100))
vars <- list(
  amy_cl = pet %>% transmute(PTID, d, v = num(PHC_CENTILOIDS)),
  bmi    = vit %>% transmute(PTID, d, v = kg / m^2) %>% filter(v > 12, v < 60),
  adas13 = adas %>% transmute(PTID, d, v),
  faq    = faq %>% transmute(PTID, d, v = log1p(num(FAQTOTAL))) %>% filter(v >= 0))
dem <- ld("PTDEMOG") %>% mutate(PTID = as.character(PTID)) %>% group_by(PTID) %>%
  summarise(sex = first(na.omit(as.character(PTGENDER))), educ = first(na.omit(num(PTEDUCAT))),
            dob = first(na.omit(num(substr(as.character(PTDOBYY), 1, 4)))),
            dob_my = first(na.omit(as.character(PTDOB))), .groups = "drop")
pcs <- read.table(file.path(P, "gpc.eigenvec"))[, -1]; names(pcs) <- c("PTID", paste0("PC", 1:10))

W <- anchor
for (nm in names(vars)) W <- left_join(W, near(vars[[nm]], nm), "PTID")
W <- W %>% left_join(dem, "PTID") %>% inner_join(pcs, "PTID") %>%
  mutate(age = as.numeric(format(anchor, "%Y")) - dob,
         age_exact = as.numeric(anchor - as.Date(paste0("15/", dob_my), "%d/%m/%Y")) / 365.25)   # PTDOB is MM/YYYY
rel <- read.table(file.path(P, "rel.genome"), header = TRUE); drop <- character(0)
for (r in seq_len(nrow(rel))) { a <- rel$IID1[r]; b <- rel$IID2[r]
  if (a %in% W$PTID && b %in% W$PTID && !(a %in% drop) && !(b %in% drop)) drop <- c(drop, b) }
W <- W %>% filter(!PTID %in% drop, !is.na(sex), !is.na(educ), !is.na(age))
cat("subjects:", nrow(W), "\n")
write_csv(W, file.path(P, "pheno.csv"))
