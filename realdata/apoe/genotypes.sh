#!/usr/bin/env bash
# Section 5.2: genotypes of the APOE locus analysis (PLINK 1.9 as `plink`). Run from the
# repository root:  bash realdata/apoe/genotypes.sh
#   apoe_region   variants on chr19:44.8-46.0 Mb (GRCh37)
#   win.raw       allele counts of the variants on chr19:45.0-45.8 Mb
#   gpc           10 genome-wide PCs (LD-pruned autosomal variants, chr19:44-47 Mb excluded)
#   rel.genome    pairs with pi-hat > 0.2 (LD-pruned common SNPs)
set -euo pipefail
G=data/raw/filtered_ADNI_808_plink2
O=data/processed/apoe
mkdir -p $O
echo "19 44000000 47000000 apoe" > $O/apoe.range
plink --bfile $G --chr 19 --from-bp 44800000 --to-bp 46000000 --make-bed --out $O/apoe_region
plink --bfile $O/apoe_region --chr 19 --from-bp 45000000 --to-bp 45800000 --recode A --out $O/win
plink --bfile $G --autosome --exclude range $O/apoe.range --indep-pairwise 200 50 0.1 --threads 16 --out $O/prune
plink --bfile $G --extract $O/prune.prune.in --pca 10 --threads 16 --out $O/gpc
plink --bfile $G --autosome --snps-only just-acgt --geno 0.02 --maf 0.05 --hwe 1e-6 \
  --exclude range $O/apoe.range --indep-pairwise 200 50 0.1 --threads 16 --out $O/prune_rel
plink --bfile $G --extract $O/prune_rel.prune.in --genome --min 0.2 --threads 16 --out $O/rel
