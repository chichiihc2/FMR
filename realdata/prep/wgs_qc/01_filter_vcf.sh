#!/usr/bin/env bash
# Genotype QC, step 1 of 2 (supplement, Section S.4): ADNI WGS VCF -> PLINK binary.
#   vcftools: keep PASS sites; genotypes with GQ < 20 or DP < 5 set to missing; drop sites with
#             mean depth < 10.
#   plink2:   biallelic variants, PLINK 1 binary.
# Input:  data/raw/ADNI.808_indiv.minGQ_21.pass_combined.vcf[.gz]
# Output: data/raw/filtered_ADNI_808_b4plink2.{bed,bim,fam}
# Requires vcftools and plink2. Paths can be changed with DATA_DIR, VCF_IN, OUT_PREFIX, BFILE_PREFIX.
#   bash realdata/prep/wgs_qc/01_filter_vcf.sh

set -euo pipefail

# --- Paths (override via environment) ---------------------------------------
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data/raw}"
LOGS_DIR="${LOGS_DIR:-$PROJECT_ROOT/data/processed/logs}"

VCF_IN="${VCF_IN:-$DATA_DIR/ADNI.808_indiv.minGQ_21.pass_combined.vcf}"
OUT_PREFIX="${OUT_PREFIX:-$DATA_DIR/ADNI_808_vcf_final}"
BFILE_PREFIX="${BFILE_PREFIX:-$DATA_DIR/filtered_ADNI_808_b4plink2}"

mkdir -p "$LOGS_DIR"
LOG="$LOGS_DIR/01_filter_vcf_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[01_filter_vcf] start $(date)"
echo "  VCF_IN       = $VCF_IN"
echo "  OUT_PREFIX   = $OUT_PREFIX"
echo "  BFILE_PREFIX = $BFILE_PREFIX"

# the VCF may be compressed or not
if [[ ! -f "$VCF_IN" && -f "${VCF_IN}.gz" ]]; then
    VCF_IN="${VCF_IN}.gz"
fi
if [[ ! -f "$VCF_IN" ]]; then
    echo "ERROR: input VCF not found: $VCF_IN" >&2
    exit 1
fi

if [[ "$VCF_IN" == *.gz ]]; then
    VCF_FLAG="--gzvcf"
else
    VCF_FLAG="--vcf"
fi
echo "  vcftools input flag = $VCF_FLAG"

command -v vcftools >/dev/null || { echo "ERROR: vcftools not on PATH" >&2; exit 1; }
command -v plink2  >/dev/null || { echo "ERROR: plink2 not on PATH"  >&2; exit 1; }
vcftools --version
plink2 --version

# --- 1. Variant-level filtering (vcftools) ----------------------------------
vcftools "$VCF_FLAG" "$VCF_IN" \
    --remove-filtered-all \
    --minGQ 20 \
    --minDP 5 \
    --min-meanDP 10 \
    --recode \
    --out "$OUT_PREFIX"

# --- 2. Conversion to PLINK1 binary, biallelic SNPs only --------------------
plink2 --vcf "${OUT_PREFIX}.recode.vcf" \
    --make-bed \
    --max-alleles 2 \
    --out "$BFILE_PREFIX"

echo "[01_filter_vcf] done $(date)"
echo "  -> ${OUT_PREFIX}.recode.vcf"
echo "  -> ${BFILE_PREFIX}.{bed,bim,fam}"
