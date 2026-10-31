#!/usr/bin/env bash
# Genotype QC, step 2 of 2 (supplement, Section S.4): drop subjects with > 10% missing
# genotypes and variants with > 10% missing genotypes, MAF < 1% or HWE p < 1e-7 (plink2).
# Input:  data/raw/filtered_ADNI_808_b4plink2.{bed,bim,fam}
# Output: data/raw/filtered_ADNI_808_plink2.{bed,bim,fam}
#   bash realdata/prep/wgs_qc/02_qc_plink.sh

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data/raw}"
LOGS_DIR="${LOGS_DIR:-$PROJECT_ROOT/data/processed/logs}"

BFILE_IN="${BFILE_IN:-$DATA_DIR/filtered_ADNI_808_b4plink2}"
BFILE_OUT="${BFILE_OUT:-$DATA_DIR/filtered_ADNI_808_plink2}"

mkdir -p "$LOGS_DIR"
LOG="$LOGS_DIR/02_qc_plink_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[02_qc_plink] start $(date)"
echo "  BFILE_IN  = $BFILE_IN"
echo "  BFILE_OUT = $BFILE_OUT"

command -v plink2 >/dev/null || { echo "ERROR: plink2 not on PATH" >&2; exit 1; }
plink2 --version

for ext in bed bim fam; do
    [[ -f "${BFILE_IN}.${ext}" ]] || { echo "ERROR: missing ${BFILE_IN}.${ext}" >&2; exit 1; }
done

plink2 --bfile "$BFILE_IN" \
    --mind 0.1 \
    --geno 0.1 \
    --maf 0.01 \
    --hwe 1e-7 \
    --make-bed \
    --out "$BFILE_OUT"

echo "[02_qc_plink] done $(date)"
echo "  -> ${BFILE_OUT}.{bed,bim,fam}"
