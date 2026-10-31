#!/usr/bin/env bash
# Genotype principal components (plink2 --pca 50); PC1-PC6 adjust the instruments of the gene
# expression analysis (realdata/prep/R/05_adjust_FIV_with_PCA.R). plink2 is required: PLINK 1.9
# gives PCs that differ in the last digits.
# Input:  data/raw/filtered_ADNI_808_plink2.{bed,bim,fam}
# Output: data/processed/filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.{eigenvec,eigenval}
#   bash realdata/prep/wgs_qc/03_pca_plink.sh

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data/raw}"
LOGS_DIR="${LOGS_DIR:-$PROJECT_ROOT/data/processed/logs}"

BFILE_IN="${BFILE_IN:-$DATA_DIR/filtered_ADNI_808_plink2}"
OUT_PREFIX="${OUT_PREFIX:-$PROJECT_ROOT/data/processed/filtered_ADNI_808_plink2/pca50_ADNI_808_plink2}"

mkdir -p "$LOGS_DIR"
LOG="$LOGS_DIR/03_pca_plink_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "[03_pca_plink] start $(date)"
echo "  BFILE_IN   = $BFILE_IN"
echo "  OUT_PREFIX = $OUT_PREFIX"

PLINK=$(command -v plink2 || true)
[[ -n "$PLINK" ]] || { echo "ERROR: plink2 is not on the PATH (https://www.cog-genomics.org/plink/2.0/)" >&2; exit 1; }
mkdir -p "$(dirname "$OUT_PREFIX")"
"$PLINK" --version

for ext in bed bim fam; do
    [[ -f "${BFILE_IN}.${ext}" ]] || { echo "ERROR: missing ${BFILE_IN}.${ext}" >&2; exit 1; }
done

"$PLINK" --bfile "$BFILE_IN" \
    --pca 50 \
    --out "$OUT_PREFIX"

echo "[03_pca_plink] done $(date)"
echo "  -> ${OUT_PREFIX}.eigenvec"
echo "  -> ${OUT_PREFIX}.eigenval"
