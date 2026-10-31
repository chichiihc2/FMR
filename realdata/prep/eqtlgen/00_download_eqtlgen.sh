#!/usr/bin/env bash
# Downloads the eQTLGen cis-eQTL summary statistics (Vosa et al., 2021) used to select the
# instruments of the conventional 2SLS in Table 2: the full cis-eQTL file (all tested SNP-gene
# pairs), the allele-frequency file and the READMEs, into data/external/eqtlgen/ (about 4 GB).
# Resumes a partial download.
#   bash realdata/prep/eqtlgen/00_download_eqtlgen.sh

set -euo pipefail

BASE="https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl"
OUTDIR="data/external/eqtlgen"

# needs free disk space for the download and the next step
MIN_FREE_GB=40

if [[ ! -f run_realdata.R ]]; then
  echo "ERROR: run this from the repository root." >&2
  exit 1
fi

mkdir -p "$OUTDIR"

free_gb=$(df -Pk . | awk 'NR==2 {print int($4 / 1048576)}')
if (( free_gb < MIN_FREE_GB )); then
  echo "ERROR: only ${free_gb}GB free; need >= ${MIN_FREE_GB}GB." >&2
  exit 1
fi
echo "Disk check OK (${free_gb}GB free)."

fetch () {
  local fname="$1"
  echo ""
  echo "==> $fname"
  # -C - resumes a partial file; --retry survives transient drops on a long pull.
  curl -fL -C - --retry 5 --retry-delay 10 \
       -o "${OUTDIR}/${fname}" "${BASE}/${fname}"
}

# 1. Full cis-eQTL summary statistics (all tested SNP-gene pairs).
fetch "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

# 2. Allele frequencies for the assessed allele (AlleleB).  Required for Z -> beta.
fetch "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"

# 3. Provenance.
fetch "README_cis"
fetch "README_allele_frequency"

echo ""
echo "=== Downloaded into ${OUTDIR} ==="
ls -lh "$OUTDIR"

echo ""
echo "Integrity check (gzip):"
for f in "${OUTDIR}"/*.txt.gz; do
  if gzip -t "$f" 2>/dev/null; then
    echo "  OK       $(basename "$f")"
  else
    echo "  CORRUPT  $(basename "$f")  <-- re-run to resume/repair" >&2
    exit 1
  fi
done

echo ""
echo "Next: Rscript realdata/prep/eqtlgen/01_build_eqtlgen_gene_tables.R"
