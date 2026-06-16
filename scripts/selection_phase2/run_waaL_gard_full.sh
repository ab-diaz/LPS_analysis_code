#!/usr/bin/env bash
set -euo pipefail

HYPhy=hyphy
GENE=waaL
ALIGNMENT="${GENE}.codon.aln.fasta"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUTDIR="gard_results/${GENE}_full_${RUN_ID}"

mkdir -p "${OUTDIR}"

"${HYPhy}" gard \
  --type codon \
  --code Universal \
  --alignment "${ALIGNMENT}" \
  --output "${OUTDIR}/${GENE}.GARD.full.json" \
  --output-lf "${OUTDIR}/${GENE}.best-gard.full" \
  --mode Normal \
  > "${OUTDIR}/${GENE}.GARD.full.log" 2>&1

echo "GARD finished for ${GENE}"
echo "Output directory: ${OUTDIR}"
