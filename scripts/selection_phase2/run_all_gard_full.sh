#!/usr/bin/env bash
set -euo pipefail

HYPHY=hyphy
RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUTDIR="gard_results/full_gard_${RUN_ID}"
GENES=(rfbA rfbB rfbC rfbD waaL)

mkdir -p "${OUTDIR}"

{
  echo "Full GARD run"
  echo "Started: $(date)"
  echo "HyPhy: ${HYPHY}"
  "${HYPHY}" --version
  echo "Output directory: ${OUTDIR}"
  echo
} > "${OUTDIR}/RUN_INFO.txt"

for gene in "${GENES[@]}"; do
  alignment="${gene}.codon.aln.fasta"
  gene_out="${OUTDIR}/${gene}"
  mkdir -p "${gene_out}"

  {
    echo "Gene: ${gene}"
    echo "Alignment: ${alignment}"
    echo "Started: $(date)"
    echo
  } | tee "${gene_out}/${gene}.GARD.full.log"

  "${HYPHY}" gard \
    --type codon \
    --code Universal \
    --alignment "${alignment}" \
    --output "${gene_out}/${gene}.GARD.full.json" \
    --output-lf "${gene_out}/${gene}.best-gard.full" \
    --mode Normal \
    >> "${gene_out}/${gene}.GARD.full.log" 2>&1

  {
    echo
    echo "Finished: $(date)"
  } | tee -a "${gene_out}/${gene}.GARD.full.log"
done

echo "Completed: $(date)" >> "${OUTDIR}/RUN_INFO.txt"
echo "Full GARD run finished."
echo "Output directory: ${OUTDIR}"
