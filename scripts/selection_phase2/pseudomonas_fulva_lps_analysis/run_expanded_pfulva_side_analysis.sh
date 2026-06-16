#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="${ROOT}/genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
if [[ ! -s "${META}" ]]; then
  META="${ROOT}/genome_metadata_pfulva_expanded.tsv"
fi
HITS="${ROOT}/lps_eggnog_hits_pfulva_expanded.tsv"

cd "${ROOT}"

if [[ ! -s "${META}" ]]; then
  echo "Missing ${META}"
  echo "Run first: python3 select_pfulva_earth_references.py"
  echo "Then:      python3 download_selected_pfulva_references.py"
  exit 1
fi

META_FILE="${META}" ./run_bakta.sh
META_FILE="${META}" ./run_eggnog.sh
META_FILE="${META}" HITS_FILE="${HITS}" python3 collect_lps_eggnog_hits.py
META_FILE="${META}" HITS_FILE="${HITS}" OUT_PREFIX="pfulva_expanded" python3 summarize_lps_hits.py

echo "Expanded Pseudomonas fulva side analysis complete."
echo "Main outputs:"
echo "  ${HITS}"
echo "  ${ROOT}/pfulva_expanded_lps_category_counts.tsv"
echo "  ${ROOT}/pfulva_expanded_o_antigen_hits.tsv"
echo "  ${ROOT}/pfulva_expanded_o_antigen_gene_presence.tsv"
