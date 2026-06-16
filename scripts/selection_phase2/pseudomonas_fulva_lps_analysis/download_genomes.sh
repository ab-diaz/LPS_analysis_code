#!/usr/bin/env bash
set -euo pipefail

# Download the three Pseudomonas fulva genome assemblies for the LPS/O-antigen
# mini-analysis. This script uses ENA's FASTA API because it can retrieve WGS
# assemblies directly from the GenBank/WGS accessions.
#
# Run from:
#   results/lps_analysis/selection_phase2/pseudomonas_fulva_lps_analysis
#
# Expected outputs:
#   genomes/<accession>_<strain>/<isolate_id>.fasta

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

download_one() {
  local accession="$1"
  local isolate_id="$2"
  local folder="$3"
  local outdir="${ROOT}/genomes/${folder}"
  local out="${outdir}/${isolate_id}.fasta"
  local tmp="${out}.tmp"

  mkdir -p "${outdir}"

  if [[ -s "${out}" ]]; then
    echo "Already exists: ${out}"
    return 0
  fi

  echo "Downloading ${accession} -> ${out}"
  curl -L --fail --retry 3 --retry-delay 5 \
    "https://www.ebi.ac.uk/ena/browser/api/fasta/${accession}?download=true" \
    -o "${tmp}"

  if [[ ! -s "${tmp}" ]]; then
    echo "ERROR: Downloaded file is empty for ${accession}" >&2
    rm -f "${tmp}"
    return 1
  fi

  if ! grep -q '^>' "${tmp}"; then
    echo "ERROR: Downloaded file does not look like FASTA for ${accession}" >&2
    head -n 5 "${tmp}" >&2 || true
    rm -f "${tmp}"
    return 1
  fi

  mv "${tmp}" "${out}"
}

download_one "JBQCNJ000000000.1" "Pseudomonas_fulva_51-6" "JBQCNJ000000000.1_Pseudomonas_fulva_51-6"
download_one "JAFDQI000000000.1" "Pseudomonas_fulva_F8_1S_1P" "JAFDQI000000000.1_Pseudomonas_fulva_F8_1S_1P"
download_one "JHYU00000000.1" "Pseudomonas_fulva_NBRC_16637" "JHYU00000000.1_Pseudomonas_fulva_NBRC_16637"

echo
echo "Downloaded FASTA files:"
find "${ROOT}/genomes" -type f -name '*.fasta' -print | sort
