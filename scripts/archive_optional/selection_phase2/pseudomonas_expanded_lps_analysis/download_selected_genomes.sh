#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="${ROOT}/selected_genomes.tsv"

mkdir -p "${ROOT}/genomes"

while IFS=$'\t' read -r accession strain isolate_id group source_note bioproject sra_accession flight location medium_temp_c contigs genome_size_bp n50_bp depth gc_percent include_reason; do
  [[ "${accession}" == "accession" ]] && continue

  outdir="${ROOT}/genomes/${accession}_${isolate_id}"
  out="${outdir}/${isolate_id}.fasta"
  tmp="${out}.tmp"
  mkdir -p "${outdir}"

  if [[ -s "${out}" ]]; then
    echo "Already exists: ${out}"
    continue
  fi

  echo "Downloading ${accession} -> ${out}"
  curl -L --fail --retry 3 --retry-delay 5 \
    "https://www.ebi.ac.uk/ena/browser/api/fasta/${accession}?download=true" \
    -o "${tmp}"

  if [[ ! -s "${tmp}" ]] || ! grep -q '^>' "${tmp}"; then
    echo "ERROR: Downloaded file does not look like FASTA for ${accession}" >&2
    head -n 5 "${tmp}" >&2 || true
    rm -f "${tmp}"
    exit 1
  fi

  mv "${tmp}" "${out}"
done < "${META}"

echo
echo "Downloaded FASTA files:"
find "${ROOT}/genomes" -type f -name '*.fasta' -print | sort
