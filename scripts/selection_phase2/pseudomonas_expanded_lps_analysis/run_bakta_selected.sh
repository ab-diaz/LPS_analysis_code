#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="${ROOT}/selected_genomes.tsv"
ENV_BIN="bin"
BAKTA_DB="${BAKTA_DB:-resources/bakta_db}"
THREADS="${BAKTA_THREADS:-2}"

export PATH="${ENV_BIN}:${PATH}"
mkdir -p "${ROOT}/logs"

while IFS=$'\t' read -r accession strain isolate_id group source_note bioproject sra_accession flight location medium_temp_c contigs genome_size_bp n50_bp depth gc_percent include_reason; do
  [[ "${accession}" == "accession" ]] && continue

  fasta="$(find "${ROOT}/genomes" -type f -path "*/${accession}_${isolate_id}/${isolate_id}.fasta" | head -n 1)"
  outdir="${ROOT}/results/${isolate_id}/bakta_output"

  if [[ -z "${fasta}" || ! -s "${fasta}" ]]; then
    echo "ERROR: missing genome FASTA for ${isolate_id} (${accession})" >&2
    exit 1
  fi

  if [[ -s "${outdir}/${isolate_id}.faa" && -s "${outdir}/${isolate_id}.ffn" ]]; then
    echo "Bakta already done: ${isolate_id}"
    continue
  fi

  echo "Running Bakta: ${isolate_id}"
  rm -rf "${outdir}"
  mkdir -p "${ROOT}/results/${isolate_id}"

  MPLCONFIGDIR=/tmp bakta \
    --db "${BAKTA_DB}" \
    --threads "${THREADS}" \
    --output "${outdir}" \
    --prefix "${isolate_id}" \
    "${fasta}" \
    > "${ROOT}/logs/${isolate_id}.bakta.log" 2>&1

  cp "${outdir}/${isolate_id}.faa" "${ROOT}/results/${isolate_id}/proteins.faa"
done < "${META}"

echo "Bakta complete"
