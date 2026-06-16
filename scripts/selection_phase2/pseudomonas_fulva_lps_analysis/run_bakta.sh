#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_BIN="bin"
BAKTA_DB="resources/bakta_db"
META_FILE="${META_FILE:-${ROOT}/genome_metadata.tsv}"

export PATH="${ENV_BIN}:${PATH}"
mkdir -p "${ROOT}/logs"

while IFS=$'\t' read -r accession strain isolate_id group source_note; do
  [[ "${accession}" == "accession" ]] && continue

  fasta="$(find "${ROOT}/genomes" -type f -name "${isolate_id}.fasta" | head -n 1)"
  outdir="${ROOT}/results/${isolate_id}/bakta_output"

  if [[ -s "${outdir}/${isolate_id}.faa" && -s "${outdir}/${isolate_id}.ffn" ]]; then
    echo "Bakta already done: ${isolate_id}"
    continue
  fi

  echo "Running Bakta: ${isolate_id}"
  rm -rf "${outdir}"
  mkdir -p "${ROOT}/results/${isolate_id}"

  MPLCONFIGDIR=/tmp bakta \
    --db "${BAKTA_DB}" \
    --threads 2 \
    --output "${outdir}" \
    --prefix "${isolate_id}" \
    "${fasta}" \
    > "${ROOT}/logs/${isolate_id}.bakta.log" 2>&1

  cp "${outdir}/${isolate_id}.faa" "${ROOT}/results/${isolate_id}/proteins.faa"
done < "${META_FILE}"

echo "Bakta complete"
