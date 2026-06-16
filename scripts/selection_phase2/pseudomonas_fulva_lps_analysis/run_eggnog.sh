#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EMAPPER="emapper.py"
WINDOWS_EGGNOG_DB="resources/eggnog"
LOCAL_EGGNOG_DB="${LOCAL_EGGNOG_DB:-/tmp/eggnog_db}"
META_FILE="${META_FILE:-${ROOT}/genome_metadata.tsv}"

export PATH="${PATH}"
mkdir -p "${ROOT}/logs"

if [[ -s "${LOCAL_EGGNOG_DB}/eggnog.db" && -s "${LOCAL_EGGNOG_DB}/eggnog_proteins.dmnd" ]]; then
  DATA_DIR="${LOCAL_EGGNOG_DB}"
else
  DATA_DIR="${WINDOWS_EGGNOG_DB}"
fi

echo "Using eggNOG data dir: ${DATA_DIR}"

while IFS=$'\t' read -r accession strain isolate_id group source_note; do
  [[ "${accession}" == "accession" ]] && continue

  faa="${ROOT}/results/${isolate_id}/proteins.faa"
  outdir="${ROOT}/results/${isolate_id}/eggnog_out"
  annot="${outdir}/${isolate_id}.emapper.annotations"

  if [[ ! -s "${faa}" ]]; then
    echo "ERROR: missing Bakta protein FASTA: ${faa}" >&2
    exit 1
  fi

  if [[ -s "${annot}" ]]; then
    echo "eggNOG already done: ${isolate_id}"
    continue
  fi

  echo "Running eggNOG: ${isolate_id}"
  mkdir -p "${outdir}"
  rm -f "${outdir}/${isolate_id}".emapper.* \
    "${outdir}/${isolate_id}".seed_orthologs \
    "${outdir}/${isolate_id}".hits \
    "${outdir}/${isolate_id}".no_annotations 2>/dev/null || true

  "${EMAPPER}" \
    -i "${faa}" \
    -o "${isolate_id}" \
    --itype proteins \
    --cpu 2 \
    --output_dir "${outdir}" \
    --data_dir "${DATA_DIR}" \
    -m diamond \
    --sensmode sensitive \
    --dmnd_iterate no \
    --pfam_realign none \
    > "${ROOT}/logs/${isolate_id}.emapper.log" 2>&1
done < "${META_FILE}"

echo "eggNOG complete"
