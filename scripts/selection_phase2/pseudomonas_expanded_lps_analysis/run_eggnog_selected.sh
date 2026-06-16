#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="${ROOT}/selected_genomes.tsv"
EMAPPER="${EMAPPER:-emapper.py}"
WINDOWS_EGGNOG_DB="${WINDOWS_EGGNOG_DB:-resources/eggnog}"
LOCAL_EGGNOG_DB="${LOCAL_EGGNOG_DB:-/tmp/eggnog_db}"
THREADS="${EGGNOG_THREADS:-2}"

export PATH="${PATH}"
mkdir -p "${ROOT}/logs"

if [[ -s "${LOCAL_EGGNOG_DB}/eggnog.db" && -s "${LOCAL_EGGNOG_DB}/eggnog_proteins.dmnd" ]]; then
  DATA_DIR="${LOCAL_EGGNOG_DB}"
else
  DATA_DIR="${WINDOWS_EGGNOG_DB}"
fi

echo "Using eggNOG data dir: ${DATA_DIR}"

while IFS=$'\t' read -r accession strain isolate_id group source_note bioproject sra_accession flight location medium_temp_c contigs genome_size_bp n50_bp depth gc_percent include_reason; do
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
    --cpu "${THREADS}" \
    --output_dir "${outdir}" \
    --data_dir "${DATA_DIR}" \
    -m diamond \
    --sensmode sensitive \
    --dmnd_iterate no \
    --pfam_realign none \
    > "${ROOT}/logs/${isolate_id}.emapper.log" 2>&1
done < "${META}"

echo "eggNOG complete"
