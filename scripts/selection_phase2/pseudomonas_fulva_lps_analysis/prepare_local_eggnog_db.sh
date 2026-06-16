#!/usr/bin/env bash
set -euo pipefail

SRC="resources/eggnog"
DEST="${LOCAL_EGGNOG_DB:-/tmp/eggnog_db}"

mkdir -p "${DEST}"

for file in eggnog.db eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd; do
  if [[ -s "${DEST}/${file}" ]]; then
    echo "Already exists: ${DEST}/${file}"
  else
    echo "Copying ${file} to ${DEST}"
    cp "${SRC}/${file}" "${DEST}/${file}"
  fi
done

echo "Local eggNOG database ready: ${DEST}"
