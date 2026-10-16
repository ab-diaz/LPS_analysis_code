#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUTS="${MANUSCRIPT_INPUTS:-${ROOT}/inputs}"

required=(
  lps_strict_presence.tsv
  lps_strict_presence_long.tsv
  lps_strict_gene_fisher.tsv
  piersonii57.clean.full.aln
  pantoea92_markers_invariant_retained.faa
  pantoea92_JTT_distances.nex
)

missing=0
for name in "${required[@]}"; do
  if [[ ! -s "${INPUTS}/${name}" ]]; then
    echo "MISSING: ${INPUTS}/${name}" >&2
    missing=1
  else
    echo "OK: ${name}"
  fi
done

if [[ "${missing}" -ne 0 ]]; then
  echo "Set MANUSCRIPT_INPUTS to the frozen input directory." >&2
  exit 1
fi

python - "${INPUTS}" <<'PY'
from pathlib import Path
import sys
import pandas as pd

root = Path(sys.argv[1])
presence = pd.read_csv(root / "lps_strict_presence.tsv", sep="\t")
long = pd.read_csv(root / "lps_strict_presence_long.tsv", sep="\t")
if len(presence.columns) < 90:
    raise SystemExit(f"Expected approximately 92 genomes, found {len(presence.columns)-3}")
if set(long["group"].dropna()) != {"Earth", "ISS"}:
    raise SystemExit("Metadata must contain Earth and ISS source groups")
if long["isolate"].nunique() != 92:
    raise SystemExit(f"Expected 92 genomes, found {long['isolate'].nunique()}")
print(f"Validated {long['isolate'].nunique()} genomes and {len(presence)} feature rows")
PY
