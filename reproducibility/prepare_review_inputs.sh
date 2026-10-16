#!/usr/bin/env bash
set -euo pipefail

# Extract only the frozen inputs required by run_checks.sh and Figure 3.
# Usage: bash reproducibility/prepare_review_inputs.sh DATA_ZIP OUTPUT_DIR
DATA_ZIP="${1:?Supply the corrected data ZIP}"
OUT="${2:?Supply a new/existing output directory}"
if [[ ! -f "$DATA_ZIP" ]]; then
  echo "Data ZIP not found: $DATA_ZIP" >&2
  exit 2
fi
mkdir -p "$OUT"
python3 - "$DATA_ZIP" "$OUT" <<'PY'
from pathlib import Path
import sys
import zipfile

archive = Path(sys.argv[1])
out = Path(sys.argv[2])
mapping = {
    "processed_lps_go_ko/lps_strict_presence.tsv": "lps_strict_presence.tsv",
    "processed_lps_go_ko/lps_strict_presence_long.tsv": "lps_strict_presence_long.tsv",
    "processed_lps_go_ko/lps_strict_gene_fisher.tsv": "lps_strict_gene_fisher.tsv",
    "recombination_aware_piersonii/core/piersonii57.clean.full.aln": "piersonii57.clean.full.aln",
    "phylogenomics/pantoea92_markers_invariant_retained.faa": "pantoea92_markers_invariant_retained.faa",
    "splitstree_networks/output/pantoea92_JTT_distances.nex": "pantoea92_JTT_distances.nex",
    "splitstree_networks/output/pantoea92_JTT_NeighborNet.stree6": "pantoea92_JTT_NeighborNet.stree6",
    "splitstree_networks/output/pantoea28_taxon_representatives_JTT_NeighborNet_splits.nex": "pantoea28_taxon_representatives_JTT_NeighborNet_splits.nex",
}
with zipfile.ZipFile(archive) as z:
    missing = set(mapping) - set(z.namelist())
    if missing:
        raise SystemExit(f"Missing archived inputs: {sorted(missing)}")
    for source, target in mapping.items():
        dest = out / target
        if dest.exists():
            raise SystemExit(f"Refusing to overwrite {dest}")
        with z.open(source) as src, dest.open("wb") as dst:
            while chunk := src.read(1024 * 1024):
                dst.write(chunk)
        print(dest)
PY
