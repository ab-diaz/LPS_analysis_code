#!/usr/bin/env bash
set -euo pipefail

# Portable version of the 57-genome recombination-aware analysis.
# Required commands on PATH: snippy, snippy-core, snippy-clean_full_aln,
# run_gubbins.py and iqtree2.

usage() {
  echo "Usage: $0 --manifest manifest.tsv --reference reference.fna --outdir OUTDIR [--threads N]"
  echo "Manifest columns: genome<TAB>assembly (a header row is required)."
}

MANIFEST=""; REFERENCE=""; OUTDIR=""; THREADS=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -s "$MANIFEST" && -s "$REFERENCE" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
mkdir -p "$OUTDIR"/{snippy,core,gubbins,trees,logs}

snippy --version > "$OUTDIR/software_versions.txt" 2>&1 || true
run_gubbins.py --version >> "$OUTDIR/software_versions.txt" 2>&1 || true
iqtree2 --version >> "$OUTDIR/software_versions.txt" 2>&1 || true

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r genome assembly rest; do
  [[ -n "$genome" && -s "$assembly" ]] || { echo "Missing assembly for $genome" >&2; exit 1; }
  snippy --force --cpus 2 --outdir "$OUTDIR/snippy/$genome" \
    --prefix "$genome" --ref "$REFERENCE" --ctgs "$assembly" \
    > "$OUTDIR/logs/snippy_${genome}.stdout.log" \
    2> "$OUTDIR/logs/snippy_${genome}.stderr.log"
done

mapfile -t SNIPPY_DIRS < <(find "$OUTDIR/snippy" -mindepth 1 -maxdepth 1 -type d | sort)
snippy-core --ref "$REFERENCE" --prefix "$OUTDIR/core/piersonii57" "${SNIPPY_DIRS[@]}" \
  > "$OUTDIR/logs/snippy-core.log" 2>&1
snippy-clean_full_aln "$OUTDIR/core/piersonii57.full.aln" > "$OUTDIR/core/piersonii57.clean.full.aln"

iqtree2 -s "$OUTDIR/core/piersonii57.aln" -m MFP+ASC -B 1000 -alrt 1000 -T AUTO \
  --prefix "$OUTDIR/trees/piersonii57_unfiltered_core_snp" \
  > "$OUTDIR/logs/iqtree_unfiltered.stdout.log" 2> "$OUTDIR/logs/iqtree_unfiltered.stderr.log"

run_gubbins.py --prefix "$OUTDIR/gubbins/piersonii57" --threads "$THREADS" \
  --tree-builder iqtree "$OUTDIR/core/piersonii57.clean.full.aln" \
  > "$OUTDIR/logs/gubbins.stdout.log" 2> "$OUTDIR/logs/gubbins.stderr.log"

iqtree2 -s "$OUTDIR/gubbins/piersonii57.filtered_polymorphic_sites.fasta" \
  -m MFP+ASC -B 1000 -alrt 1000 -T AUTO \
  --prefix "$OUTDIR/trees/piersonii57_recombination_filtered" \
  > "$OUTDIR/logs/iqtree_filtered.stdout.log" 2> "$OUTDIR/logs/iqtree_filtered.stderr.log"

echo "Recombination-aware analysis completed in $OUTDIR"
