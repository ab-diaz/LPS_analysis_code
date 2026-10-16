#!/usr/bin/env bash
set -euo pipefail

# Usage: IQTREE_BIN=iqtree2 SNP_SITES_BIN=snp-sites MASK_SCRIPT=/path/mask_gubbins_aln.py \
#   bash reproduce_matched_trees.sh PROJECT_ROOT OUTPUT_DIRECTORY
# Requires IQ-TREE 2.4.0, SNP-sites 2.5.1, Gubbins mask_gubbins_aln.py 3.4.3.

PROJECT_ROOT="${1:?Supply the selection_phase2 project directory}"
OUT="${2:?Supply a new output directory}"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"
SNP_SITES_BIN="${SNP_SITES_BIN:-snp-sites}"
MASK_SCRIPT="${MASK_SCRIPT:-$(command -v mask_gubbins_aln.py || true)}"
if [[ -z "$MASK_SCRIPT" || ! -f "$MASK_SCRIPT" ]]; then
  echo "Set MASK_SCRIPT to the Gubbins 3.4.3 mask_gubbins_aln.py path" >&2
  exit 2
fi

FULL="$PROJECT_ROOT/recombination_aware_piersonii/core/piersonii57.clean.full.aln"
GFF="$PROJECT_ROOT/recombination_aware_piersonii/gubbins/piersonii57.recombination_predictions.gff"
GUBBINS_SNP="$PROJECT_ROOT/recombination_aware_piersonii/gubbins/piersonii57.filtered_polymorphic_sites.fasta"

mkdir -p "$OUT"
if [[ -e "$OUT/pre57_matched.treefile" || -e "$OUT/post57_matched.treefile" ]]; then
  echo "Refusing to overwrite existing matched trees in $OUT" >&2
  exit 1
fi

"$SNP_SITES_BIN" -c -o "$OUT/pre58_complete_snps.fasta" "$FULL"
python3 "$MASK_SCRIPT" --aln "$FULL" --gff "$GFF" \
  --out "$OUT/post58_masked_full.fasta" --missing-char N
"$SNP_SITES_BIN" -c -o "$OUT/post58_masked_complete_snps.fasta" \
  "$OUT/post58_masked_full.fasta"
"$SNP_SITES_BIN" -c -o "$OUT/gubbins_output_complete_snps.fasta" "$GUBBINS_SNP"
cmp "$OUT/post58_masked_complete_snps.fasta" "$OUT/gubbins_output_complete_snps.fasta"

"$SNP_SITES_BIN" -v -c "$FULL" | awk '!/^#/ {print $2}' > "$OUT/pre58_snp_positions.txt"
"$SNP_SITES_BIN" -v -c "$OUT/post58_masked_full.fasta" | \
  awk '!/^#/ {print $2}' > "$OUT/post58_snp_positions.txt"
if [[ $(comm -13 <(sort "$OUT/pre58_snp_positions.txt") \
                 <(sort "$OUT/post58_snp_positions.txt") | wc -l) -ne 0 ]]; then
  echo "Post-mask SNP positions were not a subset of pre-mask SNP positions" >&2
  exit 1
fi

awk '/^>/{keep=($0 != ">Reference")} keep' \
  "$OUT/pre58_complete_snps.fasta" > "$OUT/pre57_with_refonly_sites.fasta"
awk '/^>/{keep=($0 != ">Reference")} keep' \
  "$OUT/post58_masked_complete_snps.fasta" > "$OUT/post57_with_refonly_sites.fasta"
"$SNP_SITES_BIN" -c -o "$OUT/pre57_complete_snps.fasta" \
  "$OUT/pre57_with_refonly_sites.fasta"
"$SNP_SITES_BIN" -c -o "$OUT/post57_complete_snps.fasta" \
  "$OUT/post57_with_refonly_sites.fasta"

for state in pre post; do
  "$IQTREE_BIN" -s "$OUT/${state}57_complete_snps.fasta" \
    -m TVM+F+ASC+R2 -B 1000 -alrt 1000 -T AUTO -seed 20261016 \
    --prefix "$OUT/${state}57_matched"
done

python3 "$(dirname "$0")/compare_matched_trees.py" \
  "$PROJECT_ROOT/recombination_aware_piersonii/piersonii_57_manifest_qc.tsv" \
  "$OUT/pre57_matched.treefile" "$OUT/post57_matched.treefile" \
  "$OUT/matched_tree_comparison.json"
