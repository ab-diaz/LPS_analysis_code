#!/usr/bin/env bash
set -euo pipefail

cd results/lps_analysis/selection_phase2

conda run -n phylophlan_env phylophlan \
  -i phylogenomic_context_phylophlan/input_proteomes \
  -o pantoea_92_phylophlan_marker_tree \
  -d phylophlan -t a \
  -f phylogenomic_context_phylophlan/configs/phylophlan_aa_iqtree.cfg \
  --diversity low --accurate --nproc 8 \
  --databases_folder phylogenomic_context_phylophlan/databases \
  --output_folder phylogenomic_context_phylophlan/output \
  --data_folder phylogenomic_context_phylophlan/tmp \
  --proteome_extension .faa \
  --update --verbose
