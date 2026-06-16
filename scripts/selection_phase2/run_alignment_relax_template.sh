#!/usr/bin/env bash
set -euo pipefail

# Required external tools: mafft, pal2nal.pl, iqtree2 or iqtree, hyphy.
# This template starts from the DIAMOND-checked FASTAs in this directory.

for gene in rfbA rfbB rfbC rfbD waaL; do
  mafft --auto "${gene}.relax_input.faa" > "${gene}.aa.aln.faa"
  pal2nal.pl "${gene}.aa.aln.faa" "${gene}.relax_input.ffn" -output fasta > "${gene}.codon.aln.fasta"

  if command -v iqtree2 >/dev/null 2>&1; then
    iqtree2 -s "${gene}.codon.aln.fasta" -m MFP -B 1000 -T AUTO
    tree="${gene}.codon.aln.fasta.treefile"
  else
    iqtree -s "${gene}.codon.aln.fasta" -m MFP -bb 1000 -nt AUTO
    tree="${gene}.codon.aln.fasta.treefile"
  fi

  python3 label_relax_tree.py --tree "${tree}" --map sequence_id_map.tsv --gene "${gene}" --out "${gene}.relax_labeled.tree"
  hyphy relax --alignment "${gene}.codon.aln.fasta" --tree "${gene}.relax_labeled.tree" --test ISS --reference Earth
done
