# Manuscript reproduction release

This directory defines the final analysis boundary for the manuscript. It is
intended to run in a Linux, Conda/Mamba, Docker, or Code Ocean environment;
no workstation-specific paths are required.

## Analyses represented

1. Strict LPS/O-antigen feature extraction and GO/KO summaries for 92 genomes.
2. Earth-versus-ISS feature presence/absence tests with two-sided Fisher tests
   and Benjamini–Hochberg correction.
3. Broad 92-genome marker-based phylogenomic context.
4. Gubbins prediction and masking on the 57-genome *Pantoea piersonii* core
   alignment, followed by IQ-TREE reconstruction of the filtered core-SNP
   tree.
5. SplitsTree NeighborNet analyses using the 369-family, 162,914-site protein
   alignment, JTT distances, MultiWay ordering and ActiveSet split weights,
   for both the 28-taxon representative set and all 92 genomes.
6. The separate exploratory *Pseudomonas fulva* comparison used for Figure 7.
7. Founder-aware, temporal, phylogeny-aware, locus-architecture and
   RELAX/GARD sensitivity analyses represented in Figures 4–6 and the
   supplementary tables.
8. Manuscript figures and supplementary tables generated from those outputs.

Other unrelated exploratory analyses are not part of the minimal reproduction
path. The *Pseudomonas fulva* branch is retained because it is described in
the manuscript, but it remains statistically separate from the Pantoea study.

## Portable execution

From the repository root:

```bash
mamba env create -f environment.yml
mamba activate lps-oantigen-analysis
bash reproducibility/run_checks.sh
```

Alternatively, build the included container:

```bash
docker build -t lps-oantigen-analysis .
docker run --rm lps-oantigen-analysis
```

The Code Ocean default command should be:

```bash
bash reproducibility/run_manuscript_release.sh
```

All input locations are supplied through configuration or command-line
arguments. Do not edit scripts to add machine-specific paths.

Raw genomes and Bakta/eggNOG databases are intentionally not bundled in
GitHub. They are obtained from public accessions or supplied to the optional
upstream annotation workflow. A Code Ocean capsule should contain the frozen
processed inputs needed to reproduce the manuscript figures and tables.

## Required frozen inputs

The capsule or release archive should contain an `inputs/` directory with:

```text
metadata.tsv
lps_strict_presence.tsv
lps_strict_presence_long.tsv
lps_strict_gene_fisher.tsv
piersonii57.clean.full.aln
pantoea92_markers_invariant_retained.faa
pantoea92_JTT_distances.nex
pantoea92_JTT_NeighborNet.stree6
pantoea28_taxon_representatives_JTT_NeighborNet_splits.nex
```

The 57-genome manifest must identify the reference assembly used for the core
alignment. The 92-genome metadata must retain Earth/ISS labels and the taxon
names used by the network analysis.

## Analysis order

Annotation and feature matrices are prepared first. Broad marker context and
the multispecies network are generated from the protein alignment. The
recombination-aware 57-genome analysis is run separately with Snippy, Gubbins
and IQ-TREE. Figures and tables are generated only after input validation.

## Code Ocean

The Code Ocean capsule should use this repository plus the frozen `inputs/`
directory and write outputs beneath `outputs/`. The complete Gubbins run is
available as a separate command because it requires the 57 assembly files and
can take substantially longer than figure regeneration.
