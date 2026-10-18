# Manuscript reproduction release

This directory describes the manuscript reproduction workflow. It is
intended to run in a Linux, Conda/Mamba, Docker, or Code Ocean environment;
no workstation-specific paths are required.

## Analyses represented

1. Strict LPS/O-antigen feature extraction and GO/KO summaries for 92 genomes.
2. Earth-versus-ISS feature presence/absence tests with two-sided Fisher tests
   and Benjamini–Hochberg correction.
3. Broad 92-genome marker-based phylogenomic context.
4. Gubbins prediction and masking on the 57-genome *Pantoea piersonii* core
   alignment, followed by a matched before/after IQ-TREE comparison using
   the same complete-site criterion and excluding the technical reference.
5. SplitsTree NeighborNet analyses using the 369-family, 162,914-site protein
   alignment, JTT distances, MultiWay ordering and ActiveSet split weights,
   for both the 28-taxon representative set and all 92 genomes.
6. The separate exploratory *Pseudomonas fulva* comparison used for Figure 7.
7. Exact-profile, temporal, locus-architecture and RELAX/GARD sensitivity
   analyses represented in Figures 4–6 and the
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
bash reproducibility/prepare_review_inputs.sh /path/to/deposited_data_archive.zip inputs
bash reproducibility/run_checks.sh
```

Alternatively, build the included container:

```bash
docker build -t lps-oantigen-analysis .
docker run --rm lps-oantigen-analysis
```

The following command is a smoke test for frozen inputs and rebuilds Figure 3;
it is **not** a complete one-command reproduction of every figure:

```bash
bash reproducibility/run_manuscript_release.sh
```

The Figure 3 script merges `waaL`, `rfaL` and `ko:K02847` into one
O-antigen-ligase feature by genome-level logical OR, then reruns two-sided
Fisher tests and Benjamini-Hochberg adjustment over **46** nonredundant
features (48 source rows minus two duplicates). Its output table is
`outputs/figures/figure3_waaL_merged_fisher.tsv`.

The archived Supplementary Table S2 full-92 Fisher worksheet contains these
same 46 nonredundant tests. It can be regenerated from the
Figure 3 table without changing the focal-gene or within-species worksheets:

```bash
python scripts/selection_phase2/harmonize_s2_merged_figure3.py \
  --source-workbook /path/to/Supplementary_Table_S2_source.xlsx \
  --figure3-table outputs/figures/figure3_waaL_merged_fisher.tsv \
  --output outputs/Supplementary_Table_S2_association_robustness.xlsx
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
lps_strict_presence.tsv
lps_strict_presence_long.tsv
lps_strict_gene_fisher.tsv
piersonii57.clean.full.aln
pantoea92_markers_invariant_retained.faa
pantoea92_JTT_distances.nex
pantoea92_JTT_NeighborNet.stree6
pantoea28_taxon_representatives_JTT_NeighborNet_splits.nex
```

The accompanying data ZIP contains these frozen inputs. The helper above
extracts them into `inputs/` without altering the archive. Source labels are
read from `lps_strict_presence_long.tsv`; no separate metadata.tsv is required
for this smoke test.

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
available as an upstream command; it requires the 57 assembly files and can
take substantially longer than figure regeneration. For the
matched-tree comparison, use
`scripts/selection_phase2/reproduce_matched_trees.sh` with the archived full
core alignment and Gubbins GFF, SNP-sites 2.5.1, Gubbins 3.4.3, and IQ-TREE
2.4.0. Set `MASK_SCRIPT` to Gubbins' `mask_gubbins_aln.py` path if necessary.
