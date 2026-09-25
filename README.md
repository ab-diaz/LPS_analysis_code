# Comparative LPS/O-Antigen Genomics Code

This repository contains the public, sanitized code for the manuscript's
comparative genomics analysis of lipopolysaccharide (LPS) and O-antigen-
associated genes in 92 Pantoea genomes. The release is organized around the
analyses actually used in the manuscript; unrelated exploratory analyses are
not required for reproduction.

The final manuscript analysis has two reproducibility levels:

1. **Figure/table reproduction**, starting from frozen annotation, metadata,
   alignments and feature matrices.
2. **Upstream annotation**, an optional workflow requiring local Bakta and
   eggNOG-mapper databases.

The repository is organized around one primary reproduction path:

```bash
conda env create -f environment.yml
conda activate lps-oantigen-analysis
snakemake --cores 8 --configfile config/config.yaml
```

The Snakemake workflow performs the primary Bakta/eggNOG/LPS analysis.
Recombination-aware and network analyses are documented under
`reproducibility/` and `scripts/selection_phase2/`.

## Repository Layout

```text
Snakefile                       Main workflow for the primary analysis
config/config.yaml              Editable configuration template
scripts/                        Scripts called directly by the Snakefile
scripts/selection_phase2/        Downstream manuscript analyses
reproducibility/                 Manuscript release instructions and checks
MANIFEST.tsv                    Short description of key scripts
environment.yml                 Conda environment template
LICENSE                         Code license
```

The most important files for reproducing the main analysis are `Snakefile`,
`config/config.yaml`, and the scripts directly under `scripts/`.

## What the Main Workflow Does

The main Snakemake workflow automates:

- Bakta genome annotation
- eggNOG-mapper functional annotation
- LPS/O-antigen-associated gene extraction
- strict LPS/O-antigen hit extraction
- GO-term, GO-slim, and KO-category summaries
- rfb/O-antigen locus mapping
- selected LPS/O-antigen plots
- curation of rfbA, rfbB, rfbC, rfbD, and waaL orthologs for downstream
  codon-based analyses

The workflow writes outputs under:

```text
results/genome_annotations/
results/lps_analysis/
```

These paths can be changed in `config/config.yaml`.

## Expected Inputs

Edit `config/config.yaml` before running. A typical local layout is:

```text
data/genomes/                    genome FASTA files
metadata.tsv                     sample/group metadata
resources/bakta_db/              local Bakta database
resources/eggnog/                local eggNOG-mapper database
resources/go-basic.obo           Gene Ontology OBO file
resources/lps_go_terms.tsv       curated LPS-related GO terms
resources/goslim_prok_mapping.tsv
```

By default, the workflow discovers genome files in `data/genomes/` with
extensions `.fa`, `.fna`, or `.fasta`, and sample names are inferred from file
names.

Alternatively, provide a tab-delimited sample table and set `sample_table` in
`config/config.yaml`. The table must contain:

```text
sample    fasta
```

The metadata/group map should contain sample identifiers and source groups such
as `ISS` and `Earth`. Several scripts accept common column names such as
`sample`, `Sample ID`, `isolate`, `genome`, `group`, or `source`.

## Main Reproduction Command

After editing `config/config.yaml`, run:

```bash
snakemake --cores 8 --configfile config/config.yaml
```

For a dry run:

```bash
snakemake --dry-run --cores 1 --configfile config/config.yaml
```

## Final manuscript analyses

The final manuscript uses:

- LPS/O-antigen feature presence/absence and GO/KO summaries for 92 genomes;
- broad marker-based multispecies context;
- a Gubbins-filtered core-SNP phylogeny for 57 *P. piersonii* genomes;
- NeighborNet networks for 28 taxon representatives and all 92 genomes;
- founder-aware, temporal and phylogeny-aware robustness analyses;
- locus-coordinate and surface-glycan architecture summaries;
- RELAX/GARD selection-sensitivity and founder-relative mutation summaries;
- a separate exploratory *Pseudomonas fulva* comparison presented as Figure 7;
- the manuscript figures and supplementary tables.

These scripts are retained as standalone scripts because they depend on
intermediate files from phylogenomic and codon-selection tools that are not all
fully wrapped in the main Snakemake workflow.

The *Pseudomonas fulva* comparison is kept separate from the Pantoea analyses
and is not combined statistically with them.

## Software Requirements

The provided `environment.yml` is a starting point:

```bash
conda env create -f environment.yml
conda activate lps-oantigen-analysis
```

Core tools:

- Snakemake
- Bakta
- eggNOG-mapper
- DIAMOND
- Python 3 with pandas, numpy, scipy, matplotlib, and Biopython

Downstream/optional tools:

- MAFFT
- trimAl
- IQ-TREE
- HyPhy
- PhyloPhlAn
- Gubbins and Snippy (57-genome recombination-aware analysis)
- SplitsTree App v6.9.5 (NeighborNet; Java application)
- Mash
- R and packages needed by the R scripts

External tools and databases may require setup outside this repository. The
recombination analysis used Gubbins 3.4.3 and IQ-TREE 2.4.0; the network
analysis used SplitsTree App v6.9.5. A Code Ocean capsule should freeze the
solved Conda environment in addition to this template.

## What Is Not Included

This repository intentionally does not include:

- genome FASTA files
- Bakta or eggNOG databases
- generated result tables
- generated figures
- manuscript files
- scripts whose only purpose was manuscript text editing

Before publishing metadata or input tables, verify that they do not contain
restricted sample information. NCBI download helper scripts require `NCBI_EMAIL`
to be set in the environment before use.
