# Comparative LPS/O-Antigen Genomics Code

This directory contains a sanitized code bundle for the comparative genomics
analysis of lipopolysaccharide (LPS) and O-antigen-associated genes in
Pantoea/Kalamiella genomes, with an exploratory Pseudomonas fulva side analysis.

The code was recovered from the analysis workspace and cleaned for public
release. Local absolute paths, usernames, local conda paths, and manuscript
editing utilities were removed. External resources are represented by generic
relative paths or command names.

## Reproducibility Scope

The main Snakemake workflow automates the core analysis:

- Bakta genome annotation
- eggNOG-mapper functional annotation
- extraction of LPS/O-antigen-associated genes
- strict LPS/O-antigen hit extraction
- GO-term, GO-slim, and KO-category summaries
- rfb/O-antigen locus mapping
- selected LPS/O-antigen plots
- curation of rfbA, rfbB, rfbC, rfbD, and waaL orthologs for downstream
  codon-based analyses

Some manuscript analyses are provided as standalone downstream scripts rather
than fully integrated Snakemake rules. These include PhyloPhlAn/IQ-TREE
phylogenomic visualization, founder-aware and phylogeny-aware robustness tests,
locus-validation tables, RELAX/GARD sensitivity summaries, final supplementary
table assembly, and the exploratory Pseudomonas fulva side analysis.

Thus, this bundle is intended to provide transparent computational provenance
for the manuscript. It should reproduce the analysis when run with the same
input genomes, metadata, curated term lists, and external databases, but it is
not packaged as a one-command containerized workflow.

## Directory Contents

- `Snakefile`: main configurable Snakemake workflow.
- `MANIFEST.tsv`: short description of the main scripts and subdirectories.
- `environment.yml`: conda environment template for the main workflow and common
  downstream scripts.
- `CITATION.cff`: citation metadata placeholder for GitHub/Zenodo release.
- `LICENSE`: MIT license for the code.
- `scripts/`: helper scripts used by the main workflow and downstream analyses.
- `scripts/selection_phase2/`: ortholog curation, locus validation,
  phylogenomic context, robustness analyses, RELAX/GARD support, figure
  generation, and supplementary-table scripts.
- `scripts/selection_phase2/pseudomonas_fulva_lps_analysis/`: exploratory
  Pseudomonas fulva side-analysis scripts.
- `scripts/selection_phase2/pseudomonas_expanded_lps_analysis/`: expanded
  Pseudomonas reference-selection and annotation provenance scripts.

## Expected Inputs

A typical project layout is:

```text
data/genomes/                    genome FASTA files
metadata.tsv                     sample/group metadata
resources/bakta_db/              local Bakta database
resources/eggnog/                local eggNOG-mapper database
resources/go-basic.obo           Gene Ontology OBO file
resources/lps_go_terms.tsv       curated LPS-related GO terms
resources/goslim_prok_mapping.tsv
```

The default workflow discovers genome files from `data/genomes/` with extensions
`.fa`, `.fna`, or `.fasta`. Sample names are inferred from file names.

Alternatively, provide a tab-delimited sample table with columns:

```text
sample    fasta
```

The `metadata.tsv` or group map should contain sample identifiers and source
groups such as `ISS` and `Earth`. The scripts accept common column names such as
`sample`, `Sample ID`, `isolate`, `genome`, `group`, or `source`, depending on
the step.

## Main Workflow

Run from this directory:

```bash
snakemake --cores 8 \
  --config genome_dir=data/genomes \
           bakta_db=resources/bakta_db \
           eggnog_data_dir=resources/eggnog \
           go_obo=resources/go-basic.obo \
           lps_terms=resources/lps_go_terms.tsv \
           goslim_map=resources/goslim_prok_mapping.tsv \
           group_map=metadata.tsv
```

Or with an explicit sample table:

```bash
snakemake --cores 8 \
  --config sample_table=data/samples.tsv \
           bakta_db=resources/bakta_db \
           eggnog_data_dir=resources/eggnog \
           go_obo=resources/go-basic.obo \
           lps_terms=resources/lps_go_terms.tsv \
           goslim_map=resources/goslim_prok_mapping.tsv \
           group_map=metadata.tsv
```

The workflow writes outputs under:

```text
results/genome_annotations/
results/lps_analysis/
```

These output paths can be changed with `--config results_dir=... analysis_dir=...`.

## Environment Setup

Create the conda environment template with:

```bash
conda env create -f environment.yml
conda activate lps-oantigen-analysis
```

Some external tools and databases, especially Bakta, eggNOG-mapper, PhyloPhlAn,
Pathview, and HyPhy, may require additional setup outside this repository. Follow
the installation instructions for those tools and record the versions used in
the final repository release.

## Downstream Analyses

Downstream scripts in `scripts/selection_phase2/` expect the outputs produced by
the main workflow plus the corresponding phylogenomic or selection-analysis
inputs. They were used for analyses described in the manuscript, including:

- marker-gene phylogenomic context
- founder-aware strict-LPS-profile clustering
- ISS sampling-year persistence summaries
- species-restricted phylogeny-aware permutation tests
- coordinate-level O-antigen/LPS locus validation
- locus GC/flanking/mobility context
- RELAX/GARD sensitivity summaries
- final supplementary-table assembly

The Pseudomonas fulva exploratory scripts are retained separately because that
analysis was performed as a side analysis and was not combined statistically with
the primary Pantoea/Kalamiella comparison.

## Software Requirements

Install the relevant tools in the environment used to run the workflow:

- Snakemake
- Bakta
- eggNOG-mapper
- DIAMOND
- Python 3 with pandas, numpy, matplotlib, scipy, and Biopython where required
- R, plus packages needed by `plot_lps.R` and `run_pathview.R` if those scripts
  are used
- Optional downstream tools: MAFFT, PAL2NAL, IQ-TREE, HyPhy, PhyloPhlAn, Mash,
  and trimAl

Exact database versions and command-line tool versions should be reported with
the manuscript or repository release when available.

## Notes for Public Release

This bundle intentionally does not include:

- genome FASTA files
- Bakta or eggNOG databases
- generated result tables
- generated figures
- manuscript `.docx` files
- scripts whose only purpose was to edit manuscript text

Before archiving or publishing, verify that any input metadata released with the
code does not contain restricted sample information. NCBI download helper scripts
use the placeholder email `your.email@example.com`; replace it with the
appropriate contact email before use.

## Suggested GitHub/Journal Release Steps

1. Create a new GitHub repository.
2. Copy the contents of this `publication_code/` directory into the repository
   root.
3. Update `CITATION.cff` with the final repository URL, manuscript title,
   author list, DOI, and release date.
4. Add any publishable metadata tables or small curated term files if they are
   cleared for release.
5. Do not commit large databases, genome FASTA files, generated result
   directories, or manuscript files.
6. Create a versioned GitHub release and archive it with Zenodo if the journal
   requires a permanent DOI.
