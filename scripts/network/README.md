# NeighborNet analysis

The manuscript network analysis used the 369-family protein alignment with
invariant sites retained (162,914 amino-acid positions). JTT distances were
calculated in SplitsTree App v6.9.5, followed by NeighborNet with MultiWay
ordering and ActiveSet split-weight inference.

`prepare_splitstree_input.py` validates a protein alignment and writes a
portable NEXUS alignment. `make_taxon_representative_input.py` selects the
within-taxon JTT-distance medoid for *P. piersonii* and *P. septica* and keeps
singleton taxa. `make_publication_overview.py` creates a documented overview
from a native SplitsTree SVG without changing the underlying network edges.

All three scripts use explicit input/output arguments. For example:

```bash
python scripts/network/prepare_splitstree_input.py \
  inputs/markers.faa inputs/pantoea92_markers.nex \
  --metadata inputs/metadata.tsv \
  --width 80

python scripts/network/make_taxon_representative_input.py \
  --distances inputs/pantoea92_JTT_distances.nex \
  --key inputs/Supplementary_Table_NeighborNet_taxon_key.tsv \
  --out-nexus outputs/pantoea28_taxon_representatives_JTT_distances.nex \
  --out-table outputs/pantoea28_taxon_representative_selection.tsv
```

SplitsTree itself is a Java application and may be run interactively. The
exact NEXUS input, `.stree6` project, software version and settings must be
included in the Code Ocean input package.
