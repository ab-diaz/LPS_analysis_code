#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUTS="${MANUSCRIPT_INPUTS:-${ROOT}/inputs}"
OUTPUTS="${MANUSCRIPT_OUTPUTS:-${ROOT}/outputs}"

"${ROOT}/reproducibility/run_checks.sh"
python "${ROOT}/scripts/selection_phase2/make_figure3_lps_gene_enrichment.py" \
  --analysis-dir "${INPUTS}" \
  --outdir "${OUTPUTS}/figures"

echo "Input checks and corrected Figure 3 output written to ${OUTPUTS}"
echo "This command does not rebuild every submitted figure. See reproducibility/README.md."
echo "The Gubbins and SplitsTree workflows are run separately as documented in"
echo "scripts/selection_phase2/run_gubbins_piersonii.sh and scripts/network/README.md."
