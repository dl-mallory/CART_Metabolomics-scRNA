#!/usr/bin/env bash
# Regenerate every panel of the figure from the staged BD Rhapsody counts.
# Stage 0 (FASTQ -> counts) is separate; see pipeline/00_bd_rhapsody/README.md.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
export R_LIBS="$PWD/env/rlib:${R_LIBS:-}"

for s in R/01_qc_annotate.R R/02_de_volcanoes.R R/03_gsea_published.R; do
    echo "=== $s"
    Rscript "$s" 2>&1 | tee "logs/$(basename "${s%.R}").log"
done

echo
echo "Panels written to results/figures/published/"
