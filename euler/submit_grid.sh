#!/bin/bash
set -euo pipefail
if [[ $# -eq 0 ]]; then set -- ellipse sphere; fi
mkdir -p out/logs
for geom in "$@"; do
  sbatch --array=1-20%4 euler/run.sbatch "$geom" grid euler/grid.txt
done
