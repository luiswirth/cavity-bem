#!/bin/bash
set -euo pipefail
mkdir -p out/logs
for geom in "${@:-ellipse sphere}"; do
  sbatch --array=1-20 euler/run.sbatch "$geom" grid euler/grid.txt
done
