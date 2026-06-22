#!/bin/bash
set -euo pipefail
if [[ $# -eq 0 ]]; then set -- ellipse sphere; fi
mkdir -p out/logs
for geom in "$@"; do
  sbatch euler/run_ksweep.sbatch "$geom"
done
