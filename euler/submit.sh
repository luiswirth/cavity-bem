#!/bin/bash
set -euo pipefail
mkdir -p out/logs
sbatch --array=1-20 euler/run.sbatch ellipse
sbatch --array=1-20 euler/run.sbatch sphere
