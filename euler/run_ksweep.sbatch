#!/bin/bash
#SBATCH --job-name=bem-ksweep
#SBATCH --account=ls_math
#SBATCH --output=out/logs/slurm-%A.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=2G
#SBATCH --constraint=EPYC_7742
set -euo pipefail
module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
GEOM=$1
ROOT=$(pwd)
R=$ROOT/out/ksweep/$GEOM
rm -rf "$R"; mkdir -p "$R" "$ROOT/out/logs"
work=$R/work; rm -rf "$work"; mkdir -p "$work"; cd "$work"
srun --cpu-bind=cores "$ROOT/build/my_project" \
  ksweep "$ROOT/res/config_${GEOM}.txt" 3 3 1.9 2.1 400 > run.log 2>&1
cat run.log
cp out/ksweep.csv "$R/ksweep.csv"
echo "$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown),$(hostname),$(date -Is)" \
  > "$R/provenance.csv"
echo "wrote $R/ksweep.csv"
