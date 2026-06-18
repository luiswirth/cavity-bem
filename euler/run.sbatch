#!/bin/bash
#SBATCH --job-name=bem
#SBATCH --account=ls_math
#SBATCH --output=out/logs/slurm-%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=EPYC_7742
#
#   module load stack/2025-06 gcc/12.2.0 cmake/3.30.5 && cmake --build build -j8 --target my_project
#   mkdir -p out/logs && sbatch --array=1-20 euler/run.sbatch                    # grid
#   mkdir -p out/logs && sbatch --array=1    euler/run.sbatch euler/bem_ref.txt  # finer reference
set -euo pipefail
module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ROOT=$(pwd); R=$ROOT/out/ellipse; mkdir -p "$R" "$ROOT/out/logs"
TASKFILE=${1:-euler/bem_grid.txt}
read -r P M < <(grep -vE '^\s*#|^\s*$' "$ROOT/$TASKFILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
work=$R/work/p${P}_m${M}; rm -rf "$work"; mkdir -p "$work"; cd "$work"
t0=$SECONDS
srun --cpu-bind=cores /usr/bin/time -v -o time.txt "$ROOT/build/my_project" \
  "$ROOT/res/config_ellipse.txt" "$ROOT/res/sphere.dat" "$M" "$P" > run.log 2>&1
secs=$((SECONDS - t0)); log=$(cat run.log); echo "$log"
dofs=$(echo "$log" | grep -oE 'dofs:[[:space:]]*[0-9]+' | grep -oE '[0-9]+')
cond=$(echo "$log" | grep -oE 'cond\(A\) = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')
mem_kb=$(grep -i "Maximum resident set size" time.txt | grep -oE '[0-9]+' || true)
cp out/T_matrix.dat "$R/T_p${P}_m${M}.dat"
# manifest columns: P,M,dofs,secs,mem_kb,cond
echo "${P},${M},${dofs},${secs},${mem_kb},${cond}" >> "$R/manifest.csv"
echo "$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown),$(hostname),${P},${M},$(date -Is)" >> "$R/provenance.csv"
echo "done p${P}m${M}: dofs=$dofs cond=$cond secs=$secs mem_kb=$mem_kb"
