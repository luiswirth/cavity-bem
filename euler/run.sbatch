#!/bin/bash
# === SLURM site config (ETH Euler) ===========================================
# Edit the --account, --constraint, and the `module load` line below for a
# different cluster; drop them entirely to run locally via ../run_local.sh.
#SBATCH --job-name=bem
#SBATCH --account=ls_math
#SBATCH --output=out/logs/slurm-%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=2G
#SBATCH --exclusive
#SBATCH --constraint=EPYC_7742
# =============================================================================
#
#   module load stack/2025-06 gcc/12.2.0 cmake/3.30.5 && cmake --build build -j8 --target my_project
#   euler/submit_grid.sh [geom]   # convergence grid
#   euler/submit_ref.sh  [geom]   # single high-fidelity reference run
set -euo pipefail
module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
GEOM=$1; MODE=$2; TASKFILE=$3
ROOT=$(pwd); R=$ROOT/out/$MODE/$GEOM; mkdir -p "$R" "$ROOT/out/logs"
read -r P M < <(grep -vE '^\s*#|^\s*$' "$ROOT/$TASKFILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
work=$R/work/p${P}_m${M}; rm -rf "$work"; mkdir -p "$work"; cd "$work"
srun --cpu-bind=cores /usr/bin/time -f "elapsed=%e mem=%M" -o time.txt "$ROOT/build/my_project" \
  operator "$ROOT/res/config_${GEOM}.txt" "$M" "$P" > run.log 2>&1
log=$(cat run.log); echo "$log"
dofs=$(echo "$log" | grep -oE 'dofs:[[:space:]]*[0-9]+' | grep -oE '[0-9]+')
cond=$(echo "$log" | grep -oE 'cond\(A\) = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')
secs=$(grep -oE 'elapsed=[0-9.]+' time.txt | grep -oE '[0-9.]+' | cut -d. -f1 || true)
mem_kb=$(grep -oE 'mem=[0-9]+' time.txt | grep -oE '[0-9]+' || true)
cp out/T_matrix.dat "$R/T_p${P}_m${M}.dat"
# manifest columns: P,M,dofs,secs,mem_kb,cond
echo "${P},${M},${dofs},${secs},${mem_kb},${cond}" >> "$R/manifest.csv"
echo "$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown),$(hostname),${P},${M},$(date -Is)" >> "$R/provenance.csv"
echo "done p${P}m${M}: dofs=$dofs cond=$cond secs=$secs mem_kb=$mem_kb"
