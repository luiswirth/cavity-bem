#!/bin/bash
#SBATCH --job-name=bem
#SBATCH --account=ls_math
#SBATCH --output=bem_results/slurm-%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=EPYC_7742
#
# Build first on a login node:
#   module load stack/2025-06 gcc/12.2.0 cmake/3.30.5
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
#   cmake --build build -j8 --target my_project
# Submit one (p,m) per array index of a task file (default euler/bem_grid.txt):
#   sbatch --array=1-20 euler/run.sbatch                    # convergence grid
#   sbatch --array=1    euler/run.sbatch euler/bem_ref.txt  # a separate finer reference
# Then post-process locally (aggregate reads the manifest + saved T's). Mirrors
# the EP-GP side: sbatch the array, pull results, aggregate locally.
# BEM peak RSS is large (~200 GiB at a high-p corner: Bembel superspace assembly
# plus the LU copy), so 16G/cpu = 256 GiB; it fits the ~489 GiB EPYC_7742 nodes.
set -euo pipefail
module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
ROOT=$(pwd); R=$ROOT/bem_results; mkdir -p "$R"
TASKFILE=${1:-euler/bem_grid.txt}
line=$(grep -vE '^\s*#|^\s*$' "$ROOT/$TASKFILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
read -r P M <<< "$line"
work=$R/work/p${P}_m${M}; rm -rf "$work"; mkdir -p "$work"; cd "$work"
t0=$SECONDS
# /usr/bin/time -v captures peak resident set size; binary stdout goes to run.log,
# the time report to time.txt (kept separate so the log parsing is unaffected).
export OMP_PROC_BIND=close OMP_PLACES=cores
srun --cpu-bind=cores /usr/bin/time -v -o "$work/time.txt" "$ROOT/build/my_project" \
  "$ROOT/res/config.txt" "$ROOT/res/sphere.dat" "$M" "$P" > "$work/run.log" 2>&1
secs=$((SECONDS-t0)); log=$(cat "$work/run.log"); echo "$log"
dofs=$(echo "$log" | grep -oE 'dofs:[[:space:]]*[0-9]+' | grep -oE '[0-9]+')
# cond(A) = 1-norm condition estimate of the EFIE system matrix; NOT a function
# of the saved T, so it is recorded here. recip and ||T|| ARE functions of T and
# are computed in post-processing (cavity-dipoles aggregate), not recorded.
cond=$(echo "$log" | grep -oE 'cond\(A\) = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')
maxrss_kb=$(grep -i "Maximum resident set size" "$work/time.txt" | grep -oE '[0-9]+' || true)
cp out/T_matrix.dat "$R/T_bem_p${P}_m${M}.dat"
# Geometry is fixed by the output dir (out/ellipse_bem) + provenance, so it is
# not a manifest column. Columns: P, M, dofs, secs, maxrss_kb, cond.
echo "${P},${M},${dofs},${secs},${maxrss_kb},${cond}" >> "$R/manifest.csv"
commit=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
echo "${commit},$(hostname),${P},${M},$(date -Is),${SLURM_CPUS_PER_TASK}" >> "$R/provenance.csv"
echo "done p${P}m${M}: dofs=$dofs cond=$cond secs=$secs maxrss_kb=$maxrss_kb"
