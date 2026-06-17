#!/bin/bash
#SBATCH --job-name=bem
#SBATCH --account=ls_math
#SBATCH --output=bem_results/slurm-%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --constraint=EPYC_7742
#
# Build first on a login node:
#   module load stack/2025-06 gcc/12.2.0 cmake/3.30.5
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
#   cmake --build build -j8 --target my_project
# Submit one (p,m) per array index of a task file (default euler/bem_grid.txt):
#   sbatch --array=1-20 euler/run.sbatch                 # full grid
#   sbatch --array=20 --mem-per-cpu=8G euler/run.sbatch  # single heavy task
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
/usr/bin/time -v -o "$work/time.txt" "$ROOT/build/my_project" \
  "$ROOT/res/config.txt" "$ROOT/res/sphere.dat" "$M" "$P" > "$work/run.log" 2>&1
secs=$((SECONDS-t0)); log=$(cat "$work/run.log"); echo "$log"
dofs=$(echo "$log" | grep -oE 'dofs:[[:space:]]*[0-9]+' | grep -oE '[0-9]+')
recip=$(echo "$log" | grep -oE '/ \|\|T\|\|_F = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')
maxrss=$(grep -i "Maximum resident set size" "$work/time.txt" | grep -oE '[0-9]+' || true)
cp out/T_matrix.dat "$R/T_bem_p${P}_m${M}.dat"
echo "ell,${P},${M},${dofs},${recip},${secs},${maxrss}" >> "$R/manifest.csv"
echo "done p${P}m${M}: dofs=$dofs recip=$recip secs=$secs maxrss_kb=$maxrss"
