#!/bin/bash
#SBATCH --job-name=bem-grid
#SBATCH --output=bem_results/slurm-%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-8
#
# BEM convergence cross for the (4,4,6) ellipsoid, N=32 dipoles.
# One array task per (poly_deg, refinement) line in bem_tasks.txt.
# Run from ~/cavity-bem after building (see run_euler.sbatch header).
#   sbatch run_bem_grid.sbatch
#
# Each task runs in its own workdir so the binary's out/ dirs don't collide,
# then copies the operator to bem_results/ and appends one manifest row.

set -euo pipefail
module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

ROOT=$(pwd)                       # ~/cavity-bem
BIN=$ROOT/build/my_project
CONFIG=$ROOT/res/config.txt       # ellipsoid, N=32
GEOM=$ROOT/res/sphere.dat
RESULTS=$ROOT/bem_results
mkdir -p "$RESULTS"

# pick the task line (skip comments/blanks)
line=$(grep -vE '^\s*#|^\s*$' "$ROOT/bem_tasks.txt" | sed -n "${SLURM_ARRAY_TASK_ID}p")
read -r P M <<< "$line"
echo "task ${SLURM_ARRAY_TASK_ID}: poly_deg=$P refinement=$M"

work=$ROOT/bem_results/work/p${P}_m${M}
rm -rf "$work"; mkdir -p "$work"; cd "$work"

t0=$SECONDS
log=$("$BIN" "$CONFIG" "$GEOM" "$M" "$P")
secs=$((SECONDS - t0))
echo "$log"

dofs=$(echo "$log" | grep -oE 'dofs:[[:space:]]*[0-9]+' | grep -oE '[0-9]+')
recip=$(echo "$log" | grep -oE '/ \|\|T\|\|_F = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')

cp out/T_matrix.dat "$RESULTS/T_bem_ell_p${P}_m${M}.dat"
echo "ell,${P},${M},${dofs:-NA},${recip:-NA},${secs}" >> "$RESULTS/manifest.csv"
echo "done p${P}_m${M}: dofs=${dofs} recip=${recip} secs=${secs}"
