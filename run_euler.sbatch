#!/bin/bash
#SBATCH --job-name=cavity-bem
#SBATCH --output=out/slurm-%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32        # assembly is OpenMP-parallel; cores are the lever
#SBATCH --mem-per-cpu=4G          # 128 GB total (generous; matrix is small, ~MB)
#
# Converged high-refinement BEM for the (4,4,6) ellipsoidal cavity.
# Build first on a LOGIN node (FetchContent needs internet):
#   module load stack/2025-06 gcc/12.2.0 cmake/3.30.5
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j8 --target my_project
# Then submit:  sbatch run_euler.sbatch
#
# Tip: start by running refinement 3 only (edit REFINEMENTS below) to calibrate
# wall-time and memory from the logged DOF count, then add 4.

set -euo pipefail

module load stack/2025-06 gcc/12.2.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

CONFIG=res/config.txt          # must have header: 2.0 4.0 4.0 6.0 12  (ellipse)
GEOM=res/sphere.dat
POLY_DEG=2
REFINEMENTS="3 4"

mkdir -p out
echo "host=$(hostname)  threads=$OMP_NUM_THREADS  config=$CONFIG"
head -n 2 "$CONFIG"

for r in $REFINEMENTS; do
  echo "=== refinement $r ==="
  ./build/my_project "$CONFIG" "$GEOM" "$r" "$POLY_DEG"
  cp out/T_matrix.dat "out/T_r${r}.dat"
  echo "wrote out/T_r${r}.dat"
done
echo "done"
