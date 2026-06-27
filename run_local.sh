#!/bin/bash
# Run the BEM convergence grid (or reference) locally, without SLURM.
# Serial port of euler/run.sbatch for reproduction on any machine.
#
#   ./run_local.sh grid [geom...]   # 2D (p,m) convergence grid
#   ./run_local.sh ref  [geom...]   # single high-fidelity reference run
#
# geom defaults to: ellipse sphere. Build the binary first:
#   cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target my_project
#
# Note: high (p,m) need large dense complex matrices (the solver prints the
# estimate); run a subset by editing the taskfile or invoking build/my_project directly.
set -euo pipefail

MODE="${1:-grid}"; shift || true
case "$MODE" in
  grid) TASKFILE=euler/grid.txt;;
  ref)  TASKFILE=euler/ref.txt;;
  *) echo "usage: $0 grid|ref [geom...]" >&2; exit 1;;
esac
[[ $# -eq 0 ]] && set -- ellipse sphere

ROOT=$(cd "$(dirname "$0")" && pwd)
BIN="$ROOT/build/my_project"
[[ -x "$BIN" ]] || { echo "binary not built: $BIN (see README)" >&2; exit 1; }
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}"

for geom in "$@"; do
  R="$ROOT/out/$MODE/$geom"; mkdir -p "$R"
  rm -f "$R/manifest.csv" "$R/provenance.csv"
  while read -r P M; do
    work="$R/work/p${P}_m${M}"; rm -rf "$work"; mkdir -p "$work"
    (
      cd "$work"
      SECONDS=0
      "$BIN" operator "$ROOT/res/config_${geom}.txt" "$M" "$P" > run.log 2>&1 \
        || { echo "FAILED p${P}m${M}"; cat run.log; exit 1; }
      dofs=$(grep -oE 'dofs:[[:space:]]*[0-9]+' run.log | grep -oE '[0-9]+' | tail -1)
      cond=$(grep -oE 'cond\(A\) = [0-9.eE+-]+' run.log | grep -oE '[0-9.eE+-]+$')
      cp out/T_matrix.dat "$R/T_p${P}_m${M}.dat"
      # manifest columns: P,M,dofs,secs,mem_kb,cond  (mem_kb left empty locally)
      echo "${P},${M},${dofs},${SECONDS},,${cond}" >> "$R/manifest.csv"
      echo "$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown),$(hostname),${P},${M},$(date -Is)" >> "$R/provenance.csv"
      echo "done p${P}m${M}: dofs=$dofs cond=$cond secs=${SECONDS}"
    )
  done < <(grep -vE '^\s*#|^\s*$' "$ROOT/$TASKFILE")
done
