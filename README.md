# cavity-bem

Boundary element method (BEM) reference solver for the interior PEC cavity
reaction operator. It solves the indirect electric-field integral equation
(single-layer EFIE) on an ellipsoidal cavity, assembles the dipole
transmit-to-receive reaction operator `T`, and writes it to disk. Built on the
[Bembel](https://github.com/temf/bembel) isogeometric BEM library.

This is the deterministic baseline against which the EPGP solver
([cavity-maxwellgp](https://github.com/luiswirth/cavity-maxwellgp)) is cross-validated.
The cavity geometry is set entirely by the semi-axes in `res/config_{shape}.txt`;
the same base sphere NURBS mesh (`res/sphere.dat`) is scaled accordingly.

## Requirements

- A C++23 compiler and CMake >= 3.25 (Ninja recommended)
- OpenMP (on macOS: `brew install libomp`)

Eigen 3.4 and Bembel v1.0 are fetched automatically by CMake; no manual install.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j --target my_project
```

On macOS the bundled preset wires up Homebrew's libomp:

```bash
cmake --preset default && cmake --build build --target my_project
```

## Run locally

`run_local.sh` runs the convergence grid serially, no scheduler:

```bash
./run_local.sh grid            # 2D (p, m) grid, both shapes
```

The taskfile `euler/grid.txt` (columns `poly_deg refinement`) defines the runs.
The reference operator is the most refined grid run, not a separate run. To run a
single point directly:

```bash
build/my_project operator res/config_ellipse.txt <refinement> <poly_deg>
```

High `(p, m)` build large dense complex matrices (the solver prints the size
estimate); pick a subset for a laptop.

## Run on a cluster (ETH Euler)

```bash
euler/build.sh                 # module load + cmake build
euler/submit_grid.sh [shape]   # sbatch 2D (p, m) convergence grid
```

The SLURM account, node constraint, and `module load` lines in `euler/run.sbatch`
are Euler-specific and flagged at the top of that file; adjust them for another
SLURM site.

## Output

Per shape, under `out/grid/{shape}/`:

- `T_p{p}_m{m}.dat` reaction operator (whitespace `real imag` per entry)
- `manifest.csv` columns `p,m,dofs,secs,mem_kb,cond`
- `provenance.csv` git commit, host, parameters, timestamp

Collect into the benchmark harness with `cavity-benchmark/pull-euler.sh`.
