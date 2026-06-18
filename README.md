# cavity-bem

BEM reference solver for the interior PEC cavity reaction operator.
Built on [Bembel](https://github.com/temf/bembel). Runs on ETH Euler.

## Euler

    euler/build.sh                  # cmake configure + build (once)
    euler/submit_grid.sh [geom]     # sbatch 2D (p,m) convergence grid
    euler/submit_ref.sh  [geom]     # sbatch single high-fidelity reference run

Output: `out/{grid,ref}/{shape}/T_p{p}_m{m}.dat`, `manifest.csv`, `provenance.csv`.
Pull results into cavity-benchmark with `./pull-euler.sh`.
