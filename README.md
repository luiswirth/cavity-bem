# cavity-bem

BEM reference solver for the interior PEC cavity reaction operator.
Built on [Bembel](https://github.com/temf/bembel). Runs on ETH Euler.

## Euler

    euler/build.sh                  # cmake configure + build (once)
    euler/submit.sh                 # sbatch 2D (p,m) grid -> out/ellipse/

Output: `out/ellipse/T_p{p}_m{m}.dat`, `manifest.csv`, `provenance.csv`.
Pull results into cavity-benchmark with `./pull-euler.sh`.
