#!/bin/bash
set -euo pipefail
module load stack/2025-06 gcc/12.2.0 cmake/3.30.5
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8 --target my_project
