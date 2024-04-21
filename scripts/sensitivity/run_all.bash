#!/bin/bash

n_runs=$1
res_dir_base="$2"
n_proc="16"

echo "Building driver:"
driver="one_dim_chemotaxis_8_var"
cargo build --release --bin "$driver"

echo "Running simulations:"
time parallel --eta -j"$n_proc" ./scripts/sensitivity/run_one.bash {1} "$res_dir_base" ::: $(seq 1 $n_runs) > /dev/null
