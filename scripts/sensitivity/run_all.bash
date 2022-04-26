#!/bin/bash

n_runs=$1
res_dir="$2"
n_proc="16"

#echo "Building driver:"
#driver="one_dim_chemotaxis_7_var"
#cargo build --release --bin "$driver"

#echo "Running simulations:"
#time parallel -j"$n_proc" ./scripts/sensitivity/run_one.bash {1} "$res_dir" ::: $(seq 1 $n_runs) > /dev/null
