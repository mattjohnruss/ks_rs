#!/bin/bash

echo "Building driver:"
driver="one_dim_chemotaxis_7_var"
cargo build --release --bin "$driver"

n_runs=$1
n_proc="16"

echo "Running simulations:"
time parallel -j"$n_proc" ./scripts/sensitivity/run_one.bash {1} ::: $(seq 1 $n_runs) > /dev/null

# TODO we have yet to decide what quantities to extract from the simulations
# for the sensitivity analysis (the simple code commented out below is from a
# test problem), so for now don't do any output combining etc. We can do this
# later after running the simulations. See comments in the R script for other
# details.

#echo "Combining simulation outputs:"
#time for i in $(seq 0 $(expr $n_runs - 1))
#do
    #cat res/$i/output.csv
#done > all_outputs.csv
