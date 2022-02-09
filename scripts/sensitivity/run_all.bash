#!/bin/bash

n_runs=$1
res_dir="$2"
n_proc="16"

#echo "Building driver:"
#driver="one_dim_chemotaxis_7_var"
#cargo build --release --bin "$driver"

#echo "Running simulations:"
#time parallel -j"$n_proc" ./scripts/sensitivity/run_one.bash {1} "$res_dir" ::: $(seq 1 $n_runs) > /dev/null

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

echo "Combining simulation outputs:"

# Get header from first file
awk 'NR==1' "$res_dir"/0/trace.csv > "$res_dir"/all_outputs.csv

time for i in $(seq 0 $(expr $n_runs - 1))
do
    #cat res/$i/output.csv
    # get trace file at t = 1000 (line 1002)
    awk 'NR==302' "$res_dir"/$i/trace.csv
done >> "$res_dir"/all_outputs.csv
