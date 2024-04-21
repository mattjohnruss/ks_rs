#!/bin/bash

n_cell="100"
t_max="150"
output_time_interval="0.1"
res_dir_base="$2"
dir="./"
ssd_threshold="1e-4"

driver="one_dim_chemotaxis_8_var"

driver_prog="../../target/release/$driver"

# subtract 1 to make run number zero indexed
n=$(expr $1 - 1)

# Switch to simulation directory (which should have been created by the R script already)
pushd "$res_dir_base/$n"

# Run this simulation
$driver_prog \
    --n_cell $n_cell \
    --t_max $t_max \
    --output_time_interval $output_time_interval \
    --dir "$dir" \
    --config "config.json" \
    --ssd_threshold $ssd_threshold

# Go back to original directory
popd
