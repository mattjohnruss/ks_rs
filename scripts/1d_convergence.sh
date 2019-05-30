#!/bin/bash

rm -f "1d_convergence/all.csv"

for n in 50 100 150 200 400 600 800 1000
do
    echo "Running with n_interior_cell_1d = $n"
    # --output_interval is chosen so that only the last timestep is output
    ./target/release/one_dim -n $n --dt 1e-5 --t_max 1 --ics exact --exact_solve --output_interval 100000 --dir 1d_convergence/"$n"
    echo -n "$n " >> 1d_convergence/all.csv
    cut -f4- -d' ' 1d_convergence/"$n"/trace.csv | tail -n1 >> 1d_convergence/all.csv
done
