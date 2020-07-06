#/bin/bash

n_files=1000
n_cells=100
n_procs=4

while [ "$1" != "" ]; do
    case $1 in
        -f ) shift
            n_files=$1
            ;;
        -c ) shift
            n_cells=$1
            ;;
        -j ) shift
            n_procs=$1
            ;;
    esac
    shift
done

lines=$(seq -s " " 2 3 $((3 * $n_cells + 1)))

# FIXME total hack to hardcode this...
parallel -j"$n_procs" --progress /home/pmzmr/Code/rust/numerics/ks_rs/scripts/point_values_vs_time.sh -f "$n_files" -i {} ::: $lines
