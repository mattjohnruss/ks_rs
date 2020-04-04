#/bin/bash

n_files=1000
line_number=101

while [ "$1" != "" ]; do
    case $1 in
        -f ) shift
            n_files=$1
            ;;
        -i ) shift
            line_number=$1
            ;;
    esac
    shift
done

for i in $(seq -w 00000 $n_files)
do
    tail -n+"$line_number" ../output_"$i".csv | head -n1
done > "$line_number".csv
