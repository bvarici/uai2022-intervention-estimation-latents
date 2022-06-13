#!/bin/bash

modelno=1
for run in $(seq 11 1 30)
do
    for sample in 20000
    do
        python my_generate_data_cmd.py -r $run -s $sample -f "./Data_gauss/env-model-$modelno"
    done
done






