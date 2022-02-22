#!/bin/bash

pContext=1

for run in $(seq 21 1 40)
do
    for pSystem in 10 20 30 40 50
    do
        for sample in 5000
        do
            python my_generate_data.py -r $run -s $sample -pSystem $pSystem -pContext $pContext -l $pSystem -i 3 -c 2
        done
    done
done






