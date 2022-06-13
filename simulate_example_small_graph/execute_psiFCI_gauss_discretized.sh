#!/bin/bash

modelno=1
for run in $(seq 11 1 30)
do
    for sample in 20000
    do
        # discCItest on discretized_gaussian data
        Rscript test_psiFCI-commandline.R 3 "./psiFCI_outputs_gauss_discretized/graph_disc_model-$modelno-alpha0.05-$run-$sample.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-0.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-1.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-2.csv"
        python3 Graphplot-cmd.py -fn "./psiFCI_outputs_gauss_discretized/graph_disc_model-$modelno-alpha0.05-$run-$sample.csv" -fo "./psiFCI_outputs_gauss_discretized/Graphdiag_model-$modelno-alpha0.05-$run-$sample"

    done
done











