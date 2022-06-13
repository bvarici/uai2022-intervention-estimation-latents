#!/bin/bash

modelno=1
for run in $(seq 11 1 11)
do
    for sample in 20000
    do
        # discCItest on discretized_gaussian data
        Rscript test_JCI-commandline.R 3 "./JCI_outputs_gauss_discretized/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-0.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-1.csv" "./Data_gauss_discretized/env-model-$modelno-$run-$sample-2.csv"
        python3 Graphplot-cmd.py -fn "./JCI_outputs_gauss_discretized/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" -fo "./JCI_outputs_gauss_discretized/Graphdiag_model-$modelno-jci-alpha0.05-$run-$sample"

        # # discCItest on gauss_data
        # Rscript test_JCI-commandline.R 3 "./gauss_data/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" "./gauss_data/env-model-$modelno-$run-$sample-0.csv" "./Data/env-model-$modelno-$run-$sample-1.csv" "./Data/env-model-$modelno-$run-$sample-2.csv" disCItest
        # python3 Graphplot-cmd.py -fn "./gauss_data/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" -fo "./gauss_data/Graphdiag_model-$modelno-jci-alpha0.05-$run-$sample"

        # # gaussCItest on gauss_data
        # Rscript test_JCI-commandline.R 3 "./gauss_data/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" "./gauss_data/env-model-$modelno-$run-$sample-0.csv" "./Data/env-model-$modelno-$run-$sample-1.csv" "./Data/env-model-$modelno-$run-$sample-2.csv" gaussCItest
        # python3 Graphplot-cmd.py -fn "./gauss_data/graph_disc_model-$modelno-jci-alpha0.05-$run-$sample.csv" -fo "./gauss_data/Graphdiag_model-$modelno-jci-alpha0.05-$run-$sample"

    done
done











