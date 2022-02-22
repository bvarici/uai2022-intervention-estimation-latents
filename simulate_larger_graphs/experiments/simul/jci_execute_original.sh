#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# run as run_large.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>

# modelfirst and modellast: the id of the instances to run between
# run their data generating code.

pSysObs=$1
modelfirst=$2
modellast=$3

outdir=./jci-results/p$pSysObs 

mkdir -p $outdir

for run_no in $(seq $modelfirst 1 $modellast)
do
    fname=p$pSysObs-model$run_no

    # use the run_no as seed
    # alg: args[2] is fci, mode: args[3] is jci123, miniter is args[4], maxiter is args[5], alpha is args[6], jcifci_test is args[7]  

    ../Rcmd ../run.R $outdir/$fname fci jci123 0 0 0.01 gaussCItest

    ../Rcmd ../analyze.R $outdir/$fname 1

done