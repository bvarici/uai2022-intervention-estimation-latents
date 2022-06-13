#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# run as run_large.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>

# run their data generating code.

pSysObs=$1
pContext=$2
eps=$3
eta=$4
N=$5
acyclic=$6
surgical=$7
run_no_first=$8
run_no_last=$9
#seed=$8
#iters=$9

outdir=./jci-results/p$pSysObs 

mkdir -p $outdir

for run_no in $(seq $run_no_first 1 $run_no_last)
do
    fname=p$pSysObs-model$run_no

    args="$pSysObs $pContext $eps $eta $N $acyclic $surgical"
    # use the run_no as seed
    ../Rcmd run_simul.R $outdir/$fname $args $run_no 

done