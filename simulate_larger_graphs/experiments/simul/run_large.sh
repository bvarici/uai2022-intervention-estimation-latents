#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# run as run_large.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>

echo run_large.sh called with arguments "$@"

pSysObs=$1
pContext=$2
eps=$3
eta=$4
N=$5
acyclic=$6
surgical=$7
seed=$8
iters=$9

outdir=./jci-results/jci_$pSysObs 

#outdir=/dev/shm/jmooij1/jci-paper/out/linear_stochastic_"$pSysObs"_"$pContext"_"$eps"_"$eta"_"$N"_"$acyclic"_"$surgical"_gaussCItest/"$seed"
#fname=simul-$seed
mkdir -p $outdir

for run_no in $(seq 11 1 11)
do
    fname=simul-$run_no

    args="$pSysObs $pContext $eps $eta $N $acyclic $surgical"
    # use the run_no as seed
    ../Rcmd run_simul.R $outdir/$fname $args $run_no 

    # fci
    # alg: args[2] is fci, mode: args[3] is jci123, miniter is args[4], maxiter is args[5], alpha is args[6], jcifci_test is args[7]  

    ../Rcmd ../run.R $outdir/$fname fci jci123 0 $iters 0.01 gaussCItest

    ../Rcmd ../analyze.R $outdir/$fname 1

done