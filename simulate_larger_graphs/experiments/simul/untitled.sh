#!/bin/bash

# RUN FOR JCI BURAK

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

pSysObs=40
outdir=./variance/gauss_JCI/p$pSysObs

for run_no in $(seq 12 1 20)
do
    fname=p$pSysObs-model$run_no
    # FCI alg: args[1] is basefilename, args[2] is fci, mode: args[3] is jci123, miniter is args[4], maxiter is args[5], alpha is args[6], jcifci_test is args[7]  
    ../Rcmd ../run.R $outdir/$fname fci jci123 0 0 0.05 gaussCItest

    ../Rcmd ../analyze.R $outdir/$fname 1
done
