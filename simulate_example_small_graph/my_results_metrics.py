#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 18:35:39 2021

@author: Burak

my_results_analyze

"""
import numpy as np


## 0: no edge
## 1: -o
## 2: -> (arrowhead)
## 3: - (tail)

node_names = ['Z','Y','X','F_12','F_13','F_23']

A_true = np.zeros([6,6])
A_true[3,1] = 2; A_true[3,2] = 2; A_true[4,0] = 2; A_true[4,1] = 2; A_true[4,2] = 2; A_true[5,0] = 2;
A_true[1,3] = 3; A_true[2,3] = 3; A_true[0,4] = 3; A_true[1,4] = 3; A_true[2,4] = 3; A_true[0,5] = 3;
A_true[1,2] = 2; A_true[2,1] = 3; A_true[1,0] = 2; A_true[0,1] = 2

A_true_skeleton_ZYX = A_true[:3][:,:3].copy()
A_true_skeleton_ZYX[np.where(A_true_skeleton_ZYX)] = 1
A_true_ZYX = A_true[:3][:,:3].copy()
A_true_F = A_true[3:].copy()
A_true_F[np.where(A_true_F)] = 1

true_directed_edges = list(zip(np.where((A_true==2)*(A_true.T==3))[0],np.where((A_true==2)*(A_true.T==3))[1]))
true_bidirected_edges = list(zip(np.where((A_true+A_true.T)==4)[0],np.where((A_true+A_true.T)==4)[1]))
true_bidirected_edges = [edge for edge in true_bidirected_edges if edge[0]<edge[1]]
true_joker_edges = list(zip(np.where((A_true==2)*(A_true.T==1))[0],np.where((A_true==2)*(A_true.T==1))[1]))

TP_skeleton_ZYX = 0
FP_skeleton_ZYX = 0
FN_skeleton_ZYX = 0
TP_ZYX = 0
FP_ZYX = 0
FN_ZYX = 0

TP_F = 0
FP_F = 0
FN_F = 0

for run_no in range(11,31):
    for samples in [20000]:
        
        A = np.genfromtxt("./my_outputs_gauss_discretized/graph-env-model-1-"+str(run_no)+"-"+str(samples)+".csv",delimiter=',')
        A_skeleton_ZYX = A[:3][:,:3].copy()
        A_skeleton_ZYX[np.where(A_skeleton_ZYX)] = 1
        A_ZYX = A[:3][:,:3].copy()

        # check the skeleton among the observed nodes
        TP_skeleton_ZYX += np.sum(A_true_skeleton_ZYX*A_skeleton_ZYX)/2
        FP_skeleton_ZYX += np.sum(A_skeleton_ZYX)/2 - np.sum(A_true_skeleton_ZYX*A_skeleton_ZYX)/2
        FN_skeleton_ZYX += np.sum(A_true_skeleton_ZYX)/2 - np.sum(A_true_skeleton_ZYX*A_skeleton_ZYX)/2
                
        # check the edges among the observed nodes
        instant_tp = np.sum((A_true_ZYX == A_ZYX)* (A_true_ZYX.T == A_ZYX.T) * (A_true_ZYX>0))/2
        TP_ZYX += instant_tp
        FP_ZYX += np.sum(A_ZYX>0)/2 - instant_tp
        FN_ZYX += np.sum(A_true_ZYX>0)/2 - instant_tp

        # finally, check the F-edges, our main thing
        A_F = A[3:].copy()
        A_F[np.where(A_F)] = 1        
        
        TP_F += np.sum(A_true_F*A_F)
        FP_F += np.sum(A_F) - np.sum(A_true_F*A_F)
        FN_F += np.sum(A_true_F) - np.sum(A_true_F*A_F)