#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run our algorithm on discretized gaussian data from the small graph.

Reviewed on February, 2022.
"""

import numpy as np

from functions_sample import IMAG_pasp_sample_multiple
#node_names=['X','Y','Z','L1']

# create the 3-node + 1-latent node graph
# indices: L1 = 0, Z = 1, Y = 2, X = 3
# edges: L1 -> Y, L1 -> Z, Y -> X
B = np.zeros([4,4])
B[0,1] = 1
B[0,2] = 1
B[2,3] = 1
p = len(B)

O = [1,2,3]

intervention_targets = [[],[2,3],[1,2,3]]

#%%
# parameters for the algorithm
rho = 1.0
lambda_1 = 0.05 # for S_Delta estimation
lambda_2 = 0.05 # other PDEs
lambda_pasp = 0.3
th1 = 5e-3 # throw-away very small ones

K_all = []
K_pasp_all = []
A_all = []

# 0: no edge
# 1: -o
# 2: -> (arrowhead)
# 3: - (tail)
# build the small graph
A_true = np.zeros([6,6])
A_true[3,1] = 2; A_true[3,2] = 2; A_true[4,0] = 2; A_true[4,1] = 2; A_true[4,2] = 2; A_true[5,0] = 2;
A_true[1,3] = 3; A_true[2,3] = 3; A_true[0,4] = 3; A_true[1,4] = 3; A_true[2,4] = 3; A_true[0,5] = 3;
A_true[1,2] = 2; A_true[2,1] = 3; A_true[1,0] = 2; A_true[0,1] = 2

for run_no in range(11,31):
    for samples in [20000]:
        A = np.zeros([6,6]) # 3-node + 3 F-nodes:  Z, Y, X, F_12, F_13, F_23,
        # ground truth should be F_12 -> Y,X; F_13 -> Z,Y,X; F_23 -> Z, Y->X, Y<->Z
        # load the data
        X_all = []
        S_all = []
        S_Obs = []
        for i in range(3):
            X = np.genfromtxt("./Data_gauss_discretized/env-model-1-"+str(run_no)+"-"+str(samples)+"-"+str(i)+".csv",delimiter=',', skip_header = 1)
            X_all.append(X)

        X_Obs = [X_all[i][:,O] for i in range(len(X_all))]
        S_Obs = [X_Obs[i].T@X_Obs[i]/len(X) for i in range(len(X_Obs))]                 

        K_hat, K_pasp_hat, K_non_pasp_hat, J_hat, S_Delta_sizes, N_filtered_hat, t_past =  \
            IMAG_pasp_sample_multiple(S_Obs,lambda_1=lambda_1,lambda_2=lambda_2,lambda_pasp=lambda_pasp,\
                     rho=rho,th1=th1,tol=1e-9,verbose=False,only_diag=False,return_pasp=True)
                
        for i in range(3):
            # F-nodes to intervened nodes
            A[i+3,K_hat[i]] = 2
            A[K_hat[i],i+3] = 3
            for j in range(len(K_hat[i])):
                A[K_pasp_hat[i][j],K_hat[i][j]] = 2
                for k in range(len(K_pasp_hat[i][j])):
                    if A[K_hat[i][j],K_pasp_hat[i][j][k]] != 2 :
                        A[K_hat[i][j],K_pasp_hat[i][j][k]] = 1

                
        K_all.append(K_hat)
        K_pasp_all.append(K_pasp_hat)
        A_all.append(A)
        print(run_no,'K:',K_hat)
        print(run_no,'K_pasp:',K_pasp_hat)

        #np.savetxt("./my_outputs_gauss_discretized/graph-env-model-1-"+str(run_no)+"-"+str(samples)+".csv",A,delimiter=",", fmt="%d")    
