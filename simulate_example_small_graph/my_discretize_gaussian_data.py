#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
discretize data for the small graph

Reviewed on February, 2022.
"""

import numpy as np

node_names=['X','Y','Z','L1']

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

# load the data
X_all = []
S_all = []
S_Obs = []

#%%
for run_no in range(11,31):
    for samples in [20000]:
        X_all = []
        for i in range(3):
            X = np.genfromtxt("./Data_gauss/env-model-1-"+str(run_no)+"-"+str(samples)+"-"+str(i)+".csv",delimiter=',', skip_header = 1)
            X_all.append(X)
 
        bins = np.percentile(X_all,[20,40,60,80])
        X_disc = np.digitize(X_all,bins) - 2           
        
        for i in range(3):
            np.savetxt("./Data_gauss_discretized/env-model-1-"+str(run_no)+"-"+str(samples)+"-"+str(i)+".csv",X_disc[i],delimiter=",", fmt="%d",
                   header="L1, Z, Y, X", comments="")    

    print(run_no)
