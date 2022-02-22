 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate data for the small graph example

"""

import numpy as np
import numpy.linalg as LA

from my_helpers import get_precision_cov
from my_functions_sample import sample_multiple
import argparse
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--run','-r',help="Run Number",type=int)
parser.add_argument('--samples','-s',help="Number of Samples",type=int)
parser.add_argument('--file','-f',help="File Name Prefix for Data to be used",type=str)

args=parser.parse_args()
run=args.run
n_samples=args.samples
file=args.file


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

# assign random edges
weights = np.random.uniform(-1,-0.25,[p,p])* np.random.choice([-1,1],size=[p,p])
weights = np.triu(weights)
np.fill_diagonal(weights,0)
#B_weighted = B*weights    
B = B*weights

mu = 0
variance = 1.0
shift = 1.0
plus_variance = 0.5

'observational case will be setting_0'
mu1 = mu*np.ones(p)
variance1 = variance*np.ones(p)
Omega1 = mu1**2 + variance1
B1 = B.copy()
Theta1, Cov1 = get_precision_cov(B1, np.diag(Omega1))

B_all = []
Theta_all = []
Cov_all = []
mu_all = []
variance_all = []
BObs = []
CovObs = []
ThetaObs = []

for i in range(len(intervention_targets)):
    I = intervention_targets[i]
    mu2 = mu*np.ones(p)
    variance2 = variance*np.ones(p)
    # shift the noise means
    mu2[I] += shift
    # increase the noise variances
    variance2[I] += plus_variance
    Omega2 = mu2**2 + variance2
    B2 = B.copy()    
    Theta2, Cov2 = get_precision_cov(B2, np.diag(Omega2))
    B_all.append(B2.copy())
    Theta_all.append(Theta2.copy())
    Cov_all.append(Cov2.copy())
    mu_all.append(mu2.copy())
    variance_all.append(variance2.copy())

    BObs.append(B2[O][:,O].copy())
    CovObs.append(Cov2[O][:,O].copy())
    ThetaObs.append(LA.inv(CovObs[-1]).copy())

# B_all, Theta_all, Cov_all, mu_all, variance_all, BObs, ThetaObs, CovObs

samples, covariances = sample_multiple(B_all,mu_all,variance_all,n_samples)

file_name=file+"-"+str(run)+"-"+str(samples)+"-"+"1.csv"

for i in range(len(samples)):
    file_name=file+"-"+str(run)+"-"+str(n_samples)+"-"+str(i)+".csv"
    np.savetxt(file_name, samples[i], delimiter=",", fmt="%.4f",
           header="L1, Z, Y, X", comments="")   



