#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate gaussian data for larger nodes, save. also save in JCI format.

soft interventions are increased variance type as in our formulation.

Reviewed on February 2022.
"""

import numpy as np
from my_helpers import create_random_SEM, create_multiple_intervention
from my_functions_sample import sample_multiple
from my_functions_graph import reduce_2_observed_multiple
import argparse
#import pandas as pd
import pickle
import json
import os

def listToString(s): 

    # initialize an empty string
    str1 = "" 

    # traverse in the string  
    for ele in s: 
        str1 += ele + ", " 

    # return string  
    return str1[:-2]

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


parser=argparse.ArgumentParser()
parser.add_argument('--run','-r',help="Run Number",type=int)
parser.add_argument('--samples','-s',help="Number of Samples",type=int)
parser.add_argument('--pSystem','-pSystem',help="num system variables",type=int)
parser.add_argument('--pContext','-pContext',help="num context variables, i.e., num int. setting",type=int)
parser.add_argument('--l','-l',help="num latent nodes",type=int)
parser.add_argument('--i','-i',help="num intervetion targets per setting",type=int)
parser.add_argument('--density','-c',help="density",type=float)


args=parser.parse_args()
run=args.run
n_samples=args.samples
pSystem = args.pSystem
pContext = args.pContext
num_latent = args.l
num_int = args.i
c = args.density

p = pSystem + num_latent
num_total_settings = pContext + 1
print(pSystem)
print(args)

#base_gauss_JCI = './experiments/simul/jci-mydata-more-latents-more-targets/p'+str(pSystem)
base_gauss_JCI = './experiments/simul/jci-mydata/p'+str(pSystem)
base_gauss_JCI_instant = base_gauss_JCI+'/p'+str(pSystem)+'-'+'model'+str(run)

if os.path.exists(base_gauss_JCI) is False:
    os.makedirs(base_gauss_JCI)

#%%
'create soft-interventions where sigma_i^2 is increased by 1'
mu = 0
variance = 1.0
shift = 1.0
plus_variance = 0.0


B = create_random_SEM(p=p,density=c,weighted=True)
B_all, Theta_all, Cov_all, mu_all, variance_all, BObs, ThetaObs, CovObs, L, O, I_all = \
    create_multiple_intervention(B,L_size=num_latent,I_size=num_int,n_interventions=pContext,\
        shift=shift,plus_variance=plus_variance)

# ground truth IMAG.
IMAG, ThetaObs_all, CovObs_all = reduce_2_observed_multiple(B_all, Theta_all, Cov_all, L, I_all)
n_F = int(len(B_all)*(len(B_all)-1)/2)
K = [list(np.where(IMAG[i])[0]-n_F) for i in range(n_F)]
print('K:',K)
'save graph and intervention targets'
info = {'B':B,'IMAG':IMAG,'L':L,'O':O,'K':K,'mu':mu, 'variance':variance, 'shift':shift, 'plus_variance':plus_variance, 'density':c,}
f = open(base_gauss_JCI_instant+'-info.pkl','wb')
pickle.dump(info,f)
f.close()
'also save the IMAG in JCI format'
idx = list(np.arange(n_F,len(IMAG)))+list(np.arange(n_F))
IMAG_node_names = ['V'+str(i) for i in range(1,len(IMAG)+1)]
np.savetxt(base_gauss_JCI_instant+'-edge.csv',IMAG[idx][:,idx],delimiter=",",header=listToString(IMAG_node_names),comments="")


'generate Gaussian data and save'
samples, covariances = sample_multiple(B_all,mu_all,variance_all,n_samples)
# need to generate header (JCI code needs them)
# use L1, L2... for latents, V1, V2, ... for observed ones
all_node_names = [[] for i in range(p)]
L_idx = 1
V_idx = 1
for i in range(p):
    if i in L:
        all_node_names[i] = "L"+str(L_idx)
        L_idx +=1 
    else:
        all_node_names[i] = "V"+str(V_idx)
        V_idx +=1

for i in range(len(samples)):
    file_name=base_gauss_JCI_instant+'-'+str(n_samples)+"-"+str(i)+"-full-data.csv"
    np.savetxt(file_name, samples[i], delimiter=",", header=listToString(all_node_names),fmt="%.4f",comments="")   


#%%
'save in JCI format'
obs_node_names = ["V"+str(i+1) for i in range(pSystem+pContext)]

# save only observed values
X_obs = []
for i in range(num_total_settings):
    X_obs.append(samples[i][:,O])

X_JCI = np.zeros((n_samples*num_total_settings,pSystem+pContext),dtype=object)

for i in range(num_total_settings):
    X_JCI[i*n_samples:(i+1)*n_samples,:pSystem] = (X_obs[i]).astype(float)
    if i==0:
       X_JCI[i*n_samples:(i+1)*n_samples,pSystem] = int(0) 
    if i>0:
        X_JCI[i*n_samples:(i+1)*n_samples,pSystem+i-1] = int(1)
        

np.savetxt(base_gauss_JCI_instant+'-data.csv',X_JCI,delimiter=',',header=listToString(obs_node_names),comments='')    

'save metadata'
metadata = {}
metadata['SystemVars'] = list(np.arange(1,pSystem+1))
metadata['ContextVars'] = list(np.arange(pSystem+1,pSystem+pContext+1))
metadata['basefilename'] = base_gauss_JCI
metadata['PSysObs'] = int(pSystem)
metadata['PSysConf'] = int(num_latent)
metadata['pContext'] = int(pContext)
metadata['N'] = int(n_samples)
metadata['targets'] = K


with open(base_gauss_JCI_instant+'-metadata.json', 'w') as outfile:
    json.dump(metadata, outfile,cls=NpEncoder)


