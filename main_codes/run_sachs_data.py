#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run our algorithm on on Sachs et al. dataset (https://pubmed.ncbi.nlm.nih.gov/15845847/)

It consists of 11 nodes and 5 interventional setting. 

"""
import numpy as np
import pickle
import itertools as itr

from realdata.sachs.sachs_meta import SACHS_ESTIMATED_FOLDER, sachs_get_samples
from realdata.sachs.sachs_meta import true_dag_recent as sachs_gt
from functions_sample import IMAG_pasp_sample

sachs_skeleton = np.zeros(sachs_gt.shape)
sachs_skeleton[np.where(sachs_gt+sachs_gt.T)] = 1


def IMAG_sachs(S_all,max_subset_size=None, lambda_1=0.1, lambda_pasp=0.1, rho=1.0, th1=5e-3,\
        n_max_iter=500,stop_cond=1e-6,tol=1e-9,verbose=False,only_diag=True,return_pasp=True):
    '''

    Parameters
    ----------
    S_all : list of matrices
        list of sample covariance matrices.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need to change in most cases. The default is 1.0.
    th1 : float, optional
        small threshold to apply on Delta_hat. The default is 1e-3.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-6.
    tol : float, optional
        small number for thresholding. The default is 1e-9.
    verbose : Boolean
        The default is False.
    only_diag : Boolean, optional
        consider only the diagonal of Delta_Theta. The default is True.

    Returns
    -------
    est_edges : matrix
        estimated edges.
    est_skeleton : matrix
        estimated skeleton.
    K_hat_all : dict
        estimated K intervention targets.
    K_pasp_hat_all : dict
        estimated parents/spouses of intervention targets.
    S_Delta_size_all : list
        size of the affected nodes.
    time_all : list of float
        runtime.

    '''

    nnodes = S_obs.shape[0]

    K_hat_all = {}
    K_pasp_hat_all = {}
    S_Delta_size_all = {}
    time_all = {}
    
    setting_pairs = list(itr.combinations(np.arange(len(S_all)),2))
    n_F = len(setting_pairs)    

    for F in range(n_F):
        idx_1 = setting_pairs[F][0]
        idx_2 = setting_pairs[F][1]

        K_hat, K_pasp_hat, K_non_pasp_hat, J_hat, S_Delta_hat, N_filtered_hat, t_past = \
            IMAG_pasp_sample(S_all['setting_%d'%idx_1],S_all['setting_%d'%idx_2],max_subset_size,\
                lambda_1,lambda_pasp,rho,th1,n_max_iter,stop_cond,tol,verbose,only_diag,return_pasp) 


        K_hat_all['F_%d'%F] = K_hat
        K_pasp_hat_all['F_%d'%F] = K_pasp_hat
        #K_non_pasp_hat_all['setting_%d'%idx_setting] = K_non_pasp_hat
        #J_hat_all['setting_%d'%idx_setting] = J_hat
        S_Delta_size_all['F_%d'%F] = len(S_Delta_hat)
        #N_filtered_hat_all['setting_%d'%idx_setting] = N_filtered_hat
        time_all['F_%d'%F] = t_past
        print(F)

        
    # now combine the learned information for final causal structure
    est_edges = np.zeros((nnodes,nnodes))
    for F in range(n_F):
        edges_current = list(zip(K_hat_all['F_%d'%F],\
                                K_pasp_hat_all['F_%d'%F]))        
        for edge in edges_current:
            est_edges[edge[1],edge[0]] = 1

    est_skeleton = est_edges + est_edges.T
    est_skeleton[np.where(est_skeleton)] = 1

    return est_edges, est_skeleton, K_hat_all, K_pasp_hat_all, S_Delta_size_all, time_all 

#%% load the data
obs_samples, iv_samples_list, setting_list = sachs_get_samples()
    
# build the sufficient stats
S_obs = (obs_samples.T@obs_samples)/obs_samples.shape[0]
S_all = {}
# add the observational to the mix
S_all['setting_0'] = S_obs

for idx_setting in range(len(setting_list)):
    S_current = (iv_samples_list[idx_setting].T@iv_samples_list[idx_setting])/iv_samples_list[idx_setting].shape[0]
    S_all['setting_%d'%(idx_setting+1)] = S_current
    
#%%
# for S_Delta estimation
lambda_1_list = [0.3]
# other PDEs
lambda_pasp_list = [0.2]
# remove small values
th1_list = [0.01]
rho = 1.0
n_max_iter = 500
stop_cond = 1e-6
tol = 1e-9
verbose = False
only_diag = False
return_pasp = True
max_subset_size = None

parameters_list = \
    list(itr.product(lambda_1_list,lambda_pasp_list,th1_list))


results = {}
for parameters in parameters_list:
    results[parameters] = {}
    est_edges, est_skeleton, K_hat_all, K_pasp_hat_all ,S_Delta_size_all, time_all =\
        IMAG_sachs(S_all,max_subset_size,\
        parameters[0],parameters[1],rho,parameters[2],n_max_iter,stop_cond,tol,verbose,only_diag,return_pasp)


    results[parameters]['est_edges'] = est_edges
    results[parameters]['est_skeleton'] = est_skeleton
    results[parameters]['K_hat'] = K_hat_all
    results[parameters]['K_pasp_hat'] = K_pasp_hat_all
    results[parameters]['time'] = time_all
    results[parameters]['S_Delta_size'] = S_Delta_size_all
    print(parameters)
#%%
f = open(SACHS_ESTIMATED_FOLDER+'/sachs_run_2.pkl','wb')
pickle.dump(results,f)
f.close()
