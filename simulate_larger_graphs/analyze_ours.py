#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze results of our algorithm algo on larger graphs 

Reviewed on February 2022.
"""

import numpy as np
from helpers import counter
from functions_sample import IMAG_pasp_sample_multiple
#from config import SIMULATIONS_ESTIMATED_FOLDER
#import pickle

#JCI_data_folder = './experiments/simul/jci-results'
#JCI_data_folder = './experiments/simul/jci-mydata'


'following code is only for 2 total contexts.'

def run_our_algo(p,first,last,mode,n_samples=5000):
    '''
    Parameters
    ----------
    p : int
        size of the graph to read.
    first : int
        first model run_no to read.
    last : int
        last model run_no to read.
    mode : int
        mode=1: jci-original data. mode=2: jci-mydata
    n_samples : int, optional
        number of samples. The default is 5000.

    Returns
    -------
    I_tp : TYPE
        DESCRIPTION.
    I_fp : TYPE
        DESCRIPTION.
    I_fn : TYPE
        DESCRIPTION.

    '''
    if mode == 1:
        JCI_data_folder = './experiments/simul/jci-results'
    elif mode == 2:
        JCI_data_folder = './experiments/simul/jci-mydata'    
    
    pContext = 1
    I_tp = 0
    I_fp = 0
    I_fn = 0
    t_all = []
    
    K_true_all = []

    for run_no in range(first,last+1):
        # read data samples
        data_filename = JCI_data_folder + '/p'+str(p)+'/p'+str(p)+'-model'+str(run_no)+'-data.csv'
        X = np.genfromtxt(data_filename,delimiter=',',skip_header=1)
        
        X = X[:,:p]
        n_total_settings = int(X.shape[0]/n_samples)
        X_all = []
        for i in range(n_total_settings):
            X_all.append(X[i*n_samples:(i+1)*n_samples])
            
        # read the graph
        edge_filename = JCI_data_folder + '/p'+str(p)+'/p'+str(p)+'-model'+str(run_no)+'-edge.csv'
        A_true = np.genfromtxt(edge_filename,delimiter=',',skip_header=1)
        
        K_true = [list(np.where(A_true[-1])[0])]
        K_true_all.append(K_true)
        
        # hyperparameters
        rho = 1.0
        lambda_1 = 0.1 # for S_Delta estimation
        lambda_2 = 0.1 # other PDEs
        lambda_pasp = 0.05
        th1 = 1e-4 # throw-away very small ones
        
        S_Obs = [X_all[i].T@X_all[i]/len(X) for i in range(len(X_all))]                 
        
        #run the algorithm
        K_hat, K_pasp_hat, K_non_pasp_hat, J_hat, S_Delta_sizes, N_filtered_hat, t_past =  \
            IMAG_pasp_sample_multiple(S_Obs,lambda_1=lambda_1,lambda_2=lambda_2,lambda_pasp=lambda_pasp,\
                     rho=rho,th1=th1,tol=1e-9,verbose=False,only_diag=True,return_pasp=False)
                
        t_all.append(t_past)        
                        
        I_tp_run, I_fp_run, I_fn_run = counter(K_true,K_hat)
        I_tp += I_tp_run
        I_fp += I_fp_run
        I_fn += I_fn_run
        
        print(I_tp_run,I_fp_run,I_fn_run)
    
    print('TP - FP - FN: ',I_tp,I_fp,I_fn)
    return I_tp, I_fp, I_fn, t_all 

#%%
orig10 = run_our_algo(p=10, first=11, last=40, mode=1)
orig20 = run_our_algo(p=20, first=11, last=40, mode=1)
orig30 = run_our_algo(p=30, first=11, last=40, mode=1)
orig40 = run_our_algo(p=40, first=11, last=40, mode=1)


#%%
my10 = run_our_algo(p=10, first=11, last=40, mode=2)
my20 = run_our_algo(p=20, first=11, last=40, mode=2)
my30 = run_our_algo(p=30, first=11, last=40, mode=2)
my40 = run_our_algo(p=40, first=11, last=40, mode=2)
my50 = run_our_algo(p=50, first=11, last=40, mode=2)
