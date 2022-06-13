#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze results of FCI-JCI123 algo on larger graphs 

Reviewed on February 2022.
"""

import numpy as np
from helpers import counter
#from my_config import SIMULATIONS_ESTIMATED_FOLDER
#import pickle


#JCI_data_folder = './experiments/simul/jci-results'
#JCI_data_folder = './experiments/simul/jci-mydata'


def read_results(p,first,last,mode,n_samples=5000):
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
        edge_filename = JCI_data_folder + '/p'+str(p)+'/p'+str(p)+'-model'+str(run_no)+'-edge.csv'
        edge_true = np.genfromtxt(edge_filename,delimiter=',',skip_header=1)
        
        K_true = []
        for i in range(p,p+pContext):
            target_instant = list(np.where(edge_true[i]==1)[0])
            K_true.append(target_instant)
        
        K_true_all.append(K_true)
    
        edge_jci_filename = JCI_data_folder + '/p'+str(p)+'/p'+str(p)+'-model'+str(run_no)+'-fci-jci123-bs-edge.csv'
        edge_jci = np.genfromtxt(edge_jci_filename,delimiter=',',skip_header=1)
    
        K_jci = []
        for i in range(p,p+pContext):
            target_instant = list(np.where(edge_jci[i]==1)[0])
            K_jci.append(target_instant)
            
        runtime_filename = JCI_data_folder + '/p'+str(p)+'/p'+str(p)+'-model'+str(run_no)+'-fci-jci123-bs-runtime.csv'
        t = np.around(np.genfromtxt(runtime_filename),2)
        t_all.append(t)        
                        
        I_tp_run, I_fp_run, I_fn_run = counter(K_true,K_jci)
        I_tp += I_tp_run
        I_fp += I_fp_run
        I_fn += I_fn_run
        print(I_tp_run, I_fp_run, I_fn_run)
    
    print(I_tp,I_fp,I_fn)
    return I_tp, I_fp, I_fn, t_all

#%%
orig10 = read_results(p=10, first=11, last=40, mode=1)
orig20 = read_results(p=20, first=11, last=40, mode=1)
orig30 = read_results(p=30, first=11, last=30, mode=1)

#%%
my10 = read_results(p=10, first=11, last=40, mode=2)
my20 = read_results(p=20, first=11, last=40, mode=2)
my30 = read_results(p=30, first=11, last=40, mode=2)
my40 = read_results(p=40, first=11, last=40, mode=2)
my50 = read_results(p=50, first=11, last=16, mode=2)
