#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test the algorithm

run for different sizes of graphs, number of latent nodes, interventions etc.

p = number of nodes
l = number of latent nodes
i = number of interventions
c = density of the random graph
n = number of samples

Reviewed on February 2022.
"""

import numpy as np
from functions_run_simulations import run_once, run_repeat, run_sweep
from config import SIMULATIONS_ESTIMATED_FOLDER
import pickle

#%%
'''
shift interventions

shift_W_1: lambda_1 = 0.1, lambda_2 = 0.1
shift_W_2: lambda_1 = 0.1, lambda_2 = 0.2
shift_W_3: lambda_1 = 0.1, lambda_2 = 0.05

'''

p_list = [30,40,60,80,100]
num_latent_list = [5]
num_int_list = [5]
density_list = [2]
num_int_settings = 1

num_samples_list = [3000,5000,10000]
num_repeat = 100

# hyper-parameters
rho = 1.0
lambda_1 = 0.1 # for S_Delta estimation
lambda_2 = 0.05 # other PDEs
th1 = 5e-3 # throw-away very small ones
n_max_iter = 500
stop_cond = 1e-6
verbose = False
tol = 1e-9
max_subset_size = None
only_diag = True


shift = 1.0
plus_variance = 0.0

res_shift, I_tp_shift, I_fp_shift, I_fn_shift, S_Delta_sizes_shift, t_shift  = run_sweep(p_list, num_latent_list, num_int_list, density_list, \
                                               num_int_settings,num_samples_list,num_repeat,shift, plus_variance,\
                                                   max_subset_size,lambda_1,lambda_2,rho,th1,n_max_iter,stop_cond,\
                                                       tol,verbose,only_diag)

#%
f = open(SIMULATIONS_ESTIMATED_FOLDER+'/shift_W_3.pkl','wb')
pickle.dump(res_shift,f)
f.close()
   

    
#%%
'''
increased variance interventions

variance_W_1 : lambda_1 = 0.2, lambda_2 = 0.1

'''
p_list = [30,40,60,80,100]
num_latent_list = [5]
num_int_list = [5]
density_list = [2]
num_int_settings = 1

num_samples_list = [3000,5000,10000]
num_repeat = 100

# hyper-parameters
rho = 1.0
lambda_1 = 0.2 # for S_Delta estimation
lambda_2 = 0.1 # other PDEs
th1 = 5e-3 # throw-away very small ones
n_max_iter = 500
stop_cond = 1e-6
verbose = False
tol = 1e-9
max_subset_size = None
only_diag = True

shift = 0.0
plus_variance = 1.0

res_var, I_tp_var, I_fp_var, I_fn_var, S_Delta_sizes_var, t_var = run_sweep(p_list, num_latent_list, num_int_list, density_list, \
                                               num_int_settings,num_samples_list,num_repeat,shift, plus_variance,\
                                                   max_subset_size,lambda_1,lambda_2,rho,th1,n_max_iter,stop_cond,\
                                                       tol,verbose,only_diag)

    

f = open(SIMULATIONS_ESTIMATED_FOLDER+'/variance_W_1.pkl','wb')
pickle.dump(res_var,f)
f.close()


