#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wrapper functions for running simulations

run for different sizes of graphs, number of latent nodes, interventions etc.

p = number of nodes
l = number of latent nodes
i = number of interventions
c = density of the random graph
n = number of samples

Reviewed on February 2022.
"""

import numpy as np
from helpers import create_multiple_intervention, create_random_SEM, counter
from functions_graph import reduce_2_observed_multiple 
from functions_sample import sample_multiple, IMAG_naive_algorithm_sample_multiple
from config import SIMULATIONS_ESTIMATED_FOLDER
import pickle


def run_once(p,num_latent,num_int,density,num_int_settings,shift,plus_variance,num_samples_list,\
    max_subset_size=None,lambda_1=0.2,lambda_2=0.1,rho=1.0,th1=5e-3,n_max_iter=500,stop_cond=1e-5,tol=1e-9,\
        verbose=False,only_diag=True,return_pasp=False):
    '''

    Parameters
    ----------
    p : int
        size of the graph.
    num_latent : int
        num. of latent nodes.
    num_int : int
        num. of intervention targets.
    density : float
        probability of a random edge will be density/p.
    num_int_settings : int
        num. of interventional settings.
    shift : float
        amount of shift for mean-shift interventions.
    plus_variance : float
        amoung of change in the variance of noise terms.
    num_samples_list : list
        number of samples to run for.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.2.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need to change in most cases. The default is 1.0.
    th1 : float, optional
        small threshold to apply on Delta_hat. The default is 5e-3.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-5.
    tol : float, optional
        small number for thresholding. The default is 1e-9.
    verbose : Boolean
        The default is False.
    only_diag : Boolean, optional
        consider only the diagonal of Delta_Theta. The default is True.
    return_pasp: Boolean, optional
        return parents/spouses of estimated intervention targets. The default is False.


    Returns
    -------
    I_tp : list
        number of true positives for intervention target estimation.
    I_fp : list
        number of false positives for intervention target estimation.
    I_fn : TYPE
        number of false negatives for intervention target estimation.
    S_Delta_sizes_all : list
        sizes of the set of affected nodes.
    t : float
        runtime.

    '''

    # generate the graph
    B = create_random_SEM(p=p,density=density,weighted=True)
    B_all, Theta_all, Cov_all, mu_all, variance_all, BObs, ThetaObs, CovObs, L, O, I_all = \
        create_multiple_intervention(B,L_size=num_latent,I_size=num_int,n_interventions=num_int_settings,\
            shift=shift,plus_variance=plus_variance)

    # ground truth IMAG.
    IMAG, ThetaObs_all, CovObs_all = reduce_2_observed_multiple(B_all, Theta_all, Cov_all, L, I_all)
    n_F = int(len(B_all)*(len(B_all)-1)/2)    
    K = [list(np.where(IMAG[i])[0]-n_F) for i in range(n_F)]
        
    # generate data
    X_all, S_all = sample_multiple(B_all,mu_all,variance_all,max(num_samples_list))
    X_Obs = [X_all[i][:,O] for i in range(len(X_all))]


    I_tp = np.zeros(len(num_samples_list))
    I_fp = np.zeros(len(num_samples_list))
    I_fn = np.zeros(len(num_samples_list))
    
    S_Delta_sizes_all = np.zeros((len(num_samples_list),n_F))
    t = np.zeros(len(num_samples_list))

    for idx_num_samples in range(len(num_samples_list)):
        num_samples = num_samples_list[idx_num_samples]
        # need to generate S_Obs for the specified number of samples
        S_Obs = [X_Obs[i].T@X_Obs[i]/num_samples for i in range(len(X_Obs))] 


        K_hat, J_hat, S_Delta_sizes, N_filtered_hat, t_past =  \
            IMAG_naive_algorithm_sample_multiple(S_Obs,max_subset_size = max_subset_size,lambda_1=lambda_1,lambda_2=lambda_2,\
                     rho=rho,th1=th1,n_max_iter=n_max_iter,stop_cond=stop_cond,\
                    tol=tol,verbose=verbose,only_diag=only_diag)

        S_Delta_sizes_all[idx_num_samples] = S_Delta_sizes
        I_tp[idx_num_samples], I_fp[idx_num_samples], I_fn[idx_num_samples] = counter(K,K_hat)
        t[idx_num_samples] = t_past

    return I_tp, I_fp, I_fn, S_Delta_sizes_all, t

def run_repeat(p,num_latent,num_int,density,num_int_settings,shift,plus_variance,num_samples_list,\
    max_subset_size = None,lambda_1=0.2,lambda_2=0.1,rho=1.0,th1=5e-3,n_max_iter=500,stop_cond=1e-5,tol=1e-9,verbose=False,\
        only_diag=True,num_repeat=1):
    '''

    Parameters
    ----------
    p : int
        size of the graph.
    num_latent : int
        num. of latent nodes.
    num_int : int
        num. of intervention targets.
    density : float
        probability of a random edge will be density/p.
    num_int_settings : int
        num. of interventional settings.
    shift : float
        amount of shift for mean-shift interventions.
    plus_variance : float
        amoung of change in the variance of noise terms.
    num_samples_list : list
        number of samples to run for.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.2.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need to change in most cases. The default is 1.0.
    th1 : float, optional
        small threshold to apply on Delta_hat. The default is 5e-3.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-5.
    tol : float, optional
        small number for thresholding. The default is 1e-9.
    verbose : Boolean
        The default is False.
    only_diag : Boolean, optional
        consider only the diagonal of Delta_Theta. The default is True.
    return_pasp: Boolean, optional
        return parents/spouses of estimated intervention targets. The default is False.
    num_repeat: int, optional
        num. of trials to repeat the experiment

    Returns
    -------
    I_tp : 2d array: (repeat_idx, num_samples_idx)
        number of true positives for intervention target estimation.
    I_fp : 2d array: (repeat_idx, num_samples_idx)
        number of false positives for intervention target estimation.
    I_fn : 2d array (repeat_idx, num_samples_idx)
        number of false negatives for intervention target estimation.
    S_Delta_sizes_all : list
        sizes of the set of affected nodes.
    t : float
        runtime.

    '''

    I_tp = np.zeros((num_repeat,len(num_samples_list)))
    I_fp = np.zeros((num_repeat,len(num_samples_list)))
    I_fn = np.zeros((num_repeat,len(num_samples_list)))

    n_F = int(num_int_settings*(num_int_settings+1)/2)    
    S_Delta_sizes = np.zeros((num_repeat,len(num_samples_list),n_F))
    t = np.zeros((num_repeat,len(num_samples_list)))

    for i in range(num_repeat):
        res = run_once(p,num_latent,num_int,density,num_int_settings,shift,plus_variance,num_samples_list,\
                        max_subset_size,lambda_1,lambda_2,rho,th1,n_max_iter,stop_cond,tol,verbose,only_diag)

        I_tp[i] = res[0]
        I_fp[i] = res[1]
        I_fn[i] = res[2]
        S_Delta_sizes[i] = res[3]
        t[i] = res[4]

    return I_tp, I_fp, I_fn, S_Delta_sizes, t

# p = p_list[0]
# num_latent = num_latent_list[0]
# num_int = num_int_list[0]
# num_int_settings = 1
# density = density_list[0]

def run_sweep(p_list,num_latent_list,num_int_list,density_list,num_int_settings,num_samples_list,num_repeat,\
                shift=0.0,plus_variance=0.0,max_subset_size=None,lambda_1=0.1,lambda_2=0.1,\
                rho=1.0,th1=5e-3,n_max_iter=500,stop_cond=1e-6,tol=1e-9,verbose=False,only_diag=True):
    '''

    Parameters
    ----------
    p_list : list
        size of the graph. a list of p values.
    num_latent_list : list
        num. of latent nodes. a list of int values
    num_int_list : list
        num. of intervention targets. a list of int values
    density_list : list
        probability of a random edge will be density/p. a list of float values.
    num_int_settings : int
        num. of interventional settings.
    num_samples_list : list
        number of samples to run for.
    num_repeat: int, optional
        num. of trials to repeat the experiment
    shift : float
        amount of shift for mean-shift interventions.
    plus_variance : float
        amoung of change in the variance of noise terms.
    num_samples_list : list
        number of samples to run for.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.2.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need to change in most cases. The default is 1.0.
    th1 : float, optional
        small threshold to apply on Delta_hat. The default is 5e-3.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-5.
    tol : float, optional
        small number for thresholding. The default is 1e-9.
    verbose : Boolean
        The default is False.
    only_diag : Boolean, optional
        consider only the diagonal of Delta_Theta. The default is True.
    return_pasp: Boolean, optional
        return parents/spouses of estimated intervention targets. The default is False.

    Returns
    -------
    res : dict
        dictionary to store all results.
    I_tp_all : 6d array: (repeat_idx, p_idx, l_idx, i_idx, c_idx, num_samples_idx)
        number of true positives for intervention target estimation.
    I_fp_all : 6d array: (repeat_idx, p_idx, l_idx, i_idx, c_idx, num_samples_idx)
        number of false positives for intervention target estimation.
    I_fn_all : 6d array: (repeat_idx, p_idx, l_idx, i_idx, c_idx, num_samples_idx)
        number of false negatives for intervention target estimation.
    S_Delta_sizes_all : list
        sizes of the set of affected nodes.
    t_all : list
        runtimes.

    '''


    n_F = int(num_int_settings*(num_int_settings+1)/2)    

    I_tp_all = np.zeros((num_repeat,len(p_list),len(num_latent_list),len(num_int_list),len(density_list),len(num_samples_list)))
    I_fp_all = np.zeros((num_repeat,len(p_list),len(num_latent_list),len(num_int_list),len(density_list),len(num_samples_list)))
    I_fn_all = np.zeros((num_repeat,len(p_list),len(num_latent_list),len(num_int_list),len(density_list),len(num_samples_list)))
    t_all = np.zeros((num_repeat,len(p_list),len(num_latent_list),len(num_int_list),len(density_list),len(num_samples_list)))
    S_Delta_sizes_all = np.zeros((num_repeat,len(p_list),len(num_latent_list),len(num_int_list),len(density_list),len(num_samples_list),n_F))


    for r in range(num_repeat):
        for p in range(len(p_list)):
            for l in range(len(num_latent_list)):
                for i in range(len(num_int_list)):
                    for c in range(len(density_list)):


                        I_tp, I_fp, I_fn, S_Delta_sizes, t = \
                             run_repeat(p_list[p],num_latent_list[l],num_int_list[i],density_list[c],num_int_settings,\
                                shift,plus_variance,num_samples_list,max_subset_size,\
                                lambda_1,lambda_2,rho,th1,n_max_iter,stop_cond,tol,verbose,only_diag)

                        I_tp_all[r,p,l,i,c] = I_tp
                        I_fp_all[r,p,l,i,c] = I_fp
                        I_fn_all[r,p,l,i,c] = I_fn
                        t_all[r,p,l,i,c] = t
                        S_Delta_sizes_all[r,p,l,i,c] = S_Delta_sizes
                        print('r:%d, p:%d, l:%d, i:%d, c:%.1f'%(r,p_list[p],num_latent_list[l],num_int_list[i],density_list[c]))
                        
    
    res = {'num_repeat':num_repeat,'p_list':p_list,'num_latent_list':num_latent_list,'num_int_list':num_int_list,\
           'density_list':density_list,'num_int_settings':num_int_settings,'num_samples_list':num_samples_list,\
           'shift':shift,'plus_variance':plus_variance,'I_tp':I_tp_all,'I_fp':I_fp_all,'I_fn':I_fn_all,'S_Delta_sizes':S_Delta_sizes_all,\
           't':t,'lambda_1':lambda_1,'lambda_2':lambda_2,'th1':th1,'max_subset_size':max_subset_size,'only_diag':only_diag}
    

    return res, I_tp_all, I_fp_all, I_fn_all, S_Delta_sizes_all, t_all
