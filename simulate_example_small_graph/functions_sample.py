#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Algorithm to run with sample covariance matrices

Reviewed on February 2022.
"""

import numpy as np
import numpy.linalg as LA
import itertools as itr
import time
from helpers import delta_theta_estimation, marginal_theta_from_cov
from functions_concept import remove_supersets

flatten_list = lambda t: list(set([item for sublist in t for item in sublist]))

def sample(B,means,variances,n_samples):
    '''
    Parameters
    ----------
    B : 
        autoregression weight matrix. assumed to be strictly upper triangular
    means : 
        internal noise means.
    variances : 
        internal noise variances.
    n_samples : integer
        number of samples to generate

    Returns
    -------
    samples : n_samples x p matrix
        DAG samples

    '''
    # assume that nodes are given in topological order
    p = len(means)
    samples = np.zeros((n_samples,p))
    noise = np.zeros((n_samples,p))
    for ix, (mean,var) in enumerate(zip(means,variances)):
        noise[:,ix] = np.random.normal(loc=mean,scale=var ** .5, size=n_samples)
        
    for node in range(p):
        parents_node = np.where(B[:,node])[0]
        if len(parents_node)!=0:
            parents_vals = samples[:,parents_node]
            samples[:,node] = np.sum(parents_vals * B[parents_node,node],axis=1) + noise[:,node]
        else:
            samples[:,node] = noise[:,node]
            
    cov = (samples.T@samples)/n_samples

    return samples, cov

def sample_multiple(B_all,means_all,variances_all,n_samples):
    '''
    all refers to multiple settings

    Parameters
    ----------
    B_all : 
        autoregression weight matrices. assumed to be strictly upper triangular
    means_all : 
        internal noise means.
    variances_all : 
        internal noise variances.
    n_samples : integer
        number of samples to generate

    Returns
    -------
    samples : n_samples x p matrix
        DAG samples

    '''
    # assume that nodes are given in topological order
    samples_all = []
    cov_all = []
    n_settings = len(B_all)

    for i in range(n_settings):
        samples, cov = sample(B_all[i],means_all[i],variances_all[i],n_samples)
        samples_all.append(samples)
        cov_all.append(cov)
            
    return samples_all, cov_all


def compute_objective(S1,S2,Delta,lambda_l1,sym_loss=False):
    '''
    Parameters
    ----------
    S1, S2 : 2d array
        Sample covariance matrices.
    Delta : 2d array
        Parameters to compute gradient wrt.
    lambda_l1: scalar
        penalty parameter for l1 regularization
    sym_loss : Boolean, optional
        Use symmetric loss or not. The default is False.

    Returns
    -------
    scalar: loss with l1 regularization
    '''
    if sym_loss == False:
        return (Delta.T@S1@Delta@S2).trace()/2 - (Delta@(S1-S2)).trace() \
            + lambda_l1*np.sum(np.abs(Delta))
    elif sym_loss == True:
        return (Delta.T@S1@Delta@S2).trace()/4+ (Delta.T@S2@Delta@S1).trace()/4 \
            - (Delta@(S1-S2)).trace() + lambda_l1*np.sum(np.abs(Delta))        
    else:
        print('sym_loss input (False by default) should be either False or True')
        return None
    
def soft_thresholding(x,alpha):
    '''
    returns soft(x,alpha)
    '''
    return np.maximum((np.abs(x)-alpha),0) * np.sign(x)


def Delta_Theta_func(S1,S2,lambda_l1=0.1,rho=1.0,n_max_iter=500,stop_cond=1e-6,verbose=False,return_sym=True):
    '''
    Difference of inverse covariance estimation.
    A Direct Approach for Sparse Quadratic Discriminant Analysis (Jiang et al. 2018)

    Parameters
    ----------
    S1, S2 : 2d array
        Sample covariance matrices.
    lambda_l1 : float
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need change in most cases. The default is 1.0.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-6.
    verbose : Boolean
        The default is False.
    return_sym : Boolean
        Take symmetric (Delta + Delta.T)/2 in the end. The default is True.

    Returns
    -------
    Phi : 2d array
        Main output. Estimated Delta_Theta difference of inverse covariances.
    obj_hist: array
        history of objective over the iterations.

    '''
    p = len(S1)
    # initialize Delta, Phi, Lambda. Fix rho
    Delta = np.zeros([p,p])
    Phi = np.zeros([p,p])
    Lambda = np.zeros([p,p])
    
    # find the minimum and maximum eigenvalues of S1 kronecker S2, for rho heuristics
    eigen_max = LA.eigvals(S1)[0]*LA.eigvals(S2)[0]
    eigen_min = LA.eigvals(S1)[-1]*LA.eigvals(S2)[-1]
    # assign rho based on these values and penalty parameter
    if rho is None:
        if lambda_l1 <= eigen_min:
            rho = eigen_min
        elif lambda_l1 <= eigen_max:
            rho = eigen_max
        else:
            rho = lambda_l1
        
    # compute SVD's for sample covariance matrices
    [U1,D1,_] = LA.svd(S1)
    [U2,D2,_] = LA.svd(S2)
    B = 1/(D1[:,np.newaxis]*D2[np.newaxis,] + rho) 
    
    obj_hist = np.zeros(n_max_iter)
    obj_hist[0] = compute_objective(S1,S2,Delta,lambda_l1)
    # now update
    for it in range(n_max_iter):
        A = (S1-S2) - Lambda + rho*Phi
        # update Delta, notice the Hadamard product
        Delta = U1@(B*(U1.T@A@U2))@U2.T
        # update Phi
        Phi = soft_thresholding(Delta+Lambda/rho,lambda_l1/rho)
        # update Lambda
        Lambda += rho*(Delta-Phi)
    
        obj = compute_objective(S1,S2,Delta,lambda_l1)
        obj_hist[it] = obj
    
        # check stopping condition
        if np.abs(obj-obj_hist[it-1]) < stop_cond*(np.abs(obj)+1):
            sparsity = np.mean(Phi!=0)
            if verbose == True:
                print('Sparsity is %.3f, Converged in %d iterations, lambda:%.3f, rho:%.3f'%(sparsity,it,lambda_l1,rho))
            if return_sym == True:
                Phi = (Phi+Phi.T)/2
            return Phi, obj_hist[:it]
        
    if return_sym == True:
        Phi = (Phi+Phi.T)/2
        
    sparsity = np.mean(Phi!=0)
    if verbose == True:
        print('Sparsity is %.3f, Converged in %d iterations, lambda:%.3f, rho:%.3f'%(sparsity,it,lambda_l1,rho))
    return Phi, obj_hist        

    
def IMAG_naive_algorithm_sample(S1,S2,max_subset_size=None,lambda_1=0.1,lambda_2=0.1,rho=1.0,\
        th1=1e-3,n_max_iter=500,stop_cond=1e-6,tol=1e-9, \
        verbose=True,only_diag=True):
    '''
    Parameters
    ----------
    S1, S2 : matrices
        sample covariance matrices.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.1.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
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
    K : list
        set of nodes connected to F, i.e. effective intervention targets
    J : list
        set of nodes not connected to F.
    S_Delta_size : int
        size of the set of affected nodes.
    N_filtered : list of lists
        remove the supersets from N, minimal neutralizing sets (can be more than one) for J nodes.
    t_past : float
        runtime.

    '''

    t0 = time.time()
    p = len(S1)
    # set of nodes not connected to F
    J = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N = [[] for i in range(p)]
    # set of nodes connected to F
    K = []
    # estimate Delta_Theta for all nodes at first
    Delta_hat, obj_hist = Delta_Theta_func(S1,S2,lambda_1,rho,n_max_iter,stop_cond,verbose)
    # apply some small threshold to Delta_hat
    Delta_hat = np.abs(Delta_hat) > th1
    # set of all changed nodes
    if only_diag is True:
        S_Delta = np.where(np.diag(Delta_hat))[0] 
    elif only_diag is False:
        S_Delta = list(set(np.where(Delta_hat)[0]))
    S_Delta_size = len(S_Delta)

    # if no limit for S size is specified, check all subsets.
    if max_subset_size is None:
        max_subset_size = S_Delta_size
            
    for size in range(1,max_subset_size+1):
        S_size_sets = list(itr.combinations(S_Delta,size))
        S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
        for S in S_size_sets:
            # if each element of set S is already identified, continue.
            if set(S).issubset(J):
                continue

            Delta_Theta_S = Delta_Theta_func(S1[S][:,S],S2[S][:,S],lambda_2,rho,n_max_iter,stop_cond,verbose)[0]
            identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>th1))[0])
            # if there is some new j, add it to J
            for id_j in identified_j:
                N[S[id_j]].append(S)
                if S[id_j] not in J:
                    J.append(S[id_j])


    # remove the supersets from neutralizing sets
    N_filtered = N[:].copy()
    for j in J:
        N_filtered[j] = remove_supersets(N[j])

    K = [i for i in S_Delta if i not in J]
    t_past = time.time() - t0

    return K, J, S_Delta_size, N_filtered, t_past

def IMAG_pasp_sample(S1,S2,max_subset_size=None,lambda_1=0.1,lambda_pasp=0.1,rho=1.0,\
        th1=1e-3,n_max_iter=500,stop_cond=1e-6,tol=1e-9, \
        verbose=True,only_diag=True,return_pasp=True):
    '''
    Parameters
    ----------
    S1, S2 : matrices
        sample covariance matrices.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.1.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
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
    return_pasp: Boolean, optional
        return parents/spouses of estimated intervention targets. The default is True.
    
    Returns
    -------        
    K : list
        set of nodes connected to F, i.e. effective intervention targets
    J : list
        set of nodes not connected to F.
    S_Delta_size : int
        size of the set of affected nodes.
    N_filtered : list of lists
        remove the supersets from N, minimal neutralizing sets (can be more than one) for J nodes.
    t_past : float
        runtime.

    '''
    
    t0 = time.time()
    p = len(S1)
    # set of nodes not connected to F
    J = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N = [[] for i in range(p)]
    # set of nodes connected to F
    K = []
    # estimate Delta_Theta for all nodes at first
    Delta_hat, obj_hist = Delta_Theta_func(S1,S2,lambda_1,rho,n_max_iter,stop_cond,verbose)
    # apply some small threshold to Delta_hat
    Delta_hat = np.abs(Delta_hat) > th1
    # set of all changed nodes
    if only_diag is True:
        S_Delta = np.where(np.diag(Delta_hat))[0] 
    elif only_diag is False:
            S_Delta = list(set(np.where(Delta_hat)[0]))

    S_Delta_size = len(S_Delta)

    # if no limit for S size is specified, check all subsets.
    if max_subset_size is None:
        max_subset_size = S_Delta_size

            
    for size in range(1,max_subset_size+1):
        S_size_sets = list(itr.combinations(S_Delta,size))
        S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
        for S in S_size_sets:
            # if each element of set S is already identified, continue.
            if set(S).issubset(J):
                continue

            Delta_Theta_S = Delta_Theta_func(S1[S][:,S],S2[S][:,S],lambda_1,rho,n_max_iter,stop_cond,verbose)[0]
            identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>tol))[0])
            # if there is some new j, add it to J
            for id_j in identified_j:
                N[S[id_j]].append(S)
                if S[id_j] not in J:
                    J.append(S[id_j])


    K = [i for i in S_Delta if i not in J]
    K_pasp = [[] for i in range(len(K))]
    # remove the supersets from neutralizing sets
    N_filtered = N[:].copy()
    for j in J:
        N_filtered[j] = remove_supersets(N[j])
        
    # identify parent-spouses of i nodes.
    if return_pasp is True:
        K_pasp, K_non_pasp = post_pasp(K,S_Delta,S1,S2,lambda_pasp,rho,n_max_iter,stop_cond,tol,verbose)
    else:
        K_non_pasp = [[] for i in range(len(K))]


    t_past = time.time() - t0

    return K, K_pasp, K_non_pasp, J, S_Delta, N_filtered, t_past

def IMAG_pasp_sample_multiple(S_all,max_subset_size = None,lambda_1=0.1,lambda_2=0.1,lambda_pasp=0.1,rho=1.0,\
        th1=1e-3,n_max_iter=500,stop_cond=1e-6,tol=1e-9,\
        verbose=True,only_diag=True,return_pasp=True):


    t0 = time.time()
    p = len(S_all[0])
    # set of nodes not connected to F
    J_all = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N_filtered_all = []
    #N_filtered_all = [[] for i in range(p)]
    # set of nodes connected to F
    K_all = []
    K_pasp_all = []
    K_non_pasp_all = []

    # need to run the algorithm for each pair of settings basically.
    # for each pair, we have an F node, and we need to find F-k
    setting_pairs = list(itr.combinations(np.arange(len(S_all)),2))
    n_F = len(setting_pairs)
    S_Delta_sizes = np.zeros(n_F)

    for F in range(n_F):
        idx_1 = setting_pairs[F][0]
        idx_2 = setting_pairs[F][1]

        # estimate Delta_Theta for all nodes at first
        Delta_hat, obj_hist = Delta_Theta_func(S_all[idx_1],S_all[idx_2],lambda_1,rho,n_max_iter,stop_cond,verbose)
        # apply some small threshold to Delta_hat
        Delta_hat = np.abs(Delta_hat) > th1
        # set of all changed nodes
        if only_diag is True:
            S_Delta = np.where(np.diag(Delta_hat))[0] 
        elif only_diag is False:
            S_Delta = list(set(np.where(Delta_hat)[0]))
            
        S_Delta_sizes[F] = len(S_Delta)

        print(len(S_Delta))

        # if no limit for S size is specified, check all subsets.
        if max_subset_size is None:
            max_subset_size = len(S_Delta)

        J = []
        N = [[] for i in range(p)]
        K = []


        # finding J0 is not even needed for naive.
        # check every subset of S_Delta. record minimal sets which makes non-intervened j's invariant.
        for size in range(1,max_subset_size+1):
            S_size_sets = list(itr.combinations(S_Delta,size))
            S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
            for S in S_size_sets:
                # if each element of set S is already identified, continue.
                if set(S).issubset(J):
                    continue

                Delta_Theta_S = Delta_Theta_func(S_all[idx_1][S][:,S],S_all[idx_2][S][:,S],lambda_2,rho,n_max_iter,stop_cond,verbose)[0]
                identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>th1))[0])
                # if there is some new j, add it to J
                for id_j in identified_j:
                    N[S[id_j]].append(S)
                    if S[id_j] not in J:
                        J.append(S[id_j])

        K = [i for i in S_Delta if i not in J]

        # remove the supersets from neutralizing sets
        N_filtered = N[:].copy()
        for j in J:
            N_filtered[j] = remove_supersets(N[j])


        # identify parent-spouses of i nodes.
        if return_pasp is True:
            K_pasp, K_non_pasp = post_pasp(K,S_Delta,S_all[idx_1],S_all[idx_2],lambda_pasp,rho,n_max_iter,stop_cond,tol,verbose)
            K_pasp_all.append(K_pasp)
            K_non_pasp_all.append(K_non_pasp)
        else:
            K_non_pasp = [[] for i in range(len(K))]


        J_all.append(sorted(J))
        N_filtered_all.append(N_filtered)
        K_all.append(sorted(K))


    t_past = time.time() - t0

    return K_all, K_pasp_all, K_non_pasp_all, J_all, S_Delta_sizes, N_filtered_all, t_past

def IMAG_naive_algorithm_sample_multiple(S_all,max_subset_size = None,lambda_1=0.1,lambda_2=0.1,rho=1.0,\
        th1=1e-3,n_max_iter=500,stop_cond=1e-6,tol=1e-9,\
        verbose=True,only_diag=True):
    '''
    

    Parameters
    ----------
    S_all : list of matrices
        list of sample covariance matrices.
    max_subset_size : int, optional
        max. size of the subsets to compute PDE. The default is None.
    lambda_1 : float, optional
        l1 norm parameter for S_Delta estimation. The default is 0.1.
    lambda_2 : float, optional
        l1 norm parameter for Delta_Theta estimations over multiple nodes. The default is 0.1.
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
    K_all : list of lists
        set of nodes connected to F, i.e. effective intervention targets
    J_all : list of lists
        set of nodes not connected to F.
    S_Delta_sizes : list
        sizes of the set of affected nodes.
    N_filtered_all : list of lists
        remove the supersets from N, minimal neutralizing sets (can be more than one) for J nodes.
    t_past : float
        runtime.

    '''


    t0 = time.time()
    p = len(S_all[0])
    # set of nodes not connected to F
    J_all = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N_filtered_all = []
    #N_filtered_all = [[] for i in range(p)]
    # set of nodes connected to F
    K_all = []


    # need to run the algorithm for each pair of settings basically.
    # for each pair, we have an F node, and we need to find F-k
    setting_pairs = list(itr.combinations(np.arange(len(S_all)),2))
    n_F = len(setting_pairs)
    S_Delta_sizes = np.zeros(n_F)

    for F in range(n_F):
        idx_1 = setting_pairs[F][0]
        idx_2 = setting_pairs[F][1]

        # estimate Delta_Theta for all nodes at first
        Delta_hat, obj_hist = Delta_Theta_func(S_all[idx_1],S_all[idx_2],lambda_1,rho,n_max_iter,stop_cond,verbose)
        # apply some small threshold to Delta_hat
        Delta_hat = np.abs(Delta_hat) > th1
        # set of all changed nodes
        if only_diag is True:
            S_Delta = np.where(np.diag(Delta_hat))[0] 
        elif only_diag is False:
            S_Delta = list(set(np.where(Delta_hat)[0]))
            
        S_Delta_sizes[F] = len(S_Delta)

        print(len(S_Delta))

        # if no limit for S size is specified, check all subsets.
        if max_subset_size is None:
            max_subset_size = len(S_Delta)

        J = []
        N = [[] for i in range(p)]
        K = []


        # finding J0 is not even needed for naive.
        # check every subset of S_Delta. record minimal sets which makes non-intervened j's invariant.
        for size in range(1,max_subset_size+1):
            S_size_sets = list(itr.combinations(S_Delta,size))
            S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
            for S in S_size_sets:
                # if each element of set S is already identified, continue.
                if set(S).issubset(J):
                    continue

                Delta_Theta_S = Delta_Theta_func(S_all[idx_1][S][:,S],S_all[idx_2][S][:,S],lambda_2,rho,n_max_iter,stop_cond,verbose)[0]
                identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>th1))[0])
                # if there is some new j, add it to J
                for id_j in identified_j:
                    N[S[id_j]].append(S)
                    if S[id_j] not in J:
                        J.append(S[id_j])

        K = [i for i in S_Delta if i not in J]

        # remove the supersets from neutralizing sets
        N_filtered = N[:].copy()
        for j in J:
            N_filtered[j] = remove_supersets(N[j])


        J_all.append(sorted(J))
        N_filtered_all.append(N_filtered)
        K_all.append(sorted(K))

    t_past = time.time() - t0

    return K_all, J_all, S_Delta_sizes, N_filtered_all, t_past


def post_pasp(K,S_Delta,S1,S2,lambda_pasp=0.2,rho=1.0,n_max_iter=500,stop_cond=1e-6,tol=1e-9,verbose=False):
    '''
    Objective: finds parents/spouses of an intervend node K

    Parameters
    ----------
    K : list
        intervention targets.
    S_Delta : list
        affected nodes.
    S1, S2 : matrices
        sample covariance matrices.
    lambda_pasp : float, optional
        l1 norm parameter for Delta_Theta estimation. The default is 0.1.
    rho : float
        penalty parameter for ADMM. No need to change in most cases. The default is 1.0.
    n_max_iter : integer
        maximum number of iterations for ADMM. Does not need to be too large. The default is 500.
    stop_cond : float
        stopping condition for ADMM iterations. The default is 1e-6.
    tol : float, optional
        small number for thresholding. The default is 1e-9.
    verbose : Boolean
        The default is False.

    Returns
    -------
    K_pasp : list of lists
        parents/spouses for the intervened nodes.
    K_non_pasp : list of lists
        exlusion sets of K_pasp.

    '''

    J = list(np.setdiff1d(S_Delta,K))

    #K_pasp = [[] for i in range(len(K))]
    K_non_pasp = [[] for i in range(len(K))]

    for size in range(1,len(S_Delta)+1):
        S_size_sets = list(itr.combinations(S_Delta,size))
        S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
        for S in S_size_sets:
            # if S does not contain any K node, continue.
            if set(S).issubset(J):
                continue

            Delta_Theta_S = Delta_Theta_func(S1[S][:,S],S2[S][:,S],lambda_pasp,rho,n_max_iter,stop_cond,verbose)[0]
            Delta_Theta_S = np.abs(Delta_Theta_S) > tol

            K_in_S = np.intersect1d(K,S)

            # now, for K_in_S nodes, check whether there are for sure non-pasp J_in_S nodes.
            for k in K_in_S:
                k_idx = K.index(k)
                non_pasp_j_for_k = np.where(~Delta_Theta_S[S.index(k)])[0]

                for j in non_pasp_j_for_k:
                    if S[j] not in K_non_pasp[k_idx]:
                        K_non_pasp[k_idx].append(S[j])

    K_pasp = [np.setdiff1d(S_Delta,K_non_pasp[k_idx]) for k_idx in range(len(K))]
    K_pasp = [np.setdiff1d(K_pasp[k_idx],K) for k_idx in range(len(K))]

    return K_pasp, K_non_pasp




