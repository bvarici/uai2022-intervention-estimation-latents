#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for testing the algorithm with the covariance values

Reviewed on February 2022.
"""

import numpy as np
import itertools as itr
import time
from helpers import delta_theta_estimation


def remove_supersets(L):
    '''
    objective: remove the elements of L that are superset of another element.

    Parameters
    ----------
    L : list of lists (or sets)

    Returns
    -------
    L_filtered : list of lists
    '''
    L_to_filter = []
    L_filtered = L[:]
    for m in L:
        for n in L:
            if set(m).issubset(set(n)) and m != n:
                if n not in L_to_filter:
                    L_to_filter.append(n)
                
    for l in L_to_filter:
        L_filtered.remove(l)
            
    return L_filtered

def IMAG_naive_algorithm(Cov1,Cov2,tol=1e-9):
    '''
    Objective: Algorithm for population stats.
    Used to confirm the results before going to finite-sample.
    first, estimate PDE over all observed variables, to obtain S_Delta
    indexing of K is with respect to observed variables. not augmented, not full DAG.

    Parameters
    ----------
    Cov1, Cov2 : numpy matrices. 
        Covariance matrices for observed variables.
    tol : float, optional
        a small number for thresholding. The default is 1e-9.

    Returns
    -------
    K : list
        set of nodes connected to F, i.e. effective intervention targets
    J : list
        set of nodes not connected to F.
    S_Delta : list
        set of affected nodes.
    N : list of lists
        minimal neutralizing sets (can be more than one) for J nodes.
    N_filtered : list of lists
        remove the supersets from N.
    t_past : float
        runtime.

    '''
    t0 = time.time()
    p = Cov1.shape[0]
    # observed variables
    O = list(np.arange(p))

    Delta_hat = delta_theta_estimation(Cov1,Cov2,O)
    Delta_hat = np.abs(Delta_hat) > tol
    # set of affected nodes
    S_Delta = list(set(np.where(Delta_hat)[0])) 

    # set of nodes not connected to F
    J = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N = [[] for i in range(p)]
    # set of nodes connected to F
    K = []

    # check every subset of S_Delta. record minimal sets which makes non-intervened j's invariant.
    for size in range(1,len(S_Delta)+1):
        S_size_sets = list(itr.combinations(S_Delta,size))
        S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
        for S in S_size_sets:
            # if each element of set S is already identified, continue.
            if set(S).issubset(J):
                continue

            Delta_Theta_S = delta_theta_estimation(Cov1,Cov2,S)
            identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>tol))[0])
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
    
    t_past = time.time()-t0
    print(len(S_Delta))
    return K, J, S_Delta, N, N_filtered, t_past

def IMAG_naive_algorithm_multiple(Cov_all,tol=1e-9):
    '''
    Objective: Algorithm for population stats.
    Used to confirm the results before going to finite-sample.
    first, estimate PDE over all observed variables, to obtain S_Delta
    indexing of K is with respect to observed variables. not augmented, not full DAG.
    
    Parameters
    ----------
    Cov_all : list of numpy matrices
        Covariance matrices for observed variables.
    tol : float, optional
        a small number for thresholding. The default is 1e-9.

    Returns
    -------
    K_all : list of lists
        set of nodes connected to F, i.e. effective intervention targets
    J_all : list of lists
        set of nodes not connected to F.
    S_Delta_sizes : list
        size of the affected nodes in all contexts.
    N_filtered_all : list of lists
        minimal neutralizing sets (can be more than one) for J nodes, removed supersets.
    t_past : float
        runtime.

    '''

    t0 = time.time()
    p = Cov_all[0].shape[0]
    # observed variables
    O = list(np.arange(p))
    
    # need to run the algorithm for each pair of settings basically.
    # for each pair, we have an F node, and we need to find F-k
    setting_pairs = list(itr.combinations(np.arange(len(Cov_all)),2))
    n_F = len(setting_pairs)

    # set of nodes not connected to F
    J_all = []
    # minimal neutralizing sets (can be more than one) for j nodes
    N_filtered_all = []
    # set of nodes connected to F
    K_all = []      
    S_Delta_sizes = np.zeros(n_F)


    for F in range(n_F):
        idx_1 = setting_pairs[F][0]
        idx_2 = setting_pairs[F][1]

        Delta_hat = delta_theta_estimation(Cov_all[idx_1],Cov_all[idx_2],O)
        Delta_hat = np.abs(Delta_hat) > tol
        # set of changed nodes, note that indices are among O.
        S_Delta = list(set(np.where(Delta_hat)[0])) 
        S_Delta_sizes[F] = len(S_Delta)
        print(len(S_Delta))

        # set of nodes not connected to F
        J = []
        # minimal neutralizing sets (can be more than one) for j nodes
        N = [[] for i in range(p)]
        # set of nodes connected to F
        K = []        

        # check every subset of S_Delta. record minimal sets which makes non-intervened j's invariant.
        for size in range(1,len(S_Delta)+1):
            S_size_sets = list(itr.combinations(S_Delta,size))
            S_size_sets = [list(S_size_sets[i]) for i in range(len(S_size_sets))]
            for S in S_size_sets:
                # if each element of set S is already identified, continue.
                if set(S).issubset(J):
                    continue

                Delta_Theta_S = delta_theta_estimation(Cov_all[idx_1],Cov_all[idx_2],S)
                identified_j = (np.where(~np.diag(np.abs(Delta_Theta_S)>tol))[0])
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

        J_all.append(sorted(J))
        N_filtered_all.append(N_filtered)
        K_all.append(sorted(K))

    
    t_past = time.time()-t0
    return K_all, J_all, S_Delta_sizes, N_filtered_all, t_past

