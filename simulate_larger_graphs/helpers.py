#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
helper functions.

Reviewed on February 2022.
"""

import numpy as np
import numpy.linalg as LA
from config import SIMULATIONS_ESTIMATED_FOLDER
import pickle

def marginal_theta_from_cov(Cov,M):
    '''
    Objective: Given covariance matrix, take inverse of Sigma_{M,M}
    note: faster than using marginal_theta_from_theta function

    Parameters
    ----------
    Cov : matrix
        covariance matrix.
    M : list
        a subset of nodes.

    Returns
    -------
    ThetaM : matrix
        precision matrix for the subset of nodes M.
    
    '''
    
    ThetaM = LA.inv(Cov[M][:,M])
    return ThetaM

def get_precision_cov(B,Omega):
    '''
    Objective: Given B and Omega where noise is ~ N(0,Omega), return precision matrix

    Parameters
    ----------
    B : matrix
        autoregressive matrix.
    Omega : matrix
        diagonal noise matrix.

    Returns
    -------
    Theta : matrix
        precision matrix.
    Cov : matrix
        covariance matrix.

    '''

    p = len(B)
    Theta = (np.eye(p)-B)@LA.inv(Omega)@(np.eye(p)-B).T
    Cov = LA.inv(Theta)
    #Cov = LA.inv(np.eye(p)-B)@Omega@(LA.inv(np.eye(p)-B)).T
    return Theta, Cov

def delta_theta_estimation(Cov1,Cov2,M):
    '''
    Objective: compute Delta_Theta over set M, subset of [p]

    Parameters
    ----------
    Cov1, Cov2 : matrices
        covariance matrices.
    M : list
        a subset of nodes.

    Returns
    -------
    Delta_Theta_M : matrix
        precision difference matrix over set M.

    '''
    Theta1_M = marginal_theta_from_cov(Cov1,M)
    Theta2_M = marginal_theta_from_cov(Cov2,M)
    Delta_Theta_M = np.abs(Theta1_M-Theta2_M)
    return Delta_Theta_M 


def load_res(filename):
    '''
    Objective: read the pickle files from the simulations results folder.
    '''
    with open(SIMULATIONS_ESTIMATED_FOLDER+'/'+filename+'.pkl','rb') as f:
        res = pickle.load(f)
    return res

def counter(K,K_hat):
    '''
    Objective: evaluate the estimated targets, K_hat.

    Parameters
    ----------
    K : list of lists
        true intervention targets.
    K_hat : list of lists
        estimated intervention targets.

    Returns
    -------
    tp_i_num : int
        num. of true positives for intervention targets.
    fp_i_num : int
        num. of false positives for intervention targets.
    fn_i_num : int
        num. of false negatives for intervention targets.

    '''
    tp_i_num = 0
    fp_i_num = 0
    fn_i_num = 0

    for i in range(len(K)):
        tp_i = list(np.intersect1d(K[i],K_hat[i]))
        fp_i = list(np.setdiff1d(K_hat[i],K[i]))
        fn_i = list(np.setdiff1d(K[i],K_hat[i]))                
        
        tp_i_num += len(tp_i)
        fp_i_num += len(fp_i)
        fn_i_num += len(fn_i)        

    return tp_i_num, fp_i_num, fn_i_num

def counter_pasp(K,K_hat,K_pasp,K_pasp_hat):
    '''
    Objective: evaluate the estimated targets, and their parents/spouses.


    Parameters
    ----------
    K : list of lists
        true intervention targets.
    K_hat : list of lists
        estimated intervention targets.
    K_pasp : list of lists
        parents/spouses of true intervention targets.
    K_pasp_hat : TYPE
        estiamted parents/spouses of estimated intervention targets.

    Returns
    -------
    tp_k_num : int
        num. of true positives for intervention targets.
    fp_k_num : int
        num. of false positives for intervention targets.
    fn_k_num : int
        num. of false negatives for intervention targets.
    tp_k_pasp_num : int
        num. of true positives for parents/spouses of intervention targets.
    fp_k_pasp_num : int
        num. of false positives for parents/spouses of intervention targets.
    fn_k_pasp_num : int
        num. of false negatives for parents/spouses of intervention targets.

    '''
    tp_k_num = 0
    fp_k_num = 0
    fn_k_num = 0

    tp_k_pasp_num = 0
    fp_k_pasp_num = 0
    fn_k_pasp_num = 0

    for i in range(len(K)):
        tp_k = list(np.intersect1d(K[i],K_hat[i]))
        fp_k = list(np.setdiff1d(K_hat[i],K[i]))
        fn_k = list(np.setdiff1d(K[i],K_hat[i]))                
        
        tp_k_num += len(tp_k)
        fp_k_num += len(fp_k)
        fn_k_num += len(fn_k)    
        
        # all nodes appear in either K[i] or K_hat[i], for i-th setting
        K_cup = sorted(list(set(list(K[i])+list(K_hat[i]))))
        for k_idx in range(len(K_cup)):
            if K_cup[k_idx] in fp_k:
                # it is a false positive k. all PaSp are automatically fn.
                fp_k_pasp_num += len(K_pasp_hat[i][K_hat[i].index(K_cup[k_idx])])
            elif K_cup[k_idx] in fn_k:
                # it is false negatie k. all PaSp are automatically fn
                fn_k_pasp_num += len(K_pasp[i][K[i].index(K_cup[k_idx])])
            elif K_cup[k_idx] in tp_k:
                true_pasp = K_pasp[i][K[i].index(K_cup[k_idx])]
                est_pasp = K_pasp_hat[i][K_hat[i].index(K_cup[k_idx])]
                fp_k_pasp_num += len(np.setdiff1d(est_pasp,true_pasp))
                fn_k_pasp_num += len(np.setdiff1d(true_pasp,est_pasp))
                tp_k_pasp_num += len(np.intersect1d(true_pasp,est_pasp))

        
    return tp_k_num, fp_k_num, fn_k_num, tp_k_pasp_num, fp_k_pasp_num, fn_k_pasp_num


def take_union(LL):
    '''
    for list of lists LL, returns the union of lists.
    '''
    L = []
    for i in range(len(LL)):
        L += LL[i]
    return list(set(L))


def create_random_SEM(p,density,weighted=True):
    '''
    create an Erdös-Renyi linear SEM.

    Parameters
    ----------
    p : int
        size of the graph.
    density : float
        probability of a random edge will be density/p.
    weighted : boolean, optional
        if False, returns all edges with weight 1. The default is True.

    Returns
    -------
    B : matrix
        autoregressive matrix denoting the random linear SEM.

    '''

    # create base B
    B = np.random.uniform(-1,-0.25,[p,p])* np.random.choice([-1,1],size=[p,p])
    B = np.triu(B)
    np.fill_diagonal(B,0)
    # assign random edges
    edge_indices = np.triu(np.random.uniform(size=[p,p])<(density/p))
    B = B*edge_indices    
    
    if weighted is False:
        B[np.where(B)] = 1

    return B

def create_intervention(B,L_size,I_size,mu=0,shift=1.0,plus_variance=0.5,variance=1.0,tol=1e-6,\
                        B_distortion_amplitude=0):
    '''
    Objective: enables shift intervention, and changing variance of noise interventions for one DAG

    intervention should be applied before getting the SMCM and (I)MAGs. since the data
    will be created with respect to the full DAG.

    Parameters
    ----------
    B : matrix
        autoregressive matrix denoting the random linear SEM.
    L_size : int
        num. of latent nodes to be excluded
    I_size : int
        num. of intervention targets to be applied among observed nodes.
    mu : float, optional
        mean of noise terms. The default is 0.
    shift : float, optional
        amount of shift for mean-shift interventions. The default is 1.0.
    plus_variance : float, optional
        amoung of change in the variance of noise terms. The default is 0.5.
    variance : float, optional
        variance of noise terms. The default is 1.0.
    tol : float, optional
        a small threshold. The default is 1e-6.
    B_distortion_amplitude : float, optional
        change in edge weights. The default is 0.

    Returns
    -------
    B1, B2 : matrices
        post-intervention SEM matrices.
    Theta1, Theta2 : matrices
        precision matrices.
    Cov1, Cov2 : matrices
        covariance matrices.
    L : list
        list of latent nodes.
    O : list
        list of observed nodes.
    I : list
        list of intervention targets..

    '''
    # V = O + L. observed and latent nodes
    p_V = B.shape[0]
    V = list(np.arange(p_V))
    # select latent nodes at random
    L = list(np.sort(np.random.choice(p_V,L_size,replace=False)))

    O = V.copy()
    for l in L:
        O.remove(l)
        
    # select intervention set among the observed variables.
    I = list(np.sort(np.random.choice(O,I_size,replace=False)))

    # internal noises will be N(mean,variance) for G1, N(mean+intervention_side*shift,variance) for G2
    mu1 = mu*np.ones(p_V)
    mu2 = mu*np.ones(p_V)
    variance1 = variance*np.ones(p_V)
    variance2 = variance*np.ones(p_V)
    # shift the noise means
    mu2[I] += shift
    # increase the noise variances
    variance2[I] += plus_variance
    Omega1 = mu1**2 + variance1
    Omega2 = mu2**2 + variance2

    B1 = B.copy(); B2 = B.copy()    
    # apply imperfect changes to B, or keep it the same if  B_distortion_amplitude=0.0
    B_distortion = np.sign(B2) * B_distortion_amplitude
    B_distortion[:,np.delete(np.arange(p_V),I)] = 0
    B2 -= B_distortion

    # now ready to get precision (or generalized precision) matrix
    Theta1, Cov1 = get_precision_cov(B1, np.diag(Omega1))
    Theta2, Cov2 = get_precision_cov(B2, np.diag(Omega2))

    return B1, Theta1, Cov1, B2, Theta2, Cov2, L, O, I



def create_multiple_intervention(B,L_size,I_size,n_interventions=1,mu=0,shift=1.0,plus_variance=0.5,variance=1.0,tol=1e-6,\
                        B_distortion_amplitude=0,latent_must_source=False):
    '''
    
    Parameters
    ----------
    B : matrix
        autoregressive matrix denoting the random linear SEM.
    L_size : int
        num. of latent nodes to be excluded
    I_size : int
        num. of intervention targets to be applied among observed nodes.
    n_interventions : int, optional
        num. of interventional settings. The default is 1.
    mu : float, optional
        mean of noise terms. The default is 0.
    shift : float, optional
        amount of shift for mean-shift interventions. The default is 1.0.
    plus_variance : float, optional
        amoung of change in the variance of noise terms. The default is 0.5.
    variance : float, optional
        variance of noise terms. The default is 1.0.
    tol : float, optional
        a small threshold. The default is 1e-6.
    B_distortion_amplitude : float, optional
        change in edge weights. The default is 0.
    latent_must_source : boolean, optional
        if True, latent nodes must be sources and cannot have parents. The default is False.

    Returns
    -------
    B_all : list of matrices
        post-intervention SEM matrices.
    Theta_all : list of matrices
        precision matrices.
    Cov_all : list of matrices
        covariance matrices.
    mu_all : TYPE
        DESCRIPTION.
    variance_all : TYPE
        DESCRIPTION.
    BObs : list of matrices
        post-intervention SEM matrices of observed variables.
    ThetaObs : list of matrices
        precision matrices of observed variables
    CovObs : list of matrices
        covariance matrices of observed variables
    L : list
        list of latent nodes.
    O : list
        list of observed nodes.
    I_all : list
        list of intervention targets..
    
    '''
    # V = O + L. observed and latent nodes
    p_V = B.shape[0]
    V = list(np.arange(p_V))
    
    if latent_must_source is False:    
        # select latent nodes at random
        L = list(np.sort(np.random.choice(p_V,L_size,replace=False)))
    elif latent_must_source is True:
        # indices of root nodes in full DAG.
        root_nodes = np.where(np.sum(B!=0,0)==0)[0]   
        if L_size > len(root_nodes):
            print('num_latents is reduced to num_roots')
            L_size = len(root_nodes)
        
        # select latent nodes at random among root nodes
        L = list(np.sort(np.random.choice(root_nodes,L_size,replace=False)))    

    O = [i for i in V if i not in L]
        
    'observational case will be setting_0'
    mu1 = mu*np.ones(p_V)
    variance1 = variance*np.ones(p_V)
    Omega1 = mu1**2 + variance1
    B1 = B.copy()
    Theta1, Cov1 = get_precision_cov(B1, np.diag(Omega1))

    B_all = [B1]
    Theta_all = [Theta1]
    Cov_all = [Cov1]
    I_all = [[]]

    mu_all = [mu1]
    variance_all = [variance1]

    #I_all = [[] for range(n_interventions)]

    BObs = [B1[O][:,O]]
    CovObs = [Cov1[O][:,O]]
    ThetaObs = [LA.inv(CovObs[-1])]

    for idx_intervention in range(1,n_interventions+1):
        I = list(np.sort(np.random.choice(O,I_size,replace=False)))
        I_all.append(I.copy())
        mu2 = mu*np.ones(p_V)
        variance2 = variance*np.ones(p_V)
        # shift the noise means
        mu2[I] += shift
        # increase the noise variances
        variance2[I] += plus_variance
        Omega2 = mu2**2 + variance2
        B2 = B.copy()    
        # apply imperfect changes to B, or keep it the same if  B_distortion_amplitude=0.0
        B_distortion = np.sign(B2) * B_distortion_amplitude
        B_distortion[:,np.delete(np.arange(p_V),I)] = 0    
        B2 -= B_distortion
        Theta2, Cov2 = get_precision_cov(B2, np.diag(Omega2))
        B_all.append(B2.copy())
        Theta_all.append(Theta2.copy())
        Cov_all.append(Cov2.copy())
        mu_all.append(mu2.copy())
        variance_all.append(variance2.copy())

        BObs.append(B2[O][:,O].copy())
        CovObs.append(Cov2[O][:,O].copy())
        ThetaObs.append(LA.inv(CovObs[-1]).copy())

    return B_all, Theta_all, Cov_all, mu_all, variance_all, BObs, ThetaObs, CovObs, L, O, I_all

