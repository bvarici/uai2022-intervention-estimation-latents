#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for generating augmented graphs, IMAGs etc. from DAGs.

Reviewed on February 2022.
"""

import numpy as np
import itertools as itr
import causaldag as cd 
from helpers import marginal_theta_from_cov

def DiffPair(pair):
    '''
    Parameters
    ----------
    pair : list 
        list of 2 lists.

    Returns
    -------
    sdiff: list
        difference of the 2 input lists.

    '''
    s1 = pair[0]
    s2 = pair[1]
    sdiff = list(set(s1) - set(s2)) + list(set(s2) - set(s1))
    return sdiff

def DAG2MAG(A,L):
    '''
    Objective: takes a full DAG and the latent variables. Returns MAG.    
    requires causaldag package. 
    
    considers only directed and bi-directed edges (no undirected edge)

    Parameters
    ----------
    A : matrix.
        can be either weights or directed adjacency matrix for the full DAG.
    L : list
        latent nodes associated with the full DAG.

    Returns
    -------
    M : matrix
        MAG.
    '''

    D = cd.DAG.from_amat(A)
    M = (D.marginal_mag(latent_nodes=set(L),relabel='default')).to_amat()
    M[np.where(M==3)] = 0
    M[np.where(M==2)] = 1
    return M

def DAG_I_2_Aug(A,L,I_all):
    '''
    Objective: creates a large augmented graph without removing latent nodes, adds F nodes.
    this will be used to obtain IMAG

    Parameters
    ----------
    A : matrix.
        can be either weights or directed adjacency matrix for the full DAG.
    L : list
        latent nodes associated with the full DAG.
    I_all : list of lists
        set of interventions, for instance, I_all = [[],[1,2],[3]]
        empty set denotes observational model.

    Returns
    -------
    AugA : matrix
        Augmented DAG.
    setting_pairs : list
        pairs of settings. each one correspond to an F node.
    '''
    
    p = A.shape[0]
    #V = list(np.arange(p))
    #O = [i for i in V if i not in L]
    # number of setting pairs (or F nodes)
    setting_pairs = list(itr.combinations(I_all,2))
    n_F = len(setting_pairs)

    AugA = np.zeros([p+n_F,p+n_F])
    AugA[n_F:][:,n_F:] = A.copy()

    for F in range(n_F):
        I = DiffPair(setting_pairs[F])
        # assign F -> I edges.
        for i_idx in range(len(I)):
            # I consists of indices in A. add n_F
            AugA[F,n_F+I[i_idx]] = 1

    return AugA, setting_pairs

def DAG_I_2_IMAG(A,L,I_all):
    '''
    Objective: takes a full DAG, the latent variables, and intervention targets. 
    Returns I-MAG. requires causaldag package. 

    Parameters
    ----------
    A : matrix.
        can be either weights or directed adjacency matrix for the full DAG.
    L : list
        latent nodes associated with the full DAG.
    I_all : list of lists
        set of interventions, for instance, I_all = [[],[1,2],[3]]
        empty set denotes observational model.
        
    Returns
    -------
    M : matrix
        I-MAG.

    '''
    
    AugA, setting_pairs = DAG_I_2_Aug(A,L,I_all)
    # note that latent node indices are increased by n_F
    n_F = len(setting_pairs)

    D = cd.DAG.from_amat(AugA)
    L_Aug = [l+n_F for l in L]
    M = (D.marginal_mag(latent_nodes=set(L_Aug),relabel='default')).to_amat()
    M[np.where(M==3)] = 0
    M[np.where(M==2)] = 1
    return M    

def reduce_2_observed(B,Theta1,Cov1,Theta2,Cov2,L,I):
    '''
    Objective: returns the observed Covariance and precision matrices

    Parameters
    ----------
    B : matrix
        adjacency matrix.
    Theta1, Theta2 : matrices
        precision matrices.
    Cov1, Cov2 : matrices
        covariance matrices.
    L : list
        list of latent nodes..
    I : list
        list of intervention targets.

    Returns
    -------
    IMAG : matrix
        IMAG.
    ThetaObs1, ThetaObs2 : matrices
        precision matrices for observed variables.
    CovObs1, CovObs2 : matrices
        covariance matrices for observed variables.
    '''

    A = B.copy()
    A[np.where(A)] = 1
    p_V = B.shape[0]
    V = list(np.arange(p_V))
    O = V.copy()
    for l in L:
        O.remove(l)

    CovObs1 = Cov1[O][:,O]
    ThetaObs1 = marginal_theta_from_cov(Cov1,O)
    CovObs2 = Cov2[O][:,O]
    ThetaObs2 = marginal_theta_from_cov(Cov2,O)
    # just one setting for now
    I_all = [[],I]
    IMAG = DAG_I_2_IMAG(A,L,I_all)

    return IMAG, ThetaObs1, CovObs1, ThetaObs2, CovObs2


def reduce_2_observed_multiple(B_all,Theta_all,Cov_all,L,I_all):
    '''
    Objective: returns the observed Covariance and precision matrices

    Parameters
    ----------
    B_all : list of matrices
        list of adjacency matrices.
    Theta_all : list of matrices
        list of precision matrices.
    Cov_all : list of matrices
        list covariance matrices.
    L : list
        list of latent nodes..
    I_all : list of lists
        set of interventions, for instance, I_all = [[],[1,2],[3]]
        empty set denotes observational model.


    Returns
    -------
    IMAG : matrix
        IMAG.
    ThetaObs_all : list of matrices
        list of precision matrices for observed variables.
    CovObs_all : list of matrices
        list of covariance matrices for observed variables.

    '''
    
    A = B_all[0].copy()
    A[np.where(A)] = 1
    p_V = A.shape[0]
    V = list(np.arange(p_V))
    O = V.copy()
    for l in L:
        O.remove(l)

    CovObs_all = [Cov_all[i][O][:,O] for i in range(len(Cov_all))]
    ThetaObs_all = [marginal_theta_from_cov(Cov_all[i],O) for i in range(len(Cov_all))]
        
    IMAG = DAG_I_2_IMAG(A,L,I_all)

    return IMAG, ThetaObs_all, CovObs_all
