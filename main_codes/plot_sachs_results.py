#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Sachs protein signaling network results here

"""
import numpy as np
import graphviz

import pickle
import matplotlib.colors as mcolors

from realdata.sachs.sachs_meta import SACHS_ESTIMATED_FOLDER, SACHS_FIGURES_FOLDER, nnodes
from realdata.sachs.sachs_meta import true_dag_recent as true_dag_recent

reference = 'NessSachs2016'
utigsp_ci_test = 'gauss'

ALGS2COLORS = dict(zip(['ours','utigsp_gauss', 'utigsp_star_gauss', 'utigsp_hsic','utigsp_star_hsic'],\
                       mcolors.BASE_COLORS))
ALGS2MARKERS = {'ours':'o','utigsp_gauss': 'P', 'utigsp_star_gauss': '*', 'utigsp_hsic': 'X', 'utigsp_star_hsic': 'x'}
    
xticks_size = 14
yticks_size = 14
xlabel_size = 18
ylabel_size = 18
legend_size = 10
legend_loc = 'upper left'

def read_single_result(file):
    for vals in file.keys():
        est_edges = file[vals]['est_edges']
        est_skeleton = file[vals]['est_skeleton']
        K_hat = file[vals]['K_hat']
        K_pasp_hat = file[vals]['K_pasp_hat']
        time = file[vals]['time']
        S_Delta_size = file[vals]['S_Delta_size']


    return est_edges, est_skeleton,  K_hat, K_pasp_hat, S_Delta_size, time


B = true_dag_recent.copy()
np.fill_diagonal(B, 0)
B_skeleton = B + B.T
B_skeleton[np.where(B_skeleton)] = 1

B_skeleton_edges = list(zip(np.where(B_skeleton)[0],np.where(B_skeleton)[1]))

n_possible_skeleton = int(nnodes*(nnodes-1)/2)
n_true_skeleton = int(np.sum(B_skeleton)/2)


#%%
# load our results
with open(SACHS_ESTIMATED_FOLDER+'/sachs_run_2.pkl', 'rb') as f:
    res_ours = pickle.load(f)    
        
est_edges, est_skeleton,  K_hat, K_pasp_hat, S_Delta_size, time = read_single_result(res_ours)
    
    
#%%
'mini post-processing'

K = [K_hat['F_%d'%i] for i in range(len(K_hat.items()))]
K_pasp = [[list(K_pasp_hat['F_%d'%i][j]) for j in range(len(K_pasp_hat['F_%d'%i]))] for i \
          in range(len(K_pasp_hat.items()))]

pasp_process = [[] for i in range(nnodes)]
    
    
for idx_setting in range(len(K)):
    for idx_k in range(len(K[idx_setting])):
        pasp_process[K[idx_setting][idx_k]].extend(K_pasp[idx_setting][idx_k])
        
all_targets = []
for i in range(len(K)):
    all_targets.extend(K[i])
    
all_targets = list(set(all_targets))
             
pasp_process = [list(set(pasp_process[i])) for i in range(nnodes)]
pasp_edges = np.zeros((nnodes,nnodes))
for i in range(nnodes):
    pasp_edges[pasp_process[i],i] = 1
    
all_edges = list(zip(np.where(pasp_edges)[0],np.where(pasp_edges)[1]))

directed_edges = []
bi_directed_edges = []
#bi_directed_edges = list(zip(np.where(pasp_edges*pasp_edges.T)[0],np.where(pasp_edges*pasp_edges.T)[1]))
joker_edges = []

for edge in all_edges:
    # bi-directed edges
    if (edge[1],edge[0]) in all_edges:
        bi_directed_edges.append(edge)
    # it exists one-way. if both nodes are among targets, then we found the direction of the edge
    elif set(edge).issubset(all_targets) is True:
        directed_edges.append(edge)
    else:
        joker_edges.append(edge)

spouses = bi_directed_edges[:int(len(bi_directed_edges)/2)]

#%%
# visualization: observed nodes, i.e. proteins, only.

pros = ['Raf','Mek','Plcg','PIP2','PIP3','Erk','Akt','PKA','PKC','P38','Jnk']

'consider a graph without F nodes'
G = graphviz.Digraph()
G.attr('node', shape='circle')
G.node(pros[0])
G.node(pros[1])
G.node(pros[2])
G.node(pros[3])
G.node(pros[4])
G.node(pros[5])
G.node(pros[6])
G.node(pros[7])
G.node(pros[8])
G.node(pros[9])
G.node(pros[10])

for edge in directed_edges:
    if edge in B_skeleton_edges:
        G.edge(pros[edge[0]],pros[edge[1]],color='blue')
    else:
        G.edge(pros[edge[0]],pros[edge[1]])

    
for edge in spouses:
    if edge in B_skeleton_edges:
        G.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='normal',dir='both',color='blue')
    else:
        G.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='normal',dir='both')

for edge in joker_edges:
    if edge in B_skeleton_edges:
        G.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='odot',dir='both',color='blue')
    else:
        G.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='odot',dir='both')

G.render(SACHS_FIGURES_FOLDER+'/sachs_observed.gv', view=True)

#%%
# visualization: augmented graph with F-nodes.

pros = ['Raf','Mek','Plcg','PIP2','PIP3','Erk','Akt','PKA','PKC','P38','Jnk']
dot = graphviz.Digraph(comment='I-MAG')

dot.attr('node', shape='circle')
dot.node('F_1')
dot.node('F_2')
dot.node('F_3')
dot.node('F_4')
dot.node('F_5')
dot.node('F_6')
dot.node('F_7')
dot.node('F_8')
dot.node('F_9')
dot.node('F_10')
dot.node('F_11')
dot.node('F_12')
dot.node('F_13')
dot.node('F_14')
dot.node('F_15')

dot.node(pros[0])
dot.node(pros[1])
dot.node(pros[2])
dot.node(pros[3])
dot.node(pros[4])
dot.node(pros[5])
dot.node(pros[6])
dot.node(pros[7])
dot.node(pros[8])
dot.node(pros[9])
dot.node(pros[10])

'F - edges'
dot.edges([('F_1',pros[7]),('F_2',pros[0]),('F_2',pros[1]),('F_3',pros[3]),('F_4',pros[10]),\
            ('F_7',pros[2]),('F_7',pros[7]),('F_8',pros[7]),('F_9',pros[2]),('F_9',pros[7]),\
            ('F_10',pros[0]),('F_10',pros[1]),('F_10',pros[3]),('F_11',pros[0]),('F_11',pros[1]),\
            ('F_12',pros[0]),('F_12',pros[1]),('F_13',pros[3]),('F_13',pros[10]),\
            ('F_14',pros[3]),('F_14',pros[6]),('F_15',pros[3]),('F_15',pros[10])])

# 'other edges'
for edge in directed_edges:
    if edge in B_skeleton_edges:
        G.edge(pros[edge[0]],pros[edge[1]])
    else:
        G.edge(pros[edge[0]],pros[edge[1]])

    
for edge in spouses:
    if edge in B_skeleton_edges:
        dot.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='normal',dir='both')
    else:
        dot.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='normal',dir='both')

for edge in joker_edges:
    if edge in B_skeleton_edges:
        dot.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='odot',dir='both')
    else:
        dot.edge(pros[edge[0]],pros[edge[1]],arrowhead='normal',arrowtail='odot',dir='both')    

dot.render(SACHS_FIGURES_FOLDER+'/sachs_f.gv', view=True)

