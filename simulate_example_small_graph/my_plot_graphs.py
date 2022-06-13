#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot the graphs for results of our algorithm

Reviewed on February, 2022.
"""

import numpy as np
import graphviz
#import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors

save_graphs_folder = './my_outputs_gauss_discretized'

node_names = ['Z','Y','X','F_12','F_13','F_23']


## 0: no edge
## 1: -o
## 2: -> (arrowhead)
## 3: - (tail)
# build the small graph
A_true = np.zeros([6,6])
A_true[3,1] = 2; A_true[3,2] = 2; A_true[4,0] = 2; A_true[4,1] = 2; A_true[4,2] = 2; A_true[5,0] = 2;
A_true[1,3] = 3; A_true[2,3] = 3; A_true[0,4] = 3; A_true[1,4] = 3; A_true[2,4] = 3; A_true[0,5] = 3;
A_true[1,2] = 2; A_true[2,1] = 3; A_true[1,0] = 2; A_true[0,1] = 2

# print the true graph
G = graphviz.Digraph()
G.attr('node', shape='circle')
for i in range(len(node_names)):
    G.node(node_names[i])


A = A_true
directed_edges = list(zip(np.where((A==2)*(A.T==3))[0],np.where((A==2)*(A.T==3))[1]))
bidirected_edges = list(zip(np.where((A+A.T)==4)[0],np.where((A+A.T)==4)[1]))
bidirected_edges = [edge for edge in bidirected_edges if edge[0]<edge[1]]
joker_edges = list(zip(np.where((A==2)*(A.T==1))[0],np.where((A==2)*(A.T==1))[1]))

for edge in directed_edges:
    G.edge(node_names[edge[0]],node_names[edge[1]])

for edge in bidirected_edges:
    G.edge(node_names[edge[0]],node_names[edge[1]],arrowhead='normal',arrowtail='normal',dir='both')
    
for edge in joker_edges:
    G.edge(node_names[edge[0]],node_names[edge[1]],arrowhead='normal',arrowtail='odot',dir='both')

G.render("./my_outputs_gauss_discretized/graph-env-model-1-true.gv", view=False)


#%
for run_no in range(11,31):
    for samples in [20000]:
        G = graphviz.Digraph()
        G.attr('node', shape='circle')
        for i in range(len(node_names)):
            G.node(node_names[i])

        
        A = np.genfromtxt("./my_outputs_gauss_discretized/graph-env-model-1-"+str(run_no)+"-"+str(samples)+".csv",delimiter=',')
        directed_edges = list(zip(np.where((A==2)*(A.T==3))[0],np.where((A==2)*(A.T==3))[1]))
        bidirected_edges = list(zip(np.where((A+A.T)==4)[0],np.where((A+A.T)==4)[1]))
        bidirected_edges = [edge for edge in bidirected_edges if edge[0]<edge[1]]
        joker_edges = list(zip(np.where((A==2)*(A.T==1))[0],np.where((A==2)*(A.T==1))[1]))

        for edge in directed_edges:
            G.edge(node_names[edge[0]],node_names[edge[1]])

        for edge in bidirected_edges:
            G.edge(node_names[edge[0]],node_names[edge[1]],arrowhead='normal',arrowtail='normal',dir='both')
            
        for edge in joker_edges:
            G.edge(node_names[edge[0]],node_names[edge[1]],arrowhead='normal',arrowtail='odot',dir='both')

        G.render("./my_outputs_gauss_discretized/graph-env-model-1-"+str(run_no)+"-"+str(samples)+".gv", view=False)

