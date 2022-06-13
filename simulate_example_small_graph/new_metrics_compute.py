import sys
sys.path.append('/usr/local/lib/python3.7/site-packages')
#Adding path where graphviz is found. One needs to add local path where libraries are present.
from graphviz import Digraph
import csv
import numpy as np
import pandas as pd
import argparse
import pyAgrum as gum


n_sample = 20000
run_vec = list(range(11,31))
model_count = 0

precision_skeleton_matrix= np.zeros(20)
Recall_skeleton_matrix= np.zeros(20)
precision_edge_matrix= np.zeros(20)
#Recall_edge_matrix= np.zeros((3,4))

## 0: no edge
## 1: -o
## 2: -> (arrowhead)
## 3: - (tail)

model_arcs_list=list()
a={tuple(('Y','X')):2,tuple(('X','Y')):3, tuple(('Z','Y')):2, tuple(('Y','Z')):2 }
model_arcs_list.append(a)



for run_count in range(len(run_vec)):
    filename= "./JCI_outputs_gauss_discretized/graph_disc_model-"+str(model_count+1)+"-jci-alpha0.05-"+str(run_vec[run_count])+"-"+str(n_sample)+".csv"            
    df=pd.read_csv(filename)
    print(df.columns)
    considered_edges=list()
    True_arcs=[i for i in model_arcs_list[model_count].keys()]
    dict_arcs= model_arcs_list[model_count]
    print(True_arcs)
    count_edges=0
    count_true_edges=int(len(True_arcs)/2)
    for index_row,rows in df.iterrows():
        a=rows[str(df.columns[0])]
        for index_col,col in enumerate(list(df.columns)[1:]):
            b=col
            if(tuple((a,b)) not in considered_edges):
                considered_edges=considered_edges+[tuple((a,b)),tuple((b,a))]
                if ((str(a).split('_')[0]!='F') & (str(b).split('_')[0]!='F')):
                    if ((df.iloc[index_row,index_col+1]!= 0) | (df.iloc[index_col,index_row+1]!= 0)) :
                        count_edges=count_edges+1
                        if (tuple((a,b)) in True_arcs):
                            precision_skeleton_matrix[run_count]=precision_skeleton_matrix[run_count]+1

                            if ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 3)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==3)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 3) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==3) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 1)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==1)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==1) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 1)) & ((dict_arcs[tuple((a,b))]==1) &(dict_arcs[tuple((b,a))]==1)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1  

    temp=precision_skeleton_matrix[run_count]
    if count_edges==0:
        count_edges=1
    precision_skeleton_matrix[run_count]=temp/count_edges
    precision_edge_matrix[run_count]=precision_edge_matrix[run_count]/count_edges
    Recall_skeleton_matrix[run_count]=temp/count_true_edges
    print(count_edges)
    print(count_true_edges)


precision_skeleton_matrix_jci=precision_skeleton_matrix
precision_edge_matrix_jci=precision_edge_matrix
Recall_skeleton_matrix_jci=Recall_skeleton_matrix

#%%
PS_jci=np.mean(precision_skeleton_matrix_jci)
AO_jci=np.mean(precision_edge_matrix_jci)
RS_jci=np.mean(Recall_skeleton_matrix_jci)

#%%%%%

n_sample = 20000
run_vec = list(range(11,31))
model_count = 0

precision_skeleton_matrix= np.zeros(20)
Recall_skeleton_matrix= np.zeros(20)
precision_edge_matrix= np.zeros(20)
#Recall_edge_matrix= np.zeros((3,4))

## 0: no edge
## 1: -o
## 2: -> (arrowhead)
## 3: - (tail)

model_arcs_list=list()
a={tuple(('Y','X')):2,tuple(('X','Y')):3, tuple(('Z','Y')):2, tuple(('Y','Z')):2 }
model_arcs_list.append(a)

for run_count in range(len(run_vec)):
    filename= "./psiFCI_outputs_gauss_discretized/graph_disc_model-"+str(model_count+1)+"-alpha0.05-"+str(run_vec[run_count])+"-"+str(n_sample)+".csv"
    df=pd.read_csv(filename)
    print(df.columns)
    considered_edges=list()
    True_arcs=[i for i in model_arcs_list[model_count].keys()]
    dict_arcs= model_arcs_list[model_count]
    print(True_arcs)
    count_edges=0
    count_true_edges=int(len(True_arcs)/2)
    for index_row,rows in df.iterrows():
        a=rows[str(df.columns[0])]
        for index_col,col in enumerate(list(df.columns)[1:]):
            b=col
            if(tuple((a,b)) not in considered_edges):
                considered_edges=considered_edges+[tuple((a,b)),tuple((b,a))]
                if ((str(a).split('_')[0]!='F') & (str(b).split('_')[0]!='F')):
                    if ((df.iloc[index_row,index_col+1]!= 0) | (df.iloc[index_col,index_row+1]!= 0)) :
                        count_edges=count_edges+1
                        if (tuple((a,b)) in True_arcs):
                            precision_skeleton_matrix[run_count]=precision_skeleton_matrix[run_count]+1

                            if ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 3)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==3)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 3) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==3) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 1)) & ((dict_arcs[tuple((a,b))]==2) &(dict_arcs[tuple((b,a))]==1)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 2)) & ((dict_arcs[tuple((a,b))]==1) &(dict_arcs[tuple((b,a))]==2)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1
                            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 1)) & ((dict_arcs[tuple((a,b))]==1) &(dict_arcs[tuple((b,a))]==1)):
                                precision_edge_matrix[run_count]=precision_edge_matrix[run_count]+1  

    temp=precision_skeleton_matrix[run_count]
    if count_edges==0:
        count_edges=1
    precision_skeleton_matrix[run_count]=temp/count_edges
    precision_edge_matrix[run_count]=precision_edge_matrix[run_count]/count_edges
    Recall_skeleton_matrix[run_count]=temp/count_true_edges
    print(count_edges)
    print(count_true_edges)


PS=np.mean(precision_skeleton_matrix)
AO=np.mean(precision_edge_matrix)
RS=np.mean(Recall_skeleton_matrix)



# import matplotlib.pyplot as plt
# #"#E69F00", "#56B4E9", "#009E73"

# for model_count in range(len(model_vec)):

#     Legends=[r"Precision-Skeleton-$\Psi$-FCI","Precision-Skeleton-JCI","Recall-Skeleton-$\Psi$-FCI","Recall-Skeleton-JCI","Accuracy Orientation-$\Psi$-FCI ","Accuracy Orientation-JCI"]
#     fig, ax=plt.subplots()
#     #bars = ('200k', '100k', '50k', '10k')
#     bars = ('10k','50k','100k','200k')
#     y_pos1 = 2*np.arange(len(bars))
#     # Create bars
#     ax.bar(y_pos1, PS[model_count,:],0.2, yerr= np.std(precision_skeleton_matrix,axis=2)[model_count,:], color="#E69F00",edgecolor='black')
#     ax.set_title("Model"+str(model_vec[model_count])+"Metrics")

#     ax.bar(y_pos1+0.2, PS_jci[model_count,:],0.2, yerr= np.std(precision_skeleton_matrix_jci,axis=2)[model_count,:], color="#E69F00",hatch="/",edgecolor='black')

#     # Create bars
#     ax.bar(y_pos1+0.4, RS[model_count,:],0.2,yerr= np.std(Recall_skeleton_matrix,axis=2)[model_count,:],color="#009E73",edgecolor='black')
#     ax.bar(y_pos1+0.6, RS_jci[model_count,:],0.2,yerr= np.std(Recall_skeleton_matrix_jci,axis=2)[model_count,:],color="#009E73",hatch="/",edgecolor='black')
    
#     # Create bars
#     ax.bar(y_pos1+0.8, AO[model_count,:],0.2,yerr= np.std(precision_edge_matrix,axis=2)[model_count,:], color="#56B4E9",edgecolor='black')
#     ax.bar(y_pos1+1.0, AO_jci[model_count,:],0.2, yerr= np.std(precision_edge_matrix_jci,axis=2)[model_count,:], color="#56B4E9",hatch="/",edgecolor='black')
#     ax.set_xticks(y_pos1+0.4)
#     ax.set_xticklabels(bars)
    
#     ax.set_ylim(0,2)
#     plt.legend(Legends,loc=1)
#     plt.tight_layout()
#     plt.savefig("Model"+str(model_vec[model_count])+"-Comparisons.png")
#     #plt.show()




#G.render()