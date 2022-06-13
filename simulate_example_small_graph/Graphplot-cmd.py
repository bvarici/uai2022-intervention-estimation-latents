#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append('/usr/local/lib/python3.7/site-packages')
#Adding path where graphviz is found. One needs to add local path where libraries are present.
from graphviz import Digraph
import csv
import pandas as pd
import argparse
import pyAgrum as gum

parser=argparse.ArgumentParser()
parser.add_argument('--filename','-fn',help="Filename containing adjacencies",type=str)
parser.add_argument('--filename_out','-fo',help="Filename for output PDF",type=str)


args=parser.parse_args()
filename=args.filename
filename_out=args.filename_out
df=pd.read_csv(filename)
print(df.columns)
considered_edges=list()

G=Digraph(filename=filename_out)
for index_row,rows in df.iterrows():
    a=rows[str(df.columns[0])]
    for index_col,col in enumerate(list(df.columns)[1:]):
        b=col
        if(tuple((a,b)) not in considered_edges):
            considered_edges=considered_edges+[tuple((a,b)),tuple((b,a))]
            if ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 3)):
                G.edge(a,b)
            elif ((df.iloc[index_row,index_col+1]== 3) & (df.iloc[index_col,index_row+1]== 2)):
                G.edge(b,a)
            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 2)):
                G.edge(a,b,arrowhead='normal',arrowtail='normal',dir='both')
            elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 1)):
                G.edge(a,b,arrowhead='normal',arrowtail='odot',dir='both')
            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 2)):
                G.edge(a,b,arrowhead='odot',arrowtail='normal',dir='both')
            elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 1)):
                G.edge(a,b,arrowhead='odot',arrowtail='odot',dir='both')

G.render()


# In[2]:


# considered_edges=list()
# G1=Digraph(filename='graph_without_F')
# for index_row,rows in df.iterrows():
#     a=rows[str(df.columns[0])]
#     for index_col,col in enumerate(list(df.columns)[1:]):
#         b=col
#         if(tuple((a,b)) not in considered_edges):
#             considered_edges=considered_edges+[tuple((a,b)),tuple((b,a))]
#             if ((str(a).split('_')[0]!='F') & (str(b).split('_')[0]!='F')):
#                 if ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 3)):
#                     G1.edge(a,b)
#                 elif ((df.iloc[index_row,index_col+1]== 3) & (df.iloc[index_col,index_row+1]== 2)):
#                     G1.edge(b,a)
#                 elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 2)):
#                     G1.edge(a,b,arrowhead='normal',arrowtail='normal',dir='both')
#                 elif ((df.iloc[index_row,index_col+1]== 2) & (df.iloc[index_col,index_row+1]== 1)):
#                     G1.edge(a,b,arrowhead='normal',arrowtail='odot',dir='both')
#                 elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 2)):
#                     G1.edge(a,b,arrowhead='odot',arrowtail='normal',dir='both')
#                 elif ((df.iloc[index_row,index_col+1]== 1) & (df.iloc[index_col,index_row+1]== 1)):
#                     G1.edge(a,b,arrowhead='odot',arrowtail='odot',dir='both')

# G1


# In[ ]:




