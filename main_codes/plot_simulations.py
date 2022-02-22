#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
create figures for simulations run thru run_simulations.py

"""
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import pickle
from config import SIMULATIONS_ESTIMATED_FOLDER, SIMULATIONS_FIGURES_FOLDER

'''
simulations are run for (num_repeat x p x num_latent x num_int x density x num_saples)
when summed over the repeated trials, it becomes 5dim (p x num_latent x num_int x density x num_saples)
interpret accordingly.
'''

def load_simulation_res(filename):
    with open(SIMULATIONS_ESTIMATED_FOLDER+'/'+filename+'.pkl','rb') as f:
        res = pickle.load(f)
    
    p_list = res['p_list']
    num_latent_list = np.asarray(res['num_latent_list'])
    num_int_list = np.asarray(res['num_int_list'])
    density_list = np.asarray(res['density_list'])
    #num_int_settings = res['num_int_settings']
    num_samples_list = np.asarray(res['num_samples_list'])

    I_tp = res['I_tp']
    I_fp = res['I_fp']
    I_fn = res['I_fn']
    #I_pasp_tp = res['I_pasp_tp']
    #I_pasp_fp = res['I_pasp_fp']
    #I_pasp_fn = res['I_pasp_fn']
    S_Delta_sizes = res['S_Delta_sizes']
    t = res['t']

    return res, I_tp, I_fp, I_fn, S_Delta_sizes, t, \
            p_list, num_latent_list, num_int_list, density_list, num_samples_list


def scores(I_tp,I_fp,I_fn):
    # computes scores by summing up the results over num_repeat trials
    I_precision = np.sum(I_tp,0) / (np.sum(I_tp,0)+np.sum(I_fp,0))
    I_recall =  np.sum(I_tp,0) / (np.sum(I_tp,0)+np.sum(I_fn,0))
    I_f1 = np.sum(I_tp,0) / (np.sum(I_tp,0)+(np.sum(I_fp,0)+np.sum(I_fn,0))/2)


    return I_precision, I_recall, I_f1

def plot_res(precision,recall,f1,l,i,c,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='F1',save_fig=True):
    if stat == 'F1':
        scores = f1
    elif stat == 'recall':
        scores = recall
    elif stat == 'precision':
        scores = precision
        
    plt.figure(mode +' - density %.1f'%density_list[c])
    for p_idx in range(len(p_list)):
        plt.plot(num_samples_list.astype('str'),scores[p_idx,l,i,c],linestyle=linestyle,linewidth=linewidth,marker=markers[p_idx],markersize=markersize)

    
    legend_str = ['p=%d'%p for p in p_list] 
    #plt.title('|L|=%d, |I|=%d,  density = %d'%(num_latent_list[l],num_int_list[i],density_list[c]),size=title_size)
    plt.grid()
    plt.xlabel('Number of Samples',size=xlabel_size)
    #plt.ylabel('Precision of estimating intervention targets',size=12)
    plt.ylabel(stat,size=ylabel_size)
    # set ylim
    ylim = np.min(scores[:,l,i,c])
    plt.ylim([ylim-0.1,1])
    plt.xticks(fontsize=xticks_size)
    plt.yticks(fontsize=yticks_size)
    plt.legend(legend_str,fontsize=legend_size,loc=legend_loc)
    plt.tight_layout()
    if save_fig is True:
        plt.savefig(SIMULATIONS_FIGURES_FOLDER+'/'+mode+'_'+stat + '_latent%d'%num_latent_list[l]+ \
                    '_int%d'%num_int_list[i] +'_density%d'%density_list[c] + '.eps')

    plt.close()

#%%
markers = ['o','v','P','D','X']

xticks_size = 14
yticks_size = 14
xlabel_size = 18
ylabel_size = 18
legend_size = 12
legend_loc = 'lower right'
linewidth = 2.6
linestyle = '--'
markersize = 9
title_size = 18

#%% load results for F->K recovery, shift intervention
res_shift, I_tp_shift, I_fp_shift, I_fn_shift, S_Delta_sizes_shift, t_shift, \
        p_list, num_latent_list, num_int_list, density_list, num_samples_list\
         = load_simulation_res('shift_W_3')

I_precision_shift, I_recall_shift, I_f1_shift = scores(I_tp_shift, I_fp_shift, I_fn_shift)

#%% individual figures
save_fig = True

'f1 scores'
plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='F1',save_fig=save_fig)
#plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,1,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='F1',save_fig=save_fig)

'precision rates'
plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='precision',save_fig=save_fig)
#plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,1,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='precision',save_fig=save_fig)

'recall rates'
plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='recall',save_fig=save_fig)
#plot_res(I_precision_shift,I_recall_shift,I_f1_shift,0,1,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='shift',stat='recall',save_fig=save_fig)


#%% load results for F->K recovery, increased variance
res_var, I_tp_var, I_fp_var, I_fn_var, S_Delta_sizes_var, t_var, \
        p_list, num_latent_list, num_int_list, density_list, num_samples_list\
         = load_simulation_res('variance_W_1')

I_precision_var, I_recall_var, I_f1_var = scores(I_tp_var, I_fp_var, I_fn_var)
          
#%%
save_fig = True

'f1 scores'
plot_res(I_precision_var,I_recall_var,I_f1_var,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='variance',stat='F1',save_fig=save_fig)

'precision rates'
plot_res(I_precision_var,I_recall_var,I_f1_var,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='variance',stat='precision',save_fig=save_fig)

'recall rates'
plot_res(I_precision_var,I_recall_var,I_f1_var,0,0,0,p_list,num_latent_list,num_int_list,density_list,num_samples_list,mode='variance',stat='recall',save_fig=save_fig)

