#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import string
import matplotlib.pyplot as plt




def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
plt.rcParams["figure.figsize"] = cm2inch(16,8)
plt.rcParams.update({'font.size':10})
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)




systems=["FA_G-prot"]
begin_time_ps=["15000000"]
end_time_ps=["-1"]



lipid_dir = {"FA_G-prot":"POPC, POPS, POPA, POPG, POPE, CHOL, DPGM3, DPSM, PAP6"}


# # Plotting RDF


def load_rdf (lipid):
    '''Load occupancy data for considering only the headgroups of the lipids'''
    data_dir =  {"FA_G-prot":"",
                "FA_no_G-prot":"",
                "IA":"",
                "PA":""}
    for s in systems:
        data_loc = []
        for r in range(3):
            try:
                data_loc.append(np.loadtxt('RDF/{0:s}_{1:d}_nowat_{2:s}_rdf.xvg'.format(s,r,lipid), comments=('#', '@')))
            except OSError:
                #print ('No such file: RDF/{0:s}_{1:d}_nowat_{2:s}_rdf.xvg'.format(s,r,lipid))
                continue
        data_dir[s]=data_loc
    return data_dir



def plot_rdf (systems):
    '''For RDF'''
    x_max = 2
    color_list = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

    fig, ax = plt.subplots(nrows=2, ncols=4, sharey=True, sharex=True)

    unique_lipids = ['POPC','POPS','POPE','POPA','POPG', 'CHOL', 'PAP6', 'DPGM3', 'DPSM']

    for idx, l in enumerate(unique_lipids[:4]):
        dir_data = load_rdf(l)

        for ndx, s in enumerate(systems):
            if len(dir_data[s])==0:
                continue
            else:
                data = dir_data[s]
                #print (np.array(data).shape)
                min_shape = np.min([ i.shape[0] for i in data ] ) - 1 
                data_array = np.array([ i[:min_shape,:] for i in data ])
                
                Ava = np.average(data_array, axis=0)
                Std = np.std(data_array, axis=0)
                ax[0,idx].plot(Ava[:,0], Ava[:,1], label=s, color=color_list[ndx], alpha=0.7)
                ax[0,idx].fill_between(Ava[:,0], Ava[:,1]-Std[:,1], Ava[:,1]+Std[:,1], alpha=0.7, edgecolor=color_list[ndx], facecolor=color_list[ndx])
                ax[0,0].set_xlim(0,x_max)
                ax[0,idx].set_title(l)
                ax[0,idx].axvline(x=0.7,ymin=0, ymax=3, color='black', linewidth=1.5)
                ax[0,idx].axvline(x=1.2,ymin=0, ymax=3, color='black', linestyle='--', linewidth=1.5)
                ax[0,idx].axvline(x=1.6,ymin=0, ymax=3, color='black', linestyle=':', linewidth=1.5)
                #ax[0,0].legend(loc='best')

    for idx, l in enumerate(unique_lipids[4:][:-1]):
        dir_data = load_rdf(l)
        #print (l)

        for ndx, s in enumerate(systems):
            if len(dir_data[s])==0:
                continue
            else:
                data = dir_data[s]
                #print (np.array(data).shape)
                min_shape = np.min([ i.shape[0] for i in data ] ) - 1 
                data_array = np.array([ i[:min_shape,:] for i in data ])
                
                Ava = np.average(data_array, axis=0)
                Std = np.std(data_array, axis=0)
                ax[1,idx].plot(Ava[:,0], Ava[:,1], label=s, color=color_list[ndx], alpha=0.7)
                #ax[1,idx].plot(Ava[:,0], Ava[:,1], color=color_list[ndx], alpha=0.7)
                ax[1,idx].fill_between(Ava[:,0], Ava[:,1]-Std[:,1], Ava[:,1]+Std[:,1], alpha=0.7, edgecolor=color_list[ndx], facecolor=color_list[ndx])
                ax[1,0].set_xlim(0,x_max)
                ax[1,idx].set_title(l)
                #ax[1,0].legend(loc='best')
                #ax[1,idx].axvline(x=0.4,ymin=0, ymax=3)
                #ax[1,idx].axvline(x=0.4,ymin=0, ymax=3, color='black', linestyle='--')
                ax[1,idx].axvline(x=0.7,ymin=0, ymax=3, color='black', linewidth=1.5)
                ax[1,idx].axvline(x=1.2,ymin=0, ymax=3, color='black', linestyle='--', linewidth=1.5)
                ax[1,idx].axvline(x=1.6,ymin=0, ymax=3, color='black', linestyle=':', linewidth=1.5)

    ax[0,0].legend(loc='upper center', bbox_to_anchor=(2, -1.5),
              fancybox=True, shadow=True, ncol=3, frameon=False)

    #fig.tight_layout(pad=3.0)
    plt.subplots_adjust(hspace=0.3)
    
    plt.savefig('RDF/RDFs.png', dpi=300, bbox_inches='tight')
    return




systems = ['FA_G-prot','FA_no_G-prot', 'IA', 'PA']
plot_rdf (systems)
