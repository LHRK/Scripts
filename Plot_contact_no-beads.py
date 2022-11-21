#!/usr/bin/env python
# coding: utf-8

import MDAnalysis as md
import numpy as np
import matplotlib.pyplot as plt
import string
import os
import argparse


parser = argparse.ArgumentParser(description='Script for plotting only for the TMD')
parser.add_argument('-s', dest='system',nargs='+', action='store', type=str, help='-s for the systems, eg -s "name1" "name2"')

args = parser.parse_args()
system = args.system

print ("systems:")
print (system)


def get_lipid_prob (name):
    bulk = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_bulk.npy'.format(name), allow_pickle=True)
    prot = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_protein.npy'.format(name),allow_pickle=True)
    V_cyl = np.pi * 7**2 * 14 #Angstrom^3
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    effective_vol_per_bead = bulk_a / V_cyl
    lipids = bulk[0][0]
    prot_a = np.average(prot.T[:,1,:,:], axis=2) 
    P = []
    for l in range(len(lipids)):
        diff_loc = prot_a[l,:] * effective_vol_per_bead[l] / bulk_a[l] * effective_vol_per_bead[l]
        logs_loc = [log1p(i) for i in diff_loc]
        P.append(logs_loc)
    return np.array(P), lipids

def get_lipid_diff (name):
    bulk = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_bulk.npy'.format(name))
    prot = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_protein.npy'.format(name))
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    lipids = bulk[0][0]
    prot_a = np.average(prot.T[:,1,:,:], axis=2)
    P = []
    for l in range(len(lipids)):
        print (np.any(prot_a[l,:]==0))
        print (np.any(bulk_a[l]==0))
        diff_loc = prot_a[l,:] / bulk_a[l]
        diffs = [float(i) for i in diff_loc]
        #print (np.any(np.array(diffs)==0))
        #print (np.array(diff_loc).shape)
        #print (diffs)
        P.append(diffs)
    return np.array(P), lipids

def extract_prob_new (name):
    bulk = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_bulk.npy'.format(name),allow_pickle=True)
    prot = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_protein.npy'.format(name),allow_pickle=True)

    u = md.Universe('GRO/{0:s}_nowat.gro'.format(name))
    V_cyl = np.pi * 7**2 * 14 #Angstrom^3
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    effective_vol_per_bead = V_cyl/bulk_a
    
    lipids = bulk[0][0]
    
    bulk_sum = np.sum(bulk.T[:,1,:])
    prot_sum = np.sum(prot.T[:,1,:,:])
    norm_EV_prot = (V_cyl/bulk_sum) * prot_sum
    norm_EV_bulk = (V_cyl/bulk_sum) * bulk_sum
        
    #Normalizing for volume
    new_shape =  prot.T[:,1,:,:].shape
    prot_normalized = np.zeros(new_shape)
    new_shape_b = bulk.T[:,1,:].shape
    bulk_normalized = np.zeros(new_shape_b)
    
    for l in range(len(lipids)):
        #prot_normalized[l,:,:] = prot.T[l,1,:,:] * effective_vol_per_bead[l]
        #bulk_normalized[l,:]   = bulk.T[l,1,:] * effective_vol_per_bead[l]
        prot_normalized[l,:,:] = prot.T[l,1,:,:] * norm_EV_prot
        bulk_normalized[l,:]   = bulk.T[l,1,:] * norm_EV_bulk
    
    prot_a = np.average(prot_normalized, axis=2)
    
    
    ntypes_lipid, n_resids, nframes = prot.T[:,1,:,:].shape

    prio = 1

    my = np.zeros([len(lipids)])
    for ndx, l in enumerate(lipids):
        lipid_num = np.unique(u.select_atoms('resname {0:s}'.format(l)).resids).shape[0]
        all_num   = np.unique(u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)').resids).shape[0]
        ratio     = ( lipid_num / all_num ) 
        #print (ratio)
        my[ndx]   = ratio

    top = np.zeros([ntypes_lipid, n_resids])
    bot = np.zeros([ntypes_lipid])

    for i in range(ntypes_lipid):
        for r in range(n_resids):
            t = (prio * my[i] + nframes * prot_a[i,r] ) / (prio + nframes)
            b = (prio * my[i] + nframes * bulk_a[i] ) / (prio + nframes)
            top [i,r] = t 
            bot [i]   = b 

    P = np.zeros([ntypes_lipid, n_resids])

    for n in range(ntypes_lipid):
        p = np.log(top[n,:] / bot[n])
        P[n,:] = p
    return P, lipids

def extract_prob_new_last5 (name):
    bulk = np.load('CONTACT_NO_BEADS/LAST5/{0:s}_number_of_beads_bulk.npy'.format(name),allow_pickle=True)
    prot = np.load('CONTACT_NO_BEADS/LAST5/{0:s}_number_of_beads_protein.npy'.format(name),allow_pickle=True)

    u = md.Universe('GRO/{0:s}_nowat.gro'.format(name))
    V_cyl = np.pi * 7**2 * 14 #Angstrom^3
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    effective_vol_per_bead = bulk_a / V_cyl
    
    bulk_sum = np.sum(bulk.T[:,1,:])
    prot_sum = np.sum(prot.T[:,1,:,:])
    norm_EV_prot = (V_cyl/bulk_sum) * prot_sum
    norm_EV_bulk = (V_cyl/bulk_sum) * bulk_sum
    
    lipids = bulk[0][0]
    
    #Normalizing for volume
    new_shape =  prot.T[:,1,:,:].shape
    prot_normalized = np.zeros(new_shape)
    new_shape_b = bulk.T[:,1,:].shape
    bulk_normalized = np.zeros(new_shape_b)
    
    for l in range(len(lipids)):
        #prot_normalized[l,:,:] = prot.T[l,1,:,:] * effective_vol_per_bead[l]
        #bulk_normalized[l,:]   = bulk.T[l,1,:] * effective_vol_per_bead[l]
        prot_normalized[l,:,:] = prot.T[l,1,:,:] * norm_EV_prot
        bulk_normalized[l,:]   = bulk.T[l,1,:] * norm_EV_bulk
    
    prot_a = np.average(prot_normalized, axis=2)
    
    
    ntypes_lipid, n_resids, nframes = prot.T[:,1,:,:].shape

    prio = 1

    my = np.zeros([len(lipids)])
    for ndx, l in enumerate(lipids):
        lipid_num = np.unique(u.select_atoms('resname {0:s}'.format(l)).resids).shape[0]
        all_num   = np.unique(u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)').resids).shape[0]
        ratio     = ( lipid_num / all_num ) 
        #print (ratio)
        my[ndx]   = ratio

    top = np.zeros([ntypes_lipid, n_resids])
    bot = np.zeros([ntypes_lipid])

    for i in range(ntypes_lipid):
        for r in range(n_resids):
            t = (prio * my[i] + nframes * prot_a[i,r] ) / (prio + nframes)
            b = (prio * my[i] + nframes * bulk_a[i] ) / (prio + nframes)
            top [i,r] = t 
            bot [i]   = b 

    P = np.zeros([ntypes_lipid, n_resids])

    for n in range(ntypes_lipid):
        p = np.log(top[n,:] / bot[n])
        P[n,:] = p
    return P, lipids

def extract_prob_final (name):
    bulk = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_bulk.npy'.format(name),allow_pickle=True)
    prot = np.load('CONTACT_NO_BEADS/{0:s}_number_of_beads_protein.npy'.format(name), allow_pickle=True)

    u = md.Universe('GRO/{0:s}_nowat.gro'.format(name))
    V_cyl = np.pi * 7**2 * 14 #Angstrom^3
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    effective_vol_per_bead = V_cyl/np.sum(bulk_a)

    lipids = bulk[0][0]

    prot.T[:,1,:,:].shape
    #prot shape [nlipid_types, nres, nframes]

    n_lipidtypes, nres, nframes = prot.T[:,1,:,:].shape

    V_lipid_prot = np.zeros(nres)

    for r in range(nres):
        if np.mean(prot.T[:,1,r,:]) != 0 :
            V_lipid_prot[r] = np.average(np.sum(prot.T[:,1,r,:], axis=0).astype('float64')) * effective_vol_per_bead

    #Normalizing for volume
    new_shape =  prot.T[:,1,:,:].shape
    prot_normalized = np.zeros(new_shape)
    new_shape_b = bulk.T[:,1,:].shape
    bulk_normalized = np.zeros(new_shape_b)

    for l in range(n_lipidtypes):
        for f in range(nframes):
            prot_normalized[l,:,f] = prot.T[l,1,:,f] * (V_cyl / V_lipid_prot)

    for l in range(len(lipids)):
        bulk_normalized[l,:]   = bulk.T[l,1,:] * effective_vol_per_bead

    prot_a = np.average(prot_normalized, axis=2)


    ntypes_lipid, n_resids, nframes = prot.T[:,1,:,:].shape

    prio = 1

    my = np.zeros([len(lipids)])
    for ndx, l in enumerate(lipids):
        lipid_num = np.unique(u.select_atoms('resname {0:s}'.format(l)).resids).shape[0]
        all_num   = np.unique(u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)').resids).shape[0]
        ratio     = ( lipid_num / all_num ) 
        my[ndx]   = ratio

    top = np.zeros([ntypes_lipid, n_resids])
    bot = np.zeros([ntypes_lipid])

    for i in range(ntypes_lipid):
        for r in range(n_resids):
            t = (prio * my[i] + nframes * prot_a[i,r] ) / (prio + nframes)
            b = (prio * my[i] + nframes * bulk_a[i] ) / (prio + nframes)
            top [i,r] = t 
            bot [i]   = b 

    P = np.zeros([ntypes_lipid, n_resids])

    for n in range(ntypes_lipid):
        p = np.log(top[n,:] / bot[n])
        P[n,:] = p
    P[np.where(np.isnan(P)==True)] = 0
    return P, lipids

def extract_prob_final_last5 (name):
    bulk = np.load('CONTACT_NO_BEADS/LAST5/{0:s}_number_of_beads_bulk.npy'.format(name), allow_pickle=True)
    prot = np.load('CONTACT_NO_BEADS/LAST5/{0:s}_number_of_beads_protein.npy'.format(name), allow_pickle=True)

    u = md.Universe('GRO/{0:s}_nowat.gro'.format(name))
    V_cyl = np.pi * 7**2 * 14 #Angstrom^3
    bulk_a = np.average(bulk.T[:,1,:], axis=1)
    effective_vol_per_bead = V_cyl/np.sum(bulk_a)

    lipids = bulk[0][0]

    prot.T[:,1,:,:].shape
    #prot shape [nlipid_types, nres, nframes]

    n_lipidtypes, nres, nframes = prot.T[:,1,:,:].shape

    V_lipid_prot = np.zeros(nres)

    for r in range(nres):
        if np.mean(prot.T[:,1,r,:]) != 0 :
            V_lipid_prot[r] = np.average(np.sum(prot.T[:,1,r,:], axis=0).astype('float64')) * effective_vol_per_bead

    #Normalizing for volume
    new_shape =  prot.T[:,1,:,:].shape
    prot_normalized = np.zeros(new_shape)
    new_shape_b = bulk.T[:,1,:].shape
    bulk_normalized = np.zeros(new_shape_b)

    for l in range(n_lipidtypes):
        for f in range(nframes):
            prot_normalized[l,:,f] = prot.T[l,1,:,f] * (V_cyl / V_lipid_prot)

    for l in range(len(lipids)):
        bulk_normalized[l,:]   = bulk.T[l,1,:] * effective_vol_per_bead

    prot_a = np.average(prot_normalized, axis=2)
    #bulk_a = np.average(bulk_normalized, axis=1)

    ntypes_lipid, n_resids, nframes = prot.T[:,1,:,:].shape

    prio = 1

    my = np.zeros([len(lipids)])
    for ndx, l in enumerate(lipids):
        lipid_num = np.unique(u.select_atoms('resname {0:s}'.format(l)).resids).shape[0]
        all_num   = np.unique(u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)').resids).shape[0]
        ratio     = ( lipid_num / all_num ) 
        my[ndx]   = ratio

    top = np.zeros([ntypes_lipid, n_resids])
    bot = np.zeros([ntypes_lipid])

    for i in range(ntypes_lipid):
        for r in range(n_resids):
            t = (prio * my[i] + nframes * prot_a[i,r] ) / (prio + nframes)
            b = (prio * my[i] + nframes * bulk_a[i] ) / (prio + nframes)
            top [i,r] = t 
            bot [i]   = b 

    P = np.zeros([ntypes_lipid, n_resids])

    for n in range(ntypes_lipid):
        p = np.log(top[n,:] / bot[n])
        P[n,:] = p
    P[np.where(np.isnan(P)==True)] = 0
    return P, lipids

def plot (P, lipids):
    #resids = np.array([i+25 for i in range(len(P[0]))])
    resids  = np.arange(136, 419)
    nlipids, ncols = P.shape
    
    fig, ax = plt.subplots(nrows=nlipids, ncols=1, sharex=True, sharey=True, figsize=(10,20))
    TM_color_list = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
    alpha_value=0.5
    for idxp, p in enumerate(P):
        #if lipids[idx] == 'POPC':
        #    colorr = 'grey'
        #else:
        #    colorr = 'black'
        #plt.plot( p, label=lipids[idx], alpha=0.7, color=colorr)
        #print (idxp)
        #print (p.shape)
        ax[idxp].plot(p, 'o', color='black')
        ax[idxp].set_title(lipids[idxp])
        ax[idxp].axvspan(0, 30, alpha=alpha_value, color=TM_color_list[8], label='TM1')
        ax[idxp].axvspan(36, 63, alpha=alpha_value, color=TM_color_list[1], label='TM2') 
        ax[idxp].axvspan(83, 119, alpha=alpha_value, color=TM_color_list[2], label='TM3')
        ax[idxp].axvspan(125, 154, alpha=alpha_value, color=TM_color_list[3], label='TM4')
        ax[idxp].axvspan(165, 198, alpha=alpha_value, color=TM_color_list[4], label='TM5')
        ax[idxp].axvspan(205, 235, alpha=alpha_value, color=TM_color_list[5], label='TM6')
        ax[idxp].axvspan(240, 267, alpha=alpha_value, color=TM_color_list[6], label='TM7')
        ax[idxp].axvspan(271, 283, alpha=alpha_value, color=TM_color_list[7], label='H8')
        ax[idxp].grid(alpha=0.5)
        ax[idxp].set_ylim(-2,0.5)
        ax[idxp].hlines(0,0, 284)
        
    
    ax[0].legend(loc='center left', bbox_to_anchor=(1.1, 0.5),fancybox=True, shadow=False, ncol=2,frameon=False)
    xticks = np.arange(0,284, 44)
    xtick_labels = np.arange(136,420, 44)

    plt.xticks(xticks, xtick_labels)
    plt.xlabel('Resid')
    return


def main (systems):
	for sys in systems:
	    for r in range(3):
	        name = '{0:s}_{1:d}'.format(sys,r)
	        #print (name)
	        try:
	            P, lipids = extract_prob_final(name)
	        except FileNotFoundError:
	            print ('File not found')
	            print (name)
	            continue
	        plot(P*-1, lipids)
	        #title = '{0:s} Replica {1:d}'.format(sys, r+1)
	        #plt.title(title)
	        plt.tight_layout()
	        plt.savefig('CONTACT_NO_BEADS/'+name+'_Contact_no_beads_normalized_remove_prot_vol_final.png',dpi=300, bbox_inchex='tight')
	        plt.show()
	return


main(system)


