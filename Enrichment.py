#!/usr/bin/env python
# coding: utf-8



import numpy as np
import MDAnalysis as md 
import matplotlib.pyplot as plt
import os
from scipy.stats import ttest_1samp


# # Enrichment and deplemention 

# #### The cutoff for the first lipid shell is set to 7 angstrom
# This is based on both visual inspections and RDF curves, where the two first shells lies between ~5 - 14 angstrom



systems=["FA_G-prot","FA_no_G-prot","IA","PA"]
unique_lipids = ['POPC','POPS','POPE','POPA','POPG','CHOL','DPG3', 'DPSM', 'PAP6']
all_lipids_sele = " ".join(['{0:s}'.format(i) for i in unique_lipids ])
lipid_dir = {"FA_G-prot":"POPC, POPS, POPA, POPG, POPE, CHOL, DPG3, DPSM, PAP6",
            "FA_no_G-prot":"POPC, POPS, POPA, POPG, POPE, CHOL, DPG3, DPSM, PAP6",
            "IA":"POPC, POPS, POPA, POPG, POPE, CHOL, DPG3, DPSM, PAP6",
            "PA":"POPC, POPS, POPA, POPG, POPE, CHOL, DPG3, DPSM, PAP6"}




def select_prio(resname,u, priority=["PO4", "GL1", "AM1", "GM1"]):
    """Given a resname a priority selection is made only containing a single bead for the same lipid type"""
    selection=u.select_atoms("resname %s and name PO4 GM1 AM1 ROH GL1" %resname)
    #print (selection)
    if len(selection) != 0:
        if len(selection)!=selection.n_residues:
            beads=np.unique(selection.names)
            #print("%s contains multiple of the listed bead types:" %resname)
            #print("The bead with the highest priority is selected")
            bead_type=priority[np.amin([priority.index(v) for v in beads])]
            #print("Which is %s" %bead_type)
        else:
            beads=np.unique(selection.names)
            bead_type=beads[0]
        new_selection=u.select_atoms("resname %s and name %s" %(resname, bead_type))
        return new_selection
    else:
        return 0




def get_DE_index (s, r, stride):
    """"""
    pdb = "GRO/{0:s}_{1:d}_nowat_mda.pdb".format(s,r) #Has Chain IDs, with R=protein, P=peptide, G=G-protein, S=bilayer
    xtc = "XTC/{0:s}_{1:d}_nowat_res_center_fit.xtc".format(s,r)
    lipids = lipid_dir[s]
    lipid_list = lipids.split(',')
    u = md.Universe(pdb,xtc)
    lipid_tot = np.unique(u.select_atoms('resname {0:s}'.format(all_lipids_sele)).resids).shape[0]
    unique_lipids = 'POPC POPS POPA POPG POPE CHOL DPG3 DPSM PAP6'
    #unique_lipid = ['POPC', 'POPS', 'POPA', 'POPG', 'POPE', 'CHOL', 'DPG3', 'DPSM']

    # Select the protein zone
    
    prot = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'R']
    prot_sele = '(resname {0:s}) and global around 7 group protein'.format(unique_lipids)
    protein_zone = u.select_atoms(prot_sele, protein=prot, updating=True )

    #Select the bulk zone - not overlapping with the protein zone
    bulk_sele = '(resname {0:s}) and not group protein and not group protein_zone'.format(unique_lipids)
    bulk_selection = u.select_atoms(bulk_sele,protein=prot, protein_zone=protein_zone, updating=True)

    P     = np.zeros([len(lipid_list), 3])
    P_tot = np.zeros([len(lipid_list), 3])
    B     = np.zeros([len(lipid_list), 3])

    for ndx, lipid in enumerate(lipid_list):
        print ('lipid:', lipid)
        tot_frames = len(u.trajectory[::stride])

        loc_p     = np.zeros([tot_frames])
        loc_p_tot = np.zeros([tot_frames])
        loc_b     = np.zeros([tot_frames])
        for fdx, ts in enumerate(u.trajectory[::stride]):
            #print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))
            lipid_tot_loc = np.unique(protein_zone.resids).shape[0]
            loc_p_tot[fdx] = lipid_tot_loc
            lipid_circle = select_prio(lipid, protein_zone)
            if type(lipid_circle) != int:
                lipid_circle_num = np.unique(lipid_circle.resids).shape[0]
            else:
                lipid_circle_num = 0
            loc_p[fdx] = lipid_circle_num

            lipid_bulk = select_prio(lipid, bulk_selection)

            if type(lipid_bulk) != int:
                lipid_bulk_num = np.unique(lipid_bulk.resids).shape[0]
            else:
                lipid_bulk_num = 0 
            loc_b[fdx] = lipid_bulk_num
        P[ndx, r]     = np.average(loc_p)
        P_tot[ndx, r] = np.average(loc_p_tot)
        B[ndx, r]     = np.average(loc_b)

    R_cir  = P[:,r] / P_tot[:,r]
    R_bulk = B[:,r] / lipid_tot
    E = R_cir/R_bulk
    return E




DE_All = []
for s in systems:
    DE = []
    for r in range(3):
        print ('Running system {0:s} repeat {1:d}'.format(s,r))
        DE_loc = get_DE_index(s, r, stride=2)
        #np.save('DE_index_{0:s}_{1:d}_last10µs.npy'.format(s,r), np.array(DE_loc))
        DE.append(DE_loc)
        print (np.array(DE_loc).shape)
    np.save('DE_index_{0:s}_full_sim.npy'.format(s), np.array(DE))
    DE_All.append(DE)

DE_ava = np.average(np.array(DE_All), axis=1)
DE_std = np.std(np.array(DE_All), axis=1)

np.save('DE_index_all_full_sim.npy', np.array(DE_All))


lipid_list = ['POPC', 'POPS', 'POPA', 'POPG', 'POPE', 'CHOL', 'DPG3', 'DPSM', 'PAP6']
lipid_list_translated = ['POPC', 'POPS', 'POPA', 'POPG', 'POPE', 'CHOL', 'GM3', 'SM', 'PIP2']
color_list = ['#e66101','#fdb863','#b2abd2','#5e3c99']

for i in range(4):
    plt.errorbar([0,1,2,3,4,5,6,7,8], DE_ava[i], yerr=DE_std[i], ecolor='grey', fmt='o', label=systems[i], color=color_list[i], alpha=0.8)


plt.xticks([0,1,2,3,4,5,6,7,8], lipid_list_translated)
plt.grid(alpha=0.6)
plt.ylabel('DE index')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('DE_index_all_full_sim.png', dpi=300, bbox_inches='tight')



p_values =  np.array([ ttest_1samp(DE[i,:,j], popmean=1).pvalue for i in range(4) for j in range(9) ]).reshape([4,9])

np.savetxt('DE_p_values_full_sim.txt', p_values)
