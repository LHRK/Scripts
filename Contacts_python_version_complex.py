#!/usr/bin/env python

import MDAnalysis as md
import numpy as np
import matplotlib.pyplot as plt
import string
import os



# ## Calculating the number of contacts between the different types of lipids and each protein residue


def extract_contacts (u, lipid, stride=20):
    """Timestep = 500 ps, so striding every 20 frame,
    means taking every 10 ns"""
    GCGR = u.select_atoms('resid 26-421')
    GCGR_BB = GCGR.select_atoms('name BB')
    GCGR_resid = np.unique(GCGR.resids)
    tot_frames = len(u.trajectory[::stride])
    num_resids = GCGR_resid.shape[0]    
    GCGR_resid = np.unique(GCGR.resids)

    C = np.zeros([num_resids, tot_frames])
    for idx, ts in enumerate(u.trajectory[::stride]):
        for indx, r in enumerate(GCGR_resid):
            sel = u.select_atoms('resname {0:s} and global around 7 resid {1:d}'.format(lipid, r))
            num_contacts = np.unique(sel.resids).shape[0]
            C[indx, idx] = num_contacts
    sel = u.select_atoms('resname {0:s}'.format(lipid))
    num_lipids_total = np.unique(sel.resids).shape[0]
    ava_con = (np.average(C, axis=1) / num_lipids_total)
    tot_con = (np.sum(C, axis=1) / num_lipids_total)
    return C, ava_con, tot_con



def residence_time (u, lipid, stride=20):
    GCGR = u.select_atoms('resid 26-421')
    GCGR_BB = GCGR.select_atoms('name BB')
    GCGR_resid = np.unique(GCGR.resids)
    num_resids = GCGR_resid.shape[0]    
    GCGR_resid = np.unique(GCGR.resids)

    C = np.zeros([num_resids])
    for rdx, r in enumerate(GCGR_resid):
        C_res = []
        count = 0
        for fdx, ts in enumerate(u.trajectory[::stride]):
            sel = u.select_atoms('resname {0:s} and global around 7 resid {1:d}'.format(lipid, r))
            num_contacts = np.unique(sel.resids).shape[0]
            if num_contacts > 0 :
                count = count +1
            else:
                C_res.append(count)
                count = 0
        C[rdx] = (np.average(C_res))        
    return C


#systems=["FA_G-prot", "FA_no_G-prot", "IA", "PA"]
systems=[ "IA", "PA", "FA_G-prot", "FA_no_G-prot"]

lipid_dir = {"FA_G-prot":"POPC, POPS, POPE, POPA, POPG, PAP6, CHOL, DPSM, DPG3",
             "FA_no_G-prot":"POPC, POPS, POPE, POPA, POPG, PAP6, CHOL, DPSM, DPG3",
             "IA":"POPC, POPS, POPE, POPA, POPG, PAP6, CHOL, DPSM, DPG3",
             "PA":"POPC, POPS, POPE, POPA, POPG, PAP6, CHOL, DPSM, DPG3"}


#for idx, s in enumerate(systems):
#    lipids = lipid_dir[s]
#    for r in range(3):
#        gro = 'GRO/{0:s}_{1:d}_nowat_renumbered.gro'.format(s,r)
#        xtc = 'XTC/{0:s}_{1:d}_nowat_res_center.xtc'.format(s,r)
#        u = md.Universe(gro,xtc)
#        for lipid in lipids.split(','):
#            #print ('lipid:', lipid)
#            llipid = lipid.translate({ord(c): None for c in string.whitespace})
#            print ('lipid:', llipid)
#            output_dat = 'CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_ava_consecutive_contact_time.txt'.format(s,r,llipid)
#            if os.path.exists(output_dat):
#                print ("{0:s} exists".format(output_dat))
#                continue
#            else:
#                C = residence_time(u, llipid)
#                np.savetxt('CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_ava_consecutive_contact_time.txt'.format(s,r,llipid), C)
#



for idx, s in enumerate(systems):
    lipids = lipid_dir[s]
    for r in range(3):
        gro = 'GRO/{0:s}_{1:d}_nowat_renumbered.gro'.format(s,r)
        xtc = 'XTC/{0:s}_{1:d}_nowat_res_center_fit_last5.xtc'.format(s,r)
        u = md.Universe(gro,xtc)
        for lipid in lipids.split(','):
            #print ('lipid:', lipid)
            llipid = lipid.translate({ord(c): None for c in string.whitespace})
            print ('lipid:', llipid)
            output_dat = 'CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_number_of_contact_per_resid.txt'.format(s,r,llipid)
            if os.path.exists(output_dat):
                print ("{0:s} exists".format(output_dat))
                continue
            else:
                C, ava_con, tot_con = extract_contacts (u, llipid)
                np.savetxt('CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_number_of_contact_per_resid.txt'.format(s,r,llipid), C)
                np.savetxt('CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_Ava_number_of_contact_per_resid.txt'.format(s,r,llipid), ava_con)
                np.savetxt('CONTACT_PYTHON/{0:s}_{1:d}_{2:s}_total_number_of_contact_per_resid.txt'.format(s,r,llipid), tot_con)
