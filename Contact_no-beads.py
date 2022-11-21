#!/usr/bin/env python
# coding: utf-8



import MDAnalysis as md
import numpy as np
import matplotlib.pyplot as plt
import string
import os
import argparse


parser = argparse.ArgumentParser(description='Script for counting number of beads for each lipid type around the protein and in the bulk')
parser.add_argument('-s', dest='system',nargs='+', action='store', type=str, help='-s for the systems, eg -s "name1" "name2"')

args = parser.parse_args()
system = args.system

print ("systems:")
print (system)

def extract_no_beads_bulk (u):
    """Function to count the number of beads of the different lipid types for a random binding site in the bulk.
    The volume around the BS is approximated by a cylinder with radius = 7 Å and height = 14 Å"""
    prot = u.select_atoms('resid 1-394', updating=True)
    prot_rdf = prot.radius_of_gyration()+20
    memb = u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)')
    protein_complex = u.select_atoms('name BB SC1 SC2 SC3 SC4')
    #z_memb = 2 / (np.max(memb.positions[:,2]) - np.min(memb.positions[:,2]))

    # Select the protein zone
    protein_zone = u.select_atoms('cyzone {0:f} {1:f} {2:f} group protein or group protein_complex'
                  .format(prot_rdf, 60, -60), protein=prot, protein_complex=protein_complex, updating=True)
    #protein_zone = u.select_atoms('cyzone {0:f} {1:f} {2:f} group protein'
    #              .format(prot_rdf, z_memb, -z_memb), protein=prot, updating=True)

    #Select the bulk zone - not overlapping with the protein zone
    bulk_selection = u.select_atoms('group memb and not protein and not group protein_zone'
                                    ,memb=memb, protein=protein_complex, protein_zone=protein_zone, updating=True)
    bulk_res = bulk_selection.resnames
    u, c = np.unique(bulk_res, return_counts=True) 
    #bulk_selection = u.select_atoms('all and not protein and not group protein_zone'
    #                                .format(prot_rdf, z_memb, -z_memb), protein=prot, protein_zone=protein_zone, updating=True)

    # Select a random bead in the bulk selection, as a virtual binding site
    idx_random = np.random.randint(0,bulk_selection.positions.shape[0])
    ran_bs = bulk_selection.indices[idx_random]
    x_ran, y_ran, z_ran = bulk_selection.positions[idx_random]
    #print ('Random site selected')

    # Extract the number of beads around the site
    atom_sel = bulk_selection.select_atoms('prop x == {0:f} and prop y == {1:f} and prop z == {2:f}'.format(x_ran, y_ran, z_ran))    
    try:
        no_beads_bs_all = bulk_selection.select_atoms('cyzone {0:f} {1:f} {2:f} group atom_sel and not protein'.format(7, 7, -7), atom_sel=atom_sel, protein=protein_complex)
        res = no_beads_bs_all.resnames
        unique, counts = np.unique(res, return_counts=True) #Returns the unique lipids found around the BS along with the number of beads for each
    except ValueError:
        print ('bulk selection empty')
        return [], []
        
    return unique, counts




def assign_counts (u, unique, counts):
    """Function for assigning the counts (of beads) for each
    lipid type to an array."""
    memb = u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)').resnames
    unique_lipid, counts_lipid = np.unique(memb, return_counts=True)
    zeros = np.zeros([unique_lipid.shape[0]])
    counts_out = np.vstack((unique_lipid, zeros))
    
    if len(counts) != 0:
        idxs = [ np.where(counts_out[0,:] == i)[0] for i in unique ]
        shape = len(idxs)
        try:
            counts_out[1,idxs] = counts.reshape([shape,1])
        except IndexError:
            print ('indexerror in assign_counts function')
            print (unique, counts)
    return counts_out 



def main_func (u, stride=100, sel='resid 111-394'):
    """Main function
    Iterative over the systems, calling only this function."""
    ### MAIN ###
    
    #Set selections
    prot = u.select_atoms('resid 1-394', updating=True)
    protein_complex = u.select_atoms('name BB SC1 SC2 SC3 SC4')
    prot_rdf = prot.radius_of_gyration()
    #print ('protein rdf', prot_rdf)
    membrane = u.select_atoms('all and not (name BB SC1 SC2 SC3 SC4)', updating=True)
    #z_memb = 2 / (np.max(membrane.positions[:,2]) - np.min(membrane.positions[:,2]))

    protein_zone = u.select_atoms('cyzone {0:f} {1:f} {2:f} group protein'.format(prot_rdf, 60, -60), protein=prot, updating=True)
    #protein_zone = u.select_atoms('cyzone {0:f} {1:f} {2:f} group protein'.format(prot_rdf, z_memb, -z_memb), protein=prot, updating=True)
    memb_prot_zone = protein_zone.select_atoms('all and not group protein', protein=protein_complex, updating=True)
    bulk_selection = u.select_atoms('all and not group protein and not group protein_zone'
                                    .format(prot_rdf, 60, -60), protein=protein_complex, protein_zone=protein_zone, updating=True)
    #bulk_selection = u.select_atoms('all and not group protein and not group protein_zone'
    #                                .format(prot_rdf, z_memb, -z_memb), protein=prot, protein_zone=protein_zone, updating=True)
    TMD = u.select_atoms('group protein and {0:s}'.format(sel), protein=prot)
    TMD_resid = np.unique(TMD.resids)
    print ('Selection set')

    output_counts_all = []
    output_counts_all_bulk = []
        
    for idx, ts in enumerate(u.trajectory[::stride]):
        print (ts.frame)
        counts_TMD_local = []
        unique_bulk, counts_bulk = extract_no_beads_bulk (u) 
        if len(unique_bulk) != 0:
            out_bulk_local = assign_counts(u, unique_bulk, counts_bulk)
        else:
            continue
        output_counts_all_bulk.append(out_bulk_local)
        for indx, r in enumerate(TMD_resid):
            sel = u.select_atoms('(cyzone {0:f} {1:f} {2:f} resid {3:d}) and not group protein'
                                                .format(7, 7,-7, r), protein=protein_complex)
            sel_res = sel.resnames        
            unique, counts = np.unique(sel_res, return_counts=True)        
            #idxs = [ np.where(counts_out[0,:] == i)[0] for i in unique ]
            out_local = assign_counts(u, unique, counts)

            counts_TMD_local.append(out_local)
        output_counts_all.append(np.array(counts_TMD_local))

    COUNTS_PROT = np.array(output_counts_all) # [nframes, nresids, 2, nlipid_types]
    COUNTS_VIRTUAL_BS = np.array(output_counts_all_bulk) # [nframes, 2, nlipid_types]
    return COUNTS_PROT, COUNTS_VIRTUAL_BS



for idx, s in enumerate(system):
    for r in range(3):
        if s=='FA_G-prot':
            gro = 'GRO/{0:s}_{1:d}_nowat_mda_renumbered.pdb'.format(s,r)
        else:
            gro = 'GRO/{0:s}_{1:d}_nowat.gro'.format(s,r)
        xtc = 'XTC/{0:s}_{1:d}_nowat_res_center_fit_last5.xtc'.format(s,r)
        print ("Loaded Universe with: GRO/{0:s}_{1:d}_nowat.gro XTC/{0:s}_{1:d}_nowat_res_center_last5.xtc".format(s,r))
        u = md.Universe(gro,xtc)
        print ("Universe loaded with GRO/{0:s}_{1:d}_nowat.gro XTC/{0:s}_{1:d}_nowat_res_center_last5.xtc".format(s,r))
        if os.path.isdir("CONTACT_NO_BEADS") == False:
                os.system("mkdir CONTACT_NO_BEADS")
                print ("Output will be saved in the directory called CONTACT_NO_BEADS")
        else:
                print ("Output will be saved in the directory called CONTACT_NO_BEADS")
        output_dat1 = 'CONTACT_NO_BEADS/{0:s}_{1:d}_number_of_beads_protein.npy'.format(s,r)
        output_dat2 = 'CONTACT_NO_BEADS/{0:s}_{1:d}_number_of_beads_bulk.npy'.format(s,r)
        if os.path.exists(output_dat1):
            print ("{0:s} exists".format(output_dat1))
            continue
        else:
            print ("Performing calculation")
            COUNTS_PROT, COUNTS_VIRTUAL_BS = main_func (u)
            np.save(output_dat1, COUNTS_PROT)
            np.save(output_dat2, COUNTS_VIRTUAL_BS)

print ('Calculation Done!!')











