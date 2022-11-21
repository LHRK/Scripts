#!/usr/bin/env python

import numpy as np
import argparse
import MDAnalysis as md
from MDAnalysis.analysis.leaflet import LeafletFinder as LF


parser = argparse.ArgumentParser(description='Script for making index files')
parser.add_argument('-n', dest='name', action='store', type=str, help='-n for name of the system (without nowat, its the raw gro file)')

args = parser.parse_args()
name = args.name

u = md.Universe('GRO/{0:s}_mda.pdb'.format(name))


L = LF(u, 'name P*', cutoff=20.0)
leaflet0 = L.groups(0).residues.atoms
leaflet1 = L.groups(1).residues.atoms

memb = u.select_atoms('resname POP*')

lipids = np.unique(u.select_atoms('resname POP*').resnames)

sele = []

for l in lipids:
	sele.append(u.select_atoms('resname {0:s} and name PO4'.format(l)))



GCGR = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'R']
GCGR_BB = GCGR.select_atoms('name BB')
GCG  = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'P']
G_prot = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'G']
prot = u.select_atoms('protein')
BB   = u.select_atoms('name BB')
All_headgroups = u.select_atoms('name P*')
Po4_headgroup  = u.select_atoms('name PO4')
sel_resnames = " ".join(['{0:s}'.format(i) for i in np.unique(All_headgroups.resnames)])
all_lipids = u.select_atoms('resname {0:s} CHOL PAP6 DPSM DPG3'.format(sel_resnames))
#tmd = u.select_atoms('resid 111-394 and name BB')


with md.selections.gromacs.SelectionWriter('NDX/{0:s}_selection.ndx'.format(name), mode='w') as ndx:
	ndx.write(leaflet0, name='Top_leaflet')
	ndx.write(leaflet1, name='Bot_leaflet')
	ndx.write(memb, name='membrane')
	ndx.write(prot, name='Protein')
	ndx.write(GCGR, name='GCGR')
	ndx.write(GCG, name='GCG')
	ndx.write(G_prot, name='G-protein')
	ndx.write(BB, name = 'BB')
	ndx.write(All_headgroups, name='headgroups')
	ndx.write(Po4_headgroup, name='headgroups_PO4')
	ndx.write(all_lipids, name='All_lipids')
	#ndx.write(tmd, name='TMD')
	for idx, s in enumerate(sele):
		print ('Headgroup_'+lipids[idx])
		ndx.write(s, name='Headgroup_'+lipids[idx])


