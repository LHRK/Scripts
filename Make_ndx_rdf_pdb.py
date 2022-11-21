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

#prot = u.select_atoms('protein')
lipids = np.unique(u.select_atoms('resname POP*').resnames)

sele = []

for l in lipids:
	sele.append(u.select_atoms('resname {0:s} and name PO4'.format(l)))

GCGR = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'R']
TMD = u.select_atoms('group GCGR and resid 111-393', GCGR=GCGR)
TMD_BB = u.select_atoms('group GCGR and resid 111-393 and name BB', GCGR=GCGR)
GCGR_BB = GCGR.select_atoms('name BB')
GCG  = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'P']
G_prot = u.select_atoms('name BB SC1 SC2 SC3 SC4')[u.select_atoms('name BB SC1 SC2 SC3 SC4').chainIDs == 'G']
#prot = u.select_atoms('protein')
BB   = u.select_atoms('name BB')
All_headgroups = u.select_atoms('name P*')
sel_resnames = " ".join(['{0:s}'.format(i) for i in np.unique(All_headgroups.resnames)])
chol = u.select_atoms("resname CHOL and name ROH")
pap = u.select_atoms("resname PAP6 and name PO4")
gm3 = u.select_atoms("resname DPG3 and name AM1")
sm  = u.select_atoms("resname DPSM and name PO4")
systems = u.select_atoms('all')

with md.selections.gromacs.SelectionWriter('NDX/{0:s}_rdf_chains.ndx'.format(name), mode='w') as ndx:
	ndx.write(chol, name='CHOL')
	ndx.write(pap, name='PAP6')
	ndx.write(gm3, name='DPGM3')
	ndx.write(sm, name='DPSM')
	ndx.write(GCGR, name='GCGR')
	ndx.write(TMD, name='TMD')
	ndx.write(TMD_BB, name='TMD_BB')
	ndx.write(GCGR_BB, name='GCGR_BB')
	ndx.write(GCG, name='GCG')
	ndx.write(G_prot, name='G_protein')
	ndx.write(systems, name='System')
	for idx, s in enumerate(sele):
		print (lipids[idx])
		ndx.write(s, name=lipids[idx])


